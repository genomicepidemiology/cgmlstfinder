#!/usr/bin/env python3 
import os, sys, shutil, argparse, subprocess, pickle, re, gzip, time
from difflib import ndiff

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class SeqFile():
    ''' '''
    def __init__(self, seqfile, pe_file_reverse=None, phred_format=None):
        ''' Constructor.
        '''
        self.phred = phred_format

        seqfile = os.path.abspath(seqfile)
        if(not os.path.isfile(seqfile)):
            print("File not found: " + seqfile)
            quit(1)
        self.path = seqfile
        self.pe_file_reverse = None

        self.filename = SeqFile.get_read_filename(seqfile)
        self.filename_reverse = None

        self.trim_path = None
        self.trim_pe_file_reverse = None

        self.gzipped = SeqFile.is_gzipped(seqfile)

        if(pe_file_reverse):
            self.seq_format = "paired"

            # Add path to reverse pair file
            self.pe_file_reverse = os.path.abspath(pe_file_reverse)
            if(not os.path.isfile(pe_file_reverse)):
                print("Reverse pair file not found: \"" + pe_file_reverse
                      + "\"")
                quit(1)
            self.path_reverse = pe_file_reverse

            self.filename_reverse = SeqFile.get_read_filename(
                self.pe_file_reverse)

            # Check if pair is gzipped
            if(self.gzipped != SeqFile.is_gzipped(pe_file_reverse)):
                print("ERROR: It seems that only one of the read pair files is\
                      gzipped.")
                quit(1)
        elif(phred_format):
            self.seq_format = "single"
        else:
            self.seq_format = "assembly"

    def set_trim_files(self, path, path_reverse=None):
        self.trim_path = path
        self.trim_filename = SeqFile.get_read_filename(path)
        if(path_reverse):
            self.trim_pe_file_reverse = path_reverse
            self.trim_filename_reverse = SeqFile.get_read_filename(
                path_reverse)

    @staticmethod
    def get_read_filename(seq_path):
        ''' Removes path from given string and removes extensions:
            .fq .fastq .gz and .trim
        '''
        seq_path = os.path.basename(seq_path)
        seq_path = seq_path.replace(".fq", "")
        seq_path = seq_path.replace(".fastq", "")
        seq_path = seq_path.replace(".gz", "")
        seq_path = seq_path.replace(".trim", "")
        return seq_path.rstrip()

    @staticmethod
    def is_gzipped(file_path):
        ''' Returns True if file is gzipped and False otherwise.

            The result is inferred from the first two bits in the file read
            from the input path.
            On unix systems this should be: 1f 8b
            Theoretically there could be exceptions to this test but it is
            unlikely and impossible if the input files are otherwise expected
            to be encoded in utf-8.
        '''
        with open(file_path, mode='rb') as fh:
            bit_start = fh.read(2)
        if(bit_start == b'\x1f\x8b'):
            return True
        else:
            return False

    @staticmethod
    def group_fastqs(file_paths):
        '''
        '''

        re_filename = re.compile(r"(.+)_S\d+_L(\d+)(_.+)")
        re_seqdb_filename = re.compile(r"([E|S]RR\d+)(_\d+?\..+)")
        file_groups = {}

        for path in file_paths:

            filename = os.path.basename(path)
            match_filename = re_filename.search(filename)
            match_seqdb = re_seqdb_filename.search(filename)

            if(match_filename):
                name = match_filename.group(1)
                lane_no = match_filename.group(2)
                pair_id = match_filename.group(3)
                lane_no = int(lane_no)
                lane_nos = file_groups.get((name, pair_id), {})
                lane_nos[lane_no] = path
                file_groups[(name, pair_id)] = lane_nos
            elif(match_seqdb):
                # Assume that no lane number exist and are just provided an
                # artificial lane number.
                name = match_seqdb.group(1)
                lane_no = 1
                pair_id = match_seqdb.group(2)
                lane_nos = file_groups.get((name, pair_id), {})
                lane_nos[lane_no] = path
                file_groups[(name, pair_id)] = lane_nos
            else:
                eprint("Warning: Did not recognise filename: " + filename)

        return file_groups

    @staticmethod
    def concat_fastqs(file_paths, out_path=".", verbose=False):
        '''
        '''

        out_list = []
        file_groups = SeqFile.group_fastqs(file_paths)

        for (name, pair_id), lane_nos in file_groups.items():

            out_filename = SeqFile.get_read_filename(name + pair_id) + ".fq"
            sorted_lanes = sorted(lane_nos.keys())

            for lane_no in sorted_lanes:
                path = lane_nos[lane_no]
                if(SeqFile.is_gzipped(path)):
                    cmd = "gunzip -c "
                else:
                    cmd = "cat "
                cmd += path + " >> " + out_path + "/" + out_filename
                subprocess.run(cmd, shell=True)
            subprocess.run("gzip " + out_path + "/" + out_filename, shell=True)
            out_list.append(out_path + "/" + out_filename + ".gz")
            if(verbose):
                eprint("Wrote: " + out_path + "/" + out_filename + ".gz")

        return out_list

    @classmethod
    def parse_files(cls, file_paths, phred=None, headers2count=10, min_match=2,
                    force_neighbour=False):
        '''
        '''
        re_win_newline = re.compile(r"\r\n")
        re_mac_newline = re.compile(r"\r")

        re_fastq_header = re.compile(r"(@.+)")

        re_sra = re.compile(r"^(@SRR\d+\.\d+)")
        re_ena = re.compile(r"^(@ERR\d+\.\d+)")

        re_pair_no = re.compile(r"1|2")

        paired = {}
        single = {}
        old_read_headers = {}
        prev_read_headers = {}
        old_org_headers = {}
        prev_org_headers = {}

        for path in file_paths:
            head_count = 0
            is_fastq = False
            file_headers = []

            if(cls.is_gzipped(path)):
                fh = gzip.open(path, "rt", encoding="utf-8")
            else:
                fh = open(path, "r", encoding="utf-8")

            org_headers = {}

            for line in fh:
                line = re_win_newline.sub("\n", line)
                line = re_mac_newline.sub("\n", line)
                line = line.rstrip()

                line_1s_letter = line[0]

                # FASTA format
                if(line_1s_letter == ">"):
                    pass
                # FASTQ format
                elif(line_1s_letter == "@"):
                    # match_fastq_header = re_fastq_header.search(line)
                    # fastq_header = match_fastq_header.group(1)
                    match_sra = re_sra.search(line)
                    match_ena = re_ena.search(line)

                    is_fastq = True
                    head_count += 1

                    # If data is obtained from SRA
                    if(match_sra):
                        header = match_sra.group(1)
                        org_headers[header] = line
                    # If data is obtained from ENA
                    elif(match_ena):
                        header = match_ena.group(1)
                        org_headers[header] = line
                    else:
                        # Masking 1s and 2s so that pairs will match
                        header = re_pair_no.sub("x", line)
                        org_headers[header] = line

                    file_headers.append(header)
                    if(head_count == headers2count):
                        break

            fh.close()

            if(is_fastq):
                if(not phred):
                    # TODO: Implement find phred function
                    pass

                matches = 0
                read_file1 = ""
                read_file2 = ""
                # Check for mates in "neighbors"
                for header in file_headers:
                    if(header in prev_read_headers):
                        if(not match_sra):
                            pair_no = cls.detect_pair_no(
                                org_headers[header], prev_org_headers[header])
                        else:
                            # Correct pair numbering cannot be obtained
                            # from SRA headers.
                            # Attempt to get pair number from filenames
                            filename1 = cls.get_read_filename(path)
                            filename2 = cls.get_read_filename(
                                prev_read_headers[header])
                            pair_no = cls.detect_pair_no(filename1, filename2)

                        if(pair_no is not None):
                            matches += 1
                            if(pair_no == 2):
                                read_file1 = prev_read_headers[header]
                            else:
                                read_file1 = path
                                read_file2 = prev_read_headers[header]

                if(matches >= min_match):
                    if(read_file2):
                        paired[read_file1] = read_file2
                        del single[read_file2]
                    else:
                        paired[read_file1] = path
                elif(not force_neighbour):
                    matches = 0
                    read_file1 = ""
                    read_file2 = ""

                    for header in file_headers:
                        if(header in old_read_headers):
                            if(not match_sra):
                                pair_no = cls.detect_pair_no(
                                    org_headers[header],
                                    old_org_headers[header]
                                )
                            else:
                                # Correct pair numbering cannot be obtained
                                # from SRA headers.
                                # Attempt to get pair number from filenames
                                filename1 = cls.get_read_filename(path)
                                filename2 = cls.get_read_filename(
                                    old_read_headers[header])

                                pair_no = cls.detect_pair_no(filename1,
                                                             filename2)

                            if(pair_no is not None):
                                matches += 1
                                if(pair_no == 2):
                                    read_file1 = old_read_headers[header]
                                else:
                                    read_file1 = path
                                    read_file2 = old_read_headers[header]

                            # Check if there are more than one match
                            if(read_file1 and read_file1 in paired):
                                print("Header matches multiple read files.\n\
                                       Header: " + header + "\n\
                                       File 1: " + read_file1 + "\n\
                                       File 2: " + paired[read_file1] + "\n\
                                       File 3: " + path + "\n\
                                       DONE!")
                                quit(1)

                    if(matches >= min_match):
                        if(read_file2):
                            paired[read_file1] = read_file2
                            del single[read_file2]
                        else:
                            paired[read_file1] = path
                    else:
                        single[path] = 1

                else:
                    single[path] = 1

                # Moves neighbor headers to old headers
                for header in prev_read_headers.keys():
                    old_read_headers[header] = prev_read_headers[header]
                    old_org_headers[header] = prev_org_headers[header]

                # Resets neighbors
                prev_read_headers = {}
                prev_org_headers = {}
                for header in file_headers:
                    prev_read_headers[header] = path
                    prev_org_headers[header] = org_headers[header]

        output_seq_files = []

        for path in file_paths:
            if(path in paired):
                seq_file = SeqFile(path, pe_file_reverse=paired[path],
                                   phred_format=phred)
                output_seq_files.append(seq_file)
            elif(path in single):
                seq_file = SeqFile(path, phred_format=phred)
                print("Created seqfile: " + seq_file.filename)
                output_seq_files.append(seq_file)

        return output_seq_files

    @staticmethod
    def detect_pair_no(header1, header2):
        ''' Given two fastq headers, will output if header1 is either 1 or 2.
            If the headers are not does not match, method will return None
        '''
        head_diff = ndiff(header1, header2)
        mismatches = 0
        pair_no = None

        for s in head_diff:
            if(s[0] == "-"):
                if(s[2] == "2"):
                    pair_no = 2
                    mismatches += 1
                elif(s[2] == "1"):
                    pair_no = 1
                    mismatches += 1
                else:
                    mismatches += 1

        if(mismatches < 2):
            return pair_no
        else:
            return None

    @staticmethod
    def load_seq_files(file_paths):
        ''' Given a list of file paths, returns a list of SeqFile objects.
        '''
        file_paths_str = " ".join(file_paths)

        # Running parse_input
        # TODO: Should be rewritten as proper python class.
        parse_input_cmd = ("perl parse_input.pl " + file_paths_str
                           + " > parse_input.output.txt")
        print("PARSE CMD: " + parse_input_cmd)
        try:
            subprocess.check_call([parse_input_cmd], shell=True)
        except subprocess.CalledProcessError:
            print("ERROR: parse input call failed")
            print("CMD that failed: " + parse_input_cmd)
            quit(1)

        file_list = []

        with open("parse_input.output.txt", "r", encoding="utf-8") as input_fh:
            for line in input_fh:
                entries = line.split("\t")
                if(entries[0] == "paired"):
                    print("PAIRED")
                    file_list.append(
                        SeqFile(seqfile=entries[1].strip(),
                                pe_file_reverse=entries[3].strip(),
                                phred_format=entries[2].strip()))
                elif(entries[0] == "single"):
                    print("SINGLE")
                    file_list.append(SeqFile(seqfile=entries[1].strip(),
                                             phred_format=entries[2].strip()))
                elif(entries[0] == "assembled"):
                    print("ASSEMBLED")
                    file_list.append(SeqFile(seqfile=entries[1].strip()))

        return file_list

class KMA():

    def __init__(self, seqfile, tmp_dir, db, gene_list, kma_path, fasta = False):
        """ Constructor map reads from seqfile object using kma.
        """
        # Create kma command line list
        kma_call_list = [kma_path, "-i"]
        if fasta == False:
            filename = seqfile.filename
            kma_call_list.append(seqfile.path)
            # Add reverse reads if paired-end data
            if(seqfile.pe_file_reverse):
                kma_call_list.append(seqfile.pe_file_reverse)
        else: 
            filename = seqfile.split('/')[-1].split('.')[0]
            kma_call_list.append(seqfile)

        result_file_tmp = tmp_dir + "/kma_" + filename
        self.filename = filename
        self.result_file = result_file_tmp + ".res"
        self.seqfile = seqfile

        kma_call_list += [
            "-o", result_file_tmp,
            "-t_db", db,
            "-mem_mode", "-dense", "-boot", "-1t1", "-and"]

        # Call kma externally
        eprint("# KMA call: " + " ".join(kma_call_list))
        process = subprocess.Popen(kma_call_list, shell=False, stdout=subprocess.PIPE) #, stderr=subprocess.PIPE)
        out, err = process.communicate()
        eprint("KMA call ended")


    def best_allel_hits(self):
        """ 
        Extracts perfect matching allel hits from kma results file and returns
        a list(string) found allel profile ordered based on the gene list. 
        """

        best_alleles = {}
        
        # Create dict of locus and allel with the highest quality score
        with open(self.result_file, "r") as result_file:
            header = result_file.readline()
            header = header.strip().split("\t")
            template_id_index = header.index("Template_Identity") 
            loci_allel = re.compile(r"(\S+)_(\d+)")
            i = 0
            for line in result_file:
                i += 1
                data = line.rstrip().split("\t")
                loci_allel_object = loci_allel.search(line)
                locus = loci_allel_object.group(1)
                allel = loci_allel_object.group(2)
                template_id = float(data[template_id_index])
                # Check for perfect matches
                if template_id == 100:
                    best_alleles[locus] = allel

        # Get called alleles
        allel_str = self.filename
        for locus in gene_list:
            locus = locus.strip()
            if locus in best_alleles:
                 allel_str += "\t%s" %(best_alleles[locus])
            else:
                 allel_str += "\tNaN"
        return [allel_str]


def st_typing(pickle_path, inp):
    """
    Takes the path to a pickled dictionary, the inp list of the allel 
    number that each loci has been assigned, and an output file string
    where the found st type and similaity is written into it.  
    """

    # Find best ST type for all allel profiles
    st_output = ""

    # First line contains matrix column headers, which are the specific loci
    loci = inp[0].strip().split("\t")

    for sample_str in inp[1:]:
        sample = sample_str.strip().split("\t")
        sample_name = sample[0]
        st_hits = []
        for i in range(1, len(sample)):
            allel = sample[i].encode('utf-8')
            locus = loci[i].encode('utf-8')
            # Loci/Allel combination may not be found in the large profile file
            st_hits += loci_allel_dict[locus].get(allel, ["None"])

        # Find most frequent st_type in st_hits
        score = {}
        max_count = 0
        best_hit = b""
        for hit in st_hits:
            if hit in score:
                score[hit] += 1
                if max_count < score[hit]:
                    max_count = score[hit]
                    best_hit = hit
            elif(hit is not "None"):
                score[hit] = 1

        # Prepare output string
        similarity = round(max_count/(len(loci) - 1)*100, 2)
        st_output += "%s\t%s\t%d\t%.2f\n"%(sample_name, best_hit.decode("utf-8"), max_count, similarity)
        #st_output += "%s\t%s\t%d\t%.2f\n"%(sample_name, best_hit, max_count, similarity)

    return st_output


def file_format(input_files):
    """
    Takes all input files and checks their first character to assess
    the file format. 3 lists are return 1 list containing all fasta files
    1 containing all fastq files and 1 containing all invalid files
    """
    fasta_files = []
    fastq_files = []
    invalid_files = []
    # Open all input files and get the first character
    for infile in input_files:
        try:
            f =  gzip.open(infile, "rb")
            fst_char = f.read(1)
        except OSError:
            f = open(infile, "rb")
            fst_char = f.read(1)
        f.close()
        #fst_char = f.readline().decode("ascii")[0]
        #print(fst_char)
        # Return file format based in first char
        if fst_char == b'@':
            fastq_files.append(infile)
        elif fst_char == b'>':
            fasta_files.append(infile)
        else:
            invalid_files.append(infile)
    return (fasta_files, fastq_files, invalid_files)  


if __name__ == '__main__':
    #
    # Handling arguments
    #
    parser = argparse.ArgumentParser(description="")
    # Posotional arguments
    parser.add_argument("input",
                        help="FASTQ files to do cgMLST on.",
                        nargs="+",
                        metavar="FASTQ",
                        default=None)
    # Optional arguments
    parser.add_argument("-o", "--output",
                        help="Output file.",
                        default="AlleleMatrix.mat",
                        metavar="OUTPUT_FILE")
    parser.add_argument("-s", "--species_scheme",
                        help="species schemes to apply, e.g. ecoli_cgMLST. Must match the name of the species database",
                        #choices=['ecoli_cgMLST', 'yersinia_cgMLST', 'campy_cgMLST', 
                        #         'salmonella_cgMLST_v2', 'listeria_cgMLST'],
                        default=None,
                        metavar="SPECIES_SCHEME")
    parser.add_argument("-db", "--databases",
                        help="Directory containing the databases and gene\
                              lists for each species_scheme.",
                        default=("/home/data1/services/cgMLSTFinder/"
                                 "database_current"),
                        metavar="DB_DIR")
    parser.add_argument("-t", "--tmp_dir",
                        help="Temporary directory for storage of the results\
                              from the external software.",
                        default="cgMLST_tmp_dir")
    parser.add_argument("-k", "--kmapath",
                        help="Path to executable kma program",
                        default="kma",
                        metavar="KMA_PATH")

    args = parser.parse_args()

    # Handle tmp dir
    args.tmp_dir = os.path.abspath(args.tmp_dir)
    os.makedirs(args.tmp_dir, exist_ok=True)

    # Species scheme database
    species_scheme = args.species_scheme
    db_dir = args.databases + "/" + species_scheme
    db_species_scheme = db_dir + "/" + species_scheme

    # Check kma path
    if shutil.which(args.kmapath) is None:
        eprint("The path to kma, '%s', is not executable, append kma to"
               "$PATH"%(args.kmapath))
        quit(1)
    kma_path = args.kmapath

    # Test if database is found and indexed (works for both kma-1.0 and kma-2.0)
    db_files = [species_scheme + ".length.b"]#,
                #species_scheme + ".name", species_scheme + ".index.b", 
                #species_scheme + ".seq.b", species_scheme + ".comp.b" ]

    for db_file in db_files:
        if(not os.path.isfile(db_dir + "/" + db_file)):
            eprint("ERROR: A KMA index file seems to be missing from the"
                   "database directory. You may need to run kma_index.\n"
                   "Missing file: " + db_dir + "/" + db_file)
            quit(1)

    # Gene list
    gene_list_filename = (db_dir + "/list_genes.txt")
    if(not os.path.isfile(gene_list_filename)):
        eprint("Gene list not found at expected location: %s"%(gene_list_filename))
        quit(1)

    # Check file format of input files (fasta or fastq, gz or not gz)
    (fasta_files, fastq_files, invalid_files) = file_format(args.input)
    eprint("Input files: %d fasta file(s)\n"
           "             %d fastq file(s)\n"
           "             %d invalid file(s)"%(len(fasta_files), 
                                             len(fastq_files), 
                                             len(invalid_files)))

    # Load files and pair them if necessary
    eprint("Parsing files:" + str(args.input))
    fastq_files = SeqFile.parse_files(fastq_files, phred=33)

    # Get gene_list_file into list
    try:
        gene_list_file = open(gene_list_filename, "r")
    except IOError:
        eprint("NO valid genefile found")
        sys.exit(-1)
    gene_list = [locus.strip() for locus in  gene_list_file.readlines()]
    gene_list_file.close()

    # Write header to output file
    allel_output = ["Genome\t%s" %("\t".join(gene_list))]
    st_output = "Sample_Names\tcgST_Assigned\tNo_of_Allels_Found\tSimilarity\n"

    # Load ST-dict pickle
    pickle_path = db_dir + "/%s_profile.p"%(species_scheme)

    T0 = time.time()
    if os.path.isfile(pickle_path):
        try:
            loci_allel_dict = pickle.load(open(pickle_path, "rb"))
            T1 = time.time()
            eprint("pickle_loaded: %d s"%(int(T1-T0)) )
        except IOError:
            sys.stdout.write("Error, pickle not found", pickle_path)
            quit(1)

    for seqfile in fasta_files:
        # Run KMA to find alleles from fasta file
        seq_kma = KMA(seqfile, args.tmp_dir, db_species_scheme, gene_list, kma_path, fasta = True)
        # Get called allelel
        allel_output += seq_kma.best_allel_hits()

    for seqfile in fastq_files:
        # Run KMA to find alleles from fastq file
        seq_kma = KMA(seqfile, args.tmp_dir, db_species_scheme, gene_list, kma_path, fasta = False)
        # Get called allelel
        allel_output += seq_kma.best_allel_hits()

    # Create ST-type file if pickle containing profile list exist   
    st_filename = args.output + "-st.txt"
    if os.path.isfile(pickle_path):
        # Write header in output file
        st_output += st_typing(loci_allel_dict, allel_output)

        # Write ST-type output
        with open(st_filename, "w") as fh:
            fh.write(st_output)
        print(st_output)

    # Write allel matrix output
    with open(args.output + ".txt", "w") as fh:
        fh.write("\n".join(allel_output) + "\n")

    eprint("Done")
