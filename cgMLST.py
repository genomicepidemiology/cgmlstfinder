#!/usr/bin/env python3 

from __future__ import print_function
import os
import sys
import shutil
import argparse
import subprocess
import pickle
import re
import gzip
import time
from python_module_seqfilehandler.SeqFileHandler import SeqFile


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


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
            "-mem_mode", "-dense", "-boot", "-1t1"]

        # Call kma externally
        eprint("# KMA call: " + " ".join(kma_call_list))
        process = subprocess.Popen(kma_call_list, shell=False, stdout=subprocess.PIPE) #, stderr=subprocess.PIPE)
        out, err = process.communicate()
        eprint("KMA call ended")


    def best_allel_hits(self):
        """ 
        Extracts best allel hits from kma and returns a list(string) of best
        allel ordered based on the gene list. 
        """
        best_alleles = {}
        
        # Create dict of locus and allel with the highest quality score
        with open(self.result_file, "r") as result_file:
            header = result_file.readline()
            header = header.strip().split("\t")
            q_val_index = header.index("q_value")
            loci_allel = re.compile(r"(\S+)_(\d+)")
            i = 0
            for line in result_file:
                i += 1
                data = line.rstrip().split("\t")
                loci_allel_object = loci_allel.search(line)
                locus = loci_allel_object.group(1)
                allel = loci_allel_object.group(2)
                q_score = float(data[q_val_index])
                if locus in best_alleles:
                    # Check if allel has a higher q_score than the saved allel
                    if best_alleles[locus][1] < q_score:
                        best_alleles[locus] = [allel, q_score]
                else:
                    best_alleles[locus] = [allel, q_score]

        # Get called alleles
        allel_str = self.filename
        for locus in gene_list:
            locus = locus.strip()
            if locus in best_alleles:
                 allel_str += "\t%s" %(best_alleles[locus][0])
            else:
                 allel_str += "\tN"
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
            allel = sample[i]
            locus = loci[i]
            # Loci/Allel combination may not be found in the large profile file
            st_hits += loci_allel_dict[locus].get(allel, ["None"])

        # Find most frequent st_type in st_hits
        score = {}
        max_count = 1
        best_hit = ""
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
        st_output += "%s\t%s\t%d\t%.2f\n"%(sample_name, best_hit, max_count, similarity)

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
        if infile[-3:] == ".gz":
            f = gzip.open(infile, "rb")
            fst_char = f.read(1);
        else:
            f = open(infile, "rb")
            fst_char = f.read(1);
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
                        help="FASTQ files to do cgMLST or wgMLST on.",
                        nargs="+",
                        metavar="FASTQ",
                        default=None)
    # Optional arguments
    parser.add_argument("-o", "--output",
                        help="Output file.",
                        default="AlleleMatrix.mat",
                        metavar="OUTPUT_FILE")
    parser.add_argument("-s", "--species_scheme",
                        help="species_scheme scheme to apply e.g. ecoli_MLST.",
                        default=None,
                        metavar="species_scheme")
    parser.add_argument("-db", "--databases",
                        help="Directory containing the databases and gene\
                              lists for each species_scheme.",
                        default=("/home/data1/services/cgMLSTFinder/"
                                 "database_cgs"),
                        metavar="DB_DIR")
    parser.add_argument("-t", "--tmp_dir",
                        help="Temporary directory for storage of the results\
                              from the external software.",
                        default="cgMLST_tmp_dir")
    parser.add_argument("--kmapath",
                        help="Path to executable kma program",
                        default="kma",
                        metavar="kmapath")

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
        eprint(st_output)

        # Write ST-type output
        with open(st_filename, "w") as fh:
            fh.write(st_output)

    # Write allel matrix output
    with open(args.output + ".txt", "w") as fh:
        fh.write("\n".join(allel_output) + "\n")

    eprint("Done")
