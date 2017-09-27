#!/usr/local/bin/python3

from __future__ import print_function
import os
import sys
import shutil
import argparse
import subprocess
import pickle
import re

from python_module_dependencies.Dependencies import Dependencies
from python_module_seqfilehandler.SeqFileHandler import SeqFile


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class KMA():

    def __init__(self, seqfile, tmp_dir, db, gene_list, kma_path="kma"):
        """ Constructor map reads from seqfile object using kma.
        """

        result_file_tmp = tmp_dir + "/kma_" + seqfile.filename
        self.result_file = result_file_tmp + ".res"
        #self.result_file = "/home/people/s160509/kma_DTU2017_302_PRJ1085_Campylobacter_jejuni_ZTA14_00124CP_R1.res"
        self.seqfile = seqfile

        # Create kma command line list
        kma_call_list = [kma_path, "-i", seqfile.path]
        # Add reverse reads if paired-end data
        if(seqfile.pe_file_reverse):
            kma_call_list.append(seqfile.pe_file_reverse)

        kma_call_list += [
            "-t_db", args.databases + "/" + species+ "/" + species,
            "-mem_mode",
            "-delta", "1023",
            "-o", result_file_tmp]

        eprint("# KMA call: " + " ".join(kma_call_list))

        subprocess.call(kma_call_list)
        eprint("KMA call ended")


    def calc_best_hits(self):
        """ Extracts best hits from kma results into a dict that is returned
        """
        matrix = {}
        #self.cgmlst_file = matrix
        eprint("Func: Calc best hit:")
        
        # Create dict of locus and allel with the highest quality score
        with open(self.result_file, "r") as result_file:
            header = result_file.readline()
            loci_allel = re.compile("(\D+\d+)_(\d+)")
            i = 0
            for line in result_file:
                i += 1
                data = line.rstrip().split("\t")
                loci_allel_object = loci_allel.search(line)
                locus = loci_allel_object.group(1)
                allel = loci_allel_object.group(2)
                if int(allel) == 926:
                    print(locus, allel)
#                if i < 10:
#                    print(locus, allel)
                q_score = float(data[6])
                if locus in matrix:
                    if matrix[locus][1] < q_score:
                        matrix[locus] = [allel, q_score]
                else:
                    matrix[locus] = [allel, q_score]
        return matrix

def st_typing(pickle_path, inp, output_file):
    eprint("Finding ST type!!!")
    # Load pickle
    try:
        loci_allel_dict = pickle.load(open(pickle_path, "rb"))
    except IOError:
        sys.stdout.write("Error, pickle not found", pickle_path)
        quit(1)


    # Find best ST type for all allel profiles

    # First line contains matrix column headers, which are the specific loci
    loci = inp[0].strip().split("\t")
    print(len(inp))
    for sample_str in inp[1:]:
        sample = sample_str.strip().split("\t")
        sample_name = sample[0]
        st_hits = []
        print(len(sample))
        for i in range(1, len(sample)):
            allel = sample[i]
            locus = loci[i]
            if i%50 == 0: 
                print(i, allel, locus)
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
        output += (sample_name + "\t" + str(best_hit) + "\t" + str(max_count)
                   + "\t" + str(round((max_count / (len(loci) - 1)) * 100, 2))
                   + "\n")
    return output



#def st_typing(pickle_path, input_file, output_file):
#    print("Finding ST type")
#    #load pickle 
#    try:
#        loci_allel_dict = pickle.load(open( pickle_path, "rb"))
#    except:
#        sys.stdout.write("Error, pickle not found", pickle_path)
#        quit(1)
#       
#    #open output file
#    outfile = open(output_file, "w")
#    #write header in output file
#    outfile.write("Sample_Names\tcgST_Assigned\tNo_of_Found_Allels\tSimilarity\n")
#    
#    input_file = input_file.split("\n")
#    
#    #find best ST type for all allel profiles
#    loci = input_file[0].strip().split("\t")
#    for sample_str in input_file[1:]:
#        sample = sample_str.strip().split("\t")
#        sample_name = sample[0]
#        st_hits = []
#        for i in range(1, len(sample)):
#	    if i %100 ==0: print(i)
#            allel = sample[i]
#            locus = loci[i]
#            try:
#                st_hits += loci_allel_dict[locus][allel]
#            #keyError may occur if loci/allel combination not seen in the large profile file
#            except KeyError:
#                pass
#        #find most frequent st_type in st_hits
#        score = {}
#        max_count = 1
#        best_hit = ""
#        for hit in st_hits:
#            if hit in score:
#                score[hit] += 1
#                if max_count < score[hit]:
#                    max_count = score[hit]
#                    best_hit = hit
#            else:
#                score[hit] = 1
#        outfile.write(sample_name + "\t" + str(best_hit) + "\t" + str(max_count) + "\t" + str(round((max_count/(len(loci)-1))*100, 2)) + "\n")
#	#outfile.write(sample_name + "\t" + str(best_hit) + "\t" + str(max_count) + "\t" + str(round((max_count/(len(loci)-1))*100, 2)) + "\n")
#	
#    outfile.close()	

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
    parser.add_argument("-s", "--species",
                        help="Species scheme to apply.",
                        default=None,
                        metavar="SPECIES")
    parser.add_argument("-db", "--databases",
                        help="Directory containing the databases and gene\
                              lists for each species.",
                        default=("/home/data1/services/cgMLSTFinder/"
                                 "database_cgs"),
                        metavar="DB_DIR")
    parser.add_argument("-t", "--tmp_dir",
                        help="Temporary directory for storage of the results\
                              from the external software.",
                        default="cgMLST_tmp_dir")
    parser.add_argument("--config",
                        help="Configuration file containing information on\
                              software dependencies. Per default it is assumed\
                              that the file is located in the same directory\
                              as the main script. ",
                        default=None,
                        metavar="FILE")
    parser.add_argument("-st", "--st_output",
                        help="ST-typing Output file.",
                        default="cgMLST_STtypes.mat",
                        metavar="ST-TYPING_OUTPUT_FILE")

    args = parser.parse_args()

    # Handle tmp dir
    args.tmp_dir = os.path.abspath(args.tmp_dir)
    os.makedirs(args.tmp_dir, exist_ok=True)

    # Check configuration files.
    if(not args.config):
        args.config = os.path.dirname(
            os.path.realpath(__file__)) + "/cgMLST_config.txt"
    if(not os.path.isfile(args.config)):
        eprint("Configuration file not found:", args.config)
        quit(1)

    # Get external software dependencies
    prgs = Dependencies(args.config)

    # Species database
    species = args.species
    db_dir = args.databases + "/" + species
    db_species = db_dir + "/" + species

    # Test if database is found and indexed
    db_files = [species + ".b", species + ".length.b",
                species + ".name.b"]
    for db_file in db_files:
        if(not os.path.isfile(db_dir + "/" + db_file)):
            eprint("ERROR: A KMA index file seems to be missing from the"
                   "database directory. You may need to run kma_index.\n"
                   "Missing file: " + db_dir + "/" + db_file)
            quit(1)

    # Gene list
    gene_list_filename = (db_dir + "/" + "list_genes.txt")

    if(not os.path.isfile(gene_list_filename)):
        eprint("Gene list not found at expected location:", gene_list_filename)
        quit(1)

    # Load files and pair them if necessary
    eprint("Parsing files:" + str(args.input))
    files = SeqFile.parse_files(args.input, phred=33)

    # Get gene_list.txt into list
    gene_list = []
    try:
        gene_list_file = open(gene_list_filename, "r")
    except IOError:
        eprint("NO valid genefile found")
        sys.exit(-1)
    gene_list = [locus.strip() for locus in  gene_list_file.readlines()]

    gene_list_file.close()
    
    #eprint("Creating AlleleMatrix")

    # write header to output file
    output = ["Genome\t%s" %("\t".join(gene_list))]
    file_no = 0
    for seqfile in files:
        # Run KMA to find alleles
        seq_kma = KMA(seqfile=seqfile, tmp_dir=args.tmp_dir, db=db_species,
                      gene_list=gene_list_file, kma_path=prgs["kma"])
        eprint("Finished KMA for: " + seq_kma.seqfile.filename)
        print(seq_kma)
        # get called alleles
        best_alleles = seq_kma.calc_best_hits()
        # prepare output string
        output_str = args.input[file_no]
        for locus in gene_list:
            locus = locus.strip()
            if locus in best_alleles:
                 output_str += "\t%s" %(best_alleles[locus][0])
            else:
                 output_str += "\tN"
        output += [output_str]
        file_no += 1


        # create st-type file if pickle containing profile list excist
        pickle_path = db_dir + "/%s_cgMLST_profile.p"%(species)
        if os.path.isfile(pickle_path):
	    # Write header in output file
	    st_output = ("Sample_Names\tcgST_Assigned\tNo_of_Allels_Found\tSimilarity\n")
            st_output += st_typing(pickle_path, output, args.st_output)

        else:
            args.st_output = None 

    # write output file
    with open(args.output, "w") as fh:
        fh.write("\n".join(output) + "\n")
    if os.path.isfile(pickle_path):
        with open(args.st_output, "w") as fh:
            fh.write(st_output)
    eprint("Done")

