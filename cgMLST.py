#! /usr/local/bin/python3

from __future__ import print_function
import os
import sys
import shutil
import argparse
import subprocess

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
        #self.result_file = "/home/data1/services/cgMLSTFinder/test_data/cgMLST_tmp_dir/kma_SRR972394_1.res"
        self.cgmlst_file = None
        self.seqfile = seqfile

        # Create kma command line list
        kma_call_list = [kma_path, "-i", seqfile.path]
        # Add reverse reads if paired-end data
        if(seqfile.pe_file_reverse):
            kma_call_list.append(seqfile.pe_file_reverse)

        kma_call_list += [
            "-t_db", args.databases + "/" + args.species + "/" + args.species,
            "-mem_mode",
            "-delta", "1023",
            "-o", result_file_tmp]

        eprint("# KMA call: " + " ".join(kma_call_list))

        subprocess.call(kma_call_list)

    def calc_best_hits(self):
        """ Extracts best hits from kma results
        """
        eprint("Func: Calc best hit:")
        self.cgmlst_file = self.result_file[:-4] + ".cgMLST"
        cmd_best_hit = (
            "cat " + self.result_file + " | "
            "sed 's/CAMP/CAMP\\t/g' | "
            "sed 's/_/\\t/g' | "
            "sed 's/#Template/#Template\\tGene\\tAllele/g' | "
            "sort -k9 -r | "
            "gawk '{if (!($2 in taken)){print $0,$2};taken[$2]=$2}' | "
            "sort -k2 | "
            "grep -v '#' > "
            + self.cgmlst_file)
        subprocess.call(cmd_best_hit, shell=True)
        eprint("Call: " + cmd_best_hit)


class AlleleMatrix():

    def __init__(self, gene_list, output, kma_object, python2_path="python"):
        """
        """
        eprint("Creating AlleleMatrix")
        self.output = output
        self.kma = kma_object

        # Check external script file
        script = (os.path.dirname(os.path.realpath(__file__))
                  + "/scripts/Listeria_cgMLST_matrix_ver01_May17.py")
        if(not os.path.isfile(script)):
            eprint("ERROR: Something is wrong with your installation.\n"
                   "File missing: " + script)
            quit(1)

        # Calculate best hits if needed
        if(kma_object.cgmlst_file is None):
            eprint("Cond: Allele calc best hits")
            kma_object.calc_best_hits()

        # Build script command list
        cmd = [python2_path, script, gene_list, output, kma_object.cgmlst_file]
        subprocess.call(cmd)


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
                        help="Output prefix.",
                        default="./cgMLST_",
                        metavar="OUTPUT_PREFIX")
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

    args = parser.parse_args()

    # # TEST DATA # #
    # TODO This part should be replaced to handle all species.
    args.species = "listeria"

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
    db_dir = args.databases + "/" + args.species
    db_species = db_dir + "/" + args.species

    # Test if database is found and indexed
    db_files = [args.species + ".b", args.species + ".length.b",
                args.species + ".name.b"]
    for db_file in db_files:
        if(not os.path.isfile(db_dir + "/" + db_file)):
            eprint("ERROR: A KMA index file seems to be missing from the"
                   "database directory. You may need to run kma_index.\n"
                   "Missing file: " + db_dir + "/" + db_file)
            quit(1)

    # Gene list
    gene_list_file = (db_dir + "/" + "list_genes.txt")

    if(not os.path.isfile(gene_list_file)):
        eprint("Gene list not found at expected location:", gene_list_file)
        quit(1)

    # Load files and pair them if necessary
    eprint("Parsing files:" + str(args.input))
    files = SeqFile.parse_files(args.input, phred=33)

    for seqfile in files:
        # Run KMA to find alleles
        seq_kma = KMA(seqfile=seqfile, tmp_dir=args.tmp_dir, db=db_species,
                      gene_list=gene_list_file, kma_path=prgs["kma"])
        eprint("Finished KMA for: " + seq_kma.seqfile.filename)

        # Create allele matrix
        matrix = AlleleMatrix(gene_list=gene_list_file,
                              output=args.output + seqfile.filename + ".mat",
                              kma_object=seq_kma,
                              python2_path=prgs["python2.7"])
        eprint("Finished allele matrix for: " + seq_kma.seqfile.filename)

    eprint("Done")
