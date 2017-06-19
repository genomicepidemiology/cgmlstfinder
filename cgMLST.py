import os
import sys
import shutil
import argparse

from python_module_dependencies import Dependencies
from python_module_seqfilehandler import SeqFile


class KMA():

    def __init__(self, seqfile, tmp_dir):

        self.result_file = temp_dir + "/kma_" + seqfile.filename
        self.seqfile = seqfile

        kma_call_list = [sysprgs["kma"], "-i", seqfile.path]
        if(seqfile.pe_file_reverse):
            kma_call_list.append(seqfile.pe_file_reverse)
        kma_call_list.append([
            "-t_db", args.databases + "/" + args.species + "/" + args.species,
            "-mem_mode",
            "-delta 1023",
            "-o", result_file])
        subprocess.call(kma_call_list)


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
                        default="./metadata",
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

    # Handle tmp dir
    args.tmp_dir = os.path.abspath(args.tmp_dir)
    os.makedirs(args.tmp_dir, exist_ok=True)

    # Check configuration files.
    if(not args.config):
        args.config = os.path.dirname(
            os.path.realpath(__file__)) + "/cgMLST_config.txt"
    if(not os.path.isfile(args.config)):
        print("Configuration file not found:", args.config)
        quit(1)

    # Get external software dependencies
    prgs = Dependencies(args.config)

    # # TEST DATA # #
    # TODO This part should be replaced to handle all species.
    args.species = "listeria"

    # Load files and pair them if necessary
    files = SeqFile.parse_files(args.input, phred=33)

    for seqfile in files:
        seq_kma = KMA(seqfile, args.tmp_dir)
        print("Finished KMA for: " seq_kma.seqfile.filename)
