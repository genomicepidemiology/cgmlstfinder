cgMLSTFinder
===================

Plasmid Multi-Locus Sequence Typing


Documentation
=============

The cgMLSTFinder service contains one python script *cgmlstfinder.py* which is the script of the latest
version of the cgMLSTFinder service. The method assigns core genome MLST alleles based on WGS data for several isolates.
If more 3 or more isolates are run together in one go, and the '--neighbor' flag is set with a executable neighbor program,
the allele dictances between the isolates will be calculated and a neighbor joining tree will be generated.

## Content of the repository
1. cgMLST.py     - the program
2. README.md
3. Dockerfile   - dockerfile for building the pmlst docker container


## Installation

Setting up cgMLSTFinder program
```bash
# Go to wanted location for cgmlstfinder
cd /path/to/some/dir
# Clone and enter the cgmlstfinder directory
git clone https://bitbucket.org/genomicepidemiology/cgmlstfinder.git
cd cgmlstfinder
```

Build Docker container
```bash
# Build container
docker build -t pmlst .
```

#Download and install cgMLST database
```bash
# Go to the directory where you want to store the cgmlst database
cd /path/to/some/dir
# Clone install script from git repository
git clone https://bitbucket.org/genomicepidemiology/cgmlstfinder_db.git
cd cgmlstfinder_db
# Install pMLST database with executable kma_index program
python3 INSTALL.py
```

This script will install the already kma_index cgMLST schemes

## Usage

The program can be invoked with the -h option to get help and more information of the service.

```bash
usage: cgMLST.py [-h] [-o OUTPUT_FILE] [-s SPECIES_SCHEME] [-db DB_DIR]
                 [-t TMP_DIR] [-k KMA_PATH] [-n NJ_PATH]
                 FASTQ [FASTQ ...]

positional arguments:
  FASTQ                 FASTQ files to do cgMLST on.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_FILE, --outdir OUTPUT_FILE
                        Output file.
  -s SPECIES_SCHEME, --species SPECIES_SCHEME
                        species schemes to apply, e.g. ecoli_cgMLST. Must
                        match the name of the species database
  -db DB_DIR, --databases DB_DIR
                        Directory containing the databases and gene lists for
                        each species_scheme.
  -t TMP_DIR, --tmp_dir TMP_DIR
                        Temporary directory for storage of the results from
                        the external software.
  -k KMA_PATH, --kmapath KMA_PATH
                        Path to executable kma program
  -n NJ_PATH, --nj_path NJ_PATH
                        Path to executable neighbor joining program

```
 
Example of command to run cgMLSTFinder:

```bash
python3 cgMLST.py /path/to/isolate.fq.gz -s salmonella -o /path/to/outdir -db ../cgmlstfinder_db/ -k /usr/local/bin/kma -n /usr/local/bin/neighbor
``` 

## Web-server

A webserver implementing the methods is available at the [CGE website](http://www.genomicepidemiology.org/) and can be found here:
https://cge.cbs.dtu.dk/services/cgMLSTFinder/

Citation
=======

No references

License
=======

Copyright (c) 2014, Ole Lund, Technical University of Denmark
All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.