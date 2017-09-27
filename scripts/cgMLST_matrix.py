#!/usr/bin/python
# This script will make a matrix of cgMLST alleles number of Campylobacter 
# input file = cgMLST input
# output file = cgMLST matrix
#### EXAMPLE #####
# python cgMLST_matrix.py gene_list ouput_file.txt input.cgMLST

import sys

# make all loci numbers
f = open(sys.argv[1], 'r')
allLines = f.readlines()
f.close()

all_loci = []
<<<<<<< HEAD:scripts/Listeria_cgMLST_matrix_ver01_May17.py
for i in allLines:
	all_loci.append(i.split('\n')[0].replace("_", "")) #if "_" in loci name, delete "_"
=======
for line in allLines:
	all_loci.append(line.strip())) 
>>>>>>> STtyping_dev:scripts/cgMLST_matrix.py
	
#for i in range(1,int(sys.argv[1])+1):
#	all_loci.append('CAMP' + str(i).zfill(4))

# write header 
matrix = open(sys.argv[2],'a')
matrix.write('Genome')
for locus in all_loci:
	matrix.write('\t')
	matrix.write(locus)
matrix.write('\n')


# open cgMLST input
filepath2 = sys.argv[3:]
for infile in filepath2:
	print infile, '...'
	f = open(infile,'r')
	allLines = f.readlines()
	f.close()

        # for each locus get the asigned allel
        matrix.write(infile)
	alleles = {}
	for line in allLines:

                line = line.strip().split('\t')
                locus = line[0]
                allel = line[1].split("|")[0]
		alleles[locus] = allel

	for locus in all_loci:
		matrix.write('\t')
		if locus in alleles:
			matrix.write(alleles[locus])
		else:
			matrix.write(str('N'))
	matrix.write('\n')


matrix.close()

## count cgMLST difference
f = open(sys.argv[2],'r')
allLines = f.readlines()
f.close()

genomes = {}
for i in allLines[1:]:
	list = []
	for j in i.split('\t')[1:]:
		list.append(j)
	genomes[i.split('\t')[0]] = list

# calculate allele differences
a = 1
result = {}
pair = []
for i in genomes:
    for j in genomes:
        if [i,j] not in pair or [j,i] not in pair:
            same = 0
            diff = 0
            c = 0
            total = 0
            for k in genomes[i]:
                if k == genomes[j][c]:
                    same += 1
                    total += 1
                if k != genomes[j][c]:
                    diff += 1 
                    total += 1
                c += 1
            pair.append([i,j])
            result[(i,j)] = diff

## make csv file for making heatmap
pair = []
for i in genomes:
    pair.append(i)

f = open('output_alleleHomology.csv','w')
f.write('topic')
for i in pair:
    f.write(',')
    f.write(str(i))
    #f.write(',')
f.write('\n')
for i in pair:
    f.write(str(i))
    for j in pair:
        #f.write(str(i))
        f.write(',')
        score = 100 - (float(result[(i,j)])*100.0/float(total)) 
        f.write(str(score))
    f.write('\n')

f = open('output_alleleDifferent.txt','w')
f.write('strains')
for i in pair:
    f.write('\t')
    f.write(str(i))
    #f.write(',')
f.write('\n')
for i in pair:
    f.write(str(i))
    for j in pair:
        #f.write(str(i))
        f.write('\t')
        score = (float(result[(i,j)]))
        f.write(str(score))
    f.write('\n')

f.close()


    
