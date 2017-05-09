import gzip
from Bio import SeqIO
import sys
import re
import glob
import pymysql
import getpass
import os

cntGenome = 0 #global counters
cntReplicons = 0
cntCDS = 0

# dicitonaries for each table
genomes = {}
replicons = {}
genes = {}
gene_synonyms = {}
gene_functions = {}
exons = {}
gene_xrefs = {}

for a in range(1, len(sys.argv)): # args are the assembly names

    f = glob.glob('./*/*/%s*genomic.gbff.gz' % sys.argv[a])
    genomic = gzip.open(f[0], 'rb')
    gnmSize = 0 #counters for each replicon
    gnmCDS = 0
    numReplicons = 0

    cntGenome +=1
    genomes[cntGenome] = 7*[None]
    thisGenome = genomes[cntGenome] # Reference to the genome key we are working with
    for record in SeqIO.parse(genomic, "genbank"):
        cntReplicons +=1
        numReplicons += 1

        replicons[cntReplicons] = 8*[None] 
        thisReplicon = replicons[cntReplicons] # Reference to the replicon key we are working with
        # filling various parts of the tables
        thisReplicon[0] = cntGenome 

        thisReplicon[6] = record.name 
        thisReplicon[3] = record.annotations['topology']

        complete = re.search(r'Description:\s*.+, complete genome', str(record))
        if complete:
            thisReplicon[2] = 'chromosome'

        thisReplicon[7] = record.annotations['date']

        assembly = re.search(r'Assembly:\s*(\w+\.*\w+)', str(record))
        thisGenome[6] = assembly.group(1)

        thisGenome[2] = record.annotations['taxonomy'][0]

        repSize = re.search(r'\.\.(\d+)',str(record.annotations['contig'])).group(1)
        thisReplicon[5] = repSize
        gnmSize += int(repSize)

        if record.features:
            for feature in record.features:
                if feature.type == "source":
                    repCDS = 0 # counter for each gene

                    name = feature.qualifiers['organism'][0]
                    thisGenome[0] = name
                    thisReplicon[1] = name

                    if 'chromosome' in feature.qualifiers:
                        thisReplicon[2] = 'chromosome'

                    if 'plasmid' in feature.qualifiers:
                        thisReplicon[2] = 'plasmid'

                    if 'db_xref' in feature.qualifiers:
                        taxid = re.search(r'taxon:(.+)', str(feature.qualifiers['db_xref'][0]))
                        taxid = taxid.group(1)
                        thisGenome[1] = taxid
                if feature.type == 'CDS':
                    cntCDS += 1
                    gnmCDS += 1
                    repCDS += 1

                    genes[cntCDS] = 9*[None]
                    thisGene = genes[cntCDS] # Reference ot the gene we are working with
                    thisGene[0] = cntGenome
                    thisGene[1] = cntReplicons

                    gene_xrefs[cntCDS] = []
                    exons[cntCDS] = []

                    gene_length = len(feature.location)
                    thisGene[7] = gene_length
                    num_exons = len(feature.location.parts)

                    for ex in range(0, num_exons): # exon table
                        left = feature.location.parts[ex].start
                        right = feature.location.parts[ex].end
                        length = len(feature.location.parts[ex])
                        exons[cntCDS].append([ex+1, int(left), int(right), int(length)])
                    thisGene[6] = num_exons
                    strand = feature.location.strand
                    if strand == 1:
                        strand = 'F'
                    else:
                        strand = 'R'
                    thisGene[5] = strand
                    locus_tag = feature.qualifiers['locus_tag'][0]
                    thisGene[2] = locus_tag

                    # filling various part of gene, gene synonym, gene functions, xref tables
                    if 'gene' in feature.qualifiers:
                        thisGene[4] = feature.qualifiers['gene'][0]
                    else:
                        thisGene[4] = locus_tag

                    if 'gene_synonym' in feature.qualifiers or 'old_locus_tag' in feature.qualifiers:
                        gene_synonyms[cntCDS] = []

                    if 'gene_synonym' in feature.qualifiers:
                        synonyms = re.findall(r'\s*(\w+);*',str(feature.qualifiers['gene_synonym'][0]))
                        for synonym in synonyms:
                            gene_synonyms[cntCDS].append(synonym)

                    if 'old_locus_tag' in feature.qualifiers:
                        gene_synonyms[cntCDS].append([feature.qualifiers['old_locus_tag'][0]])

                    if 'protein_id' in feature.qualifiers:
                        protein_id = feature.qualifiers['protein_id'][0]
                        protein_id = re.search(r'(\w+)\.*', str(protein_id))
                        thisGene[3] = protein_id.group(1)
                        gene_xrefs[cntCDS].append(['refseq', protein_id.group(1)])
                    
                    if 'function' in feature.qualifiers:
                        gene_functions[cntCDS] = []
                        functions = feature.qualifiers['function'][0]
                        functions = re.findall(r'\s*([^;*]+)', functions)
                        for fcn in functions:
                            gene_functions[cntCDS].append(fcn)
                    
                    if 'product' in feature.qualifiers:
                        thisGene[8] = feature.qualifiers['product'][0]
                    else:
                        thisGene[8] = 'unavailable'
                    
                    if 'db_xref' in feature.qualifiers:
                        db_xref = feature.qualifiers['db_xref']
                        for i in db_xref:
                            db_xref = re.search(r'(.+):(.+)', str(i))
                            gene_xrefs[cntCDS].append([db_xref.group(1), db_xref.group(2)])

            thisReplicon[4] = repCDS

    thisGenome[3] = numReplicons
    thisGenome[4] = gnmCDS
    thisGenome[5] = gnmSize

    genomic.close()

genomes_table = open('genomes_table.txt', 'w')
replicons_table = open('replicons_table.txt', 'w')
genes_table = open('genes_table.txt', 'w')
exons_table = open('exons_table.txt', 'w')
gene_synonyms_table = open('gene_synonyms_table.txt', 'w')
gene_xrefs_table = open('gene_xrefs_table.txt', 'w')
functions_table = open('functions_table.txt', 'w')

# There are 3 kinds of outputs from the types of dictionaries I made
def WriteDicToFile(dic, out):
    for g in dic:
        out.write(str(g) + '\t')
        for i in dic[g]:
            out.write(str(i) + '\t')
        out.write('\n')
    out.close()

def WriteListToFile(dic, out):
    for g in dic:
        for i in dic[g]:
            out.write(str(g) + '\t' + str(i) + '\n')
    out.close()

def WriteDicListToFile(dic, out):
    for g in dic:
        for i in dic[g]:
            s = ''
            for x in i:
                s += '\t' + str(x)
            out.write(str(g) + s + '\n')
    out.close()

WriteDicToFile(genomes, genomes_table)
WriteDicToFile(replicons, replicons_table)
WriteDicToFile(genes, genes_table)
WriteListToFile(gene_functions, functions_table)
WriteListToFile(gene_synonyms, gene_synonyms_table)
WriteDicListToFile(exons, exons_table)
WriteDicListToFile(gene_xrefs, gene_xrefs_table)

password = getpass.getpass() 
os.system('mysql -h bm185s.ucsd.edu -u mcaplin -p %s mcaplin_db

"""password = getpass.getpass()
db = pymysql.connect('bm185s-mysql.ucsd.edu', 'mcaplin', password, 'mcaplin_db')
cursor = db.cursor()"""
