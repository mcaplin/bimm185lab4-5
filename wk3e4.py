from Bio import SeqIO
import gzip
import sys

inp = gzip.open(sys.argv[1], 'rb')
out = open('out.txt', 'w')
taxList = []

for record in SeqIO.parse(inp, 'swiss'):
    taxid = record.annotations['ncbi_taxid']
    organism = record.annotations['organism']
    taxonomy = record.annotations['taxonomy']

    if taxid not in taxList:
        taxList.append(taxid[0])

        taxonomyString = ''
        for t in taxonomy:
            taxonomyString += str(t) + ', '
        taxonomyString = taxonomyString[:-2]

        out.write(str(taxid[0]) + '\t' + str(organism) + '\t' + taxonomyString + '\n')

inp.close()
out.close()

