import os
import sys
import re

#genome1 = sys.argv[1] #E coli
#genome2 = sys.argv[2] #A tumefaciens

for i in range(1, len(sys.argv), 2):
    """cmd = 'zcat %s | blastp -query - -out blast_raw_%s -db blastdb/%s -evalue 0.01 -outfmt \'6 qseqid sseqid qlen slen bitscore evalue pident nident length qcovs qstat qend sstart send\'' % (sys.argv[i]+'.faa.gz', sys.argv[i], sys.argv[i+1])

    os.system(cmd)
    """
    f1 = open('blast_raw_'+sys.argv[i], 'r')
    f2 = open('blast_'+sys.argv[i], 'w')
    for line in f1:
        s = re.search(r'ref\|(.+)\|', line)
        a = line.split('\t')
        scov = float(a[8])/int(a[3])
        h = ''
        for j in range(2, 13):
            h += str(a[j]) + '\t'
        f2.write(a[0] + '\t' + s.group(1) + '\t' + h + str(scov) + '\n')
        #re.sub(r'\n', str(scov)+'\n', line)


