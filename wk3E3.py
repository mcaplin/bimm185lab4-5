import sys
import re
import os


os.system('wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README')
rea = open('README')
readme = rea.read()

ba = re.findall('(UP[0-9]+)',readme)

b1 = ba[42]
b2 = ba[34]
b3 = ba[32]

bacteria = [b1, b2, b3]

for b in bacteria:
    os.system('wget -P %s ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/%s_*'% (b, b) )
