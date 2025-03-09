from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML

record = SeqIO.read("contigs.fasta", format="fasta")
result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"))

# Une seule sequence query :
blast_record = NCBIXML.read(result_handle)
print(blast_record)

# Plusieurs sequences query :
blast_records = NCBIXML.parse(result_handle)
