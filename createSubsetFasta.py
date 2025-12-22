from Bio import SeqIO

N = 1000

records = []
for i, rec in enumerate(SeqIO.parse("Data/swissprot_50k.fasta", "fasta")):
    if i >= N:
        break
    records.append(rec)

SeqIO.write(records, "Data/subset_1k.fasta", "fasta")
