def microRNA_match(circRNA, microRNA):
    from Bio import SeqIO
    circRNA = SeqIO.read(circRNA, "fasta")
    microRNA = SeqIO.read(microRNA, "fasta")
    x = circRNA.seq.find(microRNA.seq)
    if x >= 0:
        print(x)
        print(circRNA.seq[x:-1])
        print(microRNA.seq)
    else:
        print("no dice...")


#microRNA_match('/Users/lachlan/Documents/Masters_Research_2019/Data/circRNA_miRNA_fasta/hsa_circ_0097827.fasta', '/Users/lachlan/Documents/Masters_Research_2019/Data/circRNA_miRNA_fasta/mir-5188.fasta')


x = "CTCACCTCATCTTTCTCATTCTTTAGTAGCTGATGCAAATTCCTGTTCTTACATTGTAACTGTCTGTCTCTTTCCTGAAGCCCAATCATTTCCCTCCTGG"
y = "GAACAACTGGAAACCTCCAGGAGGGAAATGATTGGGCTTCAGGAAAGAGACAGACAGTTACAATGTAAGAACAGGAATTTGCATCAGCTACTAAAGAATGAGAAAGATGAG"

print(x.find(y))