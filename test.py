import pysam
# infile = open('/Users/lachlan/Documents/Masters_Research_2019/Data/RBP_analysis/list.txt', 'r')
#
# out_list = []
# for i in infile:
#     protein = i.split(':')[1].split('(')[0].rstrip('\s')
#     print(protein)
#     it_list = (protein)
#     out_list.append(it_list)
#
# print(out_list)
#
# outfile = open('/Users/lachlan/Documents/Masters_Research_2019/Data/RBP_analysis/control_GO.txt', 'w')
#
# for i in out_list:
#     print(i)
#     outfile.write(str(i) + '\n')
# outfile.close()

samfile = pysam.AlignmentFile("/Users/lachlan/Documents/Masters_Research_2019/Data/IGV_alignments/enriched_sorted.bam", "rb")
counts = count(samfile, 'chr1:7833572-7839685')

print(counts)