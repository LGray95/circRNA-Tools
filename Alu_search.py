circfile = open('/Users/lachlan/Documents/Masters_Research_2019/Data/microRNA/HEK293/total_miRNA_enriched.csv').read().splitlines()
Alu_file = open('/Users/lachlan/Downloads/Alurepeats.bed').read().splitlines()

circ_list = []
Alu_list = []

for i in circfile[1:-1]:
    col = i.split(',')
    circRNA_chr = col[3].split(':')[0]
    circRNA_start = int(col[3].split(':')[1].split('-')[0])
    circRNA_end = int(col[3].split(':')[1].split('-')[1])

    list_n = (circRNA_chr, circRNA_start, circRNA_end)
    circ_list.append(list_n)

for i in Alu_file[0:-1]:
    col = i.split('\t')
    Alu_chr = col[0]
    Alu_start = int(col[1])
    Alu_end = int(col[2])
    Alu_ID = col[3]

    list_m = (Alu_chr, Alu_start, Alu_end, Alu_ID)
    Alu_list.append(list_m)

downstream_Alu = []
upstream_Alu = []

for i in circ_list:
    for j in Alu_list:
        if i[0] == j[0] and (i[1]-500) <= j[1] and i[1] >= j[2]:
            list_z = (i[0]+':'+str(i[1])+"-"+str(i[2])+"|"+j[0]+":"+str(j[1])+"-"+str(j[2])+":"+j[3])
            downstream_Alu.append(list_z)
        if i[0] == j[0] and i[2] <= j[1] and (i[2]+500) >= j[2]:
            list_x = (i[0]+':'+str(i[1])+"-"+str(i[2])+"|"+j[0]+":"+str(j[1])+"-"+str(j[2])+":"+j[3])
            upstream_Alu.append(list_x)
out_list = []
for i in downstream_Alu:
    for j in upstream_Alu:
        if i.split('|')[0] == j.split('|')[0]:
            list_c = (i, j)
            out_list.append(list_c)

outfile = open('/Users/lachlan/Documents/Masters_Research_2019/Data/ICM/HEK293/Enriched/Alu_repeats.txt', 'w')

for i in out_list:
    outfile.write(str(i) + '\t' + '\n')
outfile.close()