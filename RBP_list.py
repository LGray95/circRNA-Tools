infile = open('/Users/lachlan/Documents/Masters_Research_2019/Data/circbase/Fibroblast/CIRI2-CIRCexplorer2/control_IDs.txt', 'r')

list = []
for i in infile:
    chromosome = i.split(':')[0]
    start = int(i.split(':')[1].split('-')[0])
    end = int(i.split(':')[1].split('-')[0])
    strand = i.split("\t")[1].rstrip('\n')

    it_list = (chromosome+':'+str(start-100)+'-'+str(end+100)+':'+str(strand))
    list.append(it_list)

print(list)

outfile = open('/Users/lachlan/Documents/Masters_Research_2019/Data/RBP_analysis/control_list.txt', 'w')

for i in list:
    print(i)
    outfile.write(str(i) + '\n')
outfile.close()
