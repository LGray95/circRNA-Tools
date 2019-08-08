AS = open('/Users/lachlan/Documents/Masters_Research_2019/Data/CIRI-AS/GW13.txt.list').read().splitlines()
circBase = open('/Users/lachlan/Documents/Masters_Research_2019/Data/Filtered_circRNAs/circBase_uniqueGW13.csv').read().splitlines()

AS_list = []
circBase_list = []

for line in AS[0:-1]:
    col = line.split('\t')
    AS_chr = col[0].split(":")[0]
    AS_start = (col[0].split(":")[-1].split("|")[0])
    AS_end = col[0].split(":")[-1].split("|")[-1]

    list_n = [AS_chr, AS_start, AS_end]
    AS_list.append(list_n)

for line in circBase[0:-1]:
    col = line.split(',')
    circBase_chr = col[0]
    circBase_start = (col[1])
    circbase_end = col[2]

    list_m = [circBase_chr, circBase_start, circbase_end]
    circBase_list.append(list_m)

yes_count = 0
no_count = 0
out_list = []

for i in AS_list:
    for j in circBase_list:
        if i[0] == j[0] and i[2] == j[2]:
            yes_count = yes_count + 1
            out_list.append([i[0], i[1], i[2], j[0], j[1], j[2]])
        else:
            no_count = no_count + 1
outfile = open('/Users/lachlan/Documents/Masters_Research_2019/Data/CIRI-AS/circBase_CIRI-AS#2.csv', "w")
for j in out_list:
    for k in j:
        outfile.write(str(k) + ',')
    outfile.write('\n')
outfile.close()

print(yes_count)
print(no_count)
