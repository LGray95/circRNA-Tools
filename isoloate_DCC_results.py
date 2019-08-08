in_DCC = open('/Users/lachlan/Documents/Masters_Research_2019/Data/DCC/Epilepsy/SUDEPmRNA_103187/CircRNACount.txt').read().splitlines()

for k in range(3,51):
    out_list = []
    for i in in_DCC:
        col = i.split('\t')
        chr = col[0]
        start = col[1]
        end = col[2]
        sample = col[k]
        out_list.append(chr+":"+start+"-"+end+'\t'+sample)


    filtered_out_list = []
    sample_ID =(out_list[0].split('\t')[1])
    print(sample_ID)
    for i in out_list[1:-1]:
        if int(i.split('\t')[1]) >= 2:
            filtered_out_list.append(i)

    outfile = open('/Users/lachlan/Documents/Masters_Research_2019/Data/DCC/Epilepsy/SUDEPmRNA_103187/'+sample_ID+'.txt', 'w')
    for i in filtered_out_list:
        outfile.write(str(i) + '\n')
    outfile.close()