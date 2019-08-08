# import sys
# circBase_file = sys.argv[1]
# circRNA_file = sys.argv[2]
# circRNA_ID = sys.argv[3]
# gene = sys.argv[4]
# outfile = sys.argv[5]

def circBase_validator(circBase_file, circRNA_file, circRNA_ID, counts, depth, outfile):
    circBase = open(circBase_file).read().splitlines()
    circRNA = open(circRNA_file).read().splitlines()

    circBase_list = []
    circRNA_list = []

    for line in circBase[1:-1]:
        col = line.split("\t")
        circBase_chr = col[0]
        circBase_start = int(col[1])
        circBase_end = int(col[2])
        circBase_strand = col[3]
        circBase_ID = col[4]
        circBase_tissue = col[7]
        circBase_gene = col[11]


        list_n = [circBase_chr, circBase_start, circBase_end, circBase_strand, circBase_ID, circBase_tissue, circBase_gene]
        circBase_list.append(list_n)

    for line in circRNA[1:-1]:
        col = line.split(',')
        circRNA_chr = col[circRNA_ID].split(":")[0].strip('\"')
        circRNA_start = int(col[circRNA_ID].split(":")[1].split("-")[0])
        circRNA_end = int(col[circRNA_ID].split(":")[1].split("-")[1].split(",")[0].strip('\"'))
        circRNA_counts = int(col[counts])
        circRNA_CPM = (circRNA_counts/int(depth))*(10**6)

        list_m = [circRNA_chr, circRNA_start, circRNA_end, circRNA_CPM]
        circRNA_list.append(list_m)


    yes_count = 0
    no_count = 0
    out_list = []
    for i in circBase_list:
        for j in circRNA_list:
            if i[0] == j[0] and i[2] == j[2] and j[3] >= 0.1:
                yes_count = yes_count + 1
                out_list.append([i[0]+":"+str(i[1])+"-"+str(i[2]),i[3],i[4],i[5],i[6],(j[0]+":"+str(j[1])+"-"+str(j[2])),j[3],(j[1])-(i[1])])
            else:
                no_count = no_count

    outfile = open(outfile, 'w')
    outfile.write('circBase_ID' + ',' + 'circBase_strand' + ',' + 'circRNA_name' + ',' 'circBase_tissue' + ',' 'circBase_gene' + ','+
              'circRNA_ID' + ',' + 'CPM' + ',' + 'Difference' + '\n')

    for j in out_list:
        for k in j:
            outfile.write(str(k) + ',')
        outfile.write('\n')
    outfile.close()
    print(yes_count)
    
# calling the function with command line arguments.
# circBase_validator(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int((sys.argv[5])), (sys.argv[6]))

circBase_validator('/Users/lachlan/Documents/Masters_Research_2019/Data/circBase/hsa_hg19_circRNA.txt','/Users/lachlan/Documents/Masters_Research_2019/Data/Filtered_circRNAs/HEK293/enriched_CIRC-DCC.csv',0, -1, 42307449, '/Users/lachlan/Documents/Masters_Research_2019/Data/circBase/HEK293/CIRCexplorer2-DCC/enriched.csv')

# Add assertion test for file type. If its a csv file then ask for a tsv