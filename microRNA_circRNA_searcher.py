def miRNA_finder(in_gff3,infile, outfile):

    #reading in both .gff3 and filtered circRNA files
    gff3 = open(in_gff3).read().splitlines()
    circRNA_file = open(infile).read().splitlines()


    #creating empty variables for the lists
    microRNA_list = []
    circRNA_list = []
    #looping through .gff3 file
    for line in gff3[1:-1]:
        col = line.split('\t')
        microRNA_chr = col[0]
        microRNA_start = int(col[2])
        microRNA_end = int(col[3])
        microRNA_type = col[1]
        microRNA_ID = col[5]


        #creating a list with chr, start and end. Adding each list to the list variable
        list_n = [microRNA_chr, microRNA_start, microRNA_end, microRNA_type, microRNA_ID]
        microRNA_list.append(list_n)

    for line in circRNA_file[1:-1]:
        col = line.split('\t')
        circRNA_chr = col[0].split(":")[0].strip('"')
        circRNA_start = int(col[0].split(":")[1].split("-")[0])
        circRNA_end = (col[0].split(":")[1].split("-")[1]).strip('"')
        circRNA_strand = col[1]
        #circBase_tissue = col[2]
        circRNA_gene = col[2]
        CPM = col[3]
        exon_count = float(col[4])
        #Difference = col[6]

        #Add in circBase_tissue and Difference
        list_m = [circRNA_chr, circRNA_start, circRNA_end, circRNA_strand, circRNA_gene, CPM, exon_count]

        circRNA_list.append(list_m)


    yes_count = 0
    no_count = 0
    out_list = []
    for i in microRNA_list:
        for j in circRNA_list:
            if i[0] == j[0] and i[1] > j[1] and i[2] < int(j[2]) and i[3] == "miRNA_primary_transcript":
                yes_count = yes_count +1
                out_list.append([(i[0]+":"+str(i[1])+"-"+str(i[2])),i[3],i[4],(j[0]+":"+str(j[1])+"-"+str(j[2])),j[3],j[4],j[5], j[6]])
        else:
            no_count = no_count + 1
    print(yes_count)

    # Set column headings for outfile
    outfile = open(outfile, 'w')
    outfile.write('microRNA_coordinates' + ',' + 'microRNA_type' + ',' +
                  'microRNA_ID' + ',' + 'circRNA_ID' + ',' + 'strand' + ',' + 'circRNA_gene' + ',' + 'CPM' + ',' + 'Exon_Count' + '\n')
    for j in out_list:
        for k in j:
            outfile.write(str(k) + ',')
        outfile.write('\n')
    outfile.close()

miRNA_finder('/Users/lachlan/Documents/Masters_Research_2019/Data/microRNA/hsa.txt', '/Users/lachlan/Documents/Masters_Research_2019/Data/DiffExp/Hippocampus_foldchange.txt', '/Users/lachlan/Documents/Masters_Research_2019/Data/DiffExp/Hippocampus_miRNA.csv')

# filepath = open('/Users/lachlan/Documents/Masters_Research_2019/Data/merged_files/Epilepsy/SUDEPmRNA_103065/filepath.txt').read().splitlines()
# filenames = open('/Users/lachlan/Documents/Masters_Research_2019/Data/merged_files/Epilepsy/SUDEPmRNA_103065/filenames.txt').read().splitlines()
#
# for i, j in zip(filepath, filenames):
#     miRNA_finder('/Users/lachlan/Documents/Masters_Research_2019/Data/microRNA/hsa.txt', i, '/Users/lachlan/Documents/Masters_Research_2019/Data/microRNA/Epilepsy/SUDEPmRNA_103065/'+j+'.csv')