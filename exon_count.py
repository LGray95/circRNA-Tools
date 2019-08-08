def exon_count(annotation, circRNA_file, outfile):

    from collections import defaultdict
    # input files
    gff3 = open(annotation).read().splitlines()
    circRNA = open(circRNA_file).read().splitlines()

    # creating empty variables
    gff3_list = []
    circRNA_list = []

    # reading in exon coordinates
    for line in gff3[0:-1]:
        col = line.split('\t')
        transcript_type = col[2]
        exon_chr = col[0]
        exon_start = int(col[3])
        exon_end = int(col[4])
        if transcript_type == "exon":
            list_a = [transcript_type, exon_chr, exon_start, exon_end]
            gff3_list.append(list_a)

    # reading in circRNA coordinates
    for line in circRNA[1:-1]:
        col = line.split(',')
        circRNA_ID = col[0]
        circRNA_chr = col[0].split(":")[0]
        circRNA_start = int(col[0].split(":")[1].split("-")[0])
        circRNA_end = int(col[0].split(":")[1].split("-")[1])
        list_b = [circRNA_ID, circRNA_chr, circRNA_start, circRNA_end]
        circRNA_list.append(list_b)

    # finding 5', internal and 3' exons for each circRNA and adding each exon length to list
    exon_out_list = []
    for j in circRNA_list:
        for i in gff3_list:
            if i[1] == j[1]:
                # i[1] = exon chromosome
                # i[2] = exon start
                # i[3] = exon_end
                # j[1] = circRNA chromosome
                # j[2] = circRNA_start
                # j[3] = circRNA_end
                # matching for 5' exons
                if i[2] <= j[2] and i[3] >= j[2]:
                    exon_out_list.append(j[0]+" "+i[1]+":"+str(i[2])+"-"+str(i[3])+"="+str(i[3]-i[2]))
                # matching for 3' exons
                if i[2] <= j[3] and i[3] >= j[3]:
                    exon_out_list.append(j[0]+" "+i[1]+":"+str(i[2])+"-"+str(i[3])+"="+str(i[3]-i[2]))
                # matching for internal exons
                if i[2] >= j[2] and (j[3] - i[3] >= 0):
                    exon_out_list.append(j[0]+" "+i[1]+":"+str(i[2])+"-"+str(i[3])+"="+str(i[3]-i[2]))

    # removing duplicates with set() and returning to a list with list()
    exon_out_list = set(exon_out_list)
    exon_out_list = list(exon_out_list)


    # matches all circRNA_IDs so that multiple exons are close to each other
    matches = []
    for i in circRNA_list:
        for j in exon_out_list:
            if i[0] in j.split(" ")[0]:
                matches.append(j.split(" ")[0]+"="+j.split(" ")[1].split("=")[1])

    # Suming up exon length per circRNA
    d = defaultdict(int)
    lst = [item.split("=") for item in matches]
    for k, v in lst:
        d[k] += int(v)

    # Creating list of circRNAs to be counted below
    circRNA_count_list = []
    for i in exon_out_list:
        circRNA = i.split(" ")[0]
        circRNA_count_list.append(circRNA)

    # counting occurance of circRNA in circRNA_count_list and adding to list_c
    # list_c has the circRNA_ID and the exon count
    count_list = []
    for i in circRNA_list:
        count = circRNA_count_list.count(i[0])
        if count > 0:
            list_c = i[0], count
            count_list.append(list_c)

    # Selecting strings to write out
    lst = []
    for i in count_list:
        for key,value in d.items():
            if key == (i[0]):
                circRNA_chr = i[0].split(':')[0]
                start = int(i[0].split(':')[1].split("-")[0])
                end = int(i[0].split(':')[1].split("-")[1])
                difference = ((end - start) - value)
                list_d = (i[0], str(i[1]), str(value), str(end - start), str(difference))
                lst.append(list_d)

    # writing out to file
    outfile = open(outfile, 'w')
    outfile.write('circRNA_ID' + ',' + 'exon count' + ',' + 'total exon length' + ',' + 'total circRNA length' + ',' + 'Difference' + ',' + '\n')
    for i in lst:
        for k in i:
            outfile.write(str(k) + ',')
        outfile.write('\n')
    outfile.close()

exon_count('/Users/lachlan/Documents/UCSC_exons_genome', '/Users/lachlan/Documents/Masters_Research_2019/Data/circBase/HEK293/CIRCexplorer2-DCC/enriched.csv', '/Users/lachlan/Documents/Masters_Research_2019/Data/exon_count/HEK293/enriched_exoncount.csv')