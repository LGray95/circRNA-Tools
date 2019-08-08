from collections import Counter

def AGO_binding(infile, AGO_file, outfile):

    import os
    # Setting empty lists
    circ_list = []
    AGO_list = []
    circRNA_IDs = []

    # Reading in filenames - IN ORDER
    open_infile = open(infile).read().splitlines()
    for line in open_infile[1:-1]:
        if 'chr' in line:
            col = line.split('\t')
            #microRNA = col[0]
            #miRNA_chr = col[0].split(":")[0]
            #miRNA_start = int(col[0].split(":")[1].split("-")[0])
            #miRNA_end = int(col[0].split(":")[1].split("-")[1])
            #miRNA_length = int(miRNA_end) - int(miRNA_start)
            circRNA = col[0]
            circRNA_chr = col[0].split(":")[0]
            circRNA_start = int(col[0].split(":")[1].split("-")[0])
            circRNA_end = int(col[0].split(":")[1].split("-")[1])

            list_a = (circRNA, circRNA_chr, circRNA_start, circRNA_end)
                      #microRNA, miRNA_chr, miRNA_start, miRNA_end, miRNA_length)
            circ_list.append(list_a)
            lst = circRNA
            circRNA_IDs.append(lst)


    # Open a file of AGO binding sites
    AGO_bed = open(AGO_file).read().splitlines()
    for line in AGO_bed:
        col = line.split('\t')
        AGO_chr = col[0]
        AGO_start = int(col[1])
        AGO_end = int(col[2])
        AGO = AGO_chr+":"+str(AGO_start)+"-"+str(AGO_end)
        motif_length = AGO_end - AGO_start

        list_b = (AGO_chr, AGO_start, AGO_end, motif_length, AGO)
        AGO_list.append(list_b)

    AGO_list = set(AGO_list)
    # Setting list of matches
    circRNA_out_list = [i[0]
                        for i in circ_list
                        for j in AGO_list
                        if i[1] == j[0] and i[2] <= j[1]and i[3] >= j[2]]

    # miRNA_out_list = [(m[4] + "\t" + n[4])
    #                   for m in circ_list
    #                   for n in AGO_list
    #                   if m[5] == n[0] and m[6] <= n[1] and m[7] >= n[2]]

    sample_list = []
    score_list = []

    # removing duplicate values in list and printing how many items remain from each sample
    base = (os.path.basename(str(infile)))
    sample = os.path.splitext(base)[0]
    score = (str(len(set(circRNA_out_list)))+"/"+str(len(set(circRNA_IDs))))
    binding_to_circRNA = Counter(circRNA_out_list)
    #binding_to_miRNA = Counter(miRNA_out_list)
    sample_list.append(sample)
    score_list.append(score)

    print(sample)
    print(score)
    # print("AGO binding to complete circRNA")
    # for k, v in binding_to_circRNA.items():
    #     print(k, ":", v)
    # print("AGO binding to only miRNA")
    # for k, v in binding_to_miRNA.items():
    #     print(k, ":", v)

    outfile = open(outfile, 'a')
    for i, j in zip(sample_list, score_list):
        outfile.write(i + '\n' + j + '\n')
        outfile.write("AGO binding to complete circRNA" + "\n")
        for k, v in binding_to_circRNA.items():
            outfile.write(k + ":" + "\t" + str(v) + "\n")
        # outfile.write("AGO binding to only miRNA" + "\n")
        # for k, v in binding_to_miRNA.items():
        #     outfile.write(k + ":" + str(v) + "\n")
        # outfile.write("\t")




for i in open('/Users/lachlan/Documents/Masters_Research_2019/Data/merged_files/Epilepsy/SUDEPmRNA_103065/filepath.txt').read().splitlines():
    AGO_binding(i, '/Users/lachlan/Downloads/AGO.bed', '/Users/lachlan/Desktop/AGO_binding_total.txt')