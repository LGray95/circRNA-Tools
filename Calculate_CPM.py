import pandas

for line in open('/Users/lachlan/Desktop/cortex_filenames.txt').read().splitlines():
    col = line.split('\t')
    filename = col[0]
    print(filename)
    seq_depth = col[3]
    infile = open('/Users/lachlan/Documents/Masters_Research_2019/Data/linear_analysis/Epilepsy/'+filename).read().splitlines()

    out_list = []
    for i in infile:
        col = i.split('\t')
        gene = col[0]
        qry_gene_id = col[4]
        coverage = (col[8])
        if coverage != '':
            coverage = float(coverage)
            cpm = coverage/int(seq_depth) * (10**6)
            lst = [gene, qry_gene_id, cpm]
            out_list.append(lst)

        outfile = open('/Users/lachlan/Documents/Masters_Research_2019/Data/linear_analysis/Epilepsy/'+line.split('.')[0]+'.txt', 'w')
        outfile.write('Gene' + '\t' + 'qry_gene_id' + '\t' + 'CPM' + '\n')
        for j in out_list:
            for k in j:
                outfile.write(str(k) + '\t')
            outfile.write('\n')
        outfile.close()




