import pandas

DCC_filenames = open('/Users/lachlan/Documents/Masters_Research_2019/Data/DCC/Epilepsy/SUDEPmRNA_103187/filenames.txt').read().splitlines()
CE2_filepath = open('/Users/lachlan/Documents/Masters_Research_2019/Data/CIRCexplorer2_results/Epilepsy/SUDEPmRNA_103187/filepath.txt').read().splitlines()
read_coverage = open('/Users/lachlan/Documents/Masters_Research_2019/Data/merged_files/Epilepsy/seq_depth.txt').read().splitlines()


for i, j, k in zip(DCC_filenames, CE2_filepath, read_coverage):
        x = pandas.read_csv('/Users/lachlan/Documents/Masters_Research_2019/Data/DCC/Epilepsy/SUDEPmRNA_103187/'+i, sep='\t', header=None)
        y = pandas.read_csv(j, sep='\t', header=None)

        # Adding 1 to the circRNA start column of CIRCexplorer2
        y[1] += 1

        # Making the circRNA_ID column
        y[0] = y.apply(lambda x:'%s:%s' % (x[0],x[1]),axis=1)
        y[0] = y.apply(lambda x:'%s-%s' % (x[0],x[2]),axis=1)

        # Merging The dataframe
        merged_df = pandas.merge(x, y, on=0)

        merged_df['1_x'] = pandas.to_numeric(merged_df['1_x'])

        merged_df[12] = pandas.to_numeric(merged_df[12])

        merged_df[18] = ((merged_df['1_x'] + merged_df[12])/2)

        merged_df['CPM'] = ((merged_df[18]/int(k))*10**6)

        merged_df = merged_df[merged_df.CPM >= 0.1]

        print(merged_df[merged_df.duplicated(subset=0)])

        # merged_df = merged_df.drop_duplicates()

        merged_df.to_csv('/Users/lachlan/Documents/Masters_Research_2019/Data/merged_files/Epilepsy/'+i, sep='\t', index=False)



