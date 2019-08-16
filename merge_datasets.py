import pandas

def merge_datasets(DCC_list, DCC_dir, CE2_list, read_coverage_list, out_dir):

# Where:
# DCC_list = a .txt file of just the file names of each DCC circRNA file
# DCC_dir = The full path to the directory where those files are located
# CE2_list = .txt file with the full path to CIRCexplorer2 output files. Each line must match the samples in the DCC_list
# read_coverage = .txt file with each line matching the order of samples in the DCC and CE2 files
# out_dir = The directory where you want the results stored

        DCC_filenames = open(DCC_list).read().splitlines()
        CE2_filepath = open(CE2_list).read().splitlines()
        read_coverage = open(read_coverage_list).read().splitlines()


        for i, j, k in zip(DCC_filenames, CE2_filepath, read_coverage):
                print(i)
                x = pandas.read_csv(DCC_dir+i, sep='\t', header=None)
                y = pandas.read_csv(j, sep='\t', header=None)

                # Adding 1 to the circRNA start column of CIRCexplorer2 to compensate for different coordinate systems
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

                merged_df.to_csv(out_dir+i, sep='\t', index=False)


