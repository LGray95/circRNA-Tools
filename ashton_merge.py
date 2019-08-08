import pandas


# reading in the files
L5 = "/Users/lachlan/Downloads/CIRI2_2_CGATGT_L005.txt"
L6 = "/Users/lachlan/Downloads/CIRI2_2_CGATGT_L006.txt"

# Make a pandas dataframe
x = pandas.read_csv(L5, sep="\t", header=None)
y = pandas.read_csv(L6, sep="\t", header=None)

# Merging The dataframe
merged_df = pandas.merge(x, y, on=0)

merged_df.to_csv("/Users/lachlan/Desktop/test.txt", sep='\t', index=False)