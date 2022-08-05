'''
Statistics of mutation in corona virus sequences
Authors: Neta Solberg, Hadas Samuels
The algorithm receives a CSV file with worldwide corona virus samples and their mutations, and returns the file with extra columns of mutations frequencies and percentages, sorted in descending way
'''

# Imports
import pandas as pd
import tk_window

path=tk_window.path

# Find number of sequences
df=pd.read_csv(path)
sample_list = df["sample"]
counter = 1
sample = sample_list[0]
for x in sample_list: # Go through the samples
    if x!=sample: # A new sample
        counter = counter + 1
        sample = x

# Find mutations frequency
freq_dict = {} # Dictionary with the variants and their frequencies
for ind in df.index: # Go through the data
    full_name = df["refvar"][ind] + df["varname"][ind] + df["qvar"][ind]
    if full_name in freq_dict:
        freq_dict[full_name] = freq_dict[full_name] + 1
    else:
        freq_dict[full_name] = 1

# Edit the csv file's columns
df.drop('sample', axis=1, inplace=True)
df.drop('qpos', axis=1, inplace=True)
df.drop('qlength', axis=1, inplace=True)
freq = [0]*len(df)
perce = [0]*len(df)
df["freq"] = freq # Frequencies column
df["perce"] = perce # Percentages column

pd.options.mode.chained_assignment = None

# Add frequencies and percentages columns to the CSV file
for ind in df.index:
    full_name = df["refvar"][ind] + df["varname"][ind] + df["qvar"][ind]
    df["freq"][ind] = freq_dict[full_name]
    df["perce"][ind] = "%.3f" % ((freq_dict[full_name]/counter)*100)

df.sort_values("freq", ascending=False, inplace=True) # Sort the file according to descendent frequencies

df.drop_duplicates(inplace = True) # Remove duplicated rows

# Save the edited csv file
edited_path=path[:len(path)-4]
edited_path=edited_path+"_edited.csv"
df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
df.to_csv(edited_path, index=False)
print ("success!")




