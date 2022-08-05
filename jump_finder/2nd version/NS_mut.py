'''
2nd version of jump finding algorithm
Authors: Neta Solberg, Hadas Samuels
The algorithm receives CSV file containing mutations in worldwide corona virus' sequences, and the mutations of a specific variant, and find and returns a CSV file with jumps found in the sequences with consideration the mutation type
'''

import pandas as pd # Imports

MIN_MUT = 7 # Const obj for minimum jump size
MAX_OVERLAP = 3 # Const obj for maximum overlap allowed between samples


# Read the CSV file
file_name = "2nd_ver_file.csv"
samples = pd.read_csv(file_name)

# Save the relevant columns
nuc_sub_samples = samples["varname"].tolist() # Mutations column
sample = samples["sample"].tolist() # Samples id column
var_class = samples["varclass"].tolist() # Mutations types column
exist_on_ref = samples["C.37"].tolist() # Known variants mutations

# Create a dictionary of the mutations in each sample
dict_seq = {}
special_mut = [] # List of the mutations in each sample
s = sample[0]

# Insert each sample and his new mutations that don't exist in the variant into a dictionary, with a consideration of the mutation type
for x in range(len(sample)):
    if sample[x] != s: # Start of a new sample - insert the mutation list of the previous sample to the dictionary
        dict_seq[s] = special_mut
        s = sample[x]
        special_mut = []

    if exist_on_ref[x] == 0: # The mutation doesn't exist in the variant reference
         if var_class[x] == "extragenic" or var_class[x] == "SNP_silent": # If the mutation is extragenic or SNP silent - add 'S' to the end of it's name
            special_mut.append(nuc_sub_samples[x]+"S")
         else: # Else - this is NS mutations, add 'N' to the end of the mutation name
            special_mut.append(nuc_sub_samples[x] + "N")

dict_seq[sample[-1]] = special_mut # The last sample

# Create a dictionary of the mutations combinations between each 2 samples
combination_dict = {}
keys_list = list(dict_seq) # List of samples names

# Go through the optional pairs of samples
for i in range(len(keys_list)-1):
    for j in range(i+1, len(keys_list)):
        name = str(keys_list[i])+"+"+str(keys_list[j])
        combination_dict[name] = list(set(dict_seq[keys_list[i]]) & set(dict_seq[keys_list[j]]))


# Divide combination_dict into 2 dictionaries: LONG = longer than MIN_MUT, SHORT = shorter than MIN_MUT
long = {}
short = {}

for key in combination_dict.keys(): # Go through the combinations
    if len(combination_dict[key]) >= MIN_MUT: # Belongs to long
        long[key]=combination_dict[key]
    else: # Belongs to short
        short[key]=combination_dict[key]

# Find the jumps
jump_mut = [] # List of the mutations in each jump
jump_seq = [] # List of the samples in each jump

# For each combination in long dictionary - check if there is overlap with another sample that is longer than MAX_OVERLAP
for key in long.keys():
    names = key.split("+") # List of the samples in the combination
    flag2 = True # Boolean variable for the case of samples already in jump

    for sublist in jump_seq: # Check that the samples dont exist already in the jumps lists
        temp = [j for j in sublist if j in names]
        if len(temp) > 0: # Re-occured sample
            flag2 = False
            break

    if flag2 == False: # Not a potential jump
        break

    flag=True # Boolean variable for jump conditions

    for i in dict_seq.keys(): # Go through the samples in dict_seq
        if names[0] != i and names[1] != i: # The sample is not in the combination we check
            combine= len(list(set(long[key]) & set(dict_seq[i]))) # Combine the combination seq with the new seq (i)
            if  combine >= MAX_OVERLAP and combine < MIN_MUT: # Jump conditions don't exist
                flag = False
                break
            elif combine>= MIN_MUT: # All jump conditions exist
                names.append(i)

    if flag == True:     # If there is an optional jump
        # Check how many not extragenic or SNP_silent mutations there are in the optional jump
        s = sum(mut[-1] == "N" for mut in long[key]) # Count NS mutations
        if s >= 7: # At least 7 NS mutation in the jump - found a discent jump!
            jump_seq.append(names)
            jump_mut.append(long[key])

# Delete the alteration of mutations names (N or S)
for list in jump_mut:
    for mut in range(len(list)):
        list[mut] = list[mut][:-1]

# Print the results to an output file
output_file = open("output " + file_name[:-4] + ".txt", "w")

for mut in range(len(jump_mut)): # Go through the jumps found
    output_file.write('Jump number '+str(mut + 1)+"\n")
    output_file.write('number of sequences: ' + str(len(jump_seq[mut]))+"\n")
    output_file.write('number of mutations: ' + str(len(jump_mut[mut]))+"\n")
    output_file.write('the jump is: ' + str(jump_seq[mut])+"\n")
    output_file.write('the mutations are: ' + str(jump_mut[mut])+"\n")
    output_file.write('*****\n')
