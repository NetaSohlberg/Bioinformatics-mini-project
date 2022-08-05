'''
Jump finding algorithm
Authors: Neta Solberg, Hadas Samuels
The algorithm receives CSV files containing mutations in worldwide corona virus' sequences, and find and returns a CSV file with jumps found in the sequences, which lead to potential variants
Here includes 1st and 3rd versions. 2nd version in a different file
'''

import pandas as pd # Imports

MIN_JUMP = 5  # Const obj for minimum jump size

'''
1st case: with a reference file
'''

'''
# Read the CSV files
file_name = "1st_ver_file.csv"
reference_file_name = "reference.csv"
samples=pd.read_csv(file_name)
ref=pd.read_csv(reference_file_name)

# Save the relevant columns
nuc_sub_ref=ref["nuc sub"].tolist() # Reference mutations column
nuc_sub_samples=samples["nuc sub"].tolist() # Mutations column
sample=samples["sample"].tolist() # Samples column

# Create a dictionary of the new mutations in each sample
dict_seq = {}
special_mut=[] # List of the new mutations in each sample
s=sample[0]
intersect_list = [] # List of the intersection between the mutations of the samples

# Insert each sample and his new mutations that don't exist in the reference into a dictionary
for x in range(len(sample)):
    if sample[x] != s: # Start of a new sample - insert the mutation list of the previous sample to the dictionary
        if not dict_seq: # The first sample - insert the whole mutations list into the intersect list
            intersect_list = special_mut
        elif intersect_list: # Not the first sample - change the intersect list according to the new mutations list
            intersect_list = [i for i in intersect_list if i in special_mut]
            
        dict_seq[s]=special_mut
        s=sample[x]
        special_mut=[nuc_sub_samples[x]]
    
    else: # Continuation of the previous sample
        if not nuc_sub_samples[x] in nuc_sub_ref: # A new mutation - don't exist in the reference list
            special_mut.append(nuc_sub_samples[x])

dict_seq[sample[-1]]=special_mut # The last sample
'''

'''
************************************************
'''

'''
3rd version: without a reference file
'''

#'''
# Read the CSV file
file_name = "3rd_ver_file.csv"
samples = pd.read_csv(file_name)

# Save the relevant columns
nuc_sub_samples = samples["nuc_sub"].tolist() # Mutations column
sample = samples["accession_id"].tolist() # Samples id column

# Create a dictionary of the mutations in each sample
dict_seq = {}
mut_list = [] # List of the mutations in each sample
s = sample[0]
intersect_list = [] # List of the intersection between the mutations of the samples

# Insert each sample and his mutations into a dictionary
for x in range(len(sample)): # Go through the samples
    if sample[x]!=s: # Start of a new sample - insert the mutation list of the previous sample to the dictionary
        if not dict_seq: # The first sample - insert the whole mutations list into the intersect list
            intersect_list = mut_list
        elif intersect_list: # Not the first sample - change the intersect list according to the new mutations list
            intersect_list = [i for i in intersect_list if i in mut_list]

        dict_seq[s]=mut_list
        s=sample[x]
        mut_list=[nuc_sub_samples[x]]

    else: # Continuation of the previous sample
        mut_list.append(nuc_sub_samples[x])

dict_seq[sample[-1]]=mut_list # The last sample
#'''

'''
************************************************
From here apply to both versions
'''

# Create lists of the mutations combinations between each 2 samples and samples names
keys_list=list(dict_seq) # List of samples names
mut = [] # List for mutations combinations
seq = [] #List of samples names combinations

# Go through the optional pairs of samples
for i in range(len(keys_list)-1):
    for j in range(i+1, len(keys_list)):
        seq.append([keys_list[i], keys_list[j]])
        mut.append(list(set(dict_seq[keys_list[i]]) & set(dict_seq[keys_list[j]])))
        mut[-1] = list(set(mut[-1]) - set(intersect_list)) # Remove the mutation that exist in all samples

# Find the jumps
jump_seq = [] # List of the samples in each jump
jump_mut = [] # List of the mutations in each jump

for key in seq: # Go through the samples combinations
    flag1 = True # Boolean variable for minimum jump length
    if len(mut[seq.index(key)]) < MIN_JUMP: # The combination between the 2 samples is too small, cannot be part of a jump
        flag1 = False

#Compare the current combination to all other samples' mutations
    names = key
    flag2 = False # Boolean variable for the case of combine's alteration

    for i in dict_seq.keys(): # Go through the samples in dict_seq
        if i not in names: # The sample is not in the combination we check
            if flag2: # Combine altered in some way
                combine = list(set(combine) & set(dict_seq[i]))
            else:
                combine = list(set(mut[seq.index(names)]) & set(dict_seq[i]))  # Combine the combination seq with the new seq (i)

            if len(combine) >= MIN_JUMP:   # If combine is a jump (same to the combination between the 2 samples or smaller)
                names.append(i) # Add i to the jump samples

            elif len(combine) >0:   # If combine is smaller than MIN_JUMP but still exists (bigger than 0)
                if len(mut[seq.index(names)]) - len(combine) < MIN_JUMP:   # If the difference between combine and the combination between the 2 samples is smaller than MIN_JUMP
                    flag1 = False    # The combination between the 2 samples is not a jump
                    break
                else:   # The combination between the 2 samples is a jump without the mutations in combine
                    temp = [j for j in mut[seq.index(names)] if j not in combine] # Define new combine - with the combination without combine mutations
                    combine = temp
                    flag2 = True # Alteration of combine

            else: # If there is no combination at all (len(combine) == 0)
                combine = mut[seq.index(names)]
                flag2 = True # Alteration of combine

    if sorted(names) in jump_seq: # Reccurent jump
        flag1 = False

    if flag1 == True and len(names)>2: # Found a jump!
        jump_seq.append(sorted(names)) # Add the jump's sample's names to jump_seq
        jump_mut.append(combine) # Add the jump's mutations to jump_mut
        print('...')

# Print the results to an output file
output_file = open("output " + file_name[:-4] + ".txt", "w")

for mut in range(len(jump_mut)): # Go through the jumps found
    output_file.write('Jump number '+str(mut + 1)+"\n")
    output_file.write('number of sequences: ' + str(len(jump_seq[mut]))+"\n")
    output_file.write('number of mutations: ' + str(len(jump_mut[mut]))+"\n")
    output_file.write('the jump is: ' + str(jump_seq[mut])+"\n")
    output_file.write('the mutations are: ' + str(jump_mut[mut])+"\n")
    output_file.write('*****\n')

# Optional output - print on the screen
'''
for mut in range(len(jump_mut)):
    print ('Jump number', mut+1)
    print ('number of mutations: ', len(jump_mut[mut]))
    print ('number of sequences: ', len(jump_seq[mut]))
    print ('the jump is: ', jump_seq[mut] )
    print ('the mutations are: ', jump_mut[mut])
'''










