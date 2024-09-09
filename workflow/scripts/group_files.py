import os

def make_files(files, out1, out2):  
 
    if os.path.isfile(out1):
        os.remove(out1)
    if os.path.isfile(out2):
        os.remove(out2)

    substrings = snakemake.params.ref
    dictionary = {"ko": [], "wt": []} 
     
    for file in files:
        for substring in substrings:
            if substring in file:  
                dictionary["wt"].append(file)  
                #print(file + " " + substring + " WT")
    for file in files:
        if file not in dictionary["wt"]:
            dictionary["ko"].append(file) 
            #print(file + " KO")

    for key, value in dictionary.items():
        with open (out2 if key == "wt" else out1, 'w') as f:
                f.write(','.join(value)) 
                f.close

make_files(snakemake.input, snakemake.output[0], snakemake.output[1])
#in the rule that executes the script, define the folder with bam files as the input, and final output files as two outputsq