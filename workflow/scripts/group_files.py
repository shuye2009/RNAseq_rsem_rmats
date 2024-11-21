#!/usr/bin/env python
import pandas as pd

def make_files(files, design, comps):
    df = pd.read_csv(comps, sep=r'\s+')
    df_cond = pd.read_csv(design, sep=r'\s+', usecols=['sample','condition'])

    for [cdn1, cdn2] in df.values.tolist():
        comparison = "%s-%s" % (cdn1,cdn2)
        out1 = "config/" + comparison + "/rmats_b1.txt" 
        out2 = "config/" + comparison + "/rmats_b2.txt" 

        wtstrings = df_cond['sample'][df_cond['condition'] == cdn1].values.tolist()
        kostrings = df_cond['sample'][df_cond['condition'] == cdn2].values.tolist()
        dictionary = {"ko": [], "wt": []} 
        
        for file in files:
            for substring in wtstrings:
                if substring in file:  
                    dictionary["wt"].append(file)  
                    #print(file + " " + substring + " WT")
            for substring in kostrings:
                if substring in file:  
                    dictionary["ko"].append(file)  
                    #print(file + " " + substring + " KO")

        for key, value in dictionary.items():
            with open (out2 if key == "wt" else out1, 'w') as f:
                    f.write(','.join(value)) 
                    f.close
         
    
    

make_files(snakemake.input.files, snakemake.input.design, snakemake.input.comps)
#in the rule that executes the script, define the folder with bam files as the input, and final output files as two outputsq