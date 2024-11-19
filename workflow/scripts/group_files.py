#!/usr/bin/env python

def make_files(files, design, comparisons):
    df = pd.read_csv(comparisons, sep=r'\s+')
    df_cond = pd.read_csv(design, sep=r'\s+', usecols=['sample','condition'])

    for [cdn1, cdn2] in df.values.tolist():
        comparison = "%s-%s" % (cdn1,cdn2)
        out2 = "config/" + comparison + "/rmats_bt2.txt"
        out1 = "config/" + comparison + "/rmats_bt1.txt"
        wtstrings = df_cond['sample'][df['condition'] == cdn1].values.tolist()
        kostrings = df_cond['sample'][df['condition'] == cdn2].values.tolist()
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
         
    
    

make_files(snakemake.input.files, snakemake.input.design, snakemake.input.comparisons)
#in the rule that executes the script, define the folder with bam files as the input, and final output files as two outputsq