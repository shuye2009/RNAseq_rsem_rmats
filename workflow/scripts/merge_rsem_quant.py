#!/usr/bin/python3

import pandas as pd
from gtfparse import read_gtf
import os

outfile = snakemake.output.tpm

final_df = pd.DataFrame()
cycle = 1
sample_list = []

for q in snakemake.input.quant:

    data = pd.read_csv(q, sep='\t', usecols=['transcript_id', 'gene_id', 'TPM'])

    # get sample name
    sample = os.path.basename(os.path.dirname(q))

    # reformat dataframe
    data.set_index('transcript_id', inplace=True)
    data.rename(columns={"TPM":sample}, inplace=True)

    # Merge dataframes
    if cycle == 1:
        final_df = data
    else:
        final_df = pd.merge(final_df, data, left_index=True, right_index=True)
    
    sample_list += [sample]
    cycle+=1

# Add gene name
df_gtf = read_gtf(gtf)
id_name = df_gtf[['gene_id','gene_name']].to_pandas().drop_duplicates(ignore_index=True)

index = final_df.index.tolist()
names =[]

for i in range(len(index)):
    id = index[i]
    name = id_name[id_name['gene_id']==id].iloc[0]['gene_name']
    names.append(name)

final_df['gene_name'] = names

cols = final_df.columns.tolist()
cols = cols[-1:] + cols[:-1]
cols = cols[-1:] + cols[:-1]
final_df = final_df[cols]

# Write to file
final_df.to_csv(outfile, sep='\t')