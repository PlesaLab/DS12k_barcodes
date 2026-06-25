'''
we use a bunch of different orthogonal dna barcodes of different sizes. They all come from the same set (240K Elledge 25mers

Need to know which of our current barcodes come from which sequence in the original set. There are six subsets:
skpp20 FWD and REV - 20mers
skpp15FWD and REV - 15mers
bc384(filt_prim_12nt_Lev_3_Tm_40_42_GC_45_55_SD_2_mod_restriction_trim) - 12mers
bc1536 (filt_prim_12nt_Lev_3_Tm_38_44_GC_45_55_SD_2_trim) - 12mers

NEED: A script that parses the alignments and makes a single csv file.
INPUT: alignments
OUTPUT: single csv file

Column 1: name of sequence in 240k Elledge 25mer set
Column 2: First column should be name of sequence in 240k Elledge 25mer set. Second column should be skpp20 FWD. If an skpp20 FWD sequence matches one in the 240k Elledge 25mer set the name should be listed on the corresponding row. And the same for skpp20 REV, skpp15 FWD and REV, bc384, bc1536.

Basically want to know the origin of each sub sequences in the original set
'''
import pandas as pd

parent_df = pd.read_csv("seq_files/bc25mer.240k.fasta", header=None)
parent_df_names = parent_df.iloc[::2, :].reset_index(drop=True)  # give me everything from 0 to the end, incrementing by 4 Df = df[start:stop:step, :]
parent_df_seqs = parent_df.iloc[1::2, :].reset_index(drop=True)  # give me every sequence line
parent_df_seqs_names = pd.merge(parent_df_names, parent_df_seqs, right_index=True, left_index=True)
parent_df_seqs_names.columns = ['parent_name', 'parent_seq']
parent_df_seqs_names['parent_name'] = parent_df_seqs_names['parent_name'].str.slice(start=1)

list_of_file = ["skpp20FWD_pretty.txt", "skpp20REV_pretty.txt", "skpp15FWD_pretty.txt", "skpp15REV_pretty.txt", "bc384_pretty.txt", "bc1536_pretty.txt"]
for ff,file in enumerate(list_of_file):
    seq_name = file.split("_")[0]
    child_seq = "sequence_" + seq_name
    print(seq_name)
    df = pd.read_csv("BCs_psl/" + file, header=None)
    df_names = df.iloc[::4, :].reset_index(drop=True) # give me everything from 0 to the end, incrementing by 4 Df = df[start:stop:step, :]
    df_seqs = df.iloc[1::4, :].reset_index(drop=True) # give me every sequence line
    df_seqs_names = pd.merge(df_names,df_seqs,right_index=True, left_index=True)
    df_seqs_names.columns = [seq_name, child_seq]
    df_seqs_names[[seq_name, 'parent_name','len']] = df_seqs_names[seq_name].str.split(" of ", 3, expand=True)
    df_seqs_names[seq_name] = df_seqs_names[seq_name].str.slice(start=1)
    df_seqs_names['parent_name'] = df_seqs_names['parent_name'].str.slice(start=3, stop=-6)
    df_seqs_names.drop(columns='len', inplace=True)
    parent_df_seqs_names = parent_df_seqs_names.merge(df_seqs_names, how="outer", on='parent_name')


parent_df_seqs_names.to_csv("bc25mers_and_subsets_barcodes.csv")

