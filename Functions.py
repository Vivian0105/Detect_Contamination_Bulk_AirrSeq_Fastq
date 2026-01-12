import numpy as np
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')
from Bio import SeqIO
import gzip
from sklearn.cluster import AgglomerativeClustering
from pathlib import Path
from itertools import combinations
from presto.IO import getOutputHandle
import shutil

default_missing_chars = set(['-', '.', 'n', 'N'])

def sequence_clustering(sequence_table, mismatch_threshold=1):
    """
    Cluster sequences in sequence_table based on the sequence distance between each other

    Arguments:
      sequence_table : A dataframe of sequences
      mismatch_threshold: maximum mismatch allowed between two sequences in one cluster

    Returns:
      The input dataframe with added column 'cluster'
    """
    # Convert sequences into array of characters and calculate distance matrix
    arr = np.array([list(s) for s in sequence_table["sequence"]])
    missing_chars = np.array(['-', '.', 'n', 'N'])
    missing = np.isin(arr, missing_chars)                 # (n, L) booleans
    A = arr[:, None, :]                                   # (n, 1, L)
    B = arr[None, :, :]                                   # (1, n, L)

    diff = (A != B)                                       # (n, n, L)
    masked = missing[:, None, :] | missing[None, :, :]    # (n, n, L)  either is missing
    effective_diff = diff & (~masked)                     # ignore positions with missing

    dist = effective_diff.sum(axis=2)                     # (n, n) integer distances

    # Cluster based on distance threshold
    clust = AgglomerativeClustering(
        metric="precomputed",         # Using Precomputed Distances
        linkage="single",
        distance_threshold=mismatch_threshold + 1e-9,  
        n_clusters=None,
        compute_full_tree=True  
    )
    sequence_table['cluster']=clust.fit_predict(dist)
    return(sequence_table)

def pairwise_contamination(df1, df2, umi='umi', mismatch_threshold=1, outfile='contamination_record.tsv'):
    '''
    Detect contamination in two samples

    Arguments:
      df1 : A dataframe of sequences from sample1
      df2 : A dataframe of sequences from sample2
      umi: column name for umi
      mismatch_threshold: maximum mismatch allowed between two sequences in one cluster
      outfile: file to store the output contamination table
    '''
    
    umi_set1=set(list(df1[umi]))
    umi_set2=set(list(df2[umi]))
    common_umi=umi_set1.intersection(umi_set2)

    if 'contaminated' not in df1.columns:
        df1['contaminated']=False
    if 'contaminated' not in df2.columns:
        df2['contaminated']=False

    df1_n_contaminated=0
    df2_n_contaminated=0
    for u in common_umi:
        df1_select =  df1[df1[umi]==u]
        df2_select =  df2[df2[umi]==u]
        # put all the sequences of same umi together and group by sequence length
        df = pd.concat([df1_select, df2_select], ignore_index=True)
        df['sequence_length'] =  df["sequence"].apply(lambda x: len(x))
        
        # In each length group, cluster the sequences. 
        # If sequence from different samples are in the same cluster, label it as contaminated
        for s_length in df['sequence_length'].unique():
            df_selected = df[df['sequence_length']==s_length]
            if len(df_selected['sample_id'].unique())>1:
                contaminated_seq=[]
                df_selected=sequence_clustering(df_selected, mismatch_threshold)
                for cluster in df_selected['cluster'].unique():
                    df_cluster = df_selected[df_selected['cluster']==cluster]
                    if len(df_cluster['sample_id'].unique())>1:
                        df1_n_contaminated = df1_n_contaminated + df_cluster[df_cluster['sample_id']==df1['sample_id'].unique().item()].shape[0]
                        df2_n_contaminated = df2_n_contaminated + df_cluster[df_cluster['sample_id']==df2['sample_id'].unique().item()].shape[0]
                        contaminated_seq = contaminated_seq + list(df_cluster['seq_id'])
                df1.loc[df1['seq_id'].isin(contaminated_seq),'contaminated']=True
                df2.loc[df2['seq_id'].isin(contaminated_seq), 'contaminated']=True
                        
    # append contamination table to outfile if there is contamination
    if(df1_n_contaminated+df1_n_contaminated>0):
        contaminated_table = pd.DataFrame([{"Sample_A" : df1['sample_id'].unique().item(),
                       "Sample_B" : df2['sample_id'].unique().item(),
                       "N_Contamination_A" : df1_n_contaminated,
                        "N_Contamination_B" : df2_n_contaminated,
                        "A_Contaminated_seqs" : ",".join(df1[df1['contaminated']==True]['seq_id']),
                        "B_Contaminated_seqs" : ",".join(df2[df2['contaminated']==True]['seq_id']) 
                      }])
        contaminated_table.to_csv(
            outfile,
            sep="\t",
            mode="a",
            header=not os.path.exists(outfile),
            index=False
         )


def _open_text(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)

def fastq_to_df(fastq_path):
    # function to convert fastq file to pandas dataframe 
    rows = []
    with _open_text(fastq_path, "rt") as handle:
        for rec in SeqIO.parse(handle, "fastq"):
            # rec.description includes id + the rest of header (without leading @)
            desc = rec.description
            parts = desc.split("|")     # Need to check whether works on all fastq header format

            row = {
                "seq_id": rec.id,
                "sequence": str(rec.seq),
                "qual_str": rec.format("fastq").splitlines()[3],  # exact quality line
            }

            # Parse key=value tokens from header into columns
            for tok in parts[1:]:
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    row[k] = v

            rows.append(row)

    return pd.DataFrame(rows)


def Detect_Contamination(fastq_file_list, umi='umi', outfile='contamination_record.tsv'):
    '''
    Detect contamination in all samples

    Arguments:
      fastq_file_list : List of fastq files one for each sample
      umi: description for umi
      outfile: file to store the output contamination table
    '''
    if os.path.exists(outfile):
        os.remove(outfile)
    for a, b in combinations(fastq_file_list, 2):
        df1=fastq_to_df(a)
        df2=fastq_to_df(b)
        df1["sample_id"]=Path(a).name.split("_")[0]   # This is an arbitrary way to find sample_id, may not apply to other cases
        df2["sample_id"]=Path(b).name.split("_")[0]
        pairwise_contamination(df1, df2, umi=umi, outfile=outfile )


def save_contaminated_pass_fastq(contamination_table_path, fastq_file_list, out_label="contamination-pass"):
    contamination_table = pd.read_csv(contamination_table_path, sep='\t')
    comtamination_samples=set(list(contamination_table['Sample_A'])+list(contamination_table['Sample_B']))
    for sample, input_fastq in zip([Path(a).name.split("_")[0] for a in fastq_file_list],fastq_file_list):
        out_handle = getOutputHandle(input_fastq, out_label=out_label)
        out_fastq = out_handle.name
        if sample in comtamination_samples:
            contaminated_seq =[x.split(',') for x in contamination_table[contamination_table['Sample_A']==sample]['A_Contaminated_seqs']] + \
                                   [x.split(',') for x in contamination_table[contamination_table['Sample_B']==sample]['B_Contaminated_seqs']]
            contaminated_seq = set([x for sublist in contaminated_seq for x in sublist])
            with open(out_fastq, "w") as out:
               for record in SeqIO.parse(input_fastq, "fastq"):
                   if record.id not in contaminated_seq:
                      SeqIO.write(record, out, "fastq")
        else:
            shutil.copyfile(input_fastq, out_fastq)
