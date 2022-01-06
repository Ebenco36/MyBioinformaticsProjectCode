import os
import sys
import itertools
import gzip
import pandas as pd
from collections import OrderedDict
import math
import numpy as np
import subprocess
from src.helper.getReads import getReads

def KallistoQuant(
    read_type, 
    sra_run_file, 
    srr_column_name, 
    fastq_dir,
    kallisto_quantify_dir,
    processor,
    index_kallisto_output_dir,
    sra_accessions
    ):

    try:
        os.makedirs(kallisto_quantify_dir)
    except FileExistsError:
        # directory already exists
        pass
    if sra_run_file != "":
        fileCont = downloadSRA(sra_run_file)
        if not os.path.isfile(sra_accessions):
            fileCont[[srr_column_name]].to_csv(sra_accessions, index=False, header=False)
        else:
            pass
        # fastq.sra_bd(file="sra_accessions.txt")
        df_tab = pd.read_csv(sra_accessions, header=None)
        for i in range(0, len(df_tab)):
            SRAFile = fileCont.iloc[i]['Run']
            read_1, read_2 = getReads(read_type, SRAFile, '.fastq.gz')
            read_1_dir = "{}/{}".format(fastq_dir, read_1)
            read_2_dir = "{}/{}".format(fastq_dir, read_2)
            print("{} - {}".format(read_2, read_1))
            print("Quantification entry: {}".format(SRAFile))
            output_dir = "{}/{}_quant".format(kallisto_quantify_dir, SRAFile)
            commandSalmon= [
                'kallisto', 
                'quant', 
                '--index={}'.format(index_kallisto_output_dir), 
                '--output-dir={}'.format(output_dir), 
                '--threads={}'.format(processor), 
                '--plaintext', 
                read_1_dir,
                read_2_dir
            ]
            
            subprocess.call(commandSalmon)
        print("Done generating quantifications...")


def SalmonQuant(
    read_type, 
    sra_run_file, 
    srr_column_name, 
    fastq_dir, 
    index_output_dir,
    quantify_dir,
    processor,
    sra_accessions
    ):


    if sra_run_file != "":
        fileCont = downloadSRA(sra_run_file)
        if not os.path.isfile(sra_accessions):
            fileCont[[srr_column_name]].to_csv(sra_accessions, index=False, header=False)
        else:
            pass
        # fastq.sra_bd(file="sra_accessions.txt")
        df_tab = pd.read_csv(sra_accessions, header=None)
        for i in range(0, len(df_tab)):
            SRAFile = fileCont.iloc[i]['Run']
            read_1, read_2 = getReads(read_type, SRAFile, '.fastq.gz')
            read_1_dir = "{}/{}".format(fastq_dir, read_1)
            read_2_dir = "{}/{}".format(fastq_dir, read_2)
            print("{} - {}".format(read_2, read_1))
            print("Quantification entry: {}".format(SRAFile))
            commandSalmon= [
                'salmon', 
                'quant', '-i', 
                index_output_dir,  
                '-l', 'A',  
                '-1', read_1_dir,  
                '-2', read_2_dir,
                '-p', processor, 
                '--validateMappings', '-o', '{}/{}_quant'.format(quantify_dir, SRAFile)]
            subprocess.call(commandSalmon)
        print("Done generating quantifications...")