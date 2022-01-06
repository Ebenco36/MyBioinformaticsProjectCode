import os
import sys
import itertools
import gzip
import subprocess
from src.helper.getReads import getReads


def generateAdaptors(read_type, sra_run_file, srr_column_name, fastq_dir, sra_accessions):
    try:
        os.makedirs("adapters_gen")
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
            print("Bbmap entry: {}".format(SRAFile))
            in1 = "in1={}_1.fastq.gz".format(SRAFile)
            in2 = "in2={}_2.fastq.gz".format(SRAFile)
            outa = "outa=adapters_gen/{}.fa".format(SRAFile)
            commandBBMap = [
                'bbmerge.sh', 
                in1, 
                in2, 
                'threads=10', 
                outa
            ]
            subprocess.call(commandBBMap)
        print("Done generating adaptors")