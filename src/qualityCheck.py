import os
import sys
import itertools
import gzip
import subprocess
from src.helper.getReads import getReads


def fastQcMultiQc(sra_run_file, srr_column_name, read_type, fastqc_dir, multiqc_dir, fastq_dir, sra_accessions):
    try:
        os.makedirs(fastqc_dir)
        os.makedirs(multiqc_dir)
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
        fastqc_path_dir = "{}/".format(fastqc_dir)
        for i in range(0, len(df_tab)):
            SRAFile = fileCont.iloc[i]['Run']
            read_1, read_2 = getReads(read_type, SRAFile, '_fastqc.html')
            fastqc_1_dir = "{}/{}".format(fastqc_dir, read_1)
            fastqc_2_dir = "{}/{}".format(fastqc_dir, read_2)
            print("{} - {}".format(fastqc_1_dir, fastqc_2_dir))
            if os.path.isfile(fastqc_1_dir) and os.path.isfile(fastqc_2_dir):
                print("Done generate fastqc report for {}".format(fastqc_1_dir))
            else:
                print("Generating fastqc report for {}".format(fastqc_1_dir))
                globSearch = "{}/{}_*.fastq.gz".format(fastq_dir, SRAFile)
                print(globSearch)
                files = glob.glob(globSearch)
                print(files)
                read_1s, read_2s = files[0], files[1]
                commandFastQC = ['fastqc', '-o', "{}/".format(fastqc_dir), read_1s, read_2s]
                subprocess.call(commandFastQC)
                print("Done Analyzing your fastq files using fastqc...")
        if os.path.exists(fastqc_path_dir) and os.path.isdir(fastqc_path_dir) and len(os.listdir(fastqc_path_dir)) == 0:
            print("Combining fastqc reports using multiqc")
            commandMultiQC = ['multiqc', fastqc_path_dir, '-o', multiqc_dir]
            subprocess.call(commandMultiQC)
            print("Done combining reports.")
        else:
            print("No fastqc.html files available at the moment. Kindly try the procedure again.")
