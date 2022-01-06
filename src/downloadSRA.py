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

def downloadSRA(SRARunTable):

    try:
        X = pd.read_csv(SRARunTable)
        return X
    except Exception as e:
        print("{} File does not exist".format(SRARunTable))
        return pd.DataFrame()


def DownloadSplit(sra_run_file, srr_column_name, read_type, fastq_dir, sra_accession):
    if sra_run_file != "":
        fileCont = downloadSRA(sra_run_file)
        if fileCont.empty == True:
            print("{} not found".format(sra_run_file))
        else:
            print(fileCont)
            if not os.path.isfile(sra_accession):
                fileCont[[srr_column_name]].to_csv(sra_accession, index=False, header=False)
            else:
                pass
            df_tab = pd.read_csv(sra_accession, header=None)
            for i in range(0, len(df_tab)):
                SRAFile = fileCont.iloc[i]['Run']
                command = ['prefetch', SRAFile]
                subprocess.call(command)
                read_1, read_2 = getReads(read_type, SRAFile, '.fastq.gz')
                read_1_dir = "{}/{}".format(fastq_dir, read_1)
                read_2_dir = "{}/{}".format(fastq_dir, read_2)
                print("{} - {}".format(read_2, read_1))
                if os.path.isfile(read_1_dir) and os.path.isfile(read_2_dir):
                    print("Finish downloading {}".format(fileCont[[srr_column_name]].loc[i]))
                else:
                    SRADownloaded = "{}/{}.{}".format(SRAFile, SRAFile, 'sra')
                    if os.path.isfile(SRADownloaded):
                        print("Converting SRA file to fastq")
                        commandFastQ = ['fastq-dump', '--split-files', '--gzip', '--defline-qual', '+', SRADownloaded, '-O', fastq_dir]
                        subprocess.call(commandFastQ)
                        print('Done converting...')
                    else:
                        print("Sra file has not been downloaded.")
