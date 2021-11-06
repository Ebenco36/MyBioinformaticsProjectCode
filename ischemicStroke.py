# !/usr/bin/env base

import os
import sys
import itertools
import gzip
import argparse
import subprocess
from functools import reduce
import pandas as pd
from collections import OrderedDict
import math
import numpy as np
import warnings
import time
# from bioinfokit.analys import fastq
import glob
import re, sys
from ROperations import Load
from MachineLearningModel import machineLearning

def downloadSRA(SRARunTable):

    try:
        X = pd.read_csv(SRARunTable)
        return X
    except Exception as e:
        print("{} File does not exist".format(SRARunTable))
        return pd.DataFrame()


def getReads(read_type, SRR_, ext):
    if read_type == "paired_end_read":
        read_1 = "{}_1{}".format(SRR_, ext)
        read_2 = "{}_2{}".format(SRR_, ext)
    elif(read_type == "single_end_read"):
        read_1 = "{}_1{}".format(SRR_, ext)
        read_2 = ""
    else:
        read_1 = ""
        read_2 = ""
    return read_1, read_2

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

def applyTrimmomatic(
        trimmed_fastq_dir, 
        sra_run_file, 
        srr_column_name, 
        read_type, 
        fastq_dir, 
        fastqc_dir_trimmed, 
        multiqc_dir_trimmed,
        generate_adaptors,
        adaptors_dir,
        sra_accessions
        ):
    try:
        os.makedirs(trimmed_fastq_dir)
        os.makedirs(fastqc_dir_trimmed)
        os.makedirs(multiqc_dir_trimmed)
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
        fastqc_path_dir = "{}/".format(fastqc_dir_trimmed)
        for i in range(0, len(df_tab)):
            SRAFile = fileCont.iloc[i]['Run']
            # get fastq files if exist
            read_1, read_2 = getReads(read_type, SRAFile, '.fastq.gz')
            fastq_1_dir = "{}/{}".format(trimmed_fastq_dir, read_1)
            fastq_2_dir = "{}/{}".format(trimmed_fastq_dir, read_2)
            print("{} - {}".format(fastq_1_dir, fastq_2_dir))
            # Check if trimmed fastq file exist
            if os.path.isfile(fastq_1_dir) and os.path.isfile(fastq_2_dir):
                print("Done generate fastqc report for {}".format(fastq_1_dir))
                # fastqc analysis
                read_1_html, read_2_html = getReads(read_type, SRAFile, '_fastqc.html')
                fastqc_1_dir_html = "{}/{}".format(fastqc_dir_trimmed, read_1_html)
                fastqc_2_dir_html = "{}/{}".format(fastqc_dir_trimmed, read_2_html)
                print("{} - {}".format(fastqc_1_dir_html, fastqc_2_dir_html))
                # Check if trimmed fastq html report file exist
                if os.path.isfile(fastqc_1_dir_html) and os.path.isfile(fastqc_2_dir_html):
                    print("Finished generating report for trimmed fastq files")
                else:
                    # generate report
                    print("Generating fastqc report for {}".format(fastq_1_dir))
                    globSearch = "{}/{}_*_paired.fastq.gz".format(trimmed_fastq_dir, SRAFile)
                    files = glob.glob(globSearch)
                    read_1s, read_2s = files[0], files[1]
                    commandFastQC = ['fastqc', '-o', "{}/".format(fastqc_dir_trimmed), read_1s, read_2s]
                    subprocess.call(commandFastQC)
                    print("Done Analyzing your fastq files using fastqc...")
                if os.path.exists(fastqc_path_dir) and os.path.isdir(fastqc_path_dir) and len(os.listdir(fastqc_path_dir)) == 0:
                    print("Combining fastqc reports using multiqc")
                    commandMultiQC = ['multiqc', fastqc_path_dir, '-o', multiqc_dir_trimmed]
                    subprocess.call(commandMultiQC)
                    print("Done combining reports.")
                else:
                    print("No fastqc.html files available at the moment. Kindly try the procedure again.")

            else:
                print(" Trimmomatic entry: {} - {}".format(fastq_1_dir, fastq_2_dir))
                globSearch = "{}/{}_*.fastq.gz".format(fastq_dir, SRAFile)
                files = glob.glob(globSearch)
                read_1s, read_2s = files[0], files[1]
                x_name = lambda x: "{}/{}_{}".format(trimmed_fastq_dir, SRAFile, x)
                if generate_adaptors == "yes":
                    commandFastQC = [
                        'trimmomatic', 
                        'PE', 
                        '-threads', '14', 
                        '-phred33', 
                        read_1s, 
                        read_2s,
                        x_name("forward_paired.fastq.gz"),
                        x_name("forward_unpaired.fastq.gz"),
                        x_name("reverse_paired.fastq.gz"),
                        x_name("reverse_unpaired.fastq.gz"),
                        'ILLUMINACLIP:{}/{}.fa:2:30:10'.format(adaptors_dir, SRAFile), 
                        'LEADING:3',
                        'TRAILING:3',
                        'SLIDINGWINDOW:4:15',
                        'MINLEN:36'
                    ]
                else:
                    commandFastQC = [
                        'trimmomatic', 
                        'PE', 
                        '-threads', '14', 
                        '-phred33', 
                        read_1s, 
                        read_2s,
                        x_name("forward_paired.fastq.gz"),
                        x_name("forward_unpaired.fastq.gz"),
                        x_name("reverse_paired.fastq.gz"),
                        x_name("reverse_unpaired.fastq.gz"),
                        'ILLUMINACLIP:adaptors/TruSeq3-PE.fa:2:30:10', 
                        'LEADING:3',
                        'TRAILING:3',
                        'SLIDINGWINDOW:4:15',
                        'MINLEN:36'
                    ]
                subprocess.call(commandFastQC)
                print("Done with Trimmomatic")
        
def generateIndex(reference_genome_file, index_output_dir):
    try:
        os.makedirs(index_output_dir)
    except FileExistsError:
        # directory already exists
        pass
    commandSalmonIndex = ['salmon', 'index', '-t', reference_genome_file, '-i', index_output_dir,  '--gencode']
    subprocess.call(commandSalmonIndex)
    print("Done generating index. This can be found at {}".format(index_output_dir))

def generateKallistoIndex(reference_genome_file, index_kallisto_output_dir):
    # try:
    #     os.makedirs(index_kallisto_output_dir)
    # except FileExistsError:
    #     # directory already exists
    #     pass
    commandSalmonIndex = ['kallisto', 'index', "--index={}".format(index_kallisto_output_dir), reference_genome_file]
    subprocess.call(commandSalmonIndex)
    print("Done generating index. This can be found at {}".format(index_kallisto_output_dir))

def buildFiles(reference_genome_file):
    if os.path.exists("enst.txt"):
        os.remove("enst.txt")
    if os.path.exists("ensg.txt"):
        os.remove("ensg.txt")
    else:
        print("The file does not exist")

    f = open("gene_map.csv", "w")
    file_object = open('enst.txt', 'a')
    file_object_gene = open('ensg.txt', 'a')
    for line in open(reference_genome_file, 'r'):
        if re.search('>', line):
            transcript = re.compile(r'ENST\d{11}')
            geneTranscript = re.compile(r'ENSG\d{11}')
            file_object.writelines(transcript.search(line).group())
            file_object_gene.writelines(geneTranscript.search(line).group())
            file_object.write("\n")
            file_object_gene.write("\n")
            # transcript_list.append(transcript.search(line).group())
            if line == None:
                print('no matches found')
    file_object.close()
    file_object_gene.close()
    commandBuildFiles = ['paste', '-d', ',', 'enst.txt', 'ensg.txt']
    # print(commandBuildFiles)
    data = subprocess.call(commandBuildFiles, stdout=f)

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

def RunRAnalysis(file="SraRunTableIschemicStroke.txt", path_to_quant="quants/", gene_map="gene_map.csv", count_type="salmon", ignoreTxVersion='TRUE', sample_split="1:10, 11:15", split_group="AIS,NC", sample_count="10, 5", result_path="RData3"):
    Load(file, path_to_quant, gene_map, count_type, ignoreTxVersion, sample_split, split_group, sample_count, result_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--read_type', type = str, default="paired_end_read", choices = ['single_end_read', 'paired_end_read'])
    parser.add_argument('--sra_run_file', type = str, default="SraRunTableIschemicStroke.txt")
    parser.add_argument('--srr_column_name', type = str, default="Run")

    parser.add_argument('--fastq_dir', type = str, default="FASTQ")
    parser.add_argument('--sra_accession', type = str, default="sra_accessions.txt")
    parser.add_argument('--fastqc_dir', type = str, default="FASTQC")
    parser.add_argument('--multiqc_dir', type = str, default="MULTIQC")
    parser.add_argument('--applyTrimmomatic', type = str, default = "no", choices = ['yes', 'no'])
    parser.add_argument('--trimmed_fastq_dir', type = str, default = "trimmedData")
    parser.add_argument('--fastqc_dir_trimmed', type = str, default="FASTQC_TRIMMED")
    parser.add_argument('--multiqc_dir_trimmed', type = str, default="MULTIQC_TRIMMED")
    parser.add_argument('--generate_adaptors', type = str, default = "no", choices = ['yes', 'no'])
    parser.add_argument('--adaptors_dir', type = str, default = "adapters_gen")
    parser.add_argument('--index_genome', type = str, default = "no", choices = ['yes', 'no'])
    parser.add_argument('--reference_genome_file', type = str, default = "./gencode.v35.transcripts.fa")
    parser.add_argument('--index_output_dir', type = str, default = "gencode_v33_index")

    parser.add_argument('--quantify', type = str, default = "no", choices = ['yes', 'no'])
    parser.add_argument('--quantify_dir', type = str, default = "quants")

    parser.add_argument('--processor', type = str, default = "4")
    parser.add_argument('--index_kallisto_genome', type = str, default = "no", choices = ['yes', 'no'])
    parser.add_argument('--index_kallisto_output_dir', type = str, default = "gencode_v33_index_kallisto")
    
    parser.add_argument('--kallisto_quantify', type = str, default = "no", choices = ['yes', 'no'])
    parser.add_argument('--kallisto_quantify_dir', type = str, default = "kallisto_quants")



    parser.add_argument('--runAnalysis', type = str, default = "no", choices = ['yes', 'no'])

    parser.add_argument('--gene_map', type = str, default="gene_map.csv",)
    parser.add_argument('--count_type', type = str, default="salmon", )
    parser.add_argument('--ignoreTxVersion', type = str, default='TRUE', choices = ['TRUE', 'FALSE'])
    parser.add_argument('--sample_split', type = str, default="1:10, 11:15", )
    parser.add_argument('--split_group', type = str, default="AIS,NC", )
    parser.add_argument('--sample_count', type = str, default="10, 5")
    parser.add_argument('--result_path', type = str, default="RData3")

    parser.add_argument('--machineLearningAnalysis', type = str, default = "no", choices = ['yes', 'no'])
    parser.add_argument('--input_file', type = str, default="../../../normalizedCounts.csv")
    parser.add_argument('--graphFileName', type = str, default="result")

    parser.add_argument('--component_n', type = str, default="3")
    parser.add_argument('--path_to_save', type = str, default="machineLearningReport")
    
    


    args = parser.parse_args()

    # DownloadSplit(args.sra_run_file, args.srr_column_name, args.read_type, args.fastq_dir, args.sra_accession)
    # fastQcMultiQc(args.sra_run_file, args.srr_column_name, args.read_type, args.fastqc_dir, args.multiqc_dir, args.fastq_dir, args.sra_accession)
    if args.applyTrimmomatic == 'yes':
        applyTrimmomatic(
            args.trimmed_fastq_dir, 
            args.sra_run_file, 
            args.srr_column_name, 
            args.read_type, 
            args.fastq_dir, 
            args.fastqc_dir_trimmed, 
            args.multiqc_dir_trimmed,
            args.generate_adaptors,
            args.adaptors_dir,
            args.sra_accession
        )
    if args.generate_adaptors == 'yes':
        generateAdaptors(args.read_type, args.sra_run_file, args.srr_column_name, args.fastq_dir, args.sra_accessions)

    if args.index_genome == 'yes':
        generateIndex(args.reference_genome_file, args.index_output_dir)


    if args.index_kallisto_genome == 'yes':
        generateKallistoIndex(args.reference_genome_file, args.index_kallisto_output_dir)


    if args.quantify == 'yes':
        SalmonQuant(args.read_type, args.sra_run_file, args.srr_column_name, args.fastq_dir, args.index_output_dir,args.quantify_dir, args.processor, args.sra_accessions)

    if args.kallisto_quantify == 'yes':
        KallistoQuant(
            args.read_type, 
            args.sra_run_file, 
            args.srr_column_name, 
            args.fastq_dir, 
            args.kallisto_quantify_dir, 
            args.processor, 
            args.index_kallisto_output_dir,
            args.sra_accessions)


    if args.runAnalysis == 'yes':
        RunRAnalysis(
            file=args.sra_run_file, 
            path_to_quant=args.quantify_dir, 
            gene_map=args.gene_map, 
            count_type=args.count_type, 
            ignoreTxVersion=args.ignoreTxVersion, 
            sample_split=args.sample_split, 
            split_group=args.split_group, 
            sample_count=args.sample_count,
            result_path = args.result_path
        )
    if args.machineLearningAnalysis == 'yes':
        machineLearning(
            file=args.input_file, 
            path_to_save=args.path_to_save, 
            component_n=args.component_n, 
            condition=args.split_group, 
            sample_count=args.sample_count,
            graphFileName=args.graphFileName
        )

    # buildFiles(args.reference_genome_file)

# python3 ischemicStroke.py --sra_run_file=SraRunTableIschemicStroke.txt --kallisto_quantify=yes

# grep -e '^ENST[[:digit:]]\{11\}' ./gencode.v35.transcripts.fa > enst.txt

# python3 ischemicStroke.py --sra_run_file=/Users/mac/Downloads/SraRunTableV2.txt


#python3 ischemicStroke.py --runAnalysis="yes" --result_path="sub_acute_control" --sra_run_file=SraRunTableIschemicStroke.txt --quantify_dir="quants/" --sample_split="1:5, 11:15" --split_group="SAS,NC", --sample_count="5, 5"
#python3 ischemicStroke.py --runAnalysis="yes" --result_path="acute_control" --sra_run_file=SraRunTableIschemicStroke.txt --quantify_dir="quants/" --sample_split="6:10, 11:15" --split_group="AIS,NC", --sample_count="5, 5"

#python3 ischemicStroke.py --runAnalysis="yes" --result_path="all_samples_with_control" --sra_run_file=SraRunTableIschemicStroke.txt --quantify_dir="quants/" --sample_split="1:5, 6:10, 11:15" --split_group="SAS, AIS, NC", --sample_count="5, 5, 5"

#python3 ischemicStroke.py --machineLearningAnalysis="yes" --input_file="../../../normalizedCounts.csv" --path_to_save="machineLearningReport" --component_n="3"  --split_group="SAS, AIS, NC", --sample_count="5, 5, 5" --graphFileName="result"



#RDP Password: AN#2:KcfVv/6/|x
#student_00_c324bc431