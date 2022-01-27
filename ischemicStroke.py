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
from MachineLearningV2 import machineLearningv2

from src.downloadSRA import downloadSRA, DownloadSplit
from src.findAdaptors import generateAdaptors
from src.generateIndex import generateIndex, generateKallistoIndex
from src.qualityCheck import fastQcMultiQc
from src.quantification import KallistoQuant, SalmonQuant
from src.removeAdaptors import applyTrimmomatic
from src.buildFilesMap import buildFiles
from src.rAnalysis import RunRAnalysis



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
    parser.add_argument('--machineLearningAnalysisV2', type = str, default = "no", choices = ['yes', 'no'])
    parser.add_argument('--input_file', type = str, default="../../../normalizedCounts.csv")
    parser.add_argument('--graphFileName', type = str, default="result")
    parser.add_argument('--component_n', type = str, default="3")
    parser.add_argument('--path_to_save', type = str, default="machineLearningReport")
    parser.add_argument('--download_sra', type = str, default = "no", choices = ['yes', 'no'])
    parser.add_argument('--run_fastqc', type = str, default = "no", choices = ['yes', 'no'])
    
    


    args = parser.parse_args()

    if args.download_sra == 'yes':
        DownloadSplit(args.sra_run_file, args.srr_column_name, args.read_type, args.fastq_dir, args.sra_accession)

    if args.run_fastqc == 'yes':
        fastQcMultiQc(args.sra_run_file, args.srr_column_name, args.read_type, args.fastqc_dir, args.multiqc_dir, args.fastq_dir, args.sra_accession)
    
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

    
    if args.machineLearningAnalysisV2 == 'yes':
        machineLearningv2(
            file=args.input_file, 
            path_to_save=args.path_to_save, 
            component_n=args.component_n, 
            condition=args.split_group, 
            sample_count=args.sample_count,
            graphFileName=args.graphFileName
        )