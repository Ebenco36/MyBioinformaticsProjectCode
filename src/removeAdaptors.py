import os
import sys
import itertools
import gzip
import subprocess
from src.helper.getReads import getReads


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