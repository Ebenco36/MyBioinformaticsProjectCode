#!/usr/bin/python
#Author: Sun Lei
#Date of the lastest modification: 2014/12/4
#Usage: lncRScan-SVM-train.py -p <positive_GTF|protein_coding> -n \
#   <negative_GTF|lncRNA> -c <GONFIG> -o <output_dir>
#Function: to train a SVM model for classifying protein coding transcripts and lncRNAs 
#Third-part programs required: Biopython, LIBSVM, txCdsPredict, calcConsv and gffread
#Notes: -p is used to input the annotation of protein_coding transcripts, and -n is for lncRNAs

from __future__ import print_function

#############################################################
## parsing the command line
import sys, getopt

try:
    opts, args = getopt.getopt(sys.argv[1:],"hp:n:c:o:",["pGTF=","nGTF=","iConfig=","oDir"])
except getopt.GetoptError:
    # print (str(e))
    print("Usage: %s -p *.p.gtf -n *.n.gtf -c *.conf -o output_dir" % sys.argv[0])
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print("lncRScan-SVM-train.py -p/--pGTF <GTF> \
    -n/--nGTF <GTF> -c/--conf <configure> -o/--output <output_dir>")
        sys.exit()
    elif opt in ("-p", "--pGTF"):
        p_gtf = arg
        if p_gtf[len(p_gtf)-3:len(p_gtf)] != "gtf":
            print("Please input a positive GTF file with suffix 'gtf'.")
            sys.exit()
    elif opt in ("-n", "--nGTF"):
        n_gtf = arg
        if n_gtf[len(n_gtf)-3:len(n_gtf)] != "gtf":
            print("Please input a negative GTF file with suffix 'gtf'.")
            sys.exit()
    elif opt in ("-c", "--conf"):
        in_conf = arg
    elif opt in ("-o", "--output_dir"):
        output = arg
#############################################################

############### Functions ###############
### function 1: to generate svm features ###
def Generate_SVM_Features(gtf, label):
    flag = '0'
    if label == "p":
        flag = "1"
    elif label == "n":
        flag = "-1"
    
    # split the name of the input GTF file of PCTs and LncTs respectively
    #gtf_name = gtf.split(".")
    name = gtf[0:len(gtf)-4]
    
    # create several output files
    # test.gtf2
    gtf_2 = ".".join([name,'gtf2'])
    # FASTA files
    fa_file = ".".join([name,label,'fa'])
    # test.bed
    bed_file = ".".join([name,label,'bed'])
    # Feature files
    features_file = ".".join([name,label,'features'])
    # Feature files having SVM format
    svm_features_file = ".".join([name,label,'svm_features'])
    # test.libsvm.sub
    sub_libsvm_features = ".".join([svm_features_file,'sub'])

    # Extract lines of GTF according to chromosome names
    cmd = "extract_GTF_chr.py " + gtf + " > " + gtf_2
    os.system(cmd)

    # step 1: extract sequences of a given GTF file from the genome sequence 
    #         using gffread of Cuffflinks
    cmd = "gffread -w " + fa_file + " -g " + genome_dir + " " + gtf_2
    os.system(cmd) 
    ## Convert the input GTF to BED
    cmd = "gtf2bed.py -g " + gtf_2 + " -o " + bed_file
    os.system(cmd)

    # step 2: extract features of each transcript from the GTF and FASTA files
    # 
    cmd = "extract_features.py -c " + in_conf + " -b "+ bed_file + " -f " + fa_file + \
            " -o " + features_file
    os.system(cmd)

    # step 3: transfer features to SVM format
    cmd = "features2svm.py " + features_file + " " + flag + " > " + svm_features_file
    os.system(cmd)

    # Extract a subset of features from a standard LIBSVM feature file
    # select features
    cmd = "subsvmfeatures.py " + svm_features_file + " 1,3,6,7,8,9 > " + sub_libsvm_features
    os.system(cmd)
    
    return sub_libsvm_features
### end of generate_svm_features ###

import os

os.system("mkdir -p " + output)
'''
genome_dir = ''
phastcons_bw = ''
svm_model = ''
svm_param = ''
libsvm_mode = ''
scale_low = ''
scale_up = ''
cost = ''
gamma = ''
'''
# read the configure file
conf=open(in_conf)
for line in conf:
    line_array = line.split("\t")
    # 
    if line_array[0] == 'GENOME':
        genome_dir = line_array[1].rstrip()
    elif line_array[0] == 'PHASTCONS':
        phastcons_bw = line_array[1].rstrip()
    elif line_array[0] == 'SVM_MODEL':
        svm_model = line_array[1].rstrip()
    elif line_array[0] == 'SVM_PARAM':
        svm_param = line_array[1].rstrip()
    elif line_array[0] == 'LIBSVM_MODE':
        libsvm_mode = line_array[1].rstrip()
    elif line_array[0] == 'SCALE_L':
        scale_low = line_array[1].rstrip()
    elif line_array[0] == 'SCALE_U':
        scale_up = line_array[1].rstrip()
    elif line_array[0] == 'COST':
        cost = line_array[1].rstrip()
    elif line_array[0] == 'GAMMA':
        gamma = line_array[1].rstrip()
conf.close()

# set SVM mode
#if libsvm_mode == 0:
    ## MODE 0: no scaling && no grid searching
#    scale_low = -1
#    scale_up = 1
#elif libsvm_mode == 1:
    ## MODE 1: scaling && no grid searching
    
#elif libsvm_mode == 2:
    
    ## MODE 2: scaling && grid searching automatically by svm-easy

#### Block 1: processing positive samples sourced from the GTF file of PCTs
p_svm_features = Generate_SVM_Features(p_gtf, 'p')

#### Block 2: processing negative samples sourced from the GTF file of LncTs
n_svm_features = Generate_SVM_Features(n_gtf, 'n')

# #### cat p_svm_features and n_svm_features
# os.system("cat " + p_svm_features + " " + n_svm_features + " > training.svm_feature")

# #### scale svm_features
# os.system("svm-scale  -l " + str(scale_low) + " -u " + str(scale_up) + \
#           " -s training.scale.param training.svm_feature > training.scale")

# #### train a model
# os.system("svm-train -b 1 -c " + str(cost) + " -g " + str(gamma) + " training.scale")

# #### move files to the output directory
# os.system("mv *.p.* *.n.* training.* *.gtf2" + output)
# ## the files are named "training.scale.model" and "training.scale.param"

