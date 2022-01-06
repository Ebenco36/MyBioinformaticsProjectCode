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