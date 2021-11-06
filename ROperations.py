import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector
import rpy2.robjects as robjects
from Rpy2 import importRPack
import subprocess
def Load(file="SraRunTableIschemicStroke.txt", path_to_quant="quants/", gene_map="gene_map.csv", count_type="salmon", ignoreTxVersion='TRUE', sample_split="1:10, 11:15", split_group="AIS,NC", sample_count="10, 5"):
    base = importr('base')
    if rpackages.isinstalled("tidyverse") and rpackages.isinstalled("tximport"):
        subprocess.call(['Rscript', 'BioinformaticsAnalysisP.R', '-f', file, '-p', path_to_quant, '-g', gene_map, '-c', count_type, '-i', ignoreTxVersion, '-o', "RData3", '-S', sample_split, '-C', split_group, '-x', sample_count ])
    else:
        importRPack()

count_type="salmon"
ignoreTxVersion = True
tx2gene = ""
# Load()