

import subprocess
from src.helper.getReads import getReads

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
