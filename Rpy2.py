# import rpy2's package module
import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector

def importRPack():
    # import R's "base" package
    base = importr('base')
    # import R's "utils" package
    utils = importr('utils')
    # select a mirror for R packages
    utils.chooseCRANmirror(ind=1) # select the first mirror in the list
    package_names = ('ggplot2', 'hexbin','tidyverse')
    names_to_install = [x for x in package_names if not rpackages.isinstalled(x)]
    print(names_to_install)
    if len(names_to_install) > 0:
        utils.install_packages(StrVector(names_to_install))
    else:
        print("Packages has been installed successfully...")

    if not rpackages.isinstalled("BiocManager"):
        utils.install_packages("BiocManager")
    if rpackages.isinstalled("BiocManager"):
        if not rpackages.isinstalled("tximport"):
            BiocManager = importr('BiocManager')
            BiocManager.install('tximport', ask=False)
        else:
            print('Tximport has been installed already')
    else:
        print("Biocmanager has not been successfully installed")

    # tximport = importr('tximport')

# importRPack()