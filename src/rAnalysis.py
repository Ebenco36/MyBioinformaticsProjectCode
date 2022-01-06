from ROperations import Load

def RunRAnalysis(file="SraRunTableIschemicStroke.txt", path_to_quant="quants/", gene_map="gene_map.csv", count_type="salmon", ignoreTxVersion='TRUE', sample_split="1:10, 11:15", split_group="AIS,NC", sample_count="10, 5", result_path="RData3"):
    Load(file, path_to_quant, gene_map, count_type, ignoreTxVersion, sample_split, split_group, sample_count, result_path)
