

packages = c(
  "tidyverse", 
  "plotly",
  "pheatmap", 
  'optparse'
)

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
## NOT data.frame BUT tibble

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

my_packages <- c(
  "DESeq2", 
  "biomaRt", 
  "clusterProfiler", 
  "org.Hs.eg.db", 
  "tidyverse", 
  "tximport"
) # Specify your packages
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]    # Extract not installed packages
if(length(not_installed)) 
  BiocManager::install(not_installed, ask = FALSE)

#BiocManager::install("DESeq2", ask = FALSE)
#BiocManager::install("biomaRt", ask = FALSE)
#BiocManager::install("clusterProfiler", ask = FALSE)
#BiocManager::install("org.Hs.eg.db", ask = FALSE)
#BiocManager::install("tidyverse", ask = FALSE)
#BiocManager::install("tximport", ask = FALSE)
#install.packages("optparse", repos = "http://cran.us.r-project.org")


# Load packages
library(tximport)
library(DESeq2)
library(plotly)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(readr)

library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,help="run table file", metavar="character"),
  make_option(c("-p", "--path_to_quant"), type="character", default="quants",help="Quantification counts [default= %default]", metavar="character"),
  make_option(c("-g", "--gene_map"), type="character", default="gene_map.csv",help="gene map in csv", metavar="character"),
  make_option(c("-c", "--count_type"), type="character", default="salmon",help="method used [default= %default]", metavar="character"),
  make_option(c("-i", "--ignoreTxVersion"), type="character", default=TRUE,help="dataset file name", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default="RData",help="result", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


RAnalysis <- function(file, path_to_quant, gene_map, count_type, ignoreTxVersion=TRUE, output_dir) {
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    path <- output_dir
    newfolder <- "graphs"
    dir.create(paste0(path, '/', newfolder))
   
  } else {
    path <- output_dir
    newfolder <- "graphs"
    dir.create(paste0(path, '/', newfolder))
    print("Dir already exists!")
  }
  
  pdf(paste0(output_dir,"/graphs/", "Bioinformatics.pdf"))
  sample_table = read_csv(file)
  sample_table = dplyr::select(sample_table, `Run`, BioSample,`Assay Type`,Experiment, Platform, ReleaseDate, Condition)
  sample_table = unique(sample_table)
  view(sample_table)
  
  if (count_type == "salmon")
    sample_files = paste0(paste0(paste0(path_to_quant, pull(sample_table, `Run`)), '_quant'), '/quant.sf')
  else
    sample_files = paste0(paste0(paste0(path_to_quant, pull(sample_table, `Run`)), '_quant'), '/abundance.tsv')
  print(sample_files)
  
  names(sample_files) = pull(sample_table, `Run`)
  
  gene_map = read_csv(gene_map, col_names = c('enstid', 'ensgid'))
  
  count_data = tximport(files = sample_files, type = count_type, tx2gene = gene_map, ignoreTxVersion = ignoreTxVersion)

  print("Splitting counts into section...")
  allCounts <- count_data[['counts']]
  AcuteCounts <- count_data$counts[, 6:10]
  subAcuteCounts <- count_data$counts[, 1:5]
  controlCounts <- count_data$counts[, 11:15]
  print(controlCounts)
  
  count_data[['abundance']]
  
  sample_table = as.data.frame(sample_table)
  sample_table
  colnames(sample_table)[1] = "Sample"
  view(sample_table)
  
  conditions = c('SIS', 'AIS', 'HC')
  conditions = rep(conditions, each=5)
  conditions = factor(conditions)
  sample_table$conditions = conditions
  
  # y ~ x
  dim(count_data)
  dim(sample_table)
  view(sample_table)
  deseq_dataset = DESeqDataSetFromTximport(txi=count_data, 
                                           colData=sample_table,
                                           design=~conditions)
  count_data
  #counts allows us to get the count data from DESEQ dataset
  counts(deseq_dataset)[1:6, 1:3]
  count_data$counts[1:6, 1:3]
  dim(count_data$counts)
  deseq_dataset = estimateSizeFactors(deseq_dataset)
  normalizationFactors(deseq_dataset)
  counts(deseq_dataset, normalized=TRUE)[1:6, 1:3]
  
  boxplot(counts(deseq_dataset, normalized=TRUE))
  print("Graph 1 without Normalization is being generated in PDF...")
  vst = varianceStabilizingTransformation(deseq_dataset)
  boxplot(assay(vst))
  print("Normalized Graph 2 is being generated in PDF...")
  
  # observe the cell line effect

  plotPCA(vst, intgroup='conditions') +
    theme_bw()
  ggsave(paste(output_dir,"/graphs/AllPCA.png",sep=""))
  print("PCA Graph 3 is being generated...")
  
  d = assay(vst)
  assay(vst)
  #transpose using t metric to calculate distance between samples
  d = t(d)
  d = dist(d)
  d
  h = hclust(d)
  plot(h)
  print("Histogram Graph 4 is being generated in PDF...")
  
  #Using kmeans for clustering
  
  k = kmeans(t(assay(vst)), centers=2)
  k$cluster
  
  #Checkout what the normalized values are.
  un_normalizedVal = as.matrix(counts(deseq_dataset, normalized=FALSE))[1:6, 1:3]
  normalized_factors = normalizationFactors(deseq_dataset)[1:6, 1:3]
  
  normalizedVal = un_normalizedVal/normalized_factors
  normalizedVal
  #compare values with this below
  comparingVal = counts(deseq_dataset, normalized=TRUE)
  comparingVal[1:6, 1:3]
  
  #write to file for futher analysis
  write.csv(comparingVal, paste(output_dir,"/normalizedCounts.csv",sep=""))
  write.csv(as.matrix(counts(deseq_dataset, normalized=FALSE)), paste(output_dir,"/UnnormalizedCounts.csv",sep=""))
  
  
  # 3 steps to DESeq2 analysis
  # 1) estimate size factors (normalisation)
  # 2) estimate dispersions
  # 3) apply statistics (Wald Test)
  
  
  dispersion = estimateDispersions(deseq_dataset)
  
  
  # plotDispEsts(dispersion)
  print("Moving forward...")
  
  # step 3
  statistics_wald_test = nbinomWaldTest(dispersion)
  # plot(statistics_wald_test)
  
  ## DESeq2 shortcut
  # dds_nhbe = DESeq(dds_nhbe)
  result_table = results(statistics_wald_test)
  summary(result_table)
  # View(as.data.frame(result_table))
  #----------------------------------------------------------------------------------Continue
  # result_table is a DataFrame not a data.frame!
  result_df = as.data.frame(result_table)
  # View(result_df)
  write.csv(result_df, paste(output_dir,"/TtestValuesCounts.csv",sep=""))
  
  plotCounts(statistics_wald_test, gene='ENSG00000060709',
             intgroup='conditions')
  print("With Sample ENSG0000000060709 Graph 5 is being generated...")
  sum(!complete.cases(result_df))
  
  filter_df1 = result_df[complete.cases(result_df),]
  #Completed cases
  write.csv(filter_df1, paste(output_dir,"/TtestValuesCountsOptimizedCompleteCases.csv",sep=""))
  #uncompleted Cases
  write.csv(result_df[!complete.cases(result_df), ], paste(output_dir,"/TtestValuesCountsOptimizedUnCompleteCases.csv",sep=""))
  # View(filter_df1)
  
  # Filter results 
  # padj < 0.05
  # log2FoldChange > 1 < -1
  
  filter_df1$padj < 0.05
  
  filter_df2 = filter_df1[filter_df1$padj < 0.05,]
  write.csv(filter_df2, paste(output_dir,"/FilterpadJLessThan005.csv",sep=""))
  dim(filter_df2)
  
  abs(filter_df2$log2FoldChange) > 1
  
  filter_df3 = filter_df2[abs(filter_df2$log2FoldChange) > 1,]
  write.csv(filter_df3, paste(output_dir,"/FilterGreaterThan1withpadjlessthen005.csv",sep=""))
  dim(filter_df3)
  # View(filter_df3)
  
  plotMA(result_table)
  print("MA Plot Graph 6 is being generated in PDF...")
  
  # volcano plot
  
  filter_df1$test = filter_df1$padj < 0.05 & abs(filter_df1$log2FoldChange) > 1
  filter_df1
  
  filter_df1 = rownames_to_column(filter_df1, var='ensgene')
  filter_df1
  write.csv(filter_df1, paste(output_dir,"/rownametocolumn.csv",sep=""))
  g = ggplot(filter_df1, aes(x=log2FoldChange, y=-log10(padj), name=ensgene)) +
    geom_point(aes(colour=test), size=1, alpha=0.3) +
    scale_colour_manual(values=c('black', 'red')) +
    geom_vline(xintercept=1, colour='green', linetype=3) +
    geom_vline(xintercept=-1, colour='green', linetype=3) +
    geom_hline(yintercept=-log10(0.05), colour='blue', linetype=3) +
    theme_bw() +
    theme(legend.position = 'none')
  
  ggplotly(g)
  ggsave(paste(output_dir,"/graphs/LogFoldChangePlot.png",sep=""))
  print("Log Fold Change Graph 7 is being generated...")
  #install BioMart
  #if (!requireNamespace("BiocManager", quietly = TRUE))
  #  install.packages("BiocManager")
  
  #BiocManager::install("biomaRt")
  
  print("Gene Annotation using BioMaRT from Esemble")
  library(biomaRt)
  
  listMarts()
  ensembl99 = useEnsembl(biomart="ensembl", version=101)
  # View(listDatasets(ensembl99))
  ensembl99 = useDataset("hsapiens_gene_ensembl", 
                         mart=ensembl99)
  # View(listAttributes(ensembl99))
  # View(listFilters(ensembl99))
  getBM(attributes=c('ensembl_gene_id', 'ensembl_gene_id_version',
                     'ensembl_transcript_id', 'ensembl_transcript_id_version',
                     'external_gene_name'), filters = c('ensembl_gene_id'), 
        values = filter_df1$ensgene[1:6],
        mart = ensembl99)
  
  annotation = getBM(attributes=c('ensembl_gene_id',
                                  'chromosome_name',
                                  'start_position',
                                  'end_position',
                                  'strand',
                                  'gene_biotype',
                                  'external_gene_name',
                                  'description'),
                     filters = c('ensembl_gene_id'),
                     values = filter_df1$ensgene,
                     mart = ensembl99)
  write.csv(annotation, paste(output_dir,"/BioMartInfo.csv",sep=""))
  nrow(annotation)
  # View(annotation)
  # View(filter_df1)
  #join using dplyer
  annotated_df = left_join(filter_df1, annotation,
                           by=c('ensgene'='ensembl_gene_id'))
  write.csv(annotated_df, paste(output_dir,"/BioMartInfoJoin.csv",sep=""))
  # View(annotated_df)
  
  g = ggplot(annotated_df, aes(x=log2FoldChange, 
                               y=-log10(padj), 
                               name=external_gene_name)) +
    geom_point(aes(colour=test), size=1, alpha=0.3) +
    scale_colour_manual(values=c('black', 'red')) +
    geom_vline(xintercept=1, colour='green', linetype=3) +
    geom_vline(xintercept=-1, colour='green', linetype=3) +
    geom_hline(yintercept=-log10(0.05), colour='blue', linetype=3) +
    theme_bw() +
    theme(legend.position = 'none')
  
  
  ggplotly(g)
  ggsave(paste(output_dir,"/graphs/LogFoldChangePlotWithExternalNames.png",sep=""))
  print("Log Fold Change With names Graph 8 is being generated...")
  
  anno_df2 = annotated_df[annotated_df$padj < 0.05,]
  
  anno_df3 = anno_df2[abs(anno_df2$log2FoldChange) > 1,]
  
  write.csv(anno_df3, paste(output_dir,"/AnnotatedDFForDEG.csv",sep=""))
  degs = anno_df3$ensgene
  degs
  
  #vst_nhbe = varianceStabilizingTransformation(vst)
  vst_mat = assay(vst)
  
  #Just for the differential expressed genes
  data_for_hm = vst_mat[degs,]
  data_for_hm
  #data for Heatmap Stored
  write.csv(data_for_hm, paste(output_dir,"/HeatMapDataForDEG.csv",sep=""))
  #change the Gene ID to their respective names
  
  rownames(data_for_hm) = anno_df3$external_gene_name
  write.csv(data_for_hm, paste(output_dir,"/HeatMapDataForDEGWithNames.csv",sep=""))
  data_for_hm
  dim(data_for_hm)
  heatmap(data_for_hm)
  print("HeatMap Graph 9 is being generated in PDF...")
  
  #get Better Label
  pheatmap(data_for_hm, fontsize_row=4, scale='row')
  print("PheatMap 1 Graph 10 is being generated in PDF...")
  #To use color brewer we need to install it
  #BioManager::install('RColorBrewer')
  #see all color pelette available
  #display.brewer.all()
  
  greys = colorRampPalette(brewer.pal(9, "Greys"))(100)
  
  pheatmap(data_for_hm, fontsize_row=4, scale='row',
           color=greys)
  print("PheatMap Grey 2 Graph 11 is being generated in PDF...")
  
  #pairs = colorRampPalette(brewer.pal(12, "Paired"))(100)
  
  #pheatmap(data_for_hm, fontsize_row=4, scale='row',color=pairs)
  
  last_scheme = colorRampPalette(brewer.pal(7, "Blues"))(100)
  
  pheatmap(data_for_hm, fontsize_row=4, scale='row',
           color=last_scheme, cutree_cols = 2,
           cutree_rows = 2)
  print("PheatMap 3 Graph 12 is being generated in PDF...")
  ## GO enrichment
  
  ent_gene = getBM(attributes=c('entrezgene_id'),
                   filters = c('ensembl_gene_id'),
                   values = anno_df3$ensgene,
                   mart = ensembl99)
  ent_gene = ent_gene$entrezgene_id
  ent_gene = as.character(ent_gene)
  ent_uni = getBM(attributes=c('entrezgene_id'),
                  filters = c('ensembl_gene_id'),
                  values = annotated_df$ensgene,
                  mart = ensembl99)
  ent_uni = as.character(ent_uni$entrezgene_id)
  ent_uni = as.character(ent_uni)
  print("We are here 2nd last")
  #BiocManager::install("enrichGO")
  #BiocManager::install("enrichplot")
  #BiocManager::install("clusterProfiler")
  #BiocManager::install("org.Hs.eg.db")
  library(org.Hs.eg.db)
  library(clusterProfiler)
  ego = enrichGO(gene = ent_gene,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 universe = ent_uni)
  
  ego
  print("We are here last ")
  # barplot(ego, showCategory=20)
  # dotplot(ego, showCategory=20)
  # cnetplot(ego, showCategory=5)
  # goplot(ego)
  
  ekg = enrichKEGG(gene = ent_gene,
                   universe = ent_uni)
  
  write_tsv(anno_df3, paste(output_dir,"/filtered_nhbe_results.txt",sep=""))
  
  dev.off()
  
}
RAnalysis(file=opt$file, path_to_quant=opt$path_to_quant, gene_map=opt$gene_map, count_type=opt$count_type, ignoreTxVersion=opt$ignoreTxVersion, opt$output_dir)
