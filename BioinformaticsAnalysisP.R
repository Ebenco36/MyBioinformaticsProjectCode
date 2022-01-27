

packages = c("tidyverse", "plotly", "pheatmap", 'optparse')

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

my_packages <- c("DESeq2", "biomaRt", "clusterProfiler", "org.Hs.eg.db", "tidyverse", "tximport") # Specify your packages
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
library(stringi)
library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,help="run table file", metavar="character"),
  make_option(c("-p", "--path_to_quant"), type="character", default="quants",help="Quantification counts [default= %default]", metavar="character"),
  make_option(c("-g", "--gene_map"), type="character", default="gene_map.csv",help="gene map in csv", metavar="character"),
  make_option(c("-c", "--count_type"), type="character", default="salmon",help="method used [default= %default]", metavar="character"),
  make_option(c("-i", "--ignoreTxVersion"), type="character", default=TRUE,help="dataset file name", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default="RData",help="result", metavar="character"),
  make_option(c("-S", "--slice_group"), type="character", default="1:5, 11:15",help="select range for analysis group1 = 1:5, 11:15", metavar="character"),
  make_option(c("-C", "--condition"), type="character", default="SUB,NC",help="select range for analysis Example acute,sub", metavar="character"),
  make_option(c("-x", "--sample_count"), type="character", default="5, 5",help="sample count for each class", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


RAnalysis <- function(file, path_to_quant, gene_map, count_type, ignoreTxVersion=TRUE, output_dir, slice_group, condition, sample_count) {
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
  output_dir <- stri_trim_both(output_dir)
  pdf(paste0(output_dir,"/graphs/", "Bioinformatics.pdf"))
  sample_table = read_csv(file)
  print(sample_table)
  sample_table = dplyr::select(sample_table, `Run`, BioSample,`Assay Type`,Experiment, Platform, ReleaseDate, Condition)
  sample_table = unique(sample_table)
  #slice_group = "1:5, 6:10, 11:15"
  datat <- strsplit(slice_group, ",")[[1]]
  sample_table_sub = data.frame()
  for (splitS in datat){
    resp = strsplit(splitS, ":")[[1]]
    
    sample_table_sub <- rbind(sample_table_sub, sample_table[resp[1]:resp[2],])
  }
  view(sample_table)
  
  if (count_type == "salmon")
    sample_files_sub = paste0(paste0(paste0(path_to_quant, pull(sample_table_sub, `Run`)), '_quant'), '/quant.sf')
  else
    sample_files_sub = paste0(paste0(paste0(path_to_quant, pull(sample_table_sub, `Run`)), '_quant'), '/abundance.tsv')
  print(sample_files_sub)
  
  names(sample_files_sub) = pull(sample_table_sub, `Run`)
  
  gene_map = read_csv(gene_map, col_names = c('enstid', 'ensgid'))
  print(sample_files_sub)
  count_data_sub = tximport(files = sample_files_sub, type = count_type, tx2gene = gene_map, ignoreTxVersion = ignoreTxVersion)
  #print(count_data_sub)
  
  
  sample_table_sub = as.data.frame(sample_table_sub)
  sample_table_sub
  colnames(sample_table_sub)[1] = "Sample"
  view(sample_table_sub)
  
  condition = unlist(strsplit(condition,","))
  sample_c = unlist(strsplit(sample_count,","))
  data <- c()
  # condition = unlist(strsplit("AS, DS, FG",","))
  # sample_c = unlist(strsplit("3, 4, 3",","))
  for (x in 1:length(condition)) {
    myStr <- condition[x]
    myNum <- strtoi(sample_c[x])
    ins <- rep(myStr, myNum)
    data <- c(data,ins)
    print(ins)
  }
  print(c(data))

  conditions_sub = c(data)
  #conditions_sub = rep(conditions_sub, each=5)
  conditions_sub = factor(conditions_sub)
  sample_table_sub$conditions = conditions_sub
  sample_table$conditions = factor(c(data))
  
  split_for_header <- strsplit(slice_group, ",")[[1]]
  tail_split <- tail(split_for_header, n=1)
  header_split <-head(split_for_header, n=-1)
  header_names <- sample_table[strsplit(header_split, ":")[[1]][1]:strsplit(tail_split, ":")[[1]][2],]
  print("Header names")
  print(header_names)
  print(sample_table_sub$conditions)
  
  # y ~ x
  deseq_dataset_sub = DESeqDataSetFromTximport(txi=count_data_sub, 
                                           colData=sample_table_sub,
                                           design=~conditions)
  count_data_sub
  #counts allows us to get the count data from DESEQ dataset
  deseq_dataset_sub = estimateSizeFactors(deseq_dataset_sub)
  normalizationFactors(deseq_dataset_sub)
  counts(deseq_dataset_sub, normalized=TRUE)[1:6, 1:3]
  deseq_dataset_sub_normalized  <- counts(deseq_dataset_sub, normalized=TRUE)
  counts(deseq_dataset_sub)
  
  
  #Visuals
  # Column sun Analysis
  colSums(counts(deseq_dataset_sub))
  # Un-normalized count
  barplot(colSums(counts(deseq_dataset_sub)))
  
  # Normalized
  
  barplot(counts(deseq_dataset_sub, normalized=TRUE))
  
  counts(deseq_dataset_sub, normalized=TRUE)[1:6, 1:3]
  
  # Check how similar replicates are to each other
  Test_data = log2(counts(deseq_dataset_sub, normalized=TRUE))
  #Similar
  plot(Test_data[,7], Test_data[, 6])
  
  #Not similar
  plot(Test_data[,1], Test_data[, 6])
  
  #The same
  plot(Test_data[,1], Test_data[, 1])
  
  
  colnames(deseq_dataset_sub) <- sample_table_sub$Condition
  png(file=paste(output_dir,"/graphs/BoxPlot_for_normalizedCount.png",sep=""))
  boxplot(counts(deseq_dataset_sub, normalized=TRUE), names=header_names$conditions)
  title(main = "Boxplot for normalized gene counts")
  dev.off()
  
  print("Graph 1 without Normalization is being generated in PDF...")
  png(file=paste(output_dir,"/graphs/BoxPlot_for_VariantTransformedCount.png",sep=""))
  vst_sub = varianceStabilizingTransformation(deseq_dataset_sub)
  boxplot(assay(vst_sub), names=header_names$conditions)
  title(main = "Boxplot for variance Stabilized Transformed gene counts")
  dev.off()
  print("Normalized Graph 2 is being generated in PDF...")
  
  # observe the cell line effect
  png(file=paste(output_dir,"/graphs/PCAPlot.png",sep=""))
  plotPCA(vst_sub, intgroup='conditions') +
    theme_bw()
  dev.off()
  ggsave(paste(output_dir,"/graphs/AllPCA.png",sep=""))
  print("PCA Graph 3 is being generated...")
  
  d_sub = assay(vst_sub)
  assay(vst_sub)
  #transpose using t metric to calculate distance between samples
  d_sub = t(d_sub)
  d_sub = dist(d_sub)
  d_sub
  h_sub = hclust(d_sub)
  png(file=paste(output_dir,"/graphs/HistogramPlot.png",sep=""))
  plot(h_sub)
  dev.off()
  print("Histogram Graph 4 is being generated in PDF...")
  
  #Using kmeans for clustering
  
  k_sub = kmeans(t(assay(vst_sub)), centers=2)
  k_sub$cluster
  
  #Checkout what the normalized values are.
  un_normalizedVal = as.matrix(counts(deseq_dataset_sub, normalized=FALSE))[1:6, 1:3]
  normalized_factors = normalizationFactors(deseq_dataset_sub)[1:6, 1:3]
  
  normalizedVal = un_normalizedVal/normalized_factors
  normalizedVal
  deseq_dataset_sub_normalized  <- counts(deseq_dataset_sub, normalized=TRUE)
  
  #write to file for futher analysis
  write.csv(deseq_dataset_sub_normalized, paste(output_dir,"/normalizedCounts.csv",sep=""))
  write.csv(as.matrix(counts(deseq_dataset_sub, normalized=FALSE)), paste(output_dir,"/UnnormalizedCounts.csv",sep=""))
  
  colSums((as.matrix(counts(deseq_dataset_sub, normalized=TRUE))))
  
  colSums((as.matrix(counts(deseq_dataset_sub, normalized=FALSE))))
  
  # 3 steps to DESeq2 analysis
  # 1) estimate size factors (normalisation)
  # 2) estimate dispersions
  # 3) apply statistics (Wald Test)
  
  
  dispersion_sub = estimateDispersions(deseq_dataset_sub)
  
  
  # plotDispEsts(dispersion)
  print("Moving forward...")
  
  # step 3
  statistics_wald_test_sub = nbinomWaldTest(dispersion_sub)
  # plot(statistics_wald_test_sub)
  
  ## DESeq2 shortcut
  # dds_nhbe = DESeq(dds_nhbe)
  result_table_sub = results(statistics_wald_test_sub)
  summary(result_table_sub)
  # View(as.data.frame(result_table_sub))
  #----------------------------------------------------------------------------------Continue
  # result_table is a DataFrame not a data.frame!
  result_df_sub = as.data.frame(result_table_sub)
  # View(result_df_sub)
  write.csv(result_df_sub, paste(output_dir,"/TtestValuesCounts.csv",sep=""))
  png(file=paste(output_dir,"/graphs/SingleTestusingSingleGene.png",sep=""))
  plotCounts(statistics_wald_test_sub, gene='ENSG00000060709',
             intgroup='conditions')
  dev.off()
  print("With Sample ENSG0000000060709 Graph 5 is being generated...")
  sum(!complete.cases(result_df_sub))
  
  filter_df = result_df_sub[complete.cases(result_df_sub),]
  #Completed cases
  write.csv(filter_df, paste(output_dir,"/TtestValuesCountsOptimizedCompleteCases.csv",sep=""))
  #uncompleted Cases
  write.csv(result_df_sub[!complete.cases(result_df_sub), ], paste(output_dir,"/TtestValuesCountsOptimizedUnCompleteCases.csv",sep=""))
  # View(filter_df)
  
  # Filter results 
  # padj < 0.05
  # log2FoldChange > 1 < -1
  
  filter_df$padj < 0.05
  
  filter_df_sub = filter_df[filter_df$padj < 0.05,]
  write.csv(filter_df_sub, paste(output_dir,"/FilterpadJLessThan005.csv",sep=""))
  dim(filter_df_sub)
  
  abs(filter_df_sub$log2FoldChange) > 1
  
  filter_df_sub = filter_df_sub[abs(filter_df_sub$log2FoldChange) > 1,]
  write.csv(filter_df_sub, paste(output_dir,"/FilterGreaterThan1withpadjlessthen005.csv",sep=""))
  dim(filter_df_sub)
  # View(filter_df_sub)
  png(file=paste(output_dir,"/graphs/MAPlot.png",sep=""))
  plotMA(result_table_sub)
  dev.off()
  print("MA Plot Graph 6 is being generated in PDF...")
  
  # volcano plot
  
  filter_df_sub = rownames_to_column(filter_df_sub, var='ensgene')
  
  filter_df = rownames_to_column(filter_df, var='ensgene')
  filter_df$test = filter_df$padj < 0.05 & abs(filter_df$log2FoldChange) > 1
  abs(filter_df$log2FoldChange) 
  write.csv(filter_df_sub, 'rownametocolumnSub.csv')
  
  
  png(file=paste(output_dir,"/graphs/VolcanoPlot.png",sep=""))
  
  g_sub = ggplot(filter_df, aes(x=log2FoldChange, y=-log10(padj), name=ensgene)) +
    geom_point(aes(colour=test), size=1, alpha=0.3) +
    scale_colour_manual(values=c('black', 'red')) +
    geom_vline(xintercept=1, colour='green', linetype=3) +
    geom_vline(xintercept=-1, colour='green', linetype=3) +
    geom_hline(yintercept=-log10(0.05), colour='blue', linetype=3) +
    theme_bw() +
    theme(legend.position = 'none')
  
  #ggplotly(g)
  
  print(g_sub + ggtitle("Differential Expressed Genes"))
  ggsave(paste(output_dir,"/graphs/LogFoldChangePlot.png",sep=""))
  dev.off()
  ordered_data = filter_df[order(abs(filter_df$log2FoldChange), decreasing=TRUE),]
  # View(as.data.frame(ordered_data))
  #ggplotly(g_sub)
  
  #install BioMart
  #if (!requireNamespace("BiocManager", quietly = TRUE))
  #  install.packages("BiocManager")
  
  #BiocManager::install("biomaRt")
  print("We are here!")
  library(biomaRt)
  
  listMarts()
  ensembl104 = useEnsembl(biomart="ensembl", version=104)
  # View(listDatasets(ensembl104))
  ensembl104 = useDataset("hsapiens_gene_ensembl", 
                          mart=ensembl104)
  # View(listAttributes(ensembl104))
  # View(listFilters(ensembl104))
  getBM(attributes=c('ensembl_gene_id', 'ensembl_gene_id_version',
                     'ensembl_transcript_id', 'ensembl_transcript_id_version',
                     'external_gene_name'), 
        filters = c('ensembl_gene_id'), 
        values = filter_df_sub$ensgene[1:6],
        mart = ensembl104)
  
  dim(filter_df_sub)
  dim(filter_df_sub[filter_df_sub$log2FoldChange < 1, ])
  dim(filter_df_sub[filter_df_sub$log2FoldChange > 1, ])
  all_regulated = filter_df_sub
  up_regulated = filter_df_sub[filter_df_sub$log2FoldChange < 1, ]
  down_regulated = filter_df_sub[filter_df_sub$log2FoldChange > 1, ]
  
  
  annotation_ = getBM(attributes=c('ensembl_gene_id',
                                   'chromosome_name',
                                   'start_position',
                                   'end_position',
                                   'strand',
                                   'gene_biotype',
                                   'external_gene_name',
                                   'description'),
                      filters = c('ensembl_gene_id'),
                      values = filter_df$ensgene,
                      mart = ensembl104)
  
  # View(annotation_)
  annotation_sub = getBM(attributes=c('ensembl_gene_id',
                                      'chromosome_name',
                                      'start_position',
                                      'end_position',
                                      'strand',
                                      'gene_biotype',
                                      'external_gene_name',
                                      'description'),
                         filters = c('ensembl_gene_id'),
                         values = all_regulated$ensgene,
                         mart = ensembl104)
  
  annotation_sub_upregulated = getBM(attributes=c('ensembl_gene_id',
                                                  'chromosome_name',
                                                  'start_position',
                                                  'end_position',
                                                  'strand',
                                                  'gene_biotype',
                                                  'external_gene_name',
                                                  'description'),
                                     filters = c('ensembl_gene_id'),
                                     values = up_regulated$ensgene,
                                     mart = ensembl104)
  
  annotation_sub_downregulated = getBM(attributes=c('ensembl_gene_id',
                                                    'chromosome_name',
                                                    'start_position',
                                                    'end_position',
                                                    'strand',
                                                    'gene_biotype',
                                                    'external_gene_name',
                                                    'description'),
                                       filters = c('ensembl_gene_id'),
                                       values = down_regulated$ensgene,
                                       mart = ensembl104)
  
  
  write.csv(annotation_sub, paste(output_dir,"/BioMartInfo.csv",sep=""))
  nrow(annotation_sub)
  # View(annotation_sub)
  # View(filter_df1)
  #join using dplyer
  annotated_df_ = left_join(filter_df, annotation_,
                            by=c('ensgene'='ensembl_gene_id'))
  
  
  annotated_df_sub = left_join(all_regulated, annotation_sub,
                               by=c('ensgene'='ensembl_gene_id'))
  # View(annotated_df_)
  annotated_df_sub_up_regulated = left_join(up_regulated, annotation_sub_upregulated,
                                            by=c('ensgene'='ensembl_gene_id'))
  annotated_df_sub_down_regulated = left_join(down_regulated, annotation_sub_downregulated,
                                              by=c('ensgene'='ensembl_gene_id'))
  
  annotated_df_lncRNA_up_reg = filter(left_join(up_regulated, annotation_sub_upregulated,
                                         by=c('ensgene'='ensembl_gene_id')), gene_biotype=="lncRNA")
  annotated_df_lncRNA_down_reg = filter(left_join(down_regulated, annotation_sub_downregulated,
                                                by=c('ensgene'='ensembl_gene_id')), gene_biotype=="lncRNA")
  
  annotated_df_protein_coding_up_reg = filter(left_join(up_regulated, annotation_sub_upregulated,
                                                 by=c('ensgene'='ensembl_gene_id')), 
                                       gene_biotype=="protein_coding")
  
  annotated_df_protein_coding_down_reg = filter(left_join(down_regulated, annotation_sub_downregulated,
                                                        by=c('ensgene'='ensembl_gene_id')), 
                                              gene_biotype=="protein_coding")
  
  
  annotated_df_un_identified= filter(left_join(all_regulated, annotation_sub,
                                               by=c('ensgene'='ensembl_gene_id')), 
                                     gene_biotype!="protein_coding" & gene_biotype!="lncRNA")
  
  # View(annotated_df_sub_down_regulated)
  write.csv(annotated_df_lncRNA_up_reg, paste(output_dir,"/LncRNAUpregulated.csv",sep=""))
  write.csv(annotated_df_lncRNA_down_reg, paste(output_dir,"/LncRNADownregulated.csv",sep=""))
  write.csv(annotated_df_protein_coding_up_reg, paste(output_dir,"/ProteinCodingUpregulated.csv",sep=""))
  write.csv(annotated_df_protein_coding_down_reg, paste(output_dir,"/ProteinCodingDownregulated.csv",sep=""))
  write.csv(annotated_df_un_identified, paste(output_dir,"/unidentifiedregulated.csv",sep=""))
  
  write.csv(annotated_df_sub_up_regulated, paste(output_dir,"/Upregulated.csv",sep=""))
  write.csv(annotated_df_sub_down_regulated, paste(output_dir,"/Downregulated.csv",sep=""))
  write.csv(annotated_df_sub, paste(output_dir,"/All_regulated.csv",sep=""))
  
  resss = c(
      paste('lncRNA_up_regulated : ',sum(complete.cases(annotated_df_lncRNA_up_reg)),sep=""),
      paste('lncRNA_down_regulated : ', sum(complete.cases(annotated_df_lncRNA_down_reg)),sep=""),
      paste('Protein_coding_up_regulated : ',sum(complete.cases(annotated_df_protein_coding_up_reg)),sep=""),
      paste('Protein_coding_down_regulated : ', sum(complete.cases(annotated_df_protein_coding_down_reg)),sep=""),
      paste('Unidentified_regulated : ',sum(complete.cases(annotated_df_un_identified)),sep=""),
      paste('up_regulated : ',sum(complete.cases(annotated_df_sub_up_regulated)),sep=""),
      paste('down_regulated : ', sum(complete.cases(annotated_df_sub_down_regulated)),sep=""),
      paste('all_regulated : ', sum(complete.cases(annotated_df_sub)),sep="")
  )
  
  write.csv(resss, paste(output_dir,"/Summary_regulated.csv",sep=""))
  
  write.csv(annotated_df_sub[annotated_df_sub$test == TRUE,], paste(output_dir,"/BioMartInfoJoin.csv",sep=""))
  degs_sub = annotated_df_sub$ensgene
  
  png(file=paste(output_dir,"/graphs/VolcanoAnnotated.png",sep=""))
  g_sub = ggplot(annotated_df_, aes(x=log2FoldChange, 
                               y=-log10(padj), 
                               name=external_gene_name)) +
    geom_point(aes(colour=test), size=1, alpha=0.3) +
    scale_colour_manual(values=c('black', 'red')) +
    geom_vline(xintercept=1, colour='green', linetype=3) +
    geom_vline(xintercept=-1, colour='green', linetype=3) +
    geom_hline(yintercept=-log10(0.05), colour='blue', linetype=3) +
    theme_bw() +
    theme(legend.position = 'none')
  
  plotCounts(statistics_wald_test_sub, 'ENSG00000004142', intgroup='conditions')
  ggplotly(g_sub)
  
  ggsave(paste(output_dir,"/graphs/LogFoldChangePlotWithExternalNames.png",sep=""))
  dev.off()
  print("Log Fold Change With names Graph 8 is being generated...")
  
  
  #vst_nhbe = varianceStabilizingTransformation(vst)
  vst_mat_sub = assay(vst_sub)
  
  
  # View(vst_mat_sub)
  #Just for the differential expressed genes
  
  data_for_hm_sub = vst_mat_sub[degs_sub,]
  
  data_for_hm_sub
  
  #data for Heatmap Stored
  write.csv(data_for_hm_sub, paste(output_dir,"/HeatMapDataForDEG.csv",sep=""))
  #change the Gene ID to their respective names
  
  rownames(data_for_hm_sub) = annotated_df_sub$external_gene_name
  write.csv(data_for_hm_sub, paste(output_dir,"/HeatMapDataForDEGWithNames.csv",sep=""))
  data_for_hm_sub
  dim(data_for_hm_sub)
  png(file=paste(output_dir,"/graphs/HeatMap.png",sep=""), width=1200, height=1400)
  heatmap(data_for_hm_sub,  main = "heatmap row cluster for gene expression data", 
          scale = 'row', legend = F, fontsize = 12, cexRow = .9, cexCol = .9,)
  dev.off()
  print("HeatMap Graph 9 is being generated in PDF...")
  
  #get Better Label
  png(file=paste(output_dir,"/graphs/PheatMapPlot.png",sep=""), width=1200, height=1000)
  pheatmap(data_for_hm_sub, fontsize=12, scale='row', main = "pheatmap row cluster for gene expression data")
  dev.off()
  print("PheatMap 1 Graph 10 is being generated in PDF...")
  #To use color brewer we need to install it
  #BioManager::install('RColorBrewer')
  #see all color pelette available
  #display.brewer.all()
  
  png(file=paste(output_dir,"/graphs/ColoredPheatMapPlot.png",sep=""), width=1200, height=1000)
  greys = colorRampPalette(brewer.pal(9, "Greys"))(100)
  
  pheatmap(data_for_hm_sub, fontsize=12, scale='row',
           color=greys, main = "pheatmap grey format")
  dev.off()
  print("PheatMap Grey 2 Graph 11 is being generated in PDF...")
  
  #pairs = colorRampPalette(brewer.pal(12, "Paired"))(100)
  
  #pheatmap(data_for_hm_sub, fontsize_row=4, scale='row',color=pairs)
  
  png(file=paste(output_dir,"/graphs/PheatMapForDemercation.png",sep=""), width=1200, height=1000)
  last_scheme = colorRampPalette(brewer.pal(7, "Blues"))(100)
  
  pheatmap(data_for_hm_sub, fontsize = 12, scale='row',
           color=last_scheme, cutree_cols = 2,
           cutree_rows = 2, main = "pheatmap row cluster with boundary set")
  dev.off()
  print("PheatMap 3 Graph 12 is being generated in PDF...")
  
  
  
  ## GO enrichment
  #------------------------------------------------STart working from here
  ent_gene = getBM(attributes=c('entrezgene_id'),
                   filters = c('ensembl_gene_id'),
                   values =  annotated_df_$ensgene,
                   mart = ensembl104)
  ent_gene
  
  
  ent_gene_all = getBM(attributes=c('entrezgene_id'),
                       filters = c('ensembl_gene_id'),
                       values = annotated_df_sub$ensgene,
                       mart = ensembl104)
  print("Here 1")
  
  ent_gene_up = getBM(attributes=c('entrezgene_id'),
                      filters = c('ensembl_gene_id'),
                      values = annotated_df_sub_up_regulated$ensgene,
                      mart = ensembl104)
  # annotated_df_sub_up_regulated
  ent_gene_down = getBM(attributes=c('entrezgene_id'),
                        filters = c('ensembl_gene_id'),
                        values = annotated_df_sub_down_regulated$ensgene,
                        mart = ensembl104)
  print("Here 2")
  (ent_gene_down$entrezgene_id)
  
  ent_gene = ent_gene$entrezgene_id
  ent_gene_all = ent_gene_all$entrezgene_id
  ent_gene_up = ent_gene_up$entrezgene_id
  ent_gene_down = ent_gene_down$entrezgene_id
  
  ent_gene = as.character(ent_gene)
  ent_gene_all = as.character(ent_gene_all)
  ent_gene_up = as.character(ent_gene_up)
  ent_gene_down = as.character(ent_gene_down)
  print("Here 3")
  ent_uni = getBM(attributes=c('entrezgene_id'),
                  filters = c('ensembl_gene_id'),
                  values = annotated_df_$ensgene,
                  mart = ensembl104)
  ent_uni_all = getBM(attributes=c('entrezgene_id'),
                      filters = c('ensembl_gene_id'),
                      values = annotated_df_sub$ensgene,
                      mart = ensembl104)
  ent_uni_up = getBM(attributes=c('entrezgene_id'),
                     filters = c('ensembl_gene_id'),
                     values = annotated_df_sub_up_regulated$ensgene,
                     mart = ensembl104)
  ent_uni_down = getBM(attributes=c('entrezgene_id'),
                       filters = c('ensembl_gene_id'),
                       values = annotated_df_sub_down_regulated$ensgene,
                       mart = ensembl104)
  print("Here 4")
  # View(annotated_df_sub)
  
  
  
  
  ent_uni = as.character(as.character(ent_uni$entrezgene_id))
  ent_uni_all = as.character(as.character(ent_uni_all$entrezgene_id))
  ent_uni_up = as.character(as.character(ent_uni_up$entrezgene_id))
  ent_uni_down = as.character(as.character(ent_uni_down$entrezgene_id))
  
  ent_uni
  annotated_df_sub
  # BiocManager::install("enrichGO")
  # BiocManager::install("enrichplot")
  # BiocManager::install("clusterProfiler")
  #BiocManager::install("org.Hs.eg.db")
  print("Here 5")
  library(org.Hs.eg.db)
  library(clusterProfiler)
  
  # ego <- clusterProfiler::enrichGO(gene = ent_gene_all,
  #                OrgDb = org.Hs.eg.db,
  #                ont = "BP",
  #                keyType = "ENTREZID",
  #                universe = ent_uni, readable=TRUE)
  # as.data.frame(ego)
  
  ego_all <- enrichGO(gene = ent_gene_all,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      universe = ent_uni, readable=TRUE, pvalueCutoff=0.01)
  ego_up <- enrichGO(gene = ent_gene_up,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     universe = ent_uni, readable=TRUE)
  ego_down <- enrichGO(gene = ent_gene_down,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       universe = ent_uni, readable=TRUE, pvalueCutoff=0.01)
  info <- as.data.frame(ego_all)[1:10,]
  print("Here 6")
  info$Description
  barplot(ego_down)
  dotplot(ego_down, showCategory=100)
  
  #install.packages("ggnewscale")
  
  fold_change = annotated_df_sub_up_regulated$log2FoldChange
  names(fold_change) = annotated_df_sub_up_regulated$external_gene_name
  
  cnetplot(ego_all, showCategory=5, foldChange = fold_change)
  goplot(ego_all)
  print("Here 7")
  #BiocManager::install("GOstats")
  #BiocManager::install("GOHyperGParams")
  #library(GO.db)
  #library(GOstats)
  #print("Here 8")
  # View(as.data.frame(annotated_df_sub))
  # order_upregulated = annotated_df_sub_up_regulated[order(annotated_df_sub_up_regulated$log2FoldChange, decreasing = TRUE), ]
  # selectedUpRegGene = order_upregulated$ensgene
  # print(selectedUpRegGene)
  # View(as.data.frame(order_upregulated))
  
  # Down-regulated
  
  # order_downregulated = annotated_df_sub_down_regulated[order(annotated_df_sub_down_regulated$log2FoldChange, decreasing = TRUE), ]
  # selectedDownRegGene = unique(order_downregulated$ensgene)
  # dim(as.data.frame(order_downregulated))
  
  # AllGene = annotated_df_$ensgene
  # summary(AllGene)
  # summary(selectedUpRegGene)
  # Ontology
  
  #UpParams = new('GOHyperGParams',geneIds=ent_gene_up,universeGeneIds=ent_uni,annotation="org.Hs.eg.db",ontology="BP",pvalueCutoff=0.01,conditional=FALSE,testDirection="over")
  
  
  #DownParams = new('GOHyperGParams',geneIds=ent_gene_up,universeGeneIds=ent_uni,annotation="org.Hs.eg.db",ontology="BP",pvalueCutoff=0.01,conditional=FALSE,testDirection="over")
  
  #UpGTestBP = hyperGTest(UpParams)
  #summary(UpGTestBP)[1:10,]
  #DownGTestBP = hyperGTest(DownParams)
  
  
  # if (CatToTest == "KEGG"){
  #   params <- new("KEGGHyperGParams", geneIds = selectedEntrezIds,
  #                 universeGeneIds = entrez.universe.vec,
  #                 annotation = annotationPckg,
  #                 pvalueCutoff = hgCutoff, testDirection = testDirection)
  # } else {
  #   params <- new("GOHyperGParams", geneIds = selectedEntrezIds,
  #                 universeGeneIds = entrez.universe.vec,
  #                 annotation =  annotationPckg,
  #                 ontology = CatToTest, pvalueCutoff = hgCutoff,
  #                 conditional = TRUE, testDirection = testDirection)
  # }
  
  # Gene Ontology Analysis
  
  #ekg = enrichKEGG(gene = ent_gene_up,universe = ent_uni)
  #ekg
  
  #write_tsv(annotated_df_enrichmentAnalysis_sub, paste(output_dir,"/annotated_df_enrichmentAnalysis_sub.csv",sep=""))
  
  dev.off()
  
}


RAnalysis(
  file=opt$file, 
  path_to_quant=opt$path_to_quant, 
  gene_map=opt$gene_map, 
  count_type=opt$count_type, 
  ignoreTxVersion=opt$ignoreTxVersion, 
  opt$output_dir,
  slice_group=opt$slice_group,
  condition=opt$condition, sample_count=opt$sample_count)

