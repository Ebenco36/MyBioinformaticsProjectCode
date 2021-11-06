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
library(tibble)
library(readr)

sample_table = read_csv("Desktop/R code/SraRunTableIschemicStroke.txt")
sample_table = select(sample_table, 
                      `Sample Name`, BioSample,`Assay Type`,Experiment,
                      Platform, ReleaseDate, Condition)
sample_table = unique(sample_table)

sample_table_sub <- rbind(sample_table[1:5,], sample_table[11:15,])


sample_files_sub = paste0(paste0(paste0('Desktop/R code/quants/', 
                                        pull(sample_table_sub, `Sample Name`)), '_quant'), 
                          '/quant.sf')

names(sample_files_sub) = pull(sample_table_sub, `Sample Name`)


gene_map = read_csv("Desktop/R code/ischemoc_stroke_rnaseq/gene_map.csv", 
                    col_names = c('enstid', 'ensgid'))


count_data_sub = tximport(files = sample_files_sub,
                          type = "salmon",
                          tx2gene = gene_map,
                          ignoreTxVersion = TRUE)

# mod
summary(count_data_sub)
print(sample_table_sub)


sample_table_sub = as.data.frame(sample_table_sub)


colnames(sample_table_sub)[1] = "Sample"


print(sample_table_sub)

# mod
conditions_sub = c('SIS', 'HC')


conditions_sub = rep(conditions_sub, each=5)

conditions_sub

conditions_sub = factor(conditions_sub, levels=c('SIS', 'HC'))

conditions_sub


sample_table_sub$conditions = conditions_sub


# y ~ x


deseq_dataset_sub = DESeqDataSetFromTximport(txi=count_data_sub, 
                                             colData=sample_table_sub,
                                             design=~conditions)

deseq_dataset_sub

deseq_dataset_sub = estimateSizeFactors(deseq_dataset_sub)

#gives us the normalization factors
normalizationFactors(deseq_dataset_sub)


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










boxplot(counts(deseq_dataset_sub, normalized=TRUE))


vst_sub = varianceStabilizingTransformation(deseq_dataset_sub)

boxplot(assay(vst_sub))


# observe the cell line effect

plotPCA(vst_sub, intgroup='conditions') +
  theme_bw()



d_sub = assay(vst_sub)

d_sub = t(d_sub)


d_sub = dist(d_sub)


h_sub = hclust(d_sub)


plot(h_sub)



#Using kmeans for clustering

k_sub = kmeans(t(assay(vst_sub)), centers=2)



k_sub$cluster


#write to file for futher analysis
deseq_dataset_sub_normalized  <- counts(deseq_dataset_sub, normalized=TRUE)

write.csv(deseq_dataset_sub_normalized, 'normalizedCountsSub.csv')


write.csv(as.matrix(counts(deseq_dataset_sub, normalized=FALSE)), 'UnnormalizedCountsSub.csv')


colSums((as.matrix(counts(deseq_dataset_sub, normalized=TRUE))))

colSums((as.matrix(counts(deseq_dataset_sub, normalized=FALSE))))
#Analysis Starts from here. Previous information is used to 
#check out the relationship between data for trust in futher usage.



# 3 steps to DESeq2 analysis
# 1) estimate size factors (normalisation)
# 2) estimate dispersions
# 3) apply statistics (Wald Test)


dispersion_sub = estimateDispersions(deseq_dataset_sub)

# View(dispersion_sub)

plotDispEsts(dispersion_sub)


# step 3
statistics_wald_test_sub = nbinomWaldTest(dispersion_sub)
colData(statistics_wald_test_sub)
result_table_sub = results(statistics_wald_test_sub)
summary(result_table_sub)
# View(as.data.frame(result_table_sub))

#----------------------------------------------------------------------------------Continue
# result_table is a DataFrame not a data.frame!
result_df_sub = as.data.frame(result_table_sub)


#P-value
write.csv(result_df_sub, 'TtestValuesCountsSub.csv')



plotCounts(statistics_wald_test_sub, gene='ENSG00000060709',
           intgroup='conditions')
sum(!complete.cases(result_df_sub))

filter_df_sub = result_df_sub[complete.cases(result_df_sub),]


result_df[complete.cases(result_df_sub),]


#Completed cases
write.csv(filter_df, 'TtestValuesCountsOptimizedCompleteCasesSub.csv')





deseq_dataset_normalized_filtered_sub <- as.data.frame(as.matrix(counts(deseq_dataset_sub, normalized=TRUE)))
deseq_dataset_normalized_filtered_sub <- cbind(deseq_dataset_normalized_filtered_sub, genes=c(rownames(deseq_dataset_normalized_filtered_sub)))
format_filter_df1_sub <- cbind(filter_df_sub, genes=c(rownames(filter_df_sub)))
deseq_dataset_normalized_filtered_joined_sub = left_join(format_filter_df1_sub, deseq_dataset_normalized_filtered_sub, by=c('genes'='genes'))
deseq_dataset_normalized_filtered_joined_sub$test = deseq_dataset_normalized_filtered_joined_sub$padj < 0.05 & abs(deseq_dataset_normalized_filtered_joined_sub$log2FoldChange) > 1


write.csv(deseq_dataset_normalized_filtered_joined_sub, 'filtered_normalizedCountsSubAcute.csv')
# View(deseq_dataset_normalized_filtered_joined_sub)


#uncompleted Cases
write.csv(result_df_sub[!complete.cases(result_df_sub), ], 'TtestValuesCountsOptimizedUnCompleteCasesSub.csv')

#View(filter_df1)

# Filter results 


filter_df_sub = filter_df_sub[filter_df_sub$padj < 0.05,]
write.csv(filter_df_sub, 'FilterpadJLessThan005Sub.csv')


filter_df_sub = filter_df_sub[abs(filter_df_sub$log2FoldChange) > 1,]


write.csv(filter_df_sub, 'FilterGreaterThan1withpadjlessthen005Sub.csv')


plotMA(result_table)

# volcano plot
# Combiing the 2 conditions together
filter_df_sub$test = filter_df_sub$padj < 0.05 & abs(filter_df_sub$log2FoldChange) > 1
dim(filter_df_sub)
filter_df_sub
#filter_df1


filter_df_sub = rownames_to_column(filter_df_sub, var='ensgene')

filter_df = rownames_to_column(filter_df, var='ensgene')
filter_df$test = filter_df$padj < 0.05 & abs(filter_df$log2FoldChange) > 1
abs(filter_df$log2FoldChange) 


write.csv(filter_df_sub, 'rownametocolumnSub.csv')



g_sub = ggplot(filter_df, aes(x=log2FoldChange, 
                                   y=-log10(padj), 
                                   name=ensgene)) +
  geom_point(aes(colour=test), size=1, alpha=0.3) +
  scale_colour_manual(values=c('black', 'red')) +
  geom_vline(xintercept=1, colour='green', linetype=3) +
  geom_vline(xintercept=-1, colour='green', linetype=3) +
  geom_hline(yintercept=-log10(0.05), colour='blue', linetype=3) +
  theme_bw() +
  theme(legend.position = 'none')
print(g_sub + ggtitle("Differential Expressed Genes"))

ordered_data = filter_df[order(abs(filter_df$log2FoldChange), decreasing=TRUE),]
# View(as.data.frame(ordered_data))
#ggplotly(g_sub)

#install BioMart
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("biomaRt")

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
      values = filter_df1$ensgene[1:6],
      mart = ensembl104)

dim(filter_df_sub)
dim(filter_df_sub[filter_df_sub$log2FoldChange < 1, ])
dim(filter_df_sub[filter_df_sub$log2FoldChange > 1, ])
all_regulated = filter_df_sub
up_regulated = filter_df_sub[filter_df_sub$log2FoldChange < 1, ]
down_regulated = filter_df_sub[filter_df_sub$log2FoldChange > 1, ]

View(filter_df)

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



write.csv(annotation_sub, 'BioMartInfoSub.csv')


nrow(annotation_sub)



# View(filter_df_sub)
# View(filter_df)
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

# View(annotated_df_sub_down_regulated)

write.csv(annotated_df_sub[annotated_df_sub$test == TRUE,], 'BioMartInfoJoinSub.csv')




# View(annotated_df_sub[annotated_df_sub$test == TRUE,])


g_df_sub = ggplot(annotated_df_, aes(x=log2FoldChange, 
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
ggplotly(g_df_sub)





#vst_nhbe = varianceStabilizingTransformation(vst)
vst_mat_sub = assay(vst_sub)


# View(vst_mat_sub)
#Just for the differential expressed genes

data_for_hm_sub = vst_mat_sub[degs_sub,]

data_for_hm_sub
#data for Heatmap Stored
write.csv(data_for_hm_sub, 'HeatMapDataForDEGSUB.csv')
#change the Gene ID to their respective names

annotated_df_sub$external_gene_name
rownames(data_for_hm_sub) = annotated_df_sub$external_gene_name
write.csv(data_for_hm_sub, 'HeatMapDataForDEGWithNames.csv')
data_for_hm_sub
dim(data_for_hm_sub)
heatmap(data_for_hm_sub)

#get Better Label
pheatmap(data_for_hm_sub, fontsize_row=4, scale='row')

#To use color brewer we need to install it
#BioManager::install('RColorBrewer')
#see all color pelette available
#display.brewer.all()

greys = colorRampPalette(brewer.pal(9, "Greys"))(100)

pheatmap(data_for_hm_sub, fontsize_row=4, scale='row',
         color=greys)

#pairs = colorRampPalette(brewer.pal(12, "Paired"))(100)

last_scheme = colorRampPalette(brewer.pal(7, "Blues"))(100)

pheatmap(data_for_hm_sub, fontsize_row=4, scale='row',
         color=last_scheme, cutree_cols = 2,
         cutree_rows = 2)

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


ent_gene_up = getBM(attributes=c('entrezgene_id'),
                    filters = c('ensembl_gene_id'),
                    values = annotated_df_sub_up_regulated$ensgene,
                    mart = ensembl104)
# annotated_df_sub_up_regulated
ent_gene_down = getBM(attributes=c('entrezgene_id'),
                    filters = c('ensembl_gene_id'),
                    values = annotated_df_sub_down_regulated$ensgene,
                    mart = ensembl104)

(ent_gene_down$entrezgene_id)

ent_gene = ent_gene$entrezgene_id
ent_gene_all = ent_gene_all$entrezgene_id
ent_gene_up = ent_gene_up$entrezgene_id
ent_gene_down = ent_gene_down$entrezgene_id

ent_gene = as.character(ent_gene)
ent_gene_all = as.character(ent_gene_all)
ent_gene_up = as.character(ent_gene_up)
ent_gene_down = as.character(ent_gene_down)

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
info$Description
barplot(ego_down)
dotplot(ego_down, showCategory=100)

install.packages("ggnewscale")
fold_change = annotated_df_sub_up_regulated$log2FoldChange
names(fold_change) = annotated_df_sub_up_regulated$external_gene_name

cnetplot(ego, showCategory=5, foldChange = fold_change)
goplot(ego)

BiocManager::install("GOstats")
BiocManager::install("GOHyperGParams")
library(GO.db)
library(GOstats)
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

UpParams = new('GOHyperGParams',
               geneIds=ent_gene_up,
               universeGeneIds=ent_uni,
               annotation="org.Hs.eg.db",
               ontology="BP",
               pvalueCutoff=0.01,
               conditional=FALSE,
               testDirection="over")


DownParams = new('GOHyperGParams',
               geneIds=ent_gene_up,
               universeGeneIds=ent_uni,
               annotation="org.Hs.eg.db",
               ontology="BP",
               pvalueCutoff=0.01,
               conditional=FALSE,
               testDirection="over")

UpGTestBP = hyperGTest(UpParams)
summary(UpGTestBP)[1:10,]
DownGTestBP = hyperGTest(DownParams)


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

ekg = enrichKEGG(gene = ent_gene_up,
                 universe = ent_uni)
ekg
write_tsv(annotated_df_enrichmentAnalysis_sub, "annotated_df_enrichmentAnalysis_sub.csv")

