library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
# if (length(args)!=2) {
#   stop("please provide the path to 1) cellphoneDB output folder, 2) differentially expressed genes folder", call.=FALSE)
# }else{
#   message('cellphoneDB output folder: ', args[1])
#   message('DEGs folder: ', args[2])
# }

CPdb_folder = '~/cellphoneDB/analysis/CVID/out/' #args[1]
DEG_folder = '~/cellphoneDB/analysis/CVID/DEG_MAST_20200131/' #args[2]


# Adapt deconvoluted
my_deconvoluted = read.delim(paste0(CPdb_folder, 'deconvoluted.txt'), stringsAsFactors = F)
my_deconvoluted[, grep('celltype', colnames(my_deconvoluted)) ] = 10
colnames(my_deconvoluted) = gsub('celltype_', '', colnames(my_deconvoluted))
colnames(my_deconvoluted) = gsub('\\._', '_', colnames(my_deconvoluted))
colnames(my_deconvoluted) = gsub('\\.', '-', colnames(my_deconvoluted))


# Load DEG
DEGs_f = list.files(DEG_folder, full.names = T)
DEGs = lapply(DEGs_f, read.csv, stringsAsFactors=F)
# Filter significant DEGs
# expressed_G = lapply(DEGs, subset, percentExpr_cluster+percentExpr_rest > 0.2) %>%
  # lapply(., subset, Gene %in% my_deconvoluted$gene_name )
expressed_G = lapply(DEGs, subset, pct.1+pct.2 > 0.2) %>%
  lapply(., subset, Gene %in% my_deconvoluted$gene_name )
names(expressed_G) = sapply(strsplit(DEGs_f, '/'), tail, 1) %>% strsplit(., '_CVID_vs_') %>% sapply(., head, 1) %>% gsub('\\+', '', .)
DEGs = lapply(expressed_G, subset, adj.P.Val < 0.01)
names(DEGs) = names(expressed_G) 
DEGs = lapply(DEGs, subset, abs(logFC) >= 0.1)



# Fill deconvoluted with DEGs p-values
nrow(my_deconvoluted)
# build genes2pvalue dictionary
get_DEG_pval = function(gene)
  sapply(DEGs, function(x)  x$adj.P.Val[ x$Gene == gene] * sign(x$logFC[ x$Gene == gene]) ) %>% unlist(.)
get_DEG_foldchange = function(gene)
  sapply(DEGs, function(x)  x$logFC[ x$Gene == gene] ) %>% unlist(.)
genes2logFold = lapply(unique(my_deconvoluted$gene_name), get_DEG_foldchange)
names(genes2logFold) = unique(my_deconvoluted$gene_name)
# build genes2pvalue dictionary
# get_percent = function(gene)
#   sapply(expressed_G, function(x)  x$percentExpr_cluster[ x$Gene == gene] ) %>% unlist(.)
get_percent = function(gene)
   sapply(expressed_G, function(x)  x$pct.1[ x$Gene == gene] ) %>% unlist(.)
genes2percent = lapply(unique(my_deconvoluted$gene_name), get_percent)
names(genes2percent) = unique(my_deconvoluted$gene_name)
# Filter interactions with no DEGs
genes_in_DEGs = sapply(DEGs, function(x) x$Gene) %>% unlist(.) %>% unique(.)
my_deconvoluted = subset(my_deconvoluted, gene_name %in% genes_in_DEGs)
# Substitute value by the adj.P.Val signed according to the fold change
genes2logFold = genes2logFold[ sapply(genes2logFold, length) > 0 ]
for(gene in names(genes2logFold)){
  for ( celltype in names(genes2logFold[[gene]]) )
    my_deconvoluted[ my_deconvoluted$gene_name == gene, celltype ] = genes2logFold[[gene]][celltype]
}
# Remove genes not in the L/R collection
rows2remove = apply(my_deconvoluted[, 7:ncol(my_deconvoluted)], 1, min) != 10
my_deconvoluted = my_deconvoluted[ rows2remove, ]
# Remove celltypes with no DEGs in the L/R collection
celltype2remove = names(which(apply(my_deconvoluted[, 7:ncol(my_deconvoluted)], 2, min) == 10))
my_deconvoluted = my_deconvoluted[ , !(names(my_deconvoluted) %in% celltype2remove)]
nrow(my_deconvoluted)



# Adapt means matrix
means_file = read.delim(paste0(CPdb_folder, 'means.txt'), stringsAsFactors = F)
# Remove non-curated interactions
means_file = subset(means_file, annotation_strategy == "user_curated")
means_file = means_file[ ! duplicated(means_file$id_cp_interaction), ]
# Add genes in complexes
complexes = read.csv('~/farm/CellPhoneDB-data_smallmolecules/data/sources/complex_curated.csv', stringsAsFactors = F)
complexes$complex_name = paste0('complex:', complexes$complex_name)
genes = read.csv('~/farm/CellPhoneDB-data_smallmolecules/data/gene_input_all.csv', stringsAsFactors = F)
complexes2genes = lapply(complexes$complex_name, function(cx) subset(genes, uniprot %in% complexes[complexes$complex_name == cx, 2:5] )$gene_name )
complexes2genes = lapply(complexes2genes, unique)
names(complexes2genes) = complexes$complex_name
# Build means matrix de novo
my_means = unique(means_file[1:11])
my_means$gene_a[ my_means$partner_a %in% names(complexes2genes)] = sapply(complexes2genes[my_means$partner_a[my_means$partner_a %in% names(complexes2genes)]], paste, collapse=';')
my_means$gene_b[ my_means$partner_b %in% names(complexes2genes)] = sapply(complexes2genes[my_means$partner_b[my_means$partner_b %in% names(complexes2genes)]], paste, collapse=';')
# Add reverse partnerA -> B and vice versa
my_means_reverse = my_means
my_means_reverse$id_cp_interaction = paste0(my_means$id_cp_interaction, '_rev')
my_means_reverse$gene_a = my_means$gene_b
my_means_reverse$partner_a = my_means$partner_b
my_means_reverse$gene_b = my_means$gene_a
my_means_reverse$partner_b = my_means$partner_a
my_means_reverse$interacting_pair = paste(my_means_reverse$gene_a, my_means_reverse$gene_b, sep='_')
my_means = rbind(my_means, my_means_reverse)
my_means = my_means[ ! duplicated(my_means$interacting_pair) , ]
# We define relevant interactions as those where the partnerB have expression > 10% and any partnerA member is DEG
int_of_interest = function(int, ctA, ctB){
  partnersA = strsplit(int[1], ';') %>% unlist(.)
  partnersB = strsplit(int[2], ';') %>% unlist(.)
  A = all(ctA %in% sapply(genes2percent[partnersA], names))
  B = all(ctB %in% sapply(genes2percent[partnersB], names))
  Adeg = any(ctA %in% sapply(genes2logFold[partnersA], names))
  if(B & Adeg){
    max_fold = sapply(genes2logFold[partnersA], function(x) x[ctA] ) %>% unlist(.)
    max_fold = max_fold[ which.max(abs(max_fold)) ]
    return(max_fold)
  }else{
    return(10)
  }
}
# For each pair of interacting cell types, chek if interaction is relevant because a partner is DE and retrieve forl change
for (ctA in names(DEGs) )
  for (ctB in names(DEGs) ){
    if( ctA == ctB | length(grep('B', c(ctA, ctB))) == 0) # if any B cell there
      next()
    foldchangeA = apply(my_means[,5:6], 1, int_of_interest, ctA, ctB)
    if( all(foldchangeA == 10)  )
      next()
    df = data.frame(foldchangeA)
    names(df) = paste0(ctA, '.DEGs---', ctB)
    my_means = cbind(my_means, df)
  }
# Remove interactions that are not relevant
idx = which(apply(my_means[, 12:ncol(my_means)], 1, sum) != (10*ncol(my_means)-11) )
my_means = my_means[idx, ]
# Fix L/R names
genes_a = my_means$gene_a
genes_a[ grep('complex', my_means$partner_a) ] = grep('complex', my_means$partner_a, value = T) %>% gsub('complex:', '', .)
genes_b = my_means$gene_b
genes_b[ grep('complex', my_means$partner_b) ] = grep('complex', my_means$partner_b, value = T) %>% gsub('complex:', '', .)
rownames(my_means) = paste(genes_a, genes_b, sep = '---')
# Plot the results - as retrieved
results = as.matrix(my_means[, 12:ncol(my_means)])
results[ results == 10 ] = 0
results = results[ rowSums(results) != 0 , ]
library("RColorBrewer")
library("gplots")
col <- colorRampPalette(brewer.pal(9, "RdBu"))(256)
par(mar=c(1,1,1,1)) 
pdf('~/cellphoneDB/analysis/CVID/cellphoneDB_DEGs_significant_FDR01_heatmap_alternative.pdf', width = 22, height = 22)
heatmap.2(t(results), scale = "none", col = bluered(100), Rowv = NA, Colv = NA, 
          trace = "none", density.info = "none", 
          sepwidth=c(0.01,0.01),
          sepcolor="black",
          colsep=0:ncol(t(results)),
          rowsep=0:nrow(t(results)),
          keysize = 0.5,
          key=TRUE, symkey=FALSE, cexRow=1,cexCol=1,margins=c(12,25),srtCol=45)
graphics.off()
# Plot the results - alternative format
library(reshape2)
results = melt(as.matrix(my_means[, 12:ncol(my_means)]), factorsAsStrings = F)
results$Var1 = as.character(results$Var1)
results$Var2 = as.character(results$Var2)
results = subset(results, value != 10)
results$partnerA_DE = strsplit(results$Var1, split = '---') %>% sapply(., head, 1) 
results$partnerB = strsplit(results$Var1, split = '---') %>% sapply(., tail, 1)
results$celltypeA_DE_in_CVID = strsplit(results$Var2, split = '---') %>% sapply(., head, 1) %>% gsub('\\.DEGs', '', .)
results$celltypeB = strsplit(results$Var2, split = '---') %>% sapply(., tail, 1) 
results$logFC = results$value
results$Y = paste(results$partnerA_DE, results$celltypeA_DE_in_CVID, sep = ' --- ')
results$X = paste(results$partnerB, results$celltypeB, sep = ' --- ')
head(results)

library(ggplot2)
ggplot(results, aes(x = X, y = Y)) +
  geom_tile(aes(fill = logFC), colour = "black") +
  xlab('interacting partner / cell type') + ylab('genes differentially expressed in CVID') +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu"))) +
  ggtitle("Cell-cell communication events differentially expressed in CVID") +
  theme(#panel.background = element_blank(),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(filename = '~/cellphoneDB/analysis/CVID/cellphoneDB_DEGs_significant_FDR01_heatmap.pdf', dpi = 300, width = 30, height = 12)


# Add partner expression
partner_expression = melt(unlist(genes2percent))
results$partnerB_percentExpr = partner_expression[ gsub(' --- ', '.', results$X), ]
write.csv(results[, -c(1:3)], file = '~/cellphoneDB/analysis/CVID/cellphoneDB_DEGs_significant_FDR01.csv', quote = F, row.names = F)


