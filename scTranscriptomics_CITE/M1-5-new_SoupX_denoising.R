library(SoupX)
library(Seurat)

# 12.02.2021: reanalysis without twin data and without BCR activation

samples = list('35008_CV005_RV9039258_and_35008_CV005_RV9039282',
               '35008_CV005_RV9039259_and_35008_CV005_RV9039283',
               '35008_CV005_RV9039260_and_35008_CV005_RV9039284',
               '35008_CV005_RV9039261_and_35008_CV005_RV9039285',
               '35008_CV005_RV9039262_and_35008_CV005_RV9039286',
               '35008_CV005_RV9039263_and_35008_CV005_RV9039287',
               #'35008_CV005_RV9039264_and_35008_CV005_RV9039288',
               #'35008_CV005_RV9039265_and_35008_CV005_RV9039289',
               '35171_CV005_RV9039361',
               '35171_CV005_RV9039362',
               '35171_CV005_RV9039363'
               #'35171_CV005_RV9039364'
               )

for (sample in samples)
{
  print(sample)
  #path_raw = paste('/warehouse/team292_wh01/aa22/data/202009_CVID_revision/cellranger302_count_', sample,'_GRCh38-1_2_0/raw_feature_bc_matrix/',
  #                 sep = "")
  
  path_raw = paste('/lustre/scratch117/cellgen/team292/aa22/data/202009_CVID_revision/cellranger302_count_', sample,'_GRCh38-1_2_0/raw_feature_bc_matrix/',
                   sep = "")
  
  print(path_raw)
  emptymat <- Read10X(path_raw)[["Antibody Capture"]]
  
  # annotation of this sample's cells - but what about all the empty droplets?
  annot_path = paste('/home/jovyan/notebooks/Vento_Lab/CVID/202009_new_analysis_revision/CITE_all_samples_analysis/CVID/scTranscriptomics_CITE/annot_tables_for_SoupX/20210212_prelim_annot_all_samples_annot_sample_', sample, '_validation_cohort.csv',
                     sep = "")
  annot = read.csv(annot_path,
                   row.names = 1)
  # cluster annotation from per-sample scanpy analysis
  #clusters_path = paste('/lustre/scratch117/cellgen/team292/aa22/adata_objects/202009_CVID_revision/louvain_clusters_for_SoupX_sample_', sample, '.csv',
  #                      sep = "")
  
  # cluster annotation from per-sample scanpy analysis
  clusters_path = paste('/lustre/scratch117/cellgen/team292/aa22/adata_objects/202009_CVID_revision/louvain_clusters_for_SoupX_sample_', sample, '_new_20210211.csv',
                        sep = "")
  clusters_table = read.csv(clusters_path,
                            row.names = 1)
  mDat <- annot
  mDat$clusters = clusters_table$louvain
  
  # all the antibody channels
  genes_to_regress <- c(rownames(emptymat))
  
  cells_filter <- rownames(mDat)
  
  # cite-seq antibody counts coming from acctual cells
  Ab_counts_from_cells = emptymat[genes_to_regress, cells_filter]
  
  sc = SoupChannel(emptymat, toc = Ab_counts_from_cells,
                   # annotation table
                   metaData = mDat, calcSoupProfile=F)
  
  sc = estimateSoup(sc, 
                    # number of UMIs expected in empty droplets
                    soupRange=c(4,100) )
  
  sc = autoEstCont(sc,soupQuantile= 0.1, # was 0.15 originally 
                   tfidfMin = 0.05, # was originally 0.2
                   forceAccept=TRUE )
  
  out = adjustCounts(sc)
  
  #save_path = paste('/lustre/scratch117/cellgen/team292/aa22/adata_objects/202009_CVID_revision/20201019_counts_denoised_with_SoupX_sample_', sample, '.csv',
  #                  sep = "")
  
  save_path = paste('/lustre/scratch117/cellgen/team292/aa22/adata_objects/202009_CVID_revision/20210212_counts_denoised_with_SoupX_sample_', sample, '.csv',
                    sep = "")
  
  write.csv(out, file = save_path)
  
}




















# BELOW IS THE non-used guideline code from Louis!

#cellnames = c(cellnames, rownames(sro@meta.data)[flt])

#if (i == 1) output <- out
#else output = Matrix::cbind2(output,out)



runSoupX <- function(sro, celltype, batch, batch_folder_index= c(), cr_folders=c(), genes_to_regress = c(), 
                     
                     soupQuantile = 0.15, # --> SoupX complains about these parameters, but that's because this is cite-seq
                     tfidfMin = 0.2, # defines which features will be considered as background
                     
                     filterclustername= "filtered", check_mapping_only=F){
  #library(SoupX); library(Seurat)
  if (is.null(cr_folders)) cr_folders = paste("/lustre/scratch117/cellgen/team292/lh20/I-O-", c("1_13-N", "1_13", "11_12-N", "11_12", "2_3-N", "2_3", "4-N", "4", "5_8-N", "5_8", "6_9-N", "6_9", "7-N", "7", "10-N", "10") ,"_C-0/outs/raw_feature_bc_matrix", sep="")
  if (is.null(batch_folder_index)) batch_folder_index = c(1,2,1,2,3,4,3,4,5,6,5,6,7,8,15,16,9,10,11,12,13,14,9,10,11,12) 
  #if (is.null(batch_folder_index)) batch_folder_index = c(1,1,2,2,15,16,3,3,4,4,5,5,6,6,7,8,9,9,10,10,11,11,12,12,13,14) 
  
  if (is.null(genes_to_regress)) genes_to_regress <- 33568:33759
  cellnames = c()
  for(i in 1:length(cr_folders)){
    curbatches <- levels(sro@meta.data[[batch]])[which(batch_folder_index == i)]
    print(paste("The following samples are associated to the output folder" , cr_folders[i]))
    print(curbatches)
    if (!check_mapping_only){
      emptymat <- Read10X(cr_folders[i])[["Antibody Capture"]]
      flt <- sro@meta.data[[batch]] %in% curbatches
      flt[as.character(sro@meta.data[[celltype]]) == filterclustername] <- F
      clust = paste( as.character(sro@meta.data[[batch]][flt]),  as.character(sro@meta.data[[celltype]][flt]))
      clust <- as.factor(clust)
      mDat = data.frame(clusters=clust@.Data, row.names=rownames(sro@meta.data)[flt])
      # cite-seq counts coming from acctual cells
      toc = sro@assays$RNA@counts[genes_to_regress, flt]
      # genes_to_regress is CITE-seq row names
      rownames(emptymat) <- rownames(toc)
      print(paste("Processing soup in ", sum(flt)," cells and ", length(levels(clust)) ," clusters", sep=""))
      sc = SoupChannel(emptymat, toc = toc, 
                       # annotation table
                       metaData = mDat, calcSoupProfile=F)
      # 
      sc = estimateSoup(sc, 
                        # number of UMIs expected in empty droplets
                        soupRange=c(4,100)
      )
      sc = autoEstCont(sc,soupQuantile= soupQuantile, tfidfMin = tfidfMin,forceAccept=TRUE )
      out = adjustCounts(sc)
      cellnames = c(cellnames, rownames(sro@meta.data)[flt])
      if (i == 1) output <- out
      else output = Matrix::cbind2(output,out)
    }}
  if (check_mapping_only) return;
  fout = matrix(0, nrow(output), nrow(sro@meta.data))
  map = match(rownames(sro@meta.data), cellnames)
  print(paste("reordering available data for", ncol(output), "genes and", sum(!is.na(map)), "cells out of", length(map)))
  fout[,!is.na(map)] = as.matrix(output[,map[!is.na(map)]])
  return(fout)}