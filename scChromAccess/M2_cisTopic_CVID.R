# cisTopic installation
#devtools::install_github("aertslab/cisTopic") # requires devtools package

# loading libraries
library(cisTopic)
library(umap)
library(ggplot2)
library(plotly) 
library(fastcluster) 
library(grid) 

source("https://bioconductor.org/biocLite.R")
biocLite(c('Rsubread', 'umap', 'Rtsne', 'ComplexHeatmap', 'fastcluster', 'data.table', 'rGREAT', 'ChIPseeker', 'TxDb.Hsapiens.UCSC.hg19.knownGene', 'org.Hs.eg.db'))


pathToBams <- './picard_bam_files/'
bamFiles <- paste(pathToBams, list.files(pathToBams), sep='')
regions <- './UCSC_style_peak_aggregated_scATAC_individual.narrowPeak'

# takes about 2 minutes
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles, regions, project.name='CVID')

cellData_cvid <- read.table(file='./annotation_CVID.csv',sep='\t',header = TRUE)
rownames(cellData_cvid) <- cellData_cvid$cell
cellData_cvid <- cellData_cvid[c(2:3)]

# Renaming cells to Immunodeficiency7827658_f2q30_pmd and so on
cell.names <- cisTopicObject@cell.names
strsplit(cell.names, split = ".", fixed=TRUE)[[1]][1]
new.cell.names <- sapply(strsplit(cell.names, split = ".", fixed=TRUE), "[", 1)
cisTopicObject <- renameCells(cisTopicObject, new.cell.names)
cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = cellData_cvid)

# running models
# for CVID here takes about 20 minutes - that's for about 14 topic points
# nCores should be equal to the number of topics (? not sure)
cisTopicObject <- runModels(cisTopicObject, topic=c(10,12,16,18,20,22,24,30,40,50,60), seed=987, nCores=1, burnin = 120, iterations = 150, addModels=FALSE)
# settling on 18 topics
cisTopicObject <- runModels(cisTopicObject, topic=c(18), seed=987, nCores=1, burnin = 120, iterations = 150, addModels=FALSE)

# choosing the best model 
par(mfrow=c(1,1))

cisTopicObject <- selectModel(cisTopicObject)

logLikelihoodByIter(cisTopicObject, select=c(18,18))


# Obtaining signatures of TFs and enhancers
path_to_signatures <- './Chip_seq_signatures_CVID/for_heatmap/'
ChIP_Seq_signatures <- paste(path_to_signatures, list.files(path_to_signatures), sep='')
ChIP_Seq_signatures
labels  <- c('ATF2_GM12878', 'ATF7_GM12878', 'BATF_GM12878', 'CTCF_Bcells', 'CTCF_GM12878', 'E2F4_GM12878', 'EBF1_GM12878','ETS1_GM12878',
             'H3K27ac_GM12878', 'IRF4_GM12878', 'NFATC1_GM12878', 'PAX5_GM12878', 'SP1_GM12878', 'Naive_active', 'Naive_poised', 'Naive_primed',
             'Naive_all_3_marks', 'GC_active', 'GC_poised', 'GC_primed', 'GC_all_3_marks', 'SM_active', 'SM_poised', 'SM_primed', 'SM_all_3_marks',
             'Plasma_active', 'Plasma_poised', 'Plasma_primed', 'Plasma_all_3_marks')
cisTopicObject <- getSignaturesRegions(cisTopicObject, ChIP_Seq_signatures, labels=labels, minOverlap = 0.4)

# getting the topic-region matrix
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)

# exporting the table of topic enrichment in single cells to customise the plotting
TCmatrix <- modelMatSelection(cisTopicObject, target='cell', method='Probability') # You can also extract the region-topic matrix with target='region'
TCmatrix <- t(TCmatrix)
TCmatrix <- cbind(TCmatrix, cisTopicObject@cell.data) #You can add any annotation you need here, this is for getting the format of your table
TCmatrix <- TCmatrix[c(1:18,25,26)]
write.table(TCmatrix, './cisTopic_cell_enrichment.csv', append = FALSE, sep = "\t",
            row.names = TRUE, col.names = TRUE)

# exporting topics as bed files for HOMER analysis later
getBedFiles(cisTopicObject, path='./cisTopics_asBed')


# exporting the topic-signature matrix
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
signatureMatrix <- topicSignaturesMatrix(cisTopicObject)
write.table(signatureMatrix, './cisTopic_signatures_enrichment', append = FALSE, sep = "\t",
            row.names = TRUE, col.names = TRUE)

