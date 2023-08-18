# set directory
setwd("rstudio-files/ccbr-data/users/Jing/Taves_mTEC/ccbr1234_manuscript")

# process raw files and create seurat object 
source("Load_Process_Data.R")

# remove everything except SO_merge
rm(list=setdiff(ls(), "SO_merge"))

# plot clusters umaps
source("Functions/Plot_Metadata.R")

# plot gene expression colored umaps
source("Functions/Color_by_Gene.R")

# Figure 2A, 4A, Supplementary 4A
plotMetadata(object = SO_merge,
            metadata.to.plot = "SCT_snn_res_0_7",
            columns.to.summarize = "c()",
            summarization.cut.off = 5,
            reduction.type = "umap",
            use.cite.seq = FALSE,
            show.labels = FALSE,
            legend.text.size = 1,
            legend.position = "right",
            dot.size = 0.01)

colorByGene(object = SO_merge,
            gene = c("Cyp11b1","Aire","Hsd17b3","Cyp19a1","Star",
                     "Cyp11a1","Hsd3b1","Hsd3b6"),
            reduction.type = "umap",
            number.of.rows = 2,
            assay = "SCT")

# Supplementary 4B
source("Functions/Heatmap.R")

# cluster 4 vs all DEG - MAST
source("Functions/MAST.R")

clust4_enriched <- MAST(SO.sub = SO_merge,
                        metadata_table = SO_merge@meta.data, 
                        parameter_to_test = "SCT_snn_res.0.7",
                        contrasts = c("4-all"),
                        test_to_use = "MAST",
                        log_fc_threshold = 0,
                        assay_to_use = "SCT",
                        use_log_2 = TRUE,
                        latent_vars = c())

