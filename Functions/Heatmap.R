library(Seurat)
library(pheatmap)

Idents(SO_merge) <- SO_merge@meta.data$SCT_snn_res.0.7

# Genes to plot on heatmap - if not using MS
genes <- c("Star", "Cyp11a1", "Hsd3b1", "Hsd3b6", "Cyp11b1", "Aire", "H2-Aa",
           "Pou2f3", "Chat", "Trpm5", "Sox9", "Ccl9", "Tnfaip2", "Apoa4", 
           "Guca2a", "Foxj1", "Dynlrb2", "Dnah12", "Myog", "Actc1", "Ckm", 
           "Foxi1", "Cftr", "Atp6v1g3", "Ptf1a", "Prss2", "Clps", "Foxa3", 
           "Chgb", "Ret", "Aqp4", "Muc5b", "Spink5", "Rptn", "Gata3", "Krt10", 
           "Krtdap", "Grhl1","Ivl", "Sbsn")

# Use only genes that found in the data
genes_not_found <- genes[!genes %in% rownames(SO_merge$SCT@scale.data)]

agg_mtx <- AggregateExpression(SO_merge, assays = "SCT", 
                               group.by = "ident", slot = "scale.data")
agg_mtx <- agg_mtx[[1]]

plot_mtx <- agg_mtx[genes,]

print(pheatmap(plot_mtx, cluster_rows = FALSE, 
               cluster_cols = FALSE, scale = "row"))
