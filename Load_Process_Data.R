
GSM5831744_adult_perinate_mtx <- ReadMtx(
  mtx = "Data/Raw_Input/GSM5831744_adult_perinate_gex_matrix.mtx.gz", 
  features = "Data/Raw_Input/GSM5831744_adult_perinate_gex_features.tsv.gz",
  cells = "Data/Raw_Input/GSM5831744_adult_perinate_gex_barcodes.tsv.gz"
)

# check metadata
GSM5831744_adult_perinate_meta <- read.table(
  "Data/Raw_Input/GSM5831744_adult_perinate_metadata.txt.gz", 
  header = TRUE)

# check that rowname order of meta matches colname order of matrix
match(GSM5831744_adult_perinate_meta$barcode, 
      colnames(GSM5831744_adult_perinate_mtx))

# append barcode values to rownames of meta
rownames(GSM5831744_adult_perinate_meta) <- 
  GSM5831744_adult_perinate_meta$barcode

so.orig.nf <- list(CreateSeuratObject(
  counts = GSM5831744_adult_perinate_mtx,
  meta.data = GSM5831744_adult_perinate_meta,
  project = "adult_perinate"))

# Filter and QC
mincells = 3
mingenes = 200
organism = "Mouse"

mitoch = "^mt-"
cc.genes$g2m.genes= str_to_title(cc.genes$g2m.genes)
cc.genes$s.genes = str_to_title(cc.genes$s.genes)

seurat_object <- function(i) {
  
  so.nf <- so.orig.nf[[i]]
  so.nf <- NormalizeData(so.nf, normalization.method = "LogNormalize", scale.factor = 10000)
  so.nf[["percent.mt"]] <- PercentageFeatureSet(object = so.nf, pattern = mitoch)
  so.nf$log10GenesPerUMI <- log10(so.nf$nFeature_RNA) / log10(so.nf$nCount_RNA)
  
  so <- so.nf
  
  so.origcount = dim(so.nf)[2]

  #Start with filtering here:
  maxgenes = 2500
  complexity = 0.8
  minUMI = 500
  MAD_gene <- TRUE
  ngenestdev <- mad(so@meta.data$nFeature_RNA)
  ngenemed <- median(so@meta.data$nFeature_RNA)
  ngenemaxlim <- ngenemed+(3*ngenestdev)
  gl = format(round(ngenemaxlim,0),nsmall=0)
  
  maxmitoch = 10
  
  MAD_mitoch <- TRUE
  mitostdev <- mad(so@meta.data$percent.mt)
  mitomed <- median(so@meta.data$percent.mt)
  mitomaxlim <- mitomed+(3*mitostdev)
  ml = format(round(mitomaxlim,2),nsmall=2)
  
  so <- subset(so, cells = rownames(
    so@meta.data[which(so@meta.data$nFeature_RNA < 
                         ngenemaxlim & so@meta.data$percent.mt <= 
                         mitomaxlim & so@meta.data$log10GenesPerUMI > 
                         complexity), ]))
  
  so2.list <- list(so,so.nf)
  
  return(so2.list)
}

so.list <- lapply(seq_along(so.orig.nf), seurat_object)

SO <- lapply(so.list,function(x) x[[1]])
names(SO) <- unlist(lapply(so.list, function(x) 
  as.character(Seurat::Idents(x[[1]])[1])))

rm(so.list)

gc(full = TRUE)

# PCA & Normalization
conserve_memory = TRUE
vars_to_regress <- c()
npcs = 40

# Linearly scale data without regressing anything.
scale_so <- function(so){
  so <- CellCycleScoring(object = so, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)
  so$CC.Difference <- so$S.Score - so$G2M.Score
  so <- FindVariableFeatures(object = so, nfeatures = 2000, mean.cutoff = c(1, 8), dispersion.cutoff = c(1, 100000), selection.method = "vst")
  all.genes <- rownames(so)
  so <- ScaleData(so,features=all.genes)
  return(so)
}

# Make PCA without regressing anything, and using only SCTransform().
pca_noregress <- function(so, diet = FALSE) {
  so <- SCTransform(so,do.correct.umi = FALSE,return.only.var.genes = FALSE, conserve.memory = conserve_memory)
  so <- RunPCA(object = so, features = VariableFeatures(object = so), npcs = npcs)
  
  if (diet == TRUE){
    so <- DietSeurat(so, counts = FALSE, data = TRUE, assays = "SCT", dimreducs="pca")
  }
  
  return(so)
}

# Make PCA with SCTransform() and optional ScaleData, and do so with
# both regression (if user requests) and on all genes.
pca <- function(so) {
  # If user sets Linear Scaling toggle TRUE, also run ScaleData().
  # Use case: user has legacy project from Seurat 2 and wants to keep
  # methods consistent with pre-SCT Seurat.
  
  # Run SCTransform().
  if(is.null(vars_to_regress)){
    so <- so
  }
  else { 
    so <- SCTransform(so,do.correct.umi = TRUE, vars.to.regress = vars_to_regress, return.only.var.genes = FALSE, conserve.memory = conserve_memory)
  }
  # Make PCA using last transform run, which will always be that from
  # SCTransform().
  so <- RunPCA(object = so, npcs = npcs)
  slot(so,"commands") <- list()
  return(so)
}

# Do transformation with and without regression using SCTransform()
# and ScaleData().
so_scale <- lapply(SO, scale_so) 
SO <- lapply(so_scale, pca) 

# Combine and Renormalize

## Add commentary on how this toggle works with add.only.var.genes.
conserve_memory <- FALSE

dat = vector()

SO_merge <- SO[[1]]
allgenes <- rownames(SO_merge)

npcs = 19
Do_SCTransform = TRUE
vars_to_regress = c()

SO_merge <- SCTransform(SO_merge,do.correct.umi = TRUE, 
                        conserve.memory = conserve_memory, 
                        return.only.var.genes = FALSE)

SO_merge <- FindVariableFeatures(object = SO_merge, nfeatures = 2000, 
                                 mean.cutoff = c(0.1, 8), 
                                 dispersion.cutoff = c(1, 100000), 
                                 selection.method = "vst", verbose = FALSE)
SO_merge <- RunPCA(object = SO_merge, npcs = npcs, verbose = FALSE,
                   seed.use = 42)
SO_merge <- RunUMAP(object = SO_merge, reduction = "pca", 
                    dims = 1:npcs, seed.use=42)
SO_merge <- RunTSNE(object = SO_merge, reduction = "pca", 
                    dim.embed = 2, dims = 1:npcs, seed.use = 1)
SO_merge <- FindNeighbors(SO_merge, dims = 1:npcs)

for (i in seq(0.5,1,0.1)){
  SO_merge <- FindClusters(SO_merge, resolution = i, algorithm = 1)
}
