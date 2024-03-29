
plotMetadata <- function(
    object,
    metadata.to.plot,
    columns.to.summarize = "c()",
    summarization.cut.off = 5,
    reduction.type = "tsne",
    use.cite.seq = FALSE,
    show.labels = FALSE,
    legend.text.size = 1,
    legend.position = "right",
    dot.size = 0.01
) {

  # Helper Functions
  
  .drawMetadata <- function(m) {
    #check if there are NaNs in metadata, if there are, catch
    if (any(is.na(meta.df[[m]]))) {
      print(
        "ERROR: Metadata column appears to contain NA values.
         This is not recommended for clustering plots."
      )
      print("Please review your selected metadata column")
      print(head(meta.df[[m]]))
      print("Below are valid metadata to select for this plot:")
      print(valid.columns)
      stop("End of error message.")
    }
    #Making a plot based on tsne/umap/pca and CiteSeq settings:
    reduction = reduction.type
    
    if (!(use.cite.seq)) {
      if (reduction == "tsne") {
        p1 <- DimPlot(object,
                      reduction = "tsne",
                      group.by = "ident")
      } else if (reduction == "umap") {
        p1 <- DimPlot(object,
                      reduction = "umap",
                      group.by = "ident")
      } else {
        p1 <- DimPlot(object,
                      reduction = "pca",
                      group.by = "ident")
      }
    } else {
      if (reduction == "tsne") {
        p1 <-
          DimPlot(object,
                  reduction = "protein_tsne",
                  group.by = "ident")
      } else if (reduction == "umap") {
        p1 <-
          DimPlot(object,
                  reduction = "protein_umap",
                  group.by = "ident")
      } else {
        p1 <-
          DimPlot(object,
                  reduction = "protein_pca",
                  group.by = "ident")
      }
    }
    
    #Categorical/Qualitative Variables
    if (!is.element(class(meta.df[[m]][1]), c("numeric", "integer"))) {
      if (!(use.cite.seq)) {
        #plot RNA clusters
        if (reduction == "tsne") {
          clusmat = data.frame(
            umap1 = p1$data$tSNE_1,
            umap2 = p1$data$tSNE_2,
            clusid = as.character(object@meta.data[[m]])
          )
        } else if (reduction == "umap") {
          clusmat = data.frame(
            umap1 = p1$data$UMAP_1,
            umap2 = p1$data$UMAP_2,
            clusid = as.character(object@meta.data[[m]])
          )
        } else {
          clusmat = data.frame(
            umap1 = p1$data$PC_1,
            umap2 = p1$data$PC_2,
            clusid = as.character(object@meta.data[[m]])
          )
        }
        
      } else {
        #else plot Antibody clusters
        if (reduction == "tsne") {
          clusmat = data.frame(
            umap1 = p1$data$protein_tsne_1,
            umap2 = p1$data$protein_tsne_2,
            clusid = as.character(object@meta.data[[m]])
          )
        } else if (reduction == "umap") {
          clusmat = data.frame(
            umap1 = p1$data$protein_umap_1,
            umap2 = p1$data$protein_umap_2,
            clusid = as.character(object@meta.data[[m]])
          )
        } else {
          clusmat = data.frame(
            umap1 = p1$data$protein_pca_1,
            umap2 = p1$data$protein_pca_2,
            clusid = as.character(object@meta.data[[m]])
          )
        }
      }
      
      umap.pos <- clusmat %>% group_by(clusid) %>% summarise(umap1.mean = 
                                                               mean(umap1),
                                                 umap2.mean = mean(umap2))
      title = as.character(m)
      cols = list()
      
      lab.palette <- colorRampPalette(brewer.pal(12, "Paired"))
      n = length(unique((object@meta.data[[m]])))
      cols[[1]] = brewer.pal(8, "Set3")[-2]  #Alternative
      cols[[2]] = brewer.pal(8, "Set1")
      cols[[3]] = c(cols[[1]], brewer.pal(8, "Set2")[3:6])
      cols[[4]] = c("#F8766D",
                    "#FF9912",
                    "#a100ff",
                    "#00BA38",
                    "#619CFF",
                    "#FF1493",
                    "#010407")
      cols[[5]] = c("blue", "red", "grey")
      cols[[6]] = lab.palette(n)
      cols[[7]] = c("red", "green", "blue", "orange", "cyan", "purple")
      cols[[8]] = c(
        "#e6194B",
        "#3cb44b",
        "#4363d8",
        "#f58231",
        "#911eb4",
        "#42d4f4",
        "#f032e6",
        "#bfef45",
        "#fabebe",
        "#469990",
        "#e6beff",
        "#9A6324",
        "#800000",
        "#aaffc3",
        "#808000",
        "#000075",
        "#a9a9a9",
        "#808080",
        "#A9A9A9",
        "#8B7355"
      )
      colnum = 8
      
      n = length(unique(clusmat$clusid))
      qual.col.pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      qual.col.pals = qual.col.pals[c(7, 6, 2, 1, 8, 3, 4, 5),]
      col = unlist(mapply(
        brewer.pal,
        qual.col.pals$maxcolors,
        rownames(qual.col.pals)
      ))
      
      #Select to add labels to plot or not:
      if (show.labels) {
        g <- ggplot(clusmat) +
          theme_bw() +
          theme(legend.title = element_blank()) +
          geom_point(aes(
            x = umap1,
            y = umap2,
            colour = clusid
          ), size = dot.size) +
          scale_color_manual(values = col) +
          xlab(paste(reduction, "-1")) + ylab(paste(reduction, "-2")) +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = legend.position,
            panel.background = element_blank(),
            legend.text = element_text(size = rel(legend.text.size))
          ) +
          guides(colour = guide_legend(override.aes = list(
            size = 5, alpha = 1
          ))) +
          ggtitle(title) +
          geom_label_repel(
            data = umap.pos,
            aes(
              x = umap1.mean,
              y = umap2.mean,
              label = umap.pos$clusid
            ),
            size = 4
          )
      }else{
        g <- ggplot(clusmat, aes(x = umap1, y = umap2)) +
          theme_bw() +
          theme(legend.title = element_blank()) +
          geom_point(aes(colour = clusid), size = dot.size) +
          scale_color_manual(values = col) +
          xlab(paste(reduction, "-1")) + ylab(paste(reduction, "-2")) +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = legend.position,
            panel.background = element_blank(),
            legend.text = element_text(size = rel(legend.text.size))
          ) +
          guides(colour = guide_legend(override.aes = list(
            size = 5, alpha = 1
          ))) +
          ggtitle(title)
      }
      
      
    } else {
      ##THIS IS THE PART WE PLOT QUANTITATIVE DATA
      m = as.character(m)
      clusid = object@meta.data[[m]]
      clusid = scales::rescale(object@meta.data[[m]], to = c(0, 1))
      clus.quant = quantile(clusid[clusid > 0], probs = c(.1, .5, .9))
      midpt.1 = clus.quant[2]
      midpt.2 = clus.quant[1]
      midpt.3 = clus.quant[3]
      
      if (!(use.cite.seq)) {
        #plot RNA clusters
        if (reduction == "tsne") {
          clusmat = data.frame(
            umap1 = p1$data$tSNE_1,
            umap2 = p1$data$tSNE_2,
            clusid = as.numeric(object@meta.data[[m]])
          )
        } else if (reduction == "umap") {
          clusmat = data.frame(
            umap1 = p1$data$UMAP_1,
            umap2 = p1$data$UMAP_2,
            clusid = as.numeric(object@meta.data[[m]])
          )
        } else {
          clusmat = data.frame(
            umap1 = p1$data$PC_1,
            umap2 = p1$data$PC_2,
            clusid = as.numeric(object@meta.data[[m]])
          )
        }
        
      } else {
        #else plot Antibody clusters
        if (reduction == "tsne") {
          clusmat = data.frame(
            umap1 = p1$data$protein_tsne_1,
            umap2 = p1$data$protein_tsne_2,
            clusid = as.numeric(object@meta.data[[m]])
          )
        } else if (reduction == "umap") {
          clusmat = data.frame(
            umap1 = p1$data$protein_umap_1,
            umap2 = p1$data$protein_umap_2,
            clusid = as.numeric(object@meta.data[[m]])
          )
        } else {
          clusmat = data.frame(
            umap1 = p1$data$protein_pca_1,
            umap2 = p1$data$protein_pca_2,
            clusid = as.numeric(object@meta.data[[m]])
          )
        }
      }
      
      clusmat %>% group_by(clusid) %>% summarise(umap1.mean = mean(umap1),
                                                 umap2.mean = mean(umap2)) -> umap.pos
      title = as.character(m)
      print(environmentName(environment(arrange)))
      clusmat %>% dplyr::arrange(clusid) -> clusmat
      print(environmentName(environment(arrange)))
      g <- ggplot(clusmat, aes(x = umap1, y = umap2)) +
        theme_bw() +
        theme(legend.title = element_blank()) +
        geom_point(aes(colour = clusid), size = 1) +
        scale_color_gradientn(
          colours = c("blue4", "lightgrey", "red"),
          values = scales::rescale(c(0, midpt.2, midpt.1, midpt.3, 1),
                                   limits = c(0, 1))
        ) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
        ) +
        ggtitle(title) +
        xlab(paste0(reduction.type, "-1")) + ylab(paste0(reduction.type, "-2"))
    }
    
    return(g)
  }
  
  # Main Code
  meta.df <- object@meta.data
  summarize.cut.off <- min(summarization.cut.off, 20)
  
  print("selected object:")
  print(object)
  
  # converting dots to underscores in column names:
  colnames(object@meta.data) = gsub("\\.", "_", colnames(object@meta.data))
  
  # checking metadata for sanity
  m = metadata.to.plot
  m = m[!grepl("Barcode", m)]
  if (length(m) == 0) {
    print("No metadata columns specified.
           Plotting sample names and RNA clusters...")
    x = colnames(object@meta.data)
    x = x[grepl("RNA", x)]
    m = c("sample_name", x)
  }
  
  #ERROR CATCHING
  #collect valid names of valid columns
  valid.columns <- character()
  for (i in colnames(meta.df)) {
    if (!any(is.na(meta.df[[i]]))) {
      valid.columns <- c(valid.columns, i)
    }
  }
  
  # Checking for content of "Columns to Summarize"
  cols.to.summarize <-
    eval(parse(text = gsub('\\[\\]', 'c()', columns.to.summarize)))
  m = unique(c(m, cols.to.summarize))
  
  if (length(cols.to.summarize) > 0) {
    #Optional Summarization of Metadata
    for (i in cols.to.summarize) {
      col <- meta.df[[i]]
      val.count <- length(unique(col))
      
      if ((val.count >= summarizeCutOff) &
          (i != 'Barcode') &
          (!is.element(class(meta.df[[i]][1]), c("numeric", "integer")))) {
        freq.vals <- as.data.frame(-sort(-table(col)))$col[1:summarize.cut.off]
        print(freq.vals)
        summarized.col = list()
        count <- 0
        for (j in col) {
          print(j)
          print(paste("count is", count))
          
          if (is.na(j) || is.null(j) || (j == "None")) {
            count <- count + 1
            summarized.col[count] <- "NULLorNA"
            print("NULLorNA")
          } else if (j %in% freq.vals) {
            count <- count + 1
            summarized.col[count] <- j
            print("valid")
          } else {
            count <- count + 1
            summarized.col[count] <- "Other"
            print("Other")
          }
        }
        meta.df[[i]] <- summarized.col
      }
    }
    #assign new metadata
    colnames(meta.df) = gsub("\\.", "_", colnames(meta.df))
    object@meta.data <- meta.df
    colnames(object@meta.data) = gsub("\\.", "_", colnames(object@meta.data))
  }
  
  grobs <- lapply(m, function(x) .drawMetadata(x))
  
  print(grobs)
  
  return(NULL)
  
}

