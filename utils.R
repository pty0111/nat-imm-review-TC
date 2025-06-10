'%ni%' <- function(x,y)!('%in%'(x,y))

cols <- rev(RColorBrewer::brewer.pal(11,"Spectral"))
transitions <- c(0, 25, 30, 35, 40, 50, 60, 65, 70, 75, 100)
scaled_transitions <- scales::rescale(transitions, from = c(0, 100), to = c(0, 1))
scale.color <- scale_color_gradientn(colors = cols, values = scaled_transitions)

plot.QC.violin <- function(cell.md, ident, feature, yintercept = NULL, br, Log10 = F, pal = NULL, fontsize.xlab = 10){
  if(feature == "nCount_RNA"){ name <- "number of transcripts"
  } else if(feature == "nFeature_RNA"){ name <- "number of detected genes"
  } else if(feature == "MtFrac_RNA"){ name <- "percentage of mitochondrial transcripts"
  } else if(feature == "percent.MT"){ name <- "percentage of mitochondrial transcripts"
  } else if(feature == "S.Score"){ name <- "S-phase score"
  } else if(feature == "G2M.Score") {name <- "G2/M-phase score"
  } else if(feature == "ATAC.nFrags") {name <- "ATAC.nFrags"
  } else if(feature == "ATAC.TSSEnrichment") {name <- "ATAC.TSSEnrichment"
  } else if(feature == "ATAC.DoubletEnrichment") {name <- "ATAC.DoubletEnrichment"
  } else if(feature == "ATAC.DoubletScore") {name <- "ATAC.DoubletScore"
  } else if(feature == "nFrags") {name <- "nFrags"
  } else if(feature == "TSSEnrichment") {name <- "TSSEnrichment"
  } else if(feature == "DoubletEnrichment") {name <- "DoubletEnrichment"
  }else if(feature == "DoubletScore") {name <- "DoubletScore"
  } else name <- feature
  ft <- cell.md[, feature]
  if(Log10){
    ft <- log10(ft)
    if(!is.null(yintercept)) yintercept <- log10(yintercept)
  }
  gp <- ggplot(data.frame(id = cell.md[, ident], ft = ft),
               aes(x = id, y = ft, color = id, fill = id)) +
    geom_violin(show.legend = F) +
    geom_hline(yintercept = yintercept, linetype = 2, color = "violetred4") +
    scale_y_continuous(breaks = seq(0, max(ft), by = br)) +
    theme_classic() +
    theme(axis.text.x = element_text(size=fontsize.xlab)) +
    labs(x = ident, y = "", title = ifelse(Log10, paste0("log10(", name, ")"), name))
  if(!is.null(pal))
    gp <- gp + scale_color_manual(values = pal) + scale_fill_manual(values = pal)
  return(gp)
}

plot.all.QC <- function(obj = NULL, obj.type = "Seurat", ident = 'sample', 
                            thr = NULL, col = NULL, CC = F){
  if (obj.type == "Seurat"){
    cell.md <- obj@meta.data
    curr.assay <- "RNA"
    L <- labs(color = ident, x = "number of transcripts", y = "number of detected genes")
  } else{
    print("Invalid input")
    return(1)
  }
  if (curr.assay=='RNA'){
    gp1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = ident) + theme_classic() + L
    if(!is.null(col)){
      names(col) <- levels(cell.md[, ident])
      gp1 <- gp1 + scale_color_manual(values = col)
    }
    gp2 <- ggplot(cell.md, aes(x = nCount_RNA, y = nFeature_RNA)) + geom_density_2d_filled() + theme_classic() + L
    
    pl <- list(
      plot.QC.violin(cell.md, ident = ident, feature = "nCount_RNA", Log10 = T, yintercept = c(thr$nc.min, thr$nc.max), br = 0.2, pal = col),
      plot.QC.violin(cell.md, ident = ident, feature = "nFeature_RNA", yintercept = c(thr$nf.min, thr$nf.max), br = 500, pal = col),
      plot.QC.violin(cell.md, ident = ident, feature = "MtFrac_RNA", yintercept = thr$mp.max, br = 2, pal = col)
    )
    if(is.null(thr)){
      pl[[4]] <- gp1
      pl[[5]] <- gp2
    } else {
      hl <- geom_hline(yintercept = thr$nf.min, linetype = 2, color = "violetred4")
      vl <- geom_vline(xintercept = thr$nc.min, linetype = 2, color = "violetred4")
      pl[[4]] <- gp1 + hl + vl
      pl[[5]] <- gp2 + hl + vl
    }
    if (CC) {
      pl[[6]] <- plot.QC.violin(cell.md, ident = ident, feature = "S.Score", br = 0.2, pal = col)
      pl[[7]] <- plot.QC.violin(cell.md, ident = ident, feature = "G2M.Score", br = 0.2, pal = col)
    }
  } 
  return(pl)
}

plot.clusters <- function(
    SRO = NULL, idx = NULL,
    vis = SRO@reductions$umap@cell.embeddings, 
    axis.titles = c("UMAP1", "UMAP2"),
    groups = SRO$Clusters, clusters.col = 'Clusters', 
    col = pal[[clusters.col]], new.levels = NULL,
    pref.C = T, labels = T, line = T,
    point.size = 1.5, point.alpha = 0.7,
    label.size = point.size * 2, label.pad = 1, border.size = 0.5, 
    legend.size = 3, show_legend = T
){
  if(class(groups) == "factor"){ 
    clusters <- groups
  } else {
    clusters <- factor(groups, levels = sort(unique(groups)))
  }
  gp <- ggplot() + theme_classic() +
    labs(x = axis.titles[1], y = axis.titles[2], color = clusters.col) +
    theme(axis.ticks = element_blank(), axis.text = element_blank())
  if(!line) gp <- gp + theme(axis.line = element_blank())
  if(!is.null(idx)){
    vis <- vis[idx, ]
    if(is.null(new.levels)){
      new.levels <- as.character(sort(unique(clusters[idx])))
    }
    clusters <- factor(as.character(clusters[idx]), levels = new.levels)
    if(!is.null(col)) col <- col[new.levels]
  }
  gp <- gp + geom_point(
    mapping = aes(x = vis[, 1], y = vis[, 2], color = clusters),
    size = point.size, alpha = point.alpha, show.legend = show_legend
  )
  if(!is.null(col)) gp <- gp + scale_color_manual(values = col)
  if(labels){
    C <- levels(clusters)
    CL <- paste0(ifelse(pref.C, "C", ""), levels(clusters))
    vis.cl <- data.frame(matrix(NA, nrow = length(C), ncol = 2, dimnames = list(C, NULL)))
    for(c in C){
      vis.cl[c, 1] <- median(vis[clusters == c, 1], na.rm = T)
      vis.cl[c, 2] <- median(vis[clusters == c, 2], na.rm = T)
    }
    gp + geom_label_repel(
      aes(x = vis.cl[, 1], vis.cl[, 2], color = C, label = CL),
      color = 'black', seed = 1,
      label.size = border.size,
      size = label.size, label.padding = unit(label.pad, "mm"), show.legend = F
    ) +
      geom_label_repel(
        aes(x = vis.cl[, 1], vis.cl[, 2], color = C, label = CL),
        size = label.size, seed = 1,
        label.size = NA, label.padding = unit(label.pad, "mm"), show.legend = F
      )+ guides(colour = guide_legend(override.aes = list(size = legend.size)))
  } else gp + guides(colour = guide_legend(override.aes = list(size = legend.size)))
}

plot.clusters.highlight.one <- function(
    SRO = NULL, idx1 = NULL, idx2 = NULL,
    vis = SRO@reductions$umap@cell.embeddings, 
    axis.titles = c("UMAP1", "UMAP2"),
    groups = SRO$Clusters, clusters.col = 'Clusters', 
    col = pal[[clusters.col]], new.levels = NULL,
    pref.C = T, labels = T, line = T,
    point.size = 1.5, point.alpha = 0.7,
    label.size = point.size * 2, label.pad = 1
){
  if(class(groups) == "factor"){ 
    clusters <- groups
  } else {
    clusters <- factor(groups, levels = sort(unique(groups)))
  }
  gp <- ggplot() + theme_classic() +
    labs(x = axis.titles[1], y = axis.titles[2], color = clusters.col) +
    theme(axis.ticks = element_blank(), axis.text = element_blank())
  if(!line) gp <- gp + theme(axis.line = element_blank())
  if(!is.null(idx1)){
    vis1 <- vis[idx1, ]
    if(is.null(new.levels)){
      new.levels <- as.character(sort(unique(clusters[idx1])))
    }
    clusters1 <- factor(as.character(clusters[idx1]), levels = new.levels)
    if(!is.null(col)) col <- col[new.levels]
  }
  if(!is.null(idx2)){
    vis2 <- vis[idx2, ]
    if(is.null(new.levels)){
      new.levels <- as.character(sort(unique(clusters[idx2])))
    }
    clusters2 <- factor(as.character(clusters[idx2]), levels = new.levels)
  }
  gp <- gp + geom_point(
    mapping = aes(x = vis1[, 1], y = vis1[, 2]), color = 'lightgray',
    size = point.size, alpha = point.alpha
  ) + geom_point(
    mapping = aes(x = vis2[, 1], y = vis2[, 2], color = clusters2),
    size = point.size, alpha = point.alpha, show.legend = T
  )
  if(!is.null(col)) gp <- gp + scale_color_manual(values = col)
  if(labels){
    C <- levels(clusters2)
    CL <- paste0(ifelse(pref.C, "C", ""), levels(clusters2))
    vis.cl <- data.frame(matrix(NA, nrow = length(C), ncol = 2, dimnames = list(C, NULL)))
    for(c in C){
      vis.cl[c, 1] <- median(vis2[clusters2 == c, 1], na.rm = T)
      vis.cl[c, 2] <- median(vis2[clusters2 == c, 2], na.rm = T)
    }
    gp + geom_label_repel(
      aes(x = vis.cl[, 1], vis.cl[, 2], label = CL), 
      color = 'black',
      seed=1,
      size = label.size, 
      label.size = 0.5,
      label.padding = unit(label.pad, "mm"), show.legend = F
    )
  } else gp
}

plot.continuous.value <- function(
    sro = NULL, idx = NULL,
    vis = sro@reductions$umap@cell.embeddings,
    val, val.name, axis.blank = T, axis.titles = c("UMAP1", "UMAP2"),
    scale.color=scale_color_distiller(palette = "Spectral"), point.size = 1
){
  gp <- ggplot() +
    geom_point(
      mapping = aes(x = vis[idx, 1], y = vis[idx, 2], color = val[idx]),
      size = point.size, alpha = 0.8
    ) +
    scale.color + theme_classic() + labs(x = axis.titles[1], y = axis.titles[2], color = val.name)
  if (axis.blank){
    gp <- gp + theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
  } else{
    gp
  }
  return(gp)
}

RenameGenesSeurat <- function(obj, newnames) {
  # from Seurat.utils
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]] <- newnames
    if (length(RNA@scale.data)) rownames(RNA@scale.data) <- newnames
    if (length(RNA@meta.features)) {
      RNA@meta.features$prev.names <- rownames(RNA@meta.features)
      rownames(RNA@meta.features) <- newnames
    }
    if (length(RNA@var.features)){
      RNA@var.features <- rownames(RNA@meta.features[match(RNA@var.features, RNA@meta.features$prev.names), ])
    }
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}

get.named.vector <- function(md, colname){
  val.vec <- md[,colname]; names(val.vec) <- rownames(md)
  return(val.vec)
}

