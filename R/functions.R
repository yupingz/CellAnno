#' @import Seurat
#' @import SingleR
#' @import edgeR
#' @import Matrix
NULL
#> NULL
#> 

.clusterReorder = function(cts) {
  d = Embeddings(cts,"umap")
  d = as.data.frame(d)
  d$cluster = cts$seurat_clusters
  cls.pos = aggregate(.~ cluster, FUN=mean, data = d)
  cls = cls.pos$cluster
  
  cls.dis = lapply(1:nrow(cls.pos), function(i) {
    r = t(cls.pos[,2:3]) - as.numeric(cls.pos[i,2:3])
    colnames(r) = cls.pos$cluster
    r = r^2
    r = colSums(r)
    return(sort(r))
  })
  names(cls.dis) = cls
  
  cls.ord = as.character(cls.pos$cluster[order(cls.pos$umap_1)][1])
  
  for (i in 1:(length(cls)-1)) {
    cls.left = cls.dis[[as.character(cls.ord[i])]]
    cls.left = cls.left[setdiff(cls, cls.ord)]
    cls.left = sort(cls.left)
    cls.ord = c(cls.ord, names(cls.left)[1])
  }
  
  cls.reorder = setNames(1:length(cls),cls.ord)
  cts$cluster = as.character(cls.reorder[as.character(cts$seurat_clusters)])
  cts$cluster = as.numeric(cts$cluster)
  return(cts)
}


.scProc <- function(sc, n = 2000, pc =30, res =0.5, normalize = TRUE, hvg = TRUE, pca = TRUE, umap = TRUE, cluster = TRUE, reorder = TRUE) {
  DefaultAssay(sc) = "RNA"
  if (normalize) {
    sc <- Seurat::NormalizeData(sc, verbose = FALSE)
  } 
  if (hvg) {
    sc <- Seurat::FindVariableFeatures(sc, selection.method = "vst", nfeatures = n, verbose = FALSE)
  }
  if (pca) {
    sc <- Seurat::ScaleData(sc, verbose = FALSE)
    sc <- Seurat::RunPCA(sc, npcs = pc, verbose = FALSE)
  }
  if (umap) {
    sc <- Seurat::RunUMAP(sc, reduction = "pca", dims = 1:pc)
  }
  
  if (cluster) {
    sc <- Seurat::FindNeighbors(sc, dims=1:pc)
    sc <- Seurat::FindClusters(sc, resolution = res)
    if (reorder) {
      sc <- .clusterReorder(cts = sc)
    }
  }
  return(sc)
}

.getPseudobulk <- function(obj, cluster, keep.all.genes = FALSE, min.cell = 5, assay = "RNA") {
  
  require(edgeR)
  obj$cluster = as.character(obj@meta.data[, cluster])
  cls.no = table(obj$cluster)
  cls.no = cls.no[cls.no>=min.cell]
  obj = obj[, obj$cluster %in% names(cls.no)]
  psk = Seurat::AggregateExpression(object = obj, assays = assay, return.seurat = F, group.by = "cluster")
  psk = psk[[1]]
  colnames(psk) = names(cls.no)
  if (keep.all.genes) {
    psk = DGEList(counts = psk)
  } else {
    psk = DGEList(counts = psk[rowSums(psk)>0, ])
  }
  psk = calcNormFactors(psk)
  logtpm = log2(edgeR::cpm(psk)+1)
  
  return(list(dgelist = psk, logtpm = logtpm))
}

#' List available tissue types for reference
#'
#' @export
list.tissues <- function() {
  print(vbls$tissues)
}
#' Annotation for single cell data at cluster level
#'
#' Annotate single cell at cluster level using reference built from Tabula sapiens, human protein atlas and public datasets. The annotation algorithm used is SingleR.
#' @param query.srt query single cell data in Seurat object format
#' @param ref.tissue reference tissue type, default is "core". Other options include "immune" and specific tissue(s). To see list of available tissues, use list.tissues(). Multiple tissues separated by "|"
#' @param cluster user provided cluster assignment to use. it should be a column name in meta.data. if NULL, will use Seurat::FindCluster to find clusters
#' @param resolution resolution used to identify cell clusters, Default is 0.5 if cluster is not set.
#' @param level2 if to run an additional annotation step to get detailed cell subtypes. It only works when ref.tissue = "immune".
#' @return Seurat object with "labels","pruned.labels" and "score" added to metadata
#' @export
cellAnnotate <- function(query.srt, ref.tissue = "core", cluster = NULL, resolution = 0.5, level2 = FALSE) {
  if (ref.tissue == "core") {
    ref = all.ref
    ref.i = which(ref$meta$major !="specialized")
    ref = list(expr = ref$expr[, ref.i], meta = ref$meta[ref.i, ])
    
  } else if (ref.tissue == "immune") {
    ref = immgen.ref
  } else {
    ref = all.ref
    ref.i = which(grepl(ref.tissue, ref$meta$tissue) | ref$meta$celltype %in% c(vbls$stromal.cells, vbls$immune.cells))
    ref = list(expr = ref$expr[, ref.i], meta = ref$meta[ref.i, ])
  }
  if (is.null(cluster)) {
    query.srt = .scProc(query.srt, res = resolution)
    cls = as.character(query.srt@meta.data[, "cluster"])
    psk = .getPseudobulk(obj = query.srt,  cluster = "cluster", min.cell = 2)
  } else {
    cls =  as.character(query.srt@meta.data[, cluster])
    psk = .getPseudobulk(obj = query.srt,  cluster = cluster, min.cell = 2)
  }
  pred = SingleR::SingleR(test = psk$logtpm, ref = ref$expr, assay.type.ref = "logcounts", 
                          assay.type.test = "logcounts", labels = ref$meta$level1, fine.tune = FALSE)
  cluster.map = data.frame(cluster = colnames(psk$logtpm),labels= pred$labels,
                           pruned.labels = pred$pruned.labels, score = apply(pred$scores, 1, max),
                           stringsAsFactors = F)
  cluster.meta = cluster.map[match(cls, cluster.map$cluster),]
  rownames(cluster.meta) = colnames(query.srt)
  
  if (level2 & ref.tissue == "immune") {
    celltypes = c("b_cells", "dendritic_cells", "granulocytes",  "macrophages",  "monocytes",  "stem_cells", "stromal_cells", "t_cells")
    
    celltypes = intersect(celltypes, cluster.map$labels)
    if (length(celltypes)>0) {
      pred.level2 = lapply(celltypes, function(x) {
        query.subset = query.srt[,cluster.meta$labels ==x]
        query.subset = .scProc(sc = query.subset, res = 0.5, reorder = FALSE)
        if (length(levels(Idents(query.subset)))==1) {
          kmean.cls = kmeans(Embeddings(query.subset, "umap"), centers = 2)
          query.subset$seurat_clusters = kmean.cls$cluster
        }
        psk.subset = .getPseudobulk(obj = query.subset, cluster = "seurat_clusters",min.cell = 2)
        ref.i = which(ref$meta$level1 == x & ref$meta$level2!="")
        
        p = SingleR::SingleR(test = psk.subset$logtpm, ref = ref$expr[, ref.i], assay.type.ref = "logcounts", 
                             assay.type.test = "logcounts", labels = ref$meta$level2[ref.i], fine.tune = FALSE)
        m = data.frame(row.names = colnames(psk.subset$logtpm), level1 = x, level2 = p$labels, level2.filter = p$pruned.labels, 
                       level2.score = apply(p$scores, 1, max), stringsAsFactors = F)
        cell.meta = m[as.character(query.subset$seurat_clusters),]
        rownames(cell.meta) = colnames(query.subset)
        return(cell.meta)
      })
      pred.level2 = do.call(rbind, pred.level2)
      cluster.meta = merge(cluster.meta, pred.level2, by = "row.names", all.x = T)
      cluster.meta$level1 = cluster.meta$labels
      cluster.meta$level2[is.na(cluster.meta$level2)] = cluster.meta$level1[is.na(cluster.meta$level2)]
      cluster.meta$level2.filter[is.na(cluster.meta$level2)] = cluster.meta$level1[is.na(cluster.meta$level2)]
      rownames(cluster.meta) = cluster.meta$Row.names
      cluster.meta = cluster.meta[,-1]
      
    } else {
      print("predicted level-1 cells do not have level-2 annotation")
    }
  } 
  
  
  query.srt = AddMetaData(query.srt, cluster.meta)
  return(query.srt)
}
#' Cell-of-Origin prediction for cancer bulk RNA-seq data
#'
#' Predict COO using reference built from Tabula sapiens, human protein atlas and public datasets. Right now it only support prediction for solid tumor type
#' @param logtpm.matrix query data matrix in log2-TPM
#' @param ref reference to use, default is "core".To see list of available tissues, use list.tissues(). Multiple tissues separated by "|". Set to "core" if not sure (for samples from metastatic sites, it is possible top hit is cell from the biopsy site). 
#' @return output of SingleR. "labels" is the predicted cell type. 
#' @export
cooPredict <- function(logtpm.matrix, ref = "core", exclude.tissues = NULL) {
  refdat = all.ref
  if (ref=="core") {
    
    ref.i = which(refdat$meta$major =="epithelial" )
    refdat = list(expr = refdat$expr[, ref.i], meta = refdat$meta[ref.i, ])
    
  } else {
    
    ref.i = which(grepl(ref, refdat$meta$tissue) & refdat$meta$major =="epithelial" )
    refdat = list(expr = refdat$expr[, ref.i], meta = refdat$meta[ref.i, ])
    
  }
  if (!is.null(exclude.tissues)) {
    ref.i = which(!grepl(exclude.tissues, refdat$meta$tissue))
    refdat = list(expr = refdat$expr[, ref.i], meta = refdat$meta[ref.i, ])
  }
  pred = SingleR::SingleR(test = logtpm.matrix, ref = refdat$expr, assay.type.ref = "logcounts",
                          assay.type.test = "logcounts",
                          labels = refdat$meta$celltype)
  return(pred)
}
#' Extract Top Markers from Seurat FindAllMarkers Results
#'
#' @param mks Data frame output from Seurat FindAllMarkers, including 'gene' and 'cluster' columns.
#' @param n Number of genes to extract per cluster. Default is 5.
#' @return A vector of unique top markers from each cluster.
#' @export
getTopMarkers <- function(mks, n = 5) {
  if (!"cluster" %in% colnames(mks) || !"gene" %in% colnames(mks)) {
    stop("The input must have 'cluster' and 'gene' columns.")
  }
  
  top_genes <- lapply(unique(mks$cluster), function(cluster_id) {
    genes <- mks$gene[mks$cluster == cluster_id]
    genes[1:min(n, length(genes))]
  })
  
  unique_genes <- unique(unlist(top_genes))
  return(unique_genes)
}
