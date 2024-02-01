#' @import Seurat
#' @import SingleR
#' @import edgeR
NULL
#> NULL
#> 
.scProc <- function(sc, n = 2000, pc =30, res =0.5, normalize = TRUE, hvg = TRUE, pca = TRUE, umap = TRUE, cluster = TRUE) {
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
  }
  return(sc)
}

.getPseudobulk <- function(obj, anno, sample, keep.all.genes = FALSE, min.cell = 5, assay = "RNA") {
  
  obj$anno_sample = paste(as.character(obj@meta.data[,anno]), obj@meta.data[,sample], sep = ".")
  a = table(obj$anno_sample)
  a = a[a>min.cell]
  whole.psk = lapply(names(a), function(x) {
    if (sum(obj$anno_sample==x)>1) {
      rowSums(obj@assays[[assay]]@counts[,obj$anno_sample==x])
    } else {
      obj@assays[[assay]]@counts[,obj$anno_sample==x]
    }
    
  })
  whole.psk = do.call(cbind, whole.psk)
  colnames(whole.psk) = names(a)
  if (keep.all.genes) {
    whole.psk = edgeR::DGEList(counts = whole.psk)
  } else {
    whole.psk = edgeR::DGEList(counts = whole.psk[rowSums(whole.psk)>0, ])
  }
  
  whole.psk = edgeR::calcNormFactors(whole.psk)
  logtpm = log2(edgeR::cpm(whole.psk)+1)
  s = sapply(colnames(whole.psk), function(x) strsplit(x, "[.]")[[1]])
  s = as.data.frame(t(s))
  colnames(s) = c("anno", "sample")
  whole.psk$samples = cbind.data.frame(whole.psk$samples, s)
  return(list(dgelist = whole.psk, logtpm = logtpm))
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
#' @param ref.tissue reference tissue type, default is "all". Other options include "immune" and specific tissue(s). To see list of available tissues, use list.tissues(). Multiple tissues separated by "|"
#' @param resolution resolution used to identify cell clusters, default is 0.5
#' @return Seurat object with predicted.celltype added to metadata
#' @export
cellAnno <- function(query.srt, ref.tissue = "all", resolution = 0.5) {
  if (ref.tissue == "all") {
    ref = all.ref
  } else if (ref.tissue == "immune") {
    ref = immgen.ref
  } else {
    ref = all.ref[,grepl(ref.tissue, colnames(all.ref))| colnames(all.ref) %in% c(vbls$stromal.cells, vbls$immune.cells)]
  }
  query.srt = .scProc(query.srt, res = resolution)
  query.srt$sample = "query"
  psk = .getPseudobulk(obj = query.srt, sample = "sample", anno = "seurat_clusters", min.cell = 2)
  pred = SingleR::SingleR(test = psk$logtpm, ref = ref, assay.type.ref = "logcounts", 
                          assay.type.test = "logcounts", labels = colnames(ref))
  cluster.map = setNames(psk$dgelist$samples$anno, pred$labels)
  query.srt$predicted.celltype = cluster.map[as.character(query.srt$seurat_clusters)]
  return(query.srt)
}
#' Cell-of-Origin prediction for cancer bulk RNA-seq data
#'
#' Predict COO using reference built from Tabula sapiens, human protein atlas and public datasets. Right now it only support prediction for solid tumor type
#' @param logtpm.matrix query data matrix in log2-TPM
#' @param ref reference to use, default is "all.epithelial". To see list of available tissues, use list.tissues(). Multiple tissues separated by "|"
#' @return output of SingleR. 
#' @export
cooPredict <- function(logtpm.matrix, ref = "all.epithelial") {
  if (ref=="all.epithelial") {
    refdat = all.ref[,1:184]
  } else {
    refdat = all.ref[,grepl(ref, colnames(all.ref))]
  }
  pred = SingleR::SingleR(test = logtpm.matrix, ref = refdat, assay.type.ref = "logcounts",
                          assay.type.test = "logcounts",
                          labels = colnames(ref))
  return(pred)
}

