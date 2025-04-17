# CellAnno
CellAnno is a automatic cell annotation tool with built-in human cell reference. To make one unified reference that can be used in any dataset, we compiled, curated and integrated single cell datasets from multiple sources including Human Protein Atlas, Tabula Sapiens, ArrayExpress and GEO and created a comprehensive human cell gene expression reference database named CellMap. It covers 200 cell types/states from 42 human tissue types, including 137 epithelial cell types, 22 stromal cell types and 19 major immune cell types. CellAnno uses the correlation-based algorithm [SingleR](https://www.bioconductor.org/packages/release/bioc/html/SingleR.html) to assign the most similar cell type from the reference. To make the process efficient and increase clarity of the results, CellAnno applies annotation at cluster level rather than single cell level. For immune cells, an option to run fine grained annotation (e.g. Treg, CD8+ effector T-cell) is available. With this option, a reference generated from ImmGen datasets ([GSE15907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15907),[GSE37448](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37448)) is used.

**Installation**

`devtools::install_github("yupingz/CellAnno")`

**Functions**

Type `?cellAnnotate` and `?cooPredict` to see help information

**Tutorial**
https://rpubs.com/yupingz/1297155

