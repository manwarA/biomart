# required liberaries
# library(biomaRt) #biomaRt=2.60.1
# library(stringr)

#=========================================
# Convert gene names using BiomaRt
#=========================================

addGeneName <- function(objK,
                        martObj = eMart,
                        gene_name_col = NULL) {
   
  require("biomaRt") # or library(biomaRt)
  
  # check whether mart object has alredy been created; no need to recreat this object. 
  # TODO option to select the mart varient; currently only human gene ensembl is supported
    if(! exists("martObj", environment())){
        mart = useEnsembl("ensembl","hsapiens_gene_ensembl")
    } else {
      mart <- martObj
        }

  # where are the gene names, either as rownames or as a column, mostly they are as rownames.
    cat("\n\n setting gene name columns\n\n")
    if(is.null(gene_name_col)) {
        objK[["genes"]] <- rownames(objK)
    } else {
        objK[["genes"]] <- objK[[gene_name_col]]
    }

    print(head(objK))
    # remove "." from the ensemble gene ids, TCGA ensembl genes are coded with their version, that are difficult to convert. 
    # the version should be remove

    objK[["genes"]] = ifelse(stringr::str_detect(objK[["genes"]], "\\.") == TRUE,
                       stringr::str_split_i(objK[["genes"]], pattern = "\\.", 1),
                      objK[["genes"]])

    genes = objK[["genes"]]

    cat("\n\n this is the gene list for mart\n\n")
    print(head(genes))

    gene_IDs = getBM(filters= "ensembl_gene_id",
                     attributes= c("ensembl_gene_id","hgnc_symbol",
                                   "gene_biotype", "entrezgene_id"),
                     values = genes,
                     mart= mart)

    objK_with_names = merge(as.data.frame(objK), gene_IDs,
                            by.x = "genes",
                            by.y = "ensembl_gene_id",
                            all.x = TRUE, 
                            sort = FALSE)
  
  # TODO; assign gene names to rownames, for easier downstream usage.
    return(objK_with_names)
}

#addGeneName(df_with_ensem_ids)
