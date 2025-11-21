
#=========================================
# Convert gene names using BiomaRt
#=========================================

#biomaRt=2.60.1

addGeneName <- function(objK,
                        martObj = eMart,
                        gene_name_col = NULL) {
    require("biomaRt")
  # check whether mart has been created; no need to recreat this object
    if(! exists("martObj", environment())){
        mart = useEnsembl("ensembl","hsapiens_gene_ensembl")
    } else {
          
      mart <- martObj
        }

    cat("\n\n setting gene name columns\n\n")
    if(is.null(gene_name_col)) {
        objK[["genes"]] <- rownames(objK)
    } else {
        objK[["genes"]] <- objK[[gene_name_col]]
    }

    print(head(objK))
    # remove "." from the ensemble gene ids

    objK[["genes"]] = ifelse(stringr::str_detect(objK[["genes"]], "\\.") == TRUE,
                       stringr::str_split_i(objK[["genes"]], pattern = "\\.", 1),
                      objK[["genes"]])
    #print(head(objK))
    genes = objK[["genes"]]

    cat("\n\n this is the gene list for mart\n\n")

    #print(head(genes))

    gene_IDs = getBM(filters= "ensembl_gene_id",
                     attributes= c("ensembl_gene_id","hgnc_symbol",
                                   "gene_biotype", "entrezgene_id"),
                     values = genes,
                     mart= mart)

    objK_with_names = merge(as.data.frame(objK), gene_IDs,
                            by.x = "genes",
                            by.y = "ensembl_gene_id",
                            all.x = TRUE)
    return(objK_with_names)

}

#addGeneName(blca_test)
