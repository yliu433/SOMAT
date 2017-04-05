

require("SKAT")
source("Davis.R") # mixture chisquare approximation by Davis method
source("SOMAT.R") # SOMAT function

load("LIHC.Rdata") # data from step 1

# download hallmark gene sets (gene symbols) from gsea to the current directory
# read the gmt file of hallmark pathway by GSA package
require("GSA")
hp <- GSA.read.gmt("h.all.v5.2.symbols.gmt")

n_path <- length(hp$genesets) # 50 pathways

# pvales of the fSOMAT and oSOMAT
pvalM <- matrix(NA,n_path,2)
colnames(pvalM) <- c("fSOMAT","oSOMAT")
rownames(pvalM)<-hp$geneset.names

Y<-as.matrix(Y)
X<-as.matrix(X)

# analysis
for(i in 1:n_path){
  
  genes <- hp$genesets[[i]] # genes in the pathway
  
  genes_common <- intersect(genes,colnames(G0)) # common genes with the mutation data
  
  if(length(genes_common)>0){
    index_gene <- match(genes_common,colnames(G0))
    G1 <- as.matrix(G0[,index_gene])

    index_thres<-which(colSums(G1)>2) 
    
    if(length(index_thres)>1){
      G<-as.matrix(G1[,index_thres])
      
      p<-dim(G)[2] # number of parameters
      W<-matrix(1,1,p)

      ##########################
      # SOMAT
      ###########################
      
      pval_SOMAT<-SOMAT(Y,X,G,W)
      
      pvalM[i,"fSOMAT"]<-pval_SOMAT$pval_fSOMAT
      pvalM[i,"oSOMAT"]<-pval_SOMAT$pval_oSOMAT
    }
  }
  
  
}
  
which(pvalM<=0.05/50,arr.ind = T)  







