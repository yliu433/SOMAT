

# download the phenotypes of LIHC from TCGA
pheno <- read.table("https://tcga.xenahubs.net/download/TCGA.LIHC.sampleMap/LIHC_clinicalMatrix",
                  header=TRUE,sep='\t',fill=TRUE)

# download the somatic gene-level non-silent mutation (bcm automoated)
bcm <- read.table("https://tcga.xenahubs.net/download/TCGA.LIHC.sampleMap/mutation_bcm_gene",
                  header=TRUE,sep='\t',fill=TRUE)

# fill one missing gene name with NONE and organize the data
bcm[,1] <- as.character(bcm[,1])
bcm[is.na(bcm[,1]),1] <- "NONE"

bcm2 <- t(bcm[,-1])
rownames(bcm2) <- gsub("\\.","-",rownames(bcm2))
colnames(bcm2) <- bcm[,1]

# obtain the common sample 
samc <- intersect(rownames(bcm2),pheno$sampleID) # 373 samples

# confounders with age and gender
X0 <- pheno[match(samc,pheno$sampleID),c("age_at_initial_pathologic_diagnosis","gender")]
X0 <- cbind(1,X0)
colnames(X0) <- c("intercept","age","gender")
X0[,3] <- ifelse(X0[,3]=="MALE",1,0)
rownames(X0) <- samc

# outcomes including AFP level and PT

AFP <- pheno[match(samc,pheno$sampleID),"fetoprotein_outcome_value"]
# lower limit
AFP_lwl <- pheno[match(samc,pheno$sampleID),"fetoprotein_outcome_upper_limit"]
# dichotomize the AFP
AFP_level <- ifelse(AFP<AFP_lwl,0,1)

PT <- pheno[match(samc,pheno$sampleID),"prothrombin_time_result_value"]

Y0 <- cbind(AFP_level,PT)

# remove the missing values and obtain the final confounder matrix and the outcome matrix
XY0 <- na.omit(cbind(X0,Y0))

X <- XY0[,1:3] # confounders
Y <- XY0[,4:5] # outcomes

table(Y[,1])

# obtain the whole gene-level somatic matrix
G0 <- bcm2[match(rownames(X),rownames(bcm2)),]

save(Y,X,G0,file="LIHC.Rdata")
























