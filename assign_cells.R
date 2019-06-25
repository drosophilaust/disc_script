### This script provide the workflow that calculate Matthews correlation coefficient (MCC) scores 
### from seurat object and binary reference matrix.

# the content in reference matrix is :
#  Region      CR44334 nub CG9008 Cpr49Ah rn  ...
#Pouch&Hinge     1      1      1       1  1  ...
#Notum           0      0      0       0  0  ... 
#PM              0      0      0       0  0  ... 
#binarization    0      0      0       0  0  ...

# object is the seurat object with raw data and normalized data
# object is the binary reference matrix
mcc_calc <- function(object,bin) {
  
  library(DistMap)
  #### Distmap
  raw.data = as.matrix(object@raw.data) #raw.data
  normalized.data = as.matrix(object@data) #normlizrd data by seurat
  raw.data1 = as.data.frame(t(raw.data))
  normalized.data1 =as.data.frame(t(normalized.data))
  com_cells = intersect(rownames(raw.data1),rownames(normalized.data1))
  raw.data1 = subset(raw.data1,rownames(raw.data1) %in% com_cells)
  normalized.data1 = subset(normalized.data1,rownames(normalized.data1) %in% com_cells)
  raw.data1 = as.matrix(t(raw.data1))
  normalized.data1 = as.matrix(t(normalized.data1))
  
  ##
  gene <- colnames(bin)
  gene <- intersect(gene,rownames(raw.data1))
  bin = bin[,gene]
  quant <- bin[4,]
  
  ##
  binary.dt <- as.data.frame(matrix(nrow = ncol(object@data),ncol = ncol(bin)))
  
  for (i in 1:length(gene)) {
    dt <- object@data[gene[i],]
    dt_1 <- dt[dt > 0]
    dt_2 <- as.numeric(dt > as.numeric(quantile(dt_1,quant[1,i])))
    names(dt_2) <- names(dt)
    colnames(binary.dt)[i] <- gene[i]
    binary.dt[,i]<-dt_2
  }
  rownames(binary.dt) <- names(dt)
  
  #
  bin = bin[-4,]
  
  bin = as.matrix(bin)
  
  dm = new("DistMap",
           raw.data=raw.data1,
           data=normalized.data1,
           insitu.matrix= bin )
  
  binary.dt1 <- binary.dt
  binary.dt1 <- binary.dt1[intersect(rownames(binary.dt1),colnames(dm@data)),]
  binary.dt1 <- t(binary.dt1)
  
  dm@binarized.data <- binary.dt1
  
  dm <- mapCells(dm)
  
  # assign region
  re.mcc = t(dm@mcc.scores)
  rownames(re.mcc) = colnames(normalized.data1)
  region = as.data.frame(matrix(nrow = nrow(re.mcc),ncol = 1))
  
  for (i in 1:nrow(re.mcc)) {
    region[i,1] = paste(c(names(which(re.mcc[i,] == max(re.mcc[i,])))),collapse = "&")
  }
  re.mcc_region <- cbind(re.mcc,region)
  colnames(re.mcc_region) = c("Pouch&Hinge","Notum","PM","REGION")
  
  return(re.mcc_region)
}

### re.mcc_region will have the MCC scores of each region and the final assigned region for each cell.








