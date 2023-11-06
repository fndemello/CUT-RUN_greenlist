# Load required libraries
library(entropy)
library(doParallel)
library(foreach)
library("doRNG")
options(digits = 15)

setwd("/work2/group/prostata/metodo_normalizacao/")


# Input data
total_counts <- read.delim(file='5.Quant/reordered.clean_over2Msamp_1k-bins.tsv',header=T,row.names=1)
isAllowed <- read.delim(file='5.Quant/clean_bins_isAwayFromGene.tsv',header=F,row.names=NULL)

# Process data
target_counts <- total_counts[1:294825,]
bincount <- 2949
isAllowed <- isAllowed$V1



#####################
# Test normalization function
myfun = function(test_counts, bincount, target_counts, isAllowed){
  
  # Random sampling
  # test_counts <- total_counts[sample(1:nrow(total_counts),bincount,replace=F),]
  
  # Test sampled bins as normalizers
  smm <- summary(
    apply(base::scale(target_counts, center = F, scale = 
                                colSums(test_counts)/median(colSums(test_counts)))
                  ,1,function(x) entropy(x,method='CS'))
  )
  v <- NULL
  v <- smm[3]
  v <- append(v, smm[5]-smm[2])
  
  # Filter, then test as normalizers
  test_counts <- test_counts[rownames(test_counts)%in%isAllowed,]
  smm <- summary(
    apply(base::scale(target_counts, center = F, scale = 
                                colSums(test_counts)/median(colSums(test_counts)))
                  ,1,function(x) entropy(x,method='CS'))
  )
  v <- append(v, smm[3])
  v <- append(v, smm[5]-smm[2])
  names(v) <- NULL
  gc()
  rm()
  return(v)
}
  
# Monte-Carlo set up
results <- NULL
trials <- 100000
preload <- 2000
# To avoid passing the entire counts table to each worker, the
# script pre-samples a set number of permutations at a time. Pre-
# loading more sets at a time speeds computation, but consumes
# more memory.

# Make multithreading cluster
cl <- makeCluster(20)
registerDoParallel(cl)

RNGkind("L'Ecuyer-CMRG")
set.seed(0)
mc.reset.stream()

# MC run
start=Sys.time()

for (r in 1:(trials/preload)) {
  tests <- NULL
  intermediate_res <- NULL
  for (t in 1:preload) {
    tests <- append(tests,list(total_counts[sample(1:nrow(total_counts),bincount,replace=F),]))
  }
  intermediate_res <- foreach(i=tests, .combine = 'rbind',.packages = 'entropy',.noexport = c('total_counts','tests')) %dopar% {
    myfun(i, bincount, target_counts, isAllowed)
  }
  results <- rbind(results, intermediate_res)
  write.table(paste('Done:',r*preload,'of',trials,'||',Sys.time()-10800),file='6.Calculos/MC.progress.txt',quote=F)
}

end=Sys.time()
timing = end-start
timing
stopCluster(cl)
results <- results[1:trials,1:4]
colnames(results) <- c('Raw_median','Raw_IQR','Filtered_median','Filtered_IQR')
rownames(results) <- c(1:trials)

write.table(results,file='6.Calculos/MC.test_100kruns.tsv',quote=F,sep='\t')
write.table(timing,file='6.Calculos/MC.test_timefor_100kruns.tsv',quote=F,row.names = F,col.names = F)
