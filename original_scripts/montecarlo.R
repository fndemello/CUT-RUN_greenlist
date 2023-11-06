# Load required libraries
library(entropy)
library(doParallel)
library(foreach)
library("doRNG")
options(digits = 15)


# Input data
raw_counts <- read.delim(file='')
# Assumes bin count input as described in 'norm+entropy.R'

isAllowed <- read.delim(file='')
# Assumes binary vector of whether genomic bins are close to known genes (FALSE) or intergenic (TRUE);
# this is used for simulating the filtering step described in our manuscript, but may be discarded
# by future users if preferred. We considered a distance of at least 5kb for a bin to be intergenic.

#####
# Process data
target_counts <- raw_counts[1:294825,]
# Similarly to 'threshold_test.R' we selected the 10% of bins with the lowest entropies to test the
# normalization upon; here, they are placed as the first 294825 rows of our 'raw_counts' file to
# simplify input. 

bincount <- 2949
# This variable refers to how many bins should be selected per Monte-Carlo run. 2949 refers to 0.1%
# of non-empty bins in our pipeline.


#####################
# Main test function
myfun = function(test_counts, bincount, target_counts, isAllowed){

  # Test sampled bins as normalizers
  smm <- summary(
    apply(base::scale(target_counts, center = F, scale = 
                                colSums(test_counts)/median(colSums(test_counts)))
                  ,1,function(x) entropy(x,method='ML'))
  )
  v <- NULL
  v <- smm[3]
  v <- append(v, smm[5]-smm[2])
  
  # Filter, then test as normalizers
  test_counts <- test_counts[rownames(test_counts)%in%isAllowed,]
  smm <- summary(
    apply(base::scale(target_counts, center = F, scale = 
                                colSums(test_counts)/median(colSums(test_counts)))
                  ,1,function(x) entropy(x,method='ML'))
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
###
# To avoid passing the entire counts table to each worker thread, the script pre-samples 
# a set number of permutations at a time, to then be passed to workers. Pre-loading more 
# sets at a time speeds up computation, but consumes more memory.
###

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
    tests <- append(tests,list(raw_counts[sample(1:nrow(raw_counts),bincount,replace=F),]))
  }
  intermediate_res <- foreach(i=tests, .combine = 'rbind',.packages = 'entropy',.noexport = c('raw_counts','tests')) %dopar% {
    myfun(i, bincount, target_counts, isAllowed)
  }
  results <- rbind(results, intermediate_res)
}

end=Sys.time()
timing = end-start
timing
stopCluster(cl)
    
results <- results[1:trials,1:4]
colnames(results) <- c('Raw_median','Raw_IQR','Filtered_median','Filtered_IQR')
rownames(results) <- c(1:trials)

write.table(results,file='') #Output file
write.table(timing,file='') #Run duration output file; optional
