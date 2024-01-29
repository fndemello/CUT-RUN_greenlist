# Load required libraries
library(entropy)
library(doParallel)
library(foreach)
options(digits = 15)

# Input data
norm_counts <- read.delim(file='')
# Assumes bin count input as described in 'norm+entropy.R'


#####################
trueSet <- apply(norm_counts,1,function(x) entropy(x,method='ML'))
trueList <- names(trueSet[order(trueSet, decreasing = T)][1:ceiling(nrow(norm_counts)*0.001)])
# These are the sets against which similarity and correlation will be compared

# Main test function
myfun = function(test_counts, fraction, trueSet, trueList){
  ent = apply(test_counts,1,function(x) entropy(x,method='ML'))
  new_set = names(ent[order(ent, decreasing=T)][1:length(trueList)])
  gc()
  return(c(fraction,
           (length(intersect(new_set,trueList))/length(unique(append(new_set,trueList)))),
           cor(ent,trueSet,method='s',use='complete.obs')))
  rm()
}

# Monte-Carlo set up
clusters <- 14
trials_fun = function(fraction) {
  removed_percent <- (1-fraction)*100
  return(round(2*removed_percent))
}
# This functions controls how many permetutation will be done at each sampling level; since
# we expect large samplings (eg. 95%) to give fairly consistent results and to also cost a lot
# of memory to process, this function is inversely proportional to the sampling percentage (so
# eg. a 10% sampling will be repeated more times than a 90% sampling). The scale of how many
# tests to run based on this inverse proportion depends on computational resources available
# (we arbitrarily chose a 2x scaling).
              
steps = 0.02
# Controls the sampling intervals

###
# Make multithreading cluster
registerDoParallel(clusters)
              
set.seed(0)
mc.reset.stream()
preload <- clusters*10
# To avoid passing the entire counts table to each worker thread, the script pre-samples 
# up to a set number of permutations at a time, to then be passed to workers. Pre-loading
# more sets at a time speeds up computation, but consumes more memory.

# MC run
results <- NULL
start=Sys.time()

p = 1-steps
while (p > 0) {
  tests <- NULL
  intermediate_res <- NULL
  trials = trials_fun(p)
 # if (p == 0.50){
 #   registerDoParallel(24)
 # }
     # Optional QoL, to increase the amount of clusters halfway through; RAM usage, 
     # usually the limiting factor, tends to fall dramatically as the script goes on.
 # write(paste('Current: p =',p*100,',',trials,'trials. Start:',start,', update:',Sys.time()),file='',append = T)
     # Optional QoL, to export updates to disk when running from the CLI.
  while (trials > preload) {
    for (t in 1:preload) {
      tests <- append(tests,list(norm_counts[,sample(1:ncol(norm_counts),ceiling(ncol(norm_counts)*p),replace=F)]))
    }
    intermediate_res <- foreach(i=tests, .combine = 'rbind',.packages='entropy',.noexport = c('norm_counts','tests')) %dopar% {
      myfun(i, p, trueSet, trueList) }
    results <- rbind(results, intermediate_res)
    tests <- NULL
    intermediate_res <- NULL
    trials <- trials-preload
  }
  for (t in 1:trials) {
    tests <- append(tests,list(norm_counts[,sample(1:ncol(norm_counts),(ncol(norm_counts)*p),replace=F)]))
  }
  intermediate_res <- foreach(i=tests, .combine = 'rbind',.packages='entropy',.noexport = c('norm_counts','tests')) %dopar% {
    myfun(i, p, trueSet, trueList) }
  results <- rbind(results, intermediate_res)
  tests <- NULL
  intermediate_res <- NULL
  p <- round(p-steps,2)
}

end=Sys.time()
timing = end-start
timing

colnames(results) <- c('Fraction','Jaccard Index','Correlation')


write.table(results,file='')
  # Final output file
write.table(timing,file='')
  # Optional, output of run time, for the sake of optimization

