# Load required libraries
library(entropy)
library(doParallel)
library(foreach)
options(digits = 15)

setwd("/work2/group/prostata/metodo_normalizacao/")

# Input data
norm_counts <- read.delim(file='9.Validacao/5.CutnTag/4.Calculos/q_norm_allbins.tsv',header=T,row.names=1)
# Assumes bin count input as described in 'norm+entropy.R'




#####################
trueSet <- apply(norm_counts,1,function(x) entropy(x,method='ML'))
trueList <- names(trueSet[order(trueSet, decreasing = T)][1:ceiling(nrow(norm_counts)*0.001)])

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
steps = 0.02


# Make multithreading cluster
#cl <- makeCluster(clusters)
registerDoParallel(clusters)
#clusterExport(cl, varlist = c("trueSet","trueList"))
#clusterEvalQ(cl, {library(entropy)})

set.seed(0)
mc.reset.stream()
preload <- clusters*10
###
# To avoid passing the entire counts table to each worker thread, the script pre-samples 
# up to a set number of permutations at a time, to then be passed to workers. Pre-loading
# more sets at a time speeds up computation, but consumes more memory.
###

# MC run
results <- NULL
start=Sys.time()

p = 1-steps
while (p > 0) {
  tests <- NULL
  intermediate_res <- NULL
  trials = trials_fun(p)
  if (p == 0.50){
    #stopCluster()
    #cl <- makeCluster(12)
    registerDoParallel(24)
    #clusterExport(cl, varlist = c("trueSet","trueList"))
    #clusterEvalQ(cl, {library(entropy)})
  }
  write(paste('Current: p =',p*100,',',trials,'trials. Start:',start-10800,', update:',Sys.time()-10800),file='6.Calculos/MC_2.cnt_progress.txt',append = T)
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
  write.table(results,file='6.Calculos/MC_2.cnt_test.tsv',quote=F,sep='\t', append = F)
}

end=Sys.time()
timing = end-start
timing
#stopCluster()

colnames(results) <- c('Fraction','Jaccard Index','Correlation')


write.table(results,file='6.Calculos/MC_2.cnt_test.tsv',quote=F,sep='\t')
#write.table(timing,file='6.Calculos/MC_2.mmtest_timer.tsv',quote=F,row.names = F,col.names = F)

