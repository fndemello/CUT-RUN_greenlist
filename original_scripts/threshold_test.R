# Load required libraries
library(broman)
library(entropy)
library(ggplot2)
library(DESeq2)



# Input data
raw_counts <- read.delim(file='')
# Assumes bin count input as described in 'norm+entropy.R'

entropy_vals <- read.delim(file='')
# Assumes input as named vector of bin entropies, as calculated in 'norm+entropy.R'. Vector
# names must match row names of 'raw_counts'

meta <- read.delim(file='')
# DESeq requires a metadata table for creating the DESeq object, but the attributes present and
# design defined are not relevant for this script.

####
# Process data
target_entropy_vals = entropy_vals[entropy_vals<=quantile(entropy_vals,probs = 0.1,na.rm=T)]
target_counts = data.frame(raw_counts[names(target_entropy_vals),])
raw_target_entropy_vals = apply(target_counts,1,function(x) entropy(x,method='ML'))
dsobj = DESeqDataSetFromMatrix(countData = target_counts,
                               colData = meta,
                               design = ~ Experiment)
dsobj = estimateSizeFactors(dsobj)
dsobj = DESeq2::estimateDispersions(dsobj)
###
# For this script, 'target_counts' refers to the bins to which normalization testing is done upon;
# we selected the bottom 10% lowest entropy bins (ie. highest variability) to avoid using the entire
# dataset for testing. Additionally, 'norm_counts' refers to the counts used for normalizing in each 
# test, and 'trans_counts'/'new_entropy' refers to the results of normalizaition tests.
###
                                
tempobj = dsobj
sf = estimateSizeFactorsForMatrix(target_counts)
sizeFactors(tempobj) <- sf
trans_counts = DESeq2::counts(tempobj,normalized=T,replaced=F)
new_entropy_vals = apply(trans_counts,1,function(x) entropy(x,method='ML'))

graph_entropy_vals = rbind(data.frame(entropy=raw_target_entropy_vals,Group='Non-normalized'),
                  data.frame(entropy=new_entropy_vals,Group='Default DESeq'))
t_range = c(10,5,4,2,1,0.5,
            0.25,0.1,0.05,0.04,0.03,0.02,0.01,0.005,0.001)
# Thresholds selected for testing, in percentage points.

for (i in t_range) {
  norm_entropy_vals = entropy_vals[entropy_vals>=quantile(entropy_vals,probs = (1-(i/100)),na.rm=T)]
  norm_counts = data.frame(raw_counts[names(norm_entropy_vals),])
  test_counts = rbind(norm_counts,target_counts)
  sf = estimateSizeFactorsForMatrix(test_counts,controlGenes = 1:nrow(norm_counts))
  tempobj = dsobj
  sizeFactors(tempobj) <- sf
  trans_counts = DESeq2::counts(tempobj,normalized=T,replaced=F)
  new_entropy_vals = apply(trans_counts,1,function(x) entropy(x,method='ML'))
  graph_entropy_vals = rbind(graph_entropy_vals,data.frame(entropy=new_entropy_vals,Group=paste(i,'%',sep='')))
}


######################################################
# Example plotting (optional)
                           
graph_entropy_vals$Group = factor(graph_entropy_vals$Group,levels=c("Non-normalized","Default DESeq","10%","5%","4%","2%","1%","0.5%","0.25%","0.1%","0.05%","0.04%","0.03%","0.02%","0.01%","0.005%","0.001%"))
ggplot(graph_entropy_vals, aes(Group, entropy,fill=Group)) +
  geom_boxplot(outlier.shape = NA,show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_fill_manual(values=c('#363636','#A3A3A3',
                             '#440154FF','#481D6FFF','#453581FF','#3D4D8AFF',
                             '#34618DFF','#2B748EFF','#24878EFF','#1F998AFF',
                             '#25AC82FF','#40BC72FF','#67CC5CFF','#97D83FFF',
                             '#CBE11EFF','#FDE725FF','#FDE600FF'))

write.table(graph_entropy_vals,file='') #Output file

