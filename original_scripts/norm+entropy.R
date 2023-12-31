# Install/Load required Libraries
pacs <- c("ggplot2", "broman", "entropy", "reshape2")
             
if(sum(as.numeric(!pacs %in% installed.packages())) != 0){
  ins <- pacs[!pacs %in% installed.packages()]
  for(i in 1:length(ins)) {
    install.packages(ins, dependencies = T)
    break()}
  sapply(pacs, require, character = T) 
} else {
  sapply(pacs, require, character = T) 
}

library(broman)
library(entropy)
library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)

# Input bin counts
raw_counts <- read.delim(file=args[1], row.names=1) #Input file
###
# This script assumes that binning and sample quantification has already been previously performed;
# formatting should feature rows as independent bins, and columns as samples. Counts should be
# non-normalized, as quantile normalization is preffered. Removal of fully empty bins (ie. bins
# without counts in any sample) is highly advised. We recommend deeptools' multiBamSummary 
# for this quantification.
###

# Quantile normalize
q_norm = broman::normalize(raw_counts)
colnames(q_norm) = colnames(raw_counts)
rownames(q_norm) = rownames(raw_counts)
write.table(q_norm,
            file=paste(args[2], "-norm-entropy.txt", sep = ""), #Output file
            quote=F,sep='\t')
# Saves the quantile normalized values, if desired. Can be skipped.


# Calculate Shannon entropy for each bin across all samples
entropy_vals = apply(q_norm,1,function(x) entropy(x,method='ML'))
summary(entropy_vals)

                  
entropy_vals = entropy_vals[!is.na(entropy_vals)]
###
# The entropy function will return 'NA' if empty bins are not removed, which affects the percentage-
# based selection. Because of this, the above step removes null-entropy bins before continuing. 
###

                     
# Save files
blist = apply(q_norm,1,function(x) median(x))
bl_1 = blist[blist>=quantile(blist,probs=0.99)]
bl_0.1 = blist[blist>=quantile(blist,probs=0.999)]
gl_1=entropy_vals[entropy_vals>=quantile(entropy_vals,probs=0.99)]
gl_0.1=entropy_vals[entropy_vals>=quantile(entropy_vals,probs=0.999)]
###
# This identification of blacklist regions follows the methodology of the original ENCODE ChIPseq
# blacklist, selecting high-signal regions based on their median normalized value accross samples.
# This script saves both the top 1% and top 0.1% bins (for both greenlist and blacklist), so
# the extension step can be performed as described in the manuscript.
###


write.table(entropy_vals, file = paste (args[2], "-entropy.txt", sep = "")) #Output file
write.table(gl_1, file = paste (args[2], "-entropy-gl-1.txt", sep = "")) #Output file
write.table(gl_0.1, file = paste (args[2], "-entropy-gl-0.1.txt", sep = "")) #Output file
write.table(bl_1,file=paste (args[2], "-entropy-bl-1.txt", sep = "")) #Output file
write.table(bl_0.1,file=paste (args[2], "-entropy-bl-0.1.txt", sep = "")) #Output file


######################################################
# Example plotting (optional)

pdf (paste(args[2], "-rplot.pdf", sep = ""))

# Distribution of entropy values. Refers to figure 2a of our manuscript
tmp <- melt (q_norm)
ggplot(data.frame(tmp),aes(x = value)) +
  geom_density() + 
  theme_minimal()
#

# Calculation of sample correlation versus entropy threshold. Refers to figure 2b of our manuscript
n = 1
steps = 1000
corr_curve = NULL
while ((n/steps)<=1) {
  th = quantile(entropy_vals,probs = (1-(n/steps)),na.rm=T)
  new_ent = entropy_vals[entropy_vals>=th]
  if (length(new_ent)>1) {
    res = cor(data.frame(q_norm[names(new_ent),]),method="pearson",use="pairwise.complete.obs")
    if (is.null(corr_curve)==T){
      corr_curve = data.frame(Percentage = (n/steps)*100,Min_corr = min(res[upper.tri(res)]),Q25_corr = quantile(res[upper.tri(res)],probs=0.25),Med_corr = median(res[upper.tri(res)]),Q75_corr = quantile(res[upper.tri(res)],probs=0.75),Max_corr = max(res[upper.tri(res)]))
    } else {
      corr_curve = rbind(corr_curve,data.frame(Percentage = (n/steps)*100,Min_corr = min(res[upper.tri(res)]),Q25_corr = quantile(res[upper.tri(res)],probs=0.25),Med_corr = median(res[upper.tri(res)]),Q75_corr = quantile(res[upper.tri(res)],probs=0.75),Max_corr = max(res[upper.tri(res)])))
    }
  }
  n = n+1
}

ggplot(corr_curve,aes(Percentage,Med_corr))+
  geom_point(size=0.5) +
  scale_y_continuous(name = 'Median correlation')

dev.off()
#
