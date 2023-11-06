# Load required libraries
library(broman)
library(entropy)
#library(edgeR)
library(ggplot2)
options(mc.cores=16)


# Gene expression data
gene_data <- read.delim(file='5.CutnTag/3.Quant/nonempty_1k-bins.tsv',header=T,row.names=1)


# Normalize
q_norm = broman::normalize(gene_data)
colnames(q_norm) = colnames(gene_data)
rownames(q_norm) = rownames(gene_data)
q_norm = write.table(q_norm,file='6.Calculos/qnorm_bins.tsv',quote=F,sep='\t')


# Calculate Shannon entropy for each gene across all samples
q_entropy = apply(q_norm,1,function(x) entropy(x,method='ML'))


# Print the calculated entropy values
print(format(sort(q_entropy),scientific=F))

summary(q_entropy)
q_entropy = q_entropy[!is.na(q_entropy)]

# Save files
blist = apply(q_norm,1,function(x) median(x))
bl1 = blist[blist>=quantile(blist,probs=0.99)]
bl0.1 = blist[blist>=quantile(blist,probs=0.999)]
top1=q_entropy[q_entropy>=quantile(q_entropy,probs=0.99)]
top0.1=q_entropy[q_entropy>=quantile(q_entropy,probs=0.999)]


write.table(q_entropy,file='4.Calculos/all_qnorm_ents.tsv',quote=F,sep='\t')
write.table(q_norm,file='4.Calculos/q_norm_allbins.tsv',quote=F,sep='\t')
write.table(top1,file='4.Calculos/top_1p_ents.tsv',quote=F,sep='\t')
write.table(top0.1,file='4.Calculos/top_0.1p_ents.tsv',quote=F,sep='\t')
write.table(bl1,file='4.Calculos/bl_top_1.tsv',quote=F,sep='\t')
write.table(bl0.1,file='4.Calculos/bl_top_0.1.tsv',quote=F,sep='\t')


######################################################
# Plotting

# plot_ent = rbind(data.frame(Entropy=q_entropy,Norm='quant'),data.frame(Entropy=cpm_entropy,Norm='cpm'))
# 
# ggplot(plot_ent,aes(Entropy,Group=Norm,color=Norm)) +
#  geom_density() + 
#  scale_x_continuous(limits=c(5,6.5))
# 
ggplot(data.frame(q_norm),aes(q_norm)) +
  geom_density() +
  scale_x_continuous(limits=c(3,5.5))



##########################################################################
# Default single plot

frac = 1
n = 950
steps = 1000
corr_curve = NULL
while ((n/steps)<=frac) {
  th = quantile(q_entropy,probs = (1-(n/steps)),na.rm=T)
  new_ent = q_entropy[q_entropy>=th]
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


ggplot(corr_curve,aes(Percentage,Med_corr,color=Percentage))+
  geom_point(size=0.5) +
  scale_x_continuous(limits = c(0,0.1)) + 
  scale_y_continuous(name = 'Median correlation')

write.table(corr_curve,file='6.Calculos/ML_entropy_corr_graph.tsv',quote=F,sep='\t')
