# Load required libraries
library(broman)
library(entropy)
library(edgeR)
library(ggplot2)
library(DESeq2)
options(mc.cores=16)

setwd("/work2/group/prostata/metodo_normalizacao/")


# Input data
total_counts <- read.delim(file='5.Quant/clean_over2Msamp_1k-bins.tsv',header=T,row.names=1)
v <- read.delim(file='6.Calculos/new_qnorm_entropies.tsv',quote='')
meta <- read.delim(file='GEO_deseq_meta.tsv',header=T,quote='',row.names=1,col.names=c('Sample','Experiment','Cell line','Cell chars','Ab host','Ab target','Ab ID','CnR protocol','PI','Affiliation','Used AB','Libsize'))

# Process data
ent <- v$x
names(ent) = row.names(v)
target_ents = ent[ent<=quantile(ent,probs = 0.1,na.rm=T)]
target_counts = data.frame(total_counts[names(target_ents),])
raw_target_ents = apply(target_counts,1,function(x) entropy(x,method='ML'))
dsobj = DESeqDataSetFromMatrix(countData = target_counts,
                               colData = meta,
                               design = ~ Experiment)
dsobj = estimateSizeFactors(dsobj)
dsobj = DESeq2::estimateDispersions(dsobj)

tempobj = dsobj
sf = estimateSizeFactorsForMatrix(target_counts)
sizeFactors(tempobj) <- sf
trans_counts = DESeq2::counts(tempobj,normalized=T,replaced=F)
new_ent = apply(trans_counts,1,function(x) entropy(x,method='ML'))

graph_ent = rbind(data.frame(Entropy=raw_target_ents,Group='Non-normalized'),
                  data.frame(Entropy=new_ent,Group='Default DESeq'))
t_range = c(10,5,4,2,1,0.5,
            0.25,0.1,0.05,0.04,0.03,0.02,0.01,0.005,0.001)

for (i in t_range) {
  norm_ents = ent[ent>=quantile(ent,probs = (1-(i/100)),na.rm=T)]
  norm_counts = data.frame(total_counts[names(norm_ents),])
  test_counts = rbind(norm_counts,target_counts)
  sf = estimateSizeFactorsForMatrix(test_counts,controlGenes = 1:nrow(norm_counts))
  tempobj = dsobj
  sizeFactors(tempobj) <- sf
  trans_counts = DESeq2::counts(tempobj,normalized=T,replaced=F)
  new_ent = apply(trans_counts,1,function(x) entropy(x,method='ML'))
  graph_ent = rbind(graph_ent,data.frame(Entropy=new_ent,Group=paste(i,'%',sep='')))
}



graph_ent$Group = factor(graph_ent$Group,levels=c("Non-normalized","Default DESeq","10%","5%","4%","2%","1%","0.5%","0.25%","0.1%","0.05%","0.04%","0.03%","0.02%","0.01%","0.005%","0.001%"))
ggplot(graph_ent, aes(Group, Entropy,fill=Group)) +
  geom_boxplot(outlier.shape = NA,show.legend = F) +
  scale_y_continuous(limits=c(5,6.25)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_fill_manual(values=c('#363636','#A3A3A3',
                             '#440154FF','#481D6FFF','#453581FF','#3D4D8AFF',
                             '#34618DFF','#2B748EFF','#24878EFF','#1F998AFF',
                             '#25AC82FF','#40BC72FF','#67CC5CFF','#97D83FFF',
                             '#CBE11EFF','#FDE725FF','#FDE600FF'))

write.table(graph_ent,file='6.Calculos/threstest_qML.ent_out.tsv',quote=F,sep='\t')

