## this filed used to retrieve a new dataset as reference data set for simulating null1.1 and alt1 models
# BiocManager::install("MMUPHin")
library(MMUPHin)
data("CRC_abd", "CRC_meta")
# Find the minimum sequencing depth
N <- min(CRC_meta$number_reads)
# Create count data by multiplying relative abundances with the sequencing depth.
# Note this is the same as what SPIEC-EASI demonstrated in their package. 
counts = sweep(CRC_abd, 2, N, "*")
counts = round(counts)
counts[1:6, 1:5]
dim(counts) # p by n; 484 taxa, 551 samples
mean(counts==0) # around 0.8, high sparsity
# The count data can then be used to learn marginal distributions.

# You can create subsets with fewer species by filtering with respect to the prevalence.
## I think might make sense to random sample the species, in order to keep the sparsity level of the original data set (if only keep high prevalence, will have low sparsity)
## will leave the subsetting operation in the main file

prevdf = apply(X = CRC_abd,
               MARGIN = 1,
               FUN = function(x){sum(x > 0)})
hist(prevdf)
keepTaxa = names(prevdf)[(prevdf >= quantile(prevdf, (484-206)/484)) & (prevdf < quantile(prevdf, (484-6)/484))]


CRC_abd_sub <- CRC_abd[rownames(CRC_abd) %in% keepTaxa,]
CRC_abd_sub <- sweep(CRC_abd_sub, 2, colSums(CRC_abd_sub), "/")# Re-normalize

mean(CRC_abd_sub==0)

sub_counts = counts[rownames(counts) %in% keepTaxa,]
hist(colSums(sub_counts))
mean(sub_counts==0)
hist(sapply(1:ncol(sub_counts), function(i)mean(sub_counts[,i]==0)))

save('sub_counts', file='\\\\fs2-vip/students/yuek/Desktop/micro_net/Kun_code/reference_data.RData')
