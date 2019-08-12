# LT 06/06/2019
#
# Look at correlation between SIS (splice intolerance scores) and #of transcripts
# using downloaded https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz
#

cons.file = gzfile("../data/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz", 'rt')
cons = read.delim(cons.file, header=T)

library(dplyr)

n_trans = cons %>% group_by(gene) %>% count(gene, name="n_trans")
spl = inner_join(spl, n_trans, by=c("gene2"="gene"))
plot(gintos.cmon ~ n_trans, data = spl)
cor.test(spl$gintos.cmon,spl$n_trans, method = "spearman")
# negatively correlated, significant
# for gintos.cmon and all, lintos.cmon and all

# might be outliers
