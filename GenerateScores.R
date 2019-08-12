# LT 30/07/2019

# analyse output of spliceogen counts
#

spl = read.delim("splicecount.tsv")

# biplot for creation
plot(dg+ag ~ gain, data=spl)

#looks good for regression

# common
plot(cdg+cag ~ gain, data=spl)
# doesnt look as good

# disruption
plot(dl+al ~ loss, data=spl)
# looks really good: suspicious

plot(cdl+cal ~ loss, data=spl)
# noisier but maube better ?


####
# splice GAIN intolerance score
####

# on common
ginto.cmon = lm(cdg+cag ~ gain, data=spl)
summary(ginto.cmon)

ginto.score = studres(ginto.cmon)
spl$gintos.cmon = ginto.score

# on all
ginto.all = lm(dg+ag ~ gain, data=spl)
spl$gintos.all = studres(ginto.all)

####
# splice disrupt

linto.cmon = lm(cdl+cal ~ loss, data=spl)
summary(linto.cmon)
spl$lintos.cmon = studres(linto.cmon)

linto.all = lm(dl+al ~ loss, data=spl)
spl$lintos.all = studres(linto.all)



##########
# compare with loeuf

metrics = read.delim("../../gnomad/metrics.tsv")

res = merge(spl, metrics, by.x="gene2", by.y = "gene")
res2 = inner_join(spl, metrics, by = c("gene2"="gene"))
plot(gintos.all ~ oe_lof, data = res)
plot(gintos.all ~ oe_lof_upper, data = res)
plot(gintos.cmon ~ oe_lof, data = res)
plot(gintos.cmon ~ oe_lof_upper, data = res)
plot(lintos.all ~ oe_lof, data = res)
plot(lintos.all ~ oe_lof_upper, data = res)
plot(lintos.cmon ~ oe_lof, data = res)
plot(lintos.cmon ~ oe_lof_upper, data = res)
cor.test(res$lintos.all, res$oe_lof)

##########
# check bias with length and number of exons

plot(gintos.cmon~log(gene_length), data = res)
plot(gintos.all~log(gene_length), data = res)
plot(lintos.all~num_coding_exons, data = res)
plot(lintos.all~log(num_coding_exons), data = res)
plot(lintos.cmon~log(num_coding_exons), data = res)

#######
# load list of genes:
# 
haplo = read.delim("../data/genelists/clingen_level3_genes_2018_09_13.tsv", header=FALSE)
names(haplo) = "gene"
haplo$haplo = TRUE

res = left_join(res, haplo)
res$haplo[is.na(res$haplo)] = FALSE

plot(lintos.all ~ as.factor(haplo), data=res)
# does not look good

# olfactory
olfa = read.delim("../data/genelists/olfactory_receptors.tsv", header=FALSE)
names(olfa) = "gene"
olfa$olfa = TRUE

res = left_join(res, olfa, by=c("gene2"="gene"))
res$olfa[is.na(res$olfa)] = FALSE

plot(lintos.all ~ as.factor(olfa), data=res)
