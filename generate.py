# LT 29/07/2019
#
# Hail script to count potential splice variants per gene and sum up spliceogen probabilities

import hail as hl;

#load gnomAD 2.1.1 exome h19 data
# use filtered version that has only biallelic sites, SNVs, protein_coding transcripts and trimmed INFO
exons = hl.import_vcf("../data/gnomad.exomes.r2.1.1/filtered.vcf.bgz")
exons = exons.rows()

#load Steve's data
splice = hl.import_vcf("../data/spliceogen/spliceogen.vcf")
splice = splice.rows()

# keep only variants that are in Steve's output
# removes single exon genes
#FIXME should use join() instead (inner)
exons = exons.filter(hl.is_defined(splice[exons.key]))

# annotate gnomAD data with splice probabilities + flag indicating creation (not in splice site)
# add a flag to common variants (AF > xxx) FIXME: why 0.0001 ?
# gene is from steve's output. sometimes is an array with multiple genes
# gene2 is taken from first vep annotation (garantees unicity)
ann_expr = {
  'al': splice[exons.key].info.AL,
  'dl': splice[exons.key].info.DL,
  'dg': splice[exons.key].info.DG,
  'ag': splice[exons.key].info.AG,
  'create': hl.is_missing(splice[exons.key].info.withinSite),
  'common': exons.info.AF[0] > 0.0001,
  'gene': splice[exons.key].info.GENE,
  'gene2': exons.info.vep[0].split('\\|')[3]
}
exons = exons.annotate(**ann_expr)

# splice site creation summary
# aggregate by gene, count all variants, sum up probabilities, sum up probabilities of common variants

agg_expr = {
  'dg': hl.agg.sum(exons.dg),
  'ag': hl.agg.sum(exons.ag),
  'cdg': hl.agg.sum(hl.cond(exons.common, exons.dg, 0)),
  'cag': hl.agg.sum(hl.cond(exons.common, exons.ag, 0)),
  'dl': hl.agg.sum(exons.dl),
  'al': hl.agg.sum(exons.al),
  'cdl': hl.agg.sum(hl.cond(exons.common, exons.dl, 0)),
  'cal': hl.agg.sum(hl.cond(exons.common, exons.al, 0)),  
  'gain': hl.agg.count_where(exons.create),
  'loss': hl.agg.count_where(exons.create==False)
}

output = exons.group_by(exons.gene2).aggregate(**agg_expr)
output.export('splicecount.tsv')

 
