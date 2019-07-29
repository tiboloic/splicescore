# LT 29/07/2019
#
# Hail script to 

import hail as hl;

#load gnomAD 2.1.1 exome h19 data
exons = hl.import_vcf("../data/gnomad.exomes.r2.1.1/gnomad.exomes.r2.1.1.sites.vcf.bgz")

#tv=cv.rows()

