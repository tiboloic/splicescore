# LT 29/07/2019
# filter gnomAD exomes before looking at splice scores

cd ../data/gnomad.exomes.r2.1.1/

# filter only sites with ...
bcftools view gnomad.exomes.r2.1.1.sites.vcf.bgz -m 2 -M 2 -v snps -f PASS -i 'INFO/vep[0] ~ "protein_coding"' -Ou| bcftools annotate -x '^INFO/AC,INFO/AF,INFO/AN,INFO/dp_hist_all_bin_freq,INFO/AN_raw,INFO/vep' -o "filtered.vcf.bgz" -Oz --threads 8

# create index file (csi)
bcftools index filtered.vcf.bgz