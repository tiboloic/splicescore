# LT 29/07/2019
#
# Hail script: read gnomad constraint metrics and export as tsv

import hail as hl;

lof = hl.read_table("gnomad.v2.1.1.lof_metrics.by_transcript.ht")
lof = lof.filter(lof.canonical)
lof = lof.select(lof.gene, lof.transcript, lof.gene_length, lof.cds_length, lof.interval,
  lof.num_coding_exons, lof.oe_lof, lof.oe_lof_upper)
lof.export('metrics.tsv')  
