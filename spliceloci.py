# LT 9/08/2019

# hail script to generate splice loci covered by gnomad

#


# load gencode data
gencode = hl.experimental.import_gtf('../splicescore/data/gencode/gencode.v31lift37.annotation.gtf')

# check list of contigs in gencode
#contigs = gencode.aggregate(hl.agg.collect_as_set(gencode.interval.start.seqname))

# list of contigs:
contigs = {'GL000192.1',
 'GL000193.1',
 'GL000195.1',
 'GL000199.1',
 'GL000202.1',
 'GL000204.1',
 'GL000205.1',
 'GL000212.1',
 'GL000220.1',
 'GL000228.1',
 'GL000237.1',
 'GL000241.1',
 'chr1',
 'chr10',
 'chr11',
 'chr12',
 'chr13',
 'chr14',
 'chr15',
 'chr16',
 'chr17',
 'chr18',
 'chr19',
 'chr2',
 'chr20',
 'chr21',
 'chr22',
 'chr3',
 'chr4',
 'chr5',
 'chr6',
 'chr7',
 'chr8',
 'chr9',
 'chrM',
 'chrX',
 'chrY'}
# keep only contigs that are chromosomes (start with 'chr') and remove mitochondrial DNA
chrom = {'chr1',
 'chr10',
 'chr11',
 'chr12',
 'chr13',
 'chr14',
 'chr15',
 'chr16',
 'chr17',
 'chr18',
 'chr19',
 'chr2',
 'chr20',
 'chr21',
 'chr22',
 'chr3',
 'chr4',
 'chr5',
 'chr6',
 'chr7',
 'chr8',
 'chr9',
 'chrX',
 'chrY'}

# convert to hail expression
chrom = hl.literal(chrom)

# keep only exons on chromosomes 
gencode = gencode.filter(gencode.feature == 'exon')
gencode = gencode.annotate(contig = gencode.interval.start.seqname)
gencode = gencode.filter(chrom.contains(gencode.contig))
#gencode_sense = gencode.filter(gencode.strand == '+')
#gencode_antisense = gencode.filter(gencode.strand =='-')

# build locus intervals for donors (5' side of exon) 
gencode = gencode.annotate(
  acceptor = hl.locus_interval(gencode.contig[3:],
    hl.cond(gencode.strand=='+',
      gencode.interval.start.position - 20,
      gencode.interval.end.position - 2),
    hl.cond(gencode.strand=='+',
      gencode.interval.start.position + 2,
      gencode.interval.end.position + 20),
    includes_start = True,
    includes_end = True),
  donnor = hl.locus_interval(gencode.contig[3:],
    hl.cond(gencode.strand=='+',
      gencode.interval.end.position - 2,
      gencode.interval.start.position - 6),
    hl.cond(gencode.strand=='+',
      gencode.interval.end.position + 6,
      gencode.interval.start.position + 2),
    includes_start = True,
    includes_end = True))
  
accept = gencode.select(splice = gencode.acceptor, gene = gencode.gene_name).key_by('splice')
donor = gencode.select(splice = gencode.donnor, gene = gencode.gene_name).key_by('splice')
splicesites = accept.union(donor)

# I did all this hoping I could intersect with the exome calling regions but I can't !
# just export as a bed file and use bed tools ?
calling_regions = hl.import_locus_intervals('exome_calling_regions.v1.interval_list')

# intersect with gnomad exome calling regions

