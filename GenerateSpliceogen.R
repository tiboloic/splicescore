# LT 29/07/2019

# clean up sliceogen files by removing "." (replace by NA)
#

# load per chromosome files output by spliceogen
splice=c()
for (chr in c(as.character(1:22), 'X', 'Y')) {
  filename = paste("../data/spliceogen/gnomad", chr, "vcf_out.txt", sep=".")
  onesplice = read.delim(filename) #, na.strings=c("."))
  splice = rbind(splice, onesplice)
}
#save(splice, file="../data/spliceogen/splice.Rdata")

# save as vcf format
# first create a header

header = "##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=AL,Number=1,Type=Float,Description=\"Probability of being acceptor loss\">
##INFO=<ID=DL,Number=1,Type=Float,Description=\"Probability of being donor loss\">
##INFO=<ID=DG,Number=1,Type=Float,Description=\"Probability of being acceptor gain\">
##INFO=<ID=AG,Number=1,Type=Float,Description=\"Probability of being donor gain\">
##INFO=<ID=GENE,Number=.,Type=String,Description=\"Gene list out of spliceogen\">
##INFO=<ID=withinSite,Number=.,Type=String,Description=\"Splice site\">
##contig=<ID=1,length=249250621,assembly=gnomAD_GRCh37>
##contig=<ID=2,length=243199373,assembly=gnomAD_GRCh37>
##contig=<ID=3,length=198022430,assembly=gnomAD_GRCh37>
##contig=<ID=4,length=191154276,assembly=gnomAD_GRCh37>
##contig=<ID=5,length=180915260,assembly=gnomAD_GRCh37>
##contig=<ID=6,length=171115067,assembly=gnomAD_GRCh37>
##contig=<ID=7,length=159138663,assembly=gnomAD_GRCh37>
##contig=<ID=8,length=146364022,assembly=gnomAD_GRCh37>
##contig=<ID=9,length=141213431,assembly=gnomAD_GRCh37>
##contig=<ID=10,length=135534747,assembly=gnomAD_GRCh37>
##contig=<ID=11,length=135006516,assembly=gnomAD_GRCh37>
##contig=<ID=12,length=133851895,assembly=gnomAD_GRCh37>
##contig=<ID=13,length=115169878,assembly=gnomAD_GRCh37>
##contig=<ID=14,length=107349540,assembly=gnomAD_GRCh37>
##contig=<ID=15,length=102531392,assembly=gnomAD_GRCh37>
##contig=<ID=16,length=90354753,assembly=gnomAD_GRCh37>
##contig=<ID=17,length=81195210,assembly=gnomAD_GRCh37>
##contig=<ID=18,length=78077248,assembly=gnomAD_GRCh37>
##contig=<ID=19,length=59128983,assembly=gnomAD_GRCh37>
##contig=<ID=20,length=63025520,assembly=gnomAD_GRCh37>
##contig=<ID=21,length=48129895,assembly=gnomAD_GRCh37>
##contig=<ID=22,length=51304566,assembly=gnomAD_GRCh37>
##contig=<ID=X,length=155270560,assembly=gnomAD_GRCh37>
##contig=<ID=Y,length=59373566,assembly=gnomAD_GRCh37>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"

# create vcf file with header
fileConn = file("../data/spliceogen/spliceogen.vcf")
writeLines(header, fileConn)
close(fileConn)

# manipulations to get to vcf format
splice$CHROM = splice$X.CHR
splice$POS = splice$START
splice$ID = "."
splice$QUAL = "."
splice$FILTER="PASS"
# replace ; by , when GENE is a list
splice$GENE2 = sub(";", ",",splice$GENE)
splice$INFO = paste("AL=",splice$accLoss,
                    ";DL=", splice$donLoss,
                    ";DG=", splice$donGain,
                    ";AG=", splice$accGain,
                    ";GENE=", splice$GENE2,
                    ";withinSite=", splice$withinSite, sep ="")
vcfnames=c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
write.table(splice[,vcfnames], file="../data/spliceogen/spliceogen.vcf", na = ".", sep ="\t", quote=FALSE,
            row.names=FALSE, col.names=FALSE, append = TRUE)
