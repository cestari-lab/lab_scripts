#RNA-seq analysis: couting reads and comparison analysis
#Load required packages for featureCount analysis (package subread).
# If working in a cluster, e.g., beluga.computecanada.ca, load the modules as indicated below:

module load nixpkgs/16.09
module load gcc/7.3.0
module load StdEnv/2020
module load subread/2.0.3

#featureCounts counts mapped reads for genomic features such as genes, exons, etc. 
#The script below will generate the number of read counts for each exon. It will analyze all sorted.bam files.

featureCounts -LO -a GenomeOfReference.gtf -o "counts.txt" -F "GTF" -t "exon" -g "gene_id" --readExtension3 200 --fraction  -T8  $(ls *sorted.bam)

# The script below will combine all counts data into a matrix file for edgeR analysis
grep -v "#" counts.txt | cut -d$'\t' -f1,7- > counts.matrix

#RNA-seq analysis using edgeR:
#Initiate R in Windows, or load R in Linux. If working in a cluster, follow the steps below with module load r:

module load r
R

#Load required libraries
library (limma)
library (edgeR)

#open matrix file from featureCounts:
data <- read.table("counts.matrix", header = TRUE, sep = "", skip = 0, row.names = 1)

#change headings of data:
names(data)<- c("Control-Bio1", "Control-Bio2", "Control-Bio3", "Treatment-Bio1", "Treatment-Bio2", "Treatment-Bio3")

#Create DGEList and statistical design
y <- DGEList(counts=data, group=group)
group <- factor(c(1,1,1,2,2,2))
design <- model.matrix(~group)
y <- calcNormFactors(y)

#Generate plot MSD for sample comparisons
plotMDS.pdf <-plotMDS(y, col=c(rep("black",2), rep("red",2)))
 
#Estimate sample dispersion
y <- estimateDisp(y, design)

#Plot Biological correlate of variation
plotBCV(y)

#Generate statistical analysis using generalized linear model
fit <- glmQLFit(y, design)

# Compare groups treated vs non-treated
qlf.2vs1 <- glmQLFTest(fit, coef=2) 

