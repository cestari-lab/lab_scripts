#ChIP-seq analysis and peak calling:
#For ChIP-seq use the package DeepTools.
# If working in a cluster, e.g., beluga.computecanada.ca, load the module as indicated below:

module load python/3.8.10
source ~/ENV/bin/activate

#Create folders for the analysis output:

mkdir -p /path/to/files/bamcomp
mkdir -p /path/to/files/fingerprint
mkdir -p /path/to/files/multibams
mkdir -p /path/to/files/plotcorr
mkdir -p /path/to/files/macs

#Run bamcompare to compare ChIP vs Input samples; input files are aligned reads, i.e., sorted.bam.
#Example: reads-map_v1_sorted.bam (from minimap/samtools analysis); keep the index sorted.bam.bai files in the same folder!
#Run each biological replicate data separatelly, then run again with replicates to make sure each replicate is reproducible.
#The bigwig file can be opened in IGV for visualization: https://igv.org/
#It is also useful to repeat the analysis below outputing "bedgraph" format.

bamCompare -b1 ChIP_Bio1_sorted.bam -b2 Input_Bio1_sorted.bam \
 --outFileName /path/to/files/bamcomp/chip-inp-bpm.bw \
 --normalizeUsing BPM --extendReads 720 \
 --scaleFactorsMethod None \
 --binSize 50 --smoothLength 60 --centerReads \
 --outFileFormat bigwig --numberOfProcessors 10

#Check overall enrichment of ChIP with plotFingerprint:
 
 plotFingerprint --bamfiles ChIP_Bio1_sorted.bam -b2 Input_Bio1_sorted.bam \
 --labels chip_Bio1 Inp_Bio1 \
 --extendReads 600  --binSize=1000 \
 --plotFile /path/to/files/fingerprint/ChIPvsInput_Bio1.fingerprint.pdf \
 --outRawCounts /path/to/files/fingerprint/ChIPvsInput_Bio1_coverage.tab \
 --numberOfProcessors 10

# Compare biological replicates using multiBamSummary and plotCorrelation analysis:

multiBamSummary bins --bamfiles Input_Bio1_sorted.bam ChIP_Bio1_sorted.bam Input_Bio2_sorted.bam ChIP_Bio2_sorted.bam Input_Bio3_sorted.bam ChIP_Bio3_sorted.bam \
 --outFileName /path/to/files/multibams/multiBamArray_Inp-vs-ChIP-Bio1-3.npz --binSize=5000 \
 --extendReads=1000 --labels InputBio1 ChIP-Bio1 InputBio2 ChIP-Bio2 InputBio3 ChIP-Bio3 \
 --numberOfProcessors 10 \

plotCorrelation --corData path/files/multiBamArray_Inp-vs-ChIP-Bio1-3.npz \
 --plotFile path/to/files/plotcorr/REST_bam_correlation_bin.pdf \
 --outFileCorMatrix path/to/files/plotcorr/corr_matrix_bin.txt \
 --whatToPlot heatmap --corMethod spearman \

#MACS3 analysis for broad peak calling:
#To install macs3 in cluster  

module load python/3.8.2
virtualenv MACS3
source MACS3/bin/activate
pip install macs3 

#To load macs3:
module load python/3.8.2
virtualenv MACS3
source MACS3/bin/activate

#Use bedtools to convert bam to bed files using bedtools (do the same to all replicates):

module load bedtools

bedtools bamtobed -i ChIP_Bio1_sorted.bam > ChIP_Bio1_sorted.bed
bedtools bamtobed -i Input_Bio1_sorted.bam > Input_Bio1_sorted.bed

macs3 callpeak -t ChIP_Bio1_sorted.bed ChIP_Bio2_sorted.bed ChIP_Bio3_sorted.bed \
 -c Input_Bio1_sorted.bed Input_Bio2_sorted.bed Input_Bio3_sorted.bed \
 -f BED -g 35000000 -n chip-macs3 --nomodel --extsize 500 --fix-bimodal \
 --broad --broad-cutoff 0.05 --outdir path/to/files/macs/
 
 