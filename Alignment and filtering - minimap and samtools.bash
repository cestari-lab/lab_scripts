#Aligning data to genome using minimap2
#Files in the *.fasta.tar.gz format must be decompressed using the bash command
tar xvzf yourdata.fastq.tar.gz -C /path/to/files/fastq/

#Create folders for the analysis output:

mkdir -p /path/to/files/sam
mkdir -p /path/to/files/bam
mkdir -p /path/to/files/sortedbam
mkdir -p /path/to/files/mapstats

#Alignment of reads obtained from ONT sequencing to the genome of reference, creating a sam file.
# If working in a cluster, e.g., beluga.computecanada.ca, load the modules as indicated below:

module load minimap2/2.24

minimap2 -ax map-ont -I4G -t16 -2 GenomeOfReference.fasta \
/path/to/files/fastq/*.fastq.gz \
>/path/to/files/samreads-map_bio1.sam

#Use samtools to obtain mapping statistics, to convert files from .sam to .bam files, to filter, sort and index .bam files. 
#Use flags "2048" or "2304" to remove supplementary alignments or secondary and supplementary (keeps primary only), respectively;
#Use MAPQ to filter for mapping quality, e.g., " -q 20" keeps alignments with 99% probability of correct alignment.

module load samtools

samtools flagstat /path/to/files/sam/samreads-map_bio1.sam > /path/to/files/mapstats/samreads-map_bio1-flagstat.txt
samtools stat /path/to/files/sam/samreads-map_bio1.sam > /path/to/files/mapstats/samreads-map_bio1-stat.txt
samtools view -S -b -q 20 -F 2304 /path/to/files/sam/samreads-map_bio1.sam  > /path/to/files/bam/reads-map_bio1.bam
samtools sort /path/to/files/bam/reads-map_bio1.bam > /path/to/files/sortedbam/reads-map_bio1_sorted.bam
samtools index /path/to/files/sortedbam/reads-map_bio1_sorted.bam 

