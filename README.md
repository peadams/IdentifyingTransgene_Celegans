# Identifying transgene insertions in *Caenorhabditis elegans* genomes with Oxford Nanopore sequencing

Paula E. Adams*, John M. Sutton, Jenny L. Thies, Josh D. Millwood, Guy A. Caldwell, Kim A. Caldwell, and Janna L. Fierst

 *Corresponding Author: pea0013@auburn.edu 


## Project Summary:
Here, we sequence and assemble two transgenic strains of *Caenorhabditis elegans* used in neurodegeneration research, BY250 and UA44. We use a combination of empirically generated ONT long-read DNA sequences, publicly available Illumina short-read DNA libraries, simulated Illumina transgene libraries and reference-based scaffolding to assemble high-quality genome sequences. Our method achieves chromosome-level assemblies with high similarity to the *C. elegans* N2 reference genome. We identify the locations of the transgene insertions and confirm that all transgene sequences were inserted in intergenic regions, leaving the organismal gene content intact. Our work demonstrates that long-read sequencing is a fast, cost-effective way to assemble genome sequences and aid in the identification of transgenic insertions in model organism research

## Data
#### UA44 
BioProject:PRJNA627736
#### BY250 
BioProject:PRJNA627737 
### Illumina data used
BioProject:PRJDB2670, Accession: DRX007632


## Analyses

### Trim adaptors with porechop: 
https://github.com/rrwick/Porechop
```{}
porechop -i fastqs/UA44_alldata.fastq  -o fastqs/UA44_alldata_TRIMMED.fastq  --discard_middle 
porechop -i fastqs/BY250_alldata.fastq -o fastqs/BY250_alldata_TRIMMED.fastq --discard_middle 
```

### Correct reads with Canu: 
https://canu.readthedocs.io/en/latest/index.html
```{}
canu -p UA44  -d /UA44/canu_runs  genomeSize=100m -nanopore-raw /path/UA44/fastqs/UA44_alldata_TRIMMED.fastq
canu -p BY250 -d /BY250/canu_runs genomeSize=100m -nanopore-raw /path/BY250/fastqs/BY250_alldata_TRIMMED.fastq
```

### Assemble genomes using Canu-corrected reads with Flye: 
https://github.com/fenderglass/Flye
```{}
flye  -t 48 --genome-size 100m --nano-corr UA44.correctedReads.fasta --out-dir /path/UA44/flye/
flye  -t 48 --genome-size 100m --nano-corr BY250.correctedReads.fasta --out-dir /path/BY250/flye/
```

### Simulate short-read data for transgene insertion using ART:
https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm
```{}
art_illumina -ss HS25 -sam -i UA44_insertion.fasta -l 150 -f 100 -p -m 270 -s 30 -o insertion_lib_UA44
art_illumina -ss HS25 -sam -i Expression_vector_unc68_GFP.fasta -l 150 -f 100 -p -m 270 -s 30 -o insertion_lib_BY250
```
After simulating short-reads for the transgenes, these simulated libraries were added to our downloaded illumina data for N2: BioProject:PRJDB2670. 

### Polish genomes with illumina data using Pilon x4: 
https://github.com/broadinstitute/pilon
```{}
##ROUND1
#index genome
bwa index ./flye/assembly.fasta
#align reads to reference    
bwa mem -M -t 48 ./flye/assembly.fasta ./illumina_fastq/DR_insertion1.fastq  ./illumina_fastq/DR_insertion2.fastq > ./pilon/UA44_flye_bwa.sam
#sam to bam
samtools view -Sb ./pilon/UA44_flye_bwa.sam  > ./pilon/UA44_flye_bwa.bam

#Sort and index the BAM
samtools sort ./pilon/UA44_flye_bwa.bam -o ./pilon/UA44_flye_bwa.sort.bam
samtools index ./pilon/UA44_flye_bwa.sort.bam

#run pilon
java -Xmx300G -jar pilon-1.23-0/pilon-1.23.jar --genome ./flye/assembly.fasta --frags  ./pilon/UA44_flye_bwa.sort.bam --output ./pilon/UA44_flye_pilon1

##ROUND2
#index genome
bwa index ./pilon/UA44_flye_pilon1.fasta
#align reads to reference    
bwa mem -M -t 48  ./pilon/UA44_flye_pilon1.fasta ./illumina_fastq/DR_insertion1.fastq  ./illumina_fastq/DR_insertion2.fastq > ./pilon/UA44_flye_bwa.sam
#sam to bam
samtools view -Sb ./pilon/UA44_flye_bwa.sam  > ./pilon/UA44_flye_bwa.bam

#Sort and index the BAM
samtools sort ./pilon/UA44_flye_bwa.bam -o ./pilon/UA44_flye_bwa.sort.bam
samtools index ./pilon/UA44_flye_bwa.sort.bam

#run pilon
java -Xmx300G -jar pilon-1.23-0/pilon-1.23.jar --genome ./pilon/UA44_flye_pilon1.fasta  --frags  ./pilon/UA44_flye_bwa.sort.bam --output ./pilon/UA44_flye_pilon2

##ROUND3
#index genome
bwa index ./pilon/UA44_flye_pilon2.fasta
#align reads to reference    
bwa mem -M -t 48  ./pilon/UA44_flye_pilon2.fasta ./illumina_fastq/DR_insertion1.fastq  ./illumina_fastq/DR_insertion2.fastq > ./pilon/UA44_flye_bwa.sam
#sam to bam
samtools view -Sb ./pilon/UA44_flye_bwa.sam  > ./pilon/UA44_flye_bwa.bam

#Sort and index the BAM
samtools sort ./pilon/UA44_flye_bwa.bam -o ./pilon/UA44_flye_bwa.sort.bam
samtools index ./pilon/UA44_flye_bwa.sort.bam

#run pilon
java -Xmx300G -jar pilon-1.23-0/pilon-1.23.jar --genome ./pilon/UA44_flye_pilon2.fasta  --frags  ./pilon/UA44_flye_bwa.sort.bam --output ./pilon/UA44_flye_pilon3


##ROUND4
#index genome
bwa index ./pilon/UA44_flye_pilon3.fasta
#align reads to reference    
bwa mem -M -t 48  ./pilon/UA44_flye_pilon3.fasta ./illumina_fastq/DR_insertion1.fastq  ./illumina_fastq/DR_insertion2.fastq > ./pilon/UA44_flye_bwa.sam
#sam to bam
samtools view -Sb ./pilon/UA44_flye_bwa.sam  > ./pilon/UA44_flye_bwa.bam

#Sort and index the BAM
samtools sort ./pilon/UA44_flye_bwa.bam -o ./pilon/UA44_flye_bwa.sort.bam
samtools index ./pilon/UA44_flye_bwa.sort.bam

#run pilon
java -Xmx300G -jar pilon-1.23.jar --genome ./pilon/UA44_flye_pilon3.fasta  --frags  ./pilon/UA44_flye_bwa.sort.bam --output ./pilon/UA44_flye_pilon4
```

### Decontaminate with BLAST:
```{}
blastn \
-task megablast \
-query UA44.contigs.fasta \
-db nt \
-outfmt '6 qseqid qlen staxids bitscore std sscinames sskingdoms stitle' \
-num_threads 8 \
-evalue 1e-25 \
-max_target_seqs 2 \
-out blast_UA44total.out
```
After BLAST, remove contigs that did not align to *Caenorhabditis*

### Scaffold with RagTag: 
