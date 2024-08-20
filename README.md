# Identifying transgene insertions in *Caenorhabditis elegans* genomes with Oxford Nanopore sequencing

Paula E. Adams*, Jennifer L. Thies, John M. Sutton, Joshua D. Millwood, Guy A. Caldwell, Kim A. Caldwell, and Janna L. Fierst

 *Corresponding Author: pea0013@auburn.edu 


## Project Summary:
Here, we sequence and assemble two transgenic strains of *Caenorhabditis elegans* used in neurodegeneration research, BY250 and UA44. We use a combination of empirically generated ONT long-read DNA sequences, publicly available Illumina short-read DNA libraries, simulated Illumina transgene libraries and reference-based scaffolding to assemble high-quality genome sequences. Our method achieves chromosome-level assemblies with high similarity to the *C. elegans* N2 reference genome. We identify the locations of the transgene insertions and confirm that all transgene sequences were inserted in intergenic regions, leaving the organismal gene content intact. Our work demonstrates that long-read sequencing is a fast, cost-effective way to assemble genome sequences and aid in the identification of transgenic insertions in model organism research.

## Data Availability 
Assembled Genomes* and Oxford Nanopore Reads
||||
| --- | --- | --- | 
| UA44 | [BioProject:PRJNA627736](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA627736) | [gff file](./UA44.gff.zip) |
| BY250 | [BioProject:PRJNA627737](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA627737) | [gff file](./BY250.gff.zip) |
| [OSF: Genome Assemblies](https://osf.io/myxat/) | | |

*Genomes still processing at NCBI. Currently available at OSF link above

## Analyses

### Trim adaptors with [porechop](https://github.com/rrwick/Porechop): 
```{}
porechop -i fastqs/UA44_alldata.fastq  -o fastqs/UA44_alldata_TRIMMED.fastq  --discard_middle 
porechop -i fastqs/BY250_alldata.fastq -o fastqs/BY250_alldata_TRIMMED.fastq --discard_middle 
```

### Correct reads with [Canu](https://canu.readthedocs.io/en/latest/index.html): 
```{}
canu -p UA44  -d /UA44/canu_runs  genomeSize=100m -nanopore-raw /path/UA44/fastqs/UA44_alldata_TRIMMED.fastq
canu -p BY250 -d /BY250/canu_runs genomeSize=100m -nanopore-raw /path/BY250/fastqs/BY250_alldata_TRIMMED.fastq
```

### Assemble genomes using Canu-corrected reads with [Flye](https://github.com/fenderglass/Flye): 
```{}
flye  -t 48 --genome-size 100m --nano-corr UA44.correctedReads.fasta --out-dir /path/UA44/flye/
flye  -t 48 --genome-size 100m --nano-corr BY250.correctedReads.fasta --out-dir /path/BY250/flye/
```

### Simulate short-read data for transgene insertion using [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm):
```{}
art_illumina -ss HS25 -sam -i UA44_insertion.fasta -l 150 -f 100 -p -m 270 -s 30 -o insertion_lib_UA44
art_illumina -ss HS25 -sam -i Expression_vector_unc68_GFP.fasta -l 150 -f 100 -p -m 270 -s 30 -o insertion_lib_BY250
```
After simulating short-reads for the transgenes, these simulated libraries were added to our downloaded illumina data for N2: [BioProject:PRJDB2670](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJDB2670). 

### Polish genomes with illumina data using [Pilon](https://github.com/broadinstitute/pilon) x4: 
Example for UA44 shown here
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
```
UA44_flye_pilon1.fasta is then used for round 2. Repeat for 4 total rounds of Pilon.


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

### Scaffold with [RagTag](https://github.com/malonge/RagTag): 

RagTag scaffolds assemblies based on a reference genome and can correct missasemblies. We used the "correct" option with the -j flag to protect the contigs containing the transgene insertions from being incorrectly corrected when compared to the reference. 
For UA44 ragtagskip.txt listed contig_59.
For BY250 ragtagskip.txt listed contig_86.
```{}
ragtag.py correct -j ragtagskip.txt Caenorhabditis_elegans.WBcel235.dna.toplevel.fa  UA44_flye_keptcontigs.fasta 
ragtag.py correct -j ragtagskip.txt Caenorhabditis_elegans.WBcel235.dna.toplevel.fa  BY250_flye_keptcontigs.fasta 
```

### Identify insertion location with [minimap2](https://github.com/lh3/minimap2) and BLAST
```{}
## Pre-scaffolded assembly
minimap2  -a UA44_pilon4_keptcontigs.fasta UA44_insertion.fasta > UA44_instertion.sam 
minimap2  -a BY250_pilon4_keptcontigs.fasta Expression_vector_unc68_GFP.fasta > BY250_instertion.sam 

## Scaffolded assembly
minimap2 -L -a ragtag.scaffold.fasta  UA44_insertion.fasta > UA44_ragtag_skip.insertion.sam 
minimap2  -a ragtag.scaffold.fasta  insertion.fasta > BY250_ragtag_insertion.sam

## BLASTn
blastn -num_threads 12  -query UA44_insertion.fasta -subject UA44_genome_final.fasta
blastn -num_threads 12  -query insertion.fasta -subject BY250_genome_final.fasta
```
Figures were then created in Rstudio with [gggenomes](https://github.com/thackl/gggenomes) (an extension of [ggplot2](https://ggplot2.tidyverse.org/)). 

### Genome Quality Assessment [Quast](https://quast.sourceforge.net/) and [Busco](https://busco.ezlab.org/):
Both quality control steps were ran at many stages through analysis. Example script shown here. 
```{}
## Quast
python quast.py -t 12 --plots-format pdf -r Caenorhabditis_elegans.WBcel235.dna.toplevel.fa ragtag.scaffold.fasta -o ./quast_NEW_UA44_ragtag_correctionskip
python quast.py -t 12 --plots-format pdf -r aenorhabditis_elegans.WBcel235.dna.toplevel.fa ragtag.scaffold.fasta   -o ./quast_BY250_NEW_ragtag_correction_skip

## Busco
busco -c 12 -m genome -i ragtag.scaffold.fasta -o busco_NEW_UA44_ragtag_correction_skip --lineage_dataset nematoda_odb10 --config config.ini --update-data
busco -c 12 -m genome -i ragtag.scaffold.fasta -o busco_NEW_BY250_ragtag_correction_skip --lineage_dataset nematoda_odb10 --config config.ini --update-data
```

### Genome Annotation with [Liftoff](https://github.com/agshumate/Liftoff):
```{}
liftoff -g Caenorhabditis_elegans.WBcel235.100.gff3  -o UA44_ragtag_skip.gff  -u UA44_ragtag_skip_unmapped_features.txt   -m minimap2 ragtag.scaffold.fasta  Caenorhabditis_elegans.WBcel235.dna.toplevel.fa
liftoff -g Caenorhabditis_elegans.WBcel235.100.gff3  -o BY250_ragtag_skip.gff -u BY250_ragtag_skip_unmapped_features.txt  -m minimap2 ragtag.scaffold.fasta Caenorhabditis_elegans.WBcel235.dna.toplevel.fa
```

Only 91 (UA44) and 92 (BY250) genes were missing from our genome annotations compared to the reference. Missing gene lists were compared to WormBase Parasite BioMart [Howe et al., 2017](https://pubmed.ncbi.nlm.nih.gov/27899279/) to asess gene function and location. Missing Genes BioMart Information can be found here:  [UA44](./UA44_MissingGenes_BioMart.txt) & [BY250](./BY250_MissingGenes_BioMart.txt).


### Create Bed file for genes using BEDOPS [gff2bed](https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/gff2bed.html): 
Example shown with UA44 sequence data 
```{}
awk '$3=="mRNA" {print}'   UA44_ragtag.scaffolds.gff >   UA44_ragtag.scaffolds_gene.gff

gff2bed < UA44_ragtag.scaffolds_gene.gff > UA44_ragtag.scaffolds_gene.gff.bed

#format it 
cat UA44_ragtag.scaffolds_gene.gff.bed | awk '$8=="mRNA" {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' | sed 's/:/\t/g' | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"0"\t"$7}'  > UA44_ragtag.scaffolds_gene.gff.bed.test

mv UA44_ragtag.scaffolds_gene.gff.bed.test  UA44_scaffold.bed
```

### Extract cds sequences for genes using [AGAT extract_sequences](https://agat.readthedocs.io/en/latest/tools/agat_sp_extract_sequences.html):
```{}
# Type: conda activate agatenv
agat_sp_extract_sequences.pl -g  UA44_ragtag.scaffolds_gene.gff -f  UA44_ragtag.scaffolds.fasta  -t mRNA -o UA44_scaffold.cds.fasta

python -m jcvi.formats.fasta format UA44_scaffold.cds.fasta UA44_scaffold.cds & 

# fix contig names with sed
```


### Synteny Analysis with [MCscan (python version)](https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version)):
Using the 6 largest linkage groups (Chromosomes I, II, III, IV, V, & X), synteny mapping was performed between each genome (UA44 & BY250) with the N2 reference. UA44 and N2 example shown below. 


```{}
python -m jcvi.formats.fasta format UA44.cds.fasta UA44.cds & 
python -m jcvi.formats.fasta format N2.cds.fasta N2.cds & 


python -m jcvi.compara.catalog ortholog N2 UA44_scaffold  --no_strip_names & 

python -m jcvi.compara.synteny depth --histogram  N2.UA44_scaffold.anchors 

#make seqids and layout files 
python -m jcvi.compara.synteny screen --minspan=0 --simple N2.UA44_scaffold.anchors N2.UA44_scaffold.anchors.new

python -m jcvi.graphics.karyotype seqids layout
```



