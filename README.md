# Identifying transgene insertions in Caenorhabditis elegans genomes with Oxford Nanopore sequencing

Paula E. Adams*, John M. Sutton, Jenny L. Thies, Josh D. Millwood, Guy A. Caldwell, Kim A. Caldwell, and Janna L. Fierst

 *Corresponding Author: pea0013@auburn.edu 

## Abstract:
""

## Project Summary:
Here, we sequence and assemble two transgenic strains of Caenorhabditis elegans used in neurodegeneration research, BY250 and UA44. We use a combination of empirically generated ONT long-read DNA sequences, publicly available Illumina short-read DNA libraries, simulated Illumina transgene libraries and reference-based scaffolding to assemble high-quality genome sequences. Our method achieves chromosome-level assemblies with high similarity to the C. elegans N2 reference genome. We identify the locations of the transgene insertions and confirm that all transgene sequences were inserted in intergenic regions, leaving the organismal gene content intact. Our work demonstrates that long-read sequencing is a fast, cost-effective way to assemble genome sequences and aid in the identification of transgenic insertions in model organism research

## Data
### UA44 BioProject:PRJNA627736
#### SRA - Oxford Nanopore Reads
Accession: SRX8359988, SRX8359989, SRX8359990, SRX8359991
#### Assembled Genome
Accession: SAMN14683044
### BY250 BioProject:PRJNA627737 
#### SRA - Oxford Nanopore Reads
Accession: SRX8342585
#### Assembled Genome
Accession: SAMN14682971
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

### Polish genomes with Pilon: 

### Decontaminate with Blast

### Scaffold with RagTag: 
