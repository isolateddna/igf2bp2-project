############################################################
# IGF2BP2 PROJECT – WGCNA
############################################################

# Author: Vishnu A M
# Lab: CLOCK Lab
# Institute: Indian Institute of Technology, Gandhinagar, India
#
# This commands in this markdown files download raw sequencing files and perfomrs quantification using 
# Salmon

# Input:
# Accession IDs of the SRA files
#
#
# Output: 
# Quant.sf files for all the samples that will be used by Tximport for getting the gene-level counts
#
#





Download the raw data for GSE197726

```
**prefetch --option-file text_file_with_srr_ids --max-size 420000000000000000 -O output/path/folder**
```

Convert the SRA files to Fastq files

```
for sra in $(cat sra_ids.txt); do  
fasterq-dump $sra --split-files -O fastq/  
done
```

Quality control using FastQC

```
#!/bin/bash  
  
mkdir -p path/to/output/folder
  
for file in fastq/*.fastq; do  
fastqc "$file" -o path/to/output/folder  
done
```

Prepare Salmon Index

Reference: https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

Download the mouse Transcript and genome files

```
curl -o https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.transcripts.fa.gz

curl -o https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/GRCm39.primary_assembly.genome.fa.gz
```


Install Salmon

```
conda install --channel bioconda salmon
```

Preparing metadata

```
grep "^>" <(gunzip -c GRCm39.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i -e 's/>//g' decoys.txt
```

```
cat gencode.vM36.transcripts.fa.gz GRCm39.primary_assembly.genome.fa.gz > transcript_and_decoy.fa.gz
```

Salmon Indexing

```
salmon index -t transcript_and_decoy.fa.gz -d decoys.txt -p 30 -i salmon_index --gencode
```


Running mapping on all the samples

```
#!/bin/bash  

  
##############################  
# Paths (edit as needed)  
##############################  
salmon_index="path/to/salmon_index"  
fastq_dir="path/to/fastq_files"  
output_dir="path/to/salmon_quant"  
  
mkdir -p "$output_dir"  
  
##############################  
# Loop through samples  
##############################  
for dir in "${fastq_dir}"/SRR*; do  
  
# Find FASTQ files  
r1_file=$(find "$dir" -name "*_1.fastq")  
r2_file=$(find "$dir" -name "*_2.fastq")  
  
# Sample name  
samp=$(basename "$dir")  
  
echo "Processing sample ${samp}"  
  
salmon quant \  
-i "$salmon_index" \  
-l A \  
-1 "$r1_file" \  
-2 "$r2_file" \  
-p 8 \  
--validateMappings \  
-o "${output_dir}/${samp}_quant"  
  
done
```
