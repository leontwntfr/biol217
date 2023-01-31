# Day 7

## RNAseq - preparation

### Finding SRR numbers

Data from paper: https://www.tandfonline.com/doi/full/10.1080/15476286.2017.1306170

(search for "ncbi" or "accession")

Samples: https://www.ncbi.nlm.nih.gov/sra?term=SRP081251

SRR numbers for each sample (4 samples):

- wildtype replicate 1: `SRR4018514`
- wildtype replicate 2: `SRR4018515`
- deletion replicate 1: `SRR4018516`
- deletion replicate 2: `SRR4018517`

### Retrieving data

For retrieving the data with the SRR numbers above the following commands were used on the access node (in /Day7/fastq directory).
```
fasterq-dump SRR4018514
fasterq-dump SRR4018515
fasterq-dump SRR4018516
fasterq-dump SRR4018517
```

This retrieved the associated .fastq files.

Metadata was retrieved with the following command (does not work for `SRR4018517`):
```
grabseqs sra -l -t 4 -m metadata.csv SRP081251
```
----
The metadata for Kröger et al. (https://www.sciencedirect.com/science/article/pii/S1931312813004113?via%3Dihub) was retrieved as well.
```
grabseqs sra -l -t 4 -m metadata_kroeger.csv SRP028811
```
----

### Quality Control

It was proceeded with the four downloaded .fastq files. Qualitay control was done by using the tool FASTQC. The following script was executed (after creating /fastqc_output directory).
```
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --time=1:00:00
#SBATCH --job-name=fastqc
#SBATCH -D ./
#SBATCH --output=./fastqc.out
#SBATCH --error=./fastqc.out
#SBATCH --partition=all
#SBATCH --reservation=biol217


module load miniconda3/4.7.12.1
source activate
conda activate /home/sunam226/.conda/envs/grabseq

module load fastqc

#fastqc

fastqc -t 4 -o ../fastqc_output ../*.fastq

#this prints the required resources into your logfile
jobinfo
```
The files (reports) were merged using MULTIQC. The following command was executed in the terminal after creating a /multiqc_output directory (-d accesses folder and all its subfolders).
```
multiqc -d . -o multiqc_output
```

No sequences were flagged as poor quality (they were cleaned before). The samples do not need further trimming (by FASTP).

## RNAseq - example using Kröger et al. paper

Start conda environement "reademption":

```
conda activate /home/sunam226/.conda/envs/reademption
```
See software in current conda environment:
```
conda list
```
### Setup folders and files

Navigate to example folder. Create folder structure (change species to species worked with).
```
reademption create --project_path READemption_analysis --species salmonella="Salmonella Typhimurium"
```
Downloaded .fastq files (reads) are moved to /input/reads in new folder structure.

Download reference sequence .fasta files (first variable for long link):
```
FTP_SOURCE=ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/Bacteria/Salmonella_enterica_serovar_Typhimurium_SL1344_uid86645/
```
```
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_016810.fa $FTP_SOURCE/NC_016810.fna
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_017718.fa $FTP_SOURCE/NC_017718.fna
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_017719.fa $FTP_SOURCE/NC_017719.fna
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_017720.fa $FTP_SOURCE/NC_017720.fna
```
Download genome annotation file:
```
wget -P READemption_analysis/input/salmonella_annotations https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/210/855/GCF_000210855.2_ASM21085v2/GCF_000210855.2_ASM21085v2_genomic.gff.gz
```
Unzip the file:
```
gunzip READemption_analysis/input/salmonella_annotations/GCF_000210855.2_ASM21085v2_genomic.gff.gz
```
To be able to match annotations with reference genome sequence, a identifier is needed (NC_xxx...). To make them look the same (change header in sequence file):
```
sed -i "s/>/>NC_016810.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_016810.fa
sed -i "s/>/>NC_017718.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_017718.fa
sed -i "s/>/>NC_017719.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_017719.fa
sed -i "s/>/>NC_017720.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_017720.fa
```
Download raw reads into /reads directory.
```
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/InSPI2_R1.fa.bz2
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/InSPI2_R2.fa.bz2
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/LSP_R1.fa.bz2
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/LSP_R2.fa.bz2
```
Now the pipeline can be run.

### Pipeline

Run with bash script.
```
#!/bin/bash
#SBATCH --job-name=reademption
#SBATCH --output=reademption.out
#SBATCH --error=reademption.err
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --qos=long
#SBATCH --time=0-05:00:00
#SBATCH --partition=all
#SBATCH --export=NONE
#SBATCH --reservation=biol217

#activating conda
module load miniconda3/4.7.12.1
source activate
conda activate /home/sunam226/.conda/envs/reademption

reademption align -p 4 --poly_a_clipping --project_path READemption_analysis
reademption coverage -p 4 --project_path READemption_analysis
reademption gene_quanti -p 4 --features CDS,tRNA,rRNA --project_path READemption_analysis
reademption deseq -l InSPI2_R1.fa.bz2,InSPI2_R2.fa.bz2,LSP_R1.fa.bz2,LSP_R2.fa.bz2 -c InSPI2,InSPI2,LSP,LSP -r 1,2,1,2 --libs_by_species salmonella=InSPI2_R1,InSPI2_R2,LSP_R1,LSP_R2 --project_path READemption_analysis
reademption viz_align --project_path READemption_analysis
reademption viz_gene_quanti --project_path READemption_analysis
reademption viz_deseq --project_path READemption_analysis

conda deactivate

#this prints the required resources into your logfile
jobinfo
```

## RNAseq project

A redemption folder structure was created in a new directory by the following command.
```
reademption create --project_path READemption_analysis --species methanosarcina="Methanosarcina mazei"
```

The downloaded reads from the paper (in the morning: see RNAseq - preparation) were moved to the /reads folder.
```
cp ../fastq/*.fastq ./READemption_analysis/input/reads/
```
The reference genome and gene annotation file were retrieved from [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NC_003901.1?report=genbank) and moved into the determined folders.

The annotations can already be matched to the reference genome (same file start/identifier present).

The following bash script was executed.
```
#!/bin/bash
#SBATCH --job-name=reademption
#SBATCH --output=reademption.out
#SBATCH --error=reademption.err
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --qos=long
#SBATCH --time=0-05:00:00
#SBATCH --partition=all
#SBATCH --export=NONE
#SBATCH --reservation=biol217

#activating conda
module load miniconda3/4.7.12.1
source activate
conda activate /home/sunam226/.conda/envs/reademption

reademption align --fastq -p 4 --poly_a_clipping --project_path READemption_analysis
reademption coverage -p 4 --project_path READemption_analysis
reademption gene_quanti -p 4 --features CDS,tRNA,rRNA --project_path READemption_analysis
reademption deseq -l sRNA_R1.fastq,sRNA_R2.fastq,wt_R1.fastq,wt_R2.fastq -c sRNA,sRNA,wt,wt -r 1,2,1,2 --libs_by_species methanosarcina=sRNA_R1,sRNA_R2,wt_R1,wt_R2 --project_path READemption_analysis
reademption viz_align --project_path READemption_analysis
reademption viz_gene_quanti --project_path READemption_analysis
reademption viz_deseq --project_path READemption_analysis

conda deactivate

#this prints the required resources into your logfile
jobinfo
```



