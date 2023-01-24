# Day 2

## Metagenomics: Quality Control with FASTQC and FASTP and Assembly with MEGAHIT

The sample data folder "0_raw_reads" was copied from `/home/sunam226/Day2/` to a previously created folder (`mkdir metagenomics/`) on `/work_beegfs/sunam237/` by:

```
cp -r /home/sunam226/Day2/0_raw_reads/ /work_beegfs/sunam237/metagenomics/
```

The following bash script ("fastqc.sh") was created in the $HOME directory:

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=fastqc
#SBATCH --output=/work_beegfs/sunam237/fastqc.out
#SBATCH --error=/work_beegfs/sunam237/fastqc.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


module load miniconda3/4.7.12.1
conda activate /home/sunam226/.conda/envs/anvio
#fastqc
module load fastqc
cd /work_beegfs/sunam237/metagenomics/0_raw_reads/
for i in *.gz; do fastqc $i; done
#this prints the required resources into your logfile
jobinfo
```

The bash script was excecuted with the following command in the "caucluster":

```
sbatch fastqc.sh
```

We analyzed the graphs as output files (.html) from the analysis:
- none of the samples consisted of sequences with poor quality that needed trimming

Nevertheless, we continued to learn to trim the samples with FASTP.

\
The following bash script ("fastp.sh") was created for the FASTP process:

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=fastp
#SBATCH --output=/work_beegfs/sunam237/fastp.out
#SBATCH --error=/work_beegfs/sunam237/fastp.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


module load miniconda3/4.7.12.1
conda activate /home/sunam226/.conda/envs/anvio

#fastp
module load fastp
cd /work_beegfs/sunam237/metagenomics/0_raw_reads/
mkdir ../output_folder

for i in `ls *_R1.fastq.gz`;
do
    second=`echo ${i} | sed 's/_R1/_R2/g'`
    fastp -i ${i} -I ${second} -R _report -o ../output_folder/"${i}" -O ../output_folder/"${second}" -t 6 -q 20

done

#this prints the required resources into your logfile
jobinfo
```
The bash script was executed by:
```
sbatch fastp.sh
```

The outputs were used for assembly with MEGAHIT.

\
The following bash script ("assembly_megahit.sh") was created for the assembly (MEGAHIT requires not to create the output directory beforehand):
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=megahit
#SBATCH --output=/work_beegfs/sunam237/megahit.out
#SBATCH --error=/work_beegfs/sunam237/megahit.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


module load miniconda3/4.7.12.1
source activate
conda activate /home/sunam226/.conda/envs/anvio

#megahit
cd /work_beegfs/sunam237/metagenomics/output_folder
                                       
megahit -1 BGR_130305_R1.fastq.gz -1 BGR_130527_R1.fastq.gz -1 BGR_130708_R1.fastq.gz -2 BGR_130305_R2.fastq.gz -2 BGR_130527_R2.fastq.gz -2 BGR_130708_R2.fastq.gz --min-contig-len 1000 --presets meta-large -m 0.85 -o ../coassembly_output"$i" -t 20                      

#this prints the required resources into your logfile
jobinfo
```
The bash script was executed by:
```
sbatch assembly_megahit.sh
```