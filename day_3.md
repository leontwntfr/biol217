# Day 3

## Metagenomics: Contig Visualization, Quality Assessment of Assemblies, Genome Binning (Mapping)

Creation of a folder structure with `mkdir`.

\
Copying of initial data for today by:
```
cp -r /home/sunam226/Day3/* /work_beegfs/sunam237/Day3/
```
In "3_coassembly" folder looking for number of contigs by:
```
grep -c ">" final.contigs.fa
```
- 57414 contigs
----
Reading in terminal commands:
- head: shows beginning of file as sepcified
- tail: shows end of file as specified
- less: shows whole file in a new terminal
- cat: shows whole file in the same terminal
----

### Assembly Visualization

Conversion of fasta file "final.contigs.fa" into FASTG format:
```
megahit_toolkit contig2fastg 99 final.contigs.fa > final.contigs.fastg                   
```

This makes visualization using Bandage possible. For this, copying of FASTG file to user computer/desktop/ is necessary.

Figure of view in Bandage:

![image](./resources/contigs.png)
Observations:
- contigs are represented in different colors
- from top to bottom the assembled contigs are arranged from the largest to smallest contigs
- contigs are represented with different shapes/curvatures and sometimes loops 


### Assembly Quality Assessment
For quality assessment fo the assemblies using METAQUAST the following bash script was created and submitted via `sbatch metaquast.sh`:
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=metaquast
#SBATCH --output=/work_beegfs/sunam237/Day3/3_metaquast/metaquast.out
#SBATCH --error=/work_beegfs/sunam237/Day3/3_metaquast/metaquast.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


module load miniconda3/4.7.12.1
source activate
conda activate /home/sunam226/.conda/envs/anvio

#metaquast
                                       
metaquast -t 6 -o /work_beegfs/sunam237/Day3/3_metaquast -m 500 /work_beegfs/sunam237/Day3/3_coassembly/final.contigs.fa

#this prints the required resources into your logfile
jobinfo
```

The output "report.html" was inspected.
Questions:

What is your N50 value? Why is this value relevant?

- N50 value: 2963
- relevance: length of contigs, where equal or longer lengths account for 50% of the bases of the assembly, higher values should indicate better quality of the assembly

How many contigs are assembled

- 57.414

What is the total length of the contigs?

- 145.675.865

### Genome Binning - Mapping

Before mapping, reformatting of the fasta files is necessary. The following bash script was created and submitted:
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=reformat
#SBATCH --output=/work_beegfs/sunam237/Day3/4_mapping/reformat.out
#SBATCH --error=/work_beegfs/sunam237/Day3/4_mapping/reformat.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


module load miniconda3/4.7.12.1
source activate
conda activate /home/sunam226/.conda/envs/anvio

#reformat
                                       
anvi-script-reformat-fasta /work_beegfs/sunam237/Day3/3_coassembly/final.contigs.fa -o /work_beegfs/sunam237/Day3/4_mapping/contigs.anvio.fa --min-len 1000 --simplify-names --report-file name_conversion.txt

#this prints the required resources into your logfile
jobinfo
```

The output file "contigs.anvi.fa" was used for mapping by bowtie2 in the next step (for contig coverage). For this, the following bash script was created and submitted (for bowtie2, the paths for `reference_in` and `bt2_index_base` should be in the same folder):
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=bowtie2
#SBATCH --output=/work_beegfs/sunam237/Day3/4_mapping/bowtie2.out
#SBATCH --error=/work_beegfs/sunam237/Day3/4_mapping/bowtie2.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


module load miniconda3/4.7.12.1
#source activate
conda activate /home/sunam226/.conda/envs/anvio

#mapping
                                       
bowtie2-build /work_beegfs/sunam237/Day3/4_mapping/contigs.anvio.fa /work_beegfs/sunam237/Day3/4_mapping/contigs.anvio.fa.index

#this prints the required resources into your logfile
jobinfo
```

The output index files together with the clean reads (fastp) were used for the actual mapping using bowtie2 again. The following bash script was created and submitted (in for loop: no full patch should be between "`", use cd to navigate to corresponding path before instead):
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=mapping
#SBATCH --output=/work_beegfs/sunam237/Day3/4_mapping/mapping.out
#SBATCH --error=/work_beegfs/sunam237/Day3/4_mapping/mapping.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


module load miniconda3/4.7.12.1
source activate
conda activate /home/sunam226/.conda/envs/anvio

#mapping

bowtie2 --very-fast -x /work_beegfs/sunam237/Day3/4_mapping/contigs.anvio.fa.index -1 /work_beegfs/sunam237/Day3/2_fastp/BGR_130305_mapped_R1.fastq.gz -2 /work_beegfs/sunam237/Day3/2_fastp/BGR_130305_mapped_R2.fastq.gz -S /work_beegfs/sunam237/Day3/4_mapping/BGR_130305_sample.sam
bowtie2 --very-fast -x /work_beegfs/sunam237/Day3/4_mapping/contigs.anvio.fa.index -1 /work_beegfs/sunam237/Day3/2_fastp/BGR_130527_mapped_R1.fastq.gz -2 /work_beegfs/sunam237/Day3/2_fastp/BGR_130527_mapped_R2.fastq.gz -S /work_beegfs/sunam237/Day3/4_mapping/BGR_130527_sample.sam
bowtie2 --very-fast -x /work_beegfs/sunam237/Day3/4_mapping/contigs.anvio.fa.index -1 /work_beegfs/sunam237/Day3/2_fastp/BGR_130708_mapped_R1.fastq.gz -2 /work_beegfs/sunam237/Day3/2_fastp/BGR_130708_mapped_R2.fastq.gz -S /work_beegfs/sunam237/Day3/4_mapping/BGR_130708_sample.sam

#cd /work_beegfs/sunam237/Day3/2_fastp/
#for i in `ls *mapped_R1.fastq.gz`;
#do
#    second=`echo ${i} | sed 's/_R1/_R2/g'`
#    bowtie2 --very-fast -x /work_beegfs/sunam237/Day3/4_mapping/contigs.anvio.fa.index -1 ${i} -2 ${second} -S /work_beegfs/sunam237/Day3/4_mapping/"$i".sam 
#done

#this prints the required resources into your logfile
jobinfo
```

Transformation of the ouput files (.sam) to .bam files was executed via the following bash script:
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=transform
#SBATCH --output=/work_beegfs/sunam237/Day3/4_mapping/transform.out
#SBATCH --error=/work_beegfs/sunam237/Day3/4_mapping/transform.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


module load miniconda3/4.7.12.1
source activate
conda activate /home/sunam226/.conda/envs/anvio

#transform

module load samtools

cd /work_beegfs/sunam237/Day3/4_mapping/
for i in *.sam; do samtools view -bS $i > "$i".bam; done

#this prints the required resources into your logfile
jobinfo
```