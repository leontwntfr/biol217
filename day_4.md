# Day 4

## Metagenomics: Genome Binning, MAGs Quality Assessment, Bin Refinement

Initiation/test of ANVI'O via:
```
anvi-cluster-contigs -h
```
```
anvi-self-test --suite mini
```

### Genome Binning - Contigs data preparation

Creation of database for our contigs by, first, running the following script:
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=anvicontigsdb
#SBATCH --output=/work_beegfs/sunam237/Day3/5_anvio-profiles/scripts/anviocontigsdb.out
#SBATCH --error=/work_beegfs/sunam237/Day3/5_anvio-profiles/scripts/anviocontigsdb.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


module load miniconda3/4.7.12.1
source activate
conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

#anvicontigsdb

anvi-gen-contigs-database -f /work_beegfs/sunam237/Day3/4_mapping/contigs.anvio.fa -o ../contigs.db -n 'biol217'

#this prints the required resources into your logfile
jobinfo
```

The output .db file was used for a HMM (Hidden Markov Models) search. For that, the following script was executed.
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=anviohmmsearch
#SBATCH --output=/work_beegfs/sunam237/Day3/5_anvio-profiles/scripts/anviohmmsearch.out
#SBATCH --error=/work_beegfs/sunam237/Day3/5_anvio-profiles/scripts/anviohmmsearch.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


module load miniconda3/4.7.12.1
source activate
conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

#anvicontigsdb

anvi-run-hmms -c ../contigs.db -T 12

#this prints the required resources into your logfile
jobinfo
```

Interactive ANVI'O contigs stats were accessed and displayed by usage of two terminals and firefox. First, the following command was run in a caucluster terminal.
```
srun --reservation=biol217 --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --partition=all /bin/bash
```
This provided the following node information: 
- node010

In the same terminal the following commands were run (in folder with .db file).
```
conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1
```
```
anvi-display-contigs-stats contigs.db
```

This started a server with the address `http://0.0.0.0:8084`.

In a new terminal (not logged in into caucluster) the server was accessed with the following commands.
```
ssh -L 8060:localhost:8080 sunam237@caucluster-old.rz.uni-kiel.de
```
```
ssh -L 8080:localhost:8080 node010
```

The interactive interface of the stats were accessed in firefox with the link `http://127.0.0.1:8060`.

### Binning with ANVI'O

First, sorting and indexing was executed with the following bash script (navigation to specific folder necessary using `cd`).
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=binningsortindex
#SBATCH --output=/work_beegfs/sunam237/Day3/5_anvio-profiles/scripts/binningsortindex.out
#SBATCH --error=/work_beegfs/sunam237/Day3/5_anvio-profiles/scripts/binningsortindex.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


module load miniconda3/4.7.12.1
source activate
conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

#binningsortindex

cd /work_beegfs/sunam237/Day3/4_mapping/

for i in *.bam; do anvi-init-bam $i -o /work_beegfs/sunam237/Day3/5_anvio-profiles/"$i".sorted.bam; done

#this prints the required resources into your logfile
jobinfo
```

Now, binning with ANVI'O can be performed. For this, the following bash script for ANVI'O profiling was submitted.

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=anvioprofile
#SBATCH --output=/work_beegfs/sunam237/Day3/5_anvio-profiles/scripts/anvioprofile.out
#SBATCH --error=/work_beegfs/sunam237/Day3/5_anvio-profiles/scripts/anvioprofile.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


module load miniconda3/4.7.12.1
source activate
conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

#anvioprofile

cd /work_beegfs/sunam237/Day3/5_anvio-profiles/

mkdir /work_beegfs/sunam237/Day3/5_anvio-profiles/profiling/
for i in `ls *.sorted.bam | cut -d "." -f 1`; do anvi-profile -i "$i".bam.sorted.bam -c contigs.db -o /work_beegfs/sunam237/Day3/5_anvio-profiles/profiling/$i -T 12; done

#this prints the required resources into your logfile
jobinfo
```

ANVI'O merge was used by executing the following bash script.
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=anviomerge
#SBATCH --output=/work_beegfs/sunam237/Day3/5_anvio-profiles/scripts/anviomerge.out
#SBATCH --error=/work_beegfs/sunam237/Day3/5_anvio-profiles/scripts/anviomerge.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


module load miniconda3/4.7.12.1
source activate
conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

#anviomerge

anvi-merge /work_beegfs/sunam237/Day3/5_anvio-profiles/profiling/BGR_130305/PROFILE.db /work_beegfs/sunam237/Day3/5_anvio-profiles/profiling/BGR_130527/PROFILE.db /work_beegfs/sunam237/Day3/5_anvio-profiles/profiling/BGR_130708/PROFILE.db -o /work_beegfs/sunam237/Day3/5_anvio-profiles/merged_profiles -c /work_beegfs/sunam237/Day3/5_anvio-profiles/contigs.db --enforce-hierarchical-clustering

#this prints the required resources into your logfile
jobinfo
```

Binning with Metabat2 was performed with the following bash script.
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=metabat2
#SBATCH --output=/work_beegfs/sunam237/Day3/5_anvio-profiles/scripts/metabat2.out
#SBATCH --error=/work_beegfs/sunam237/Day3/5_anvio-profiles/scripts/metabat2.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


module load miniconda3/4.7.12.1
source activate
conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

#metabat2

anvi-cluster-contigs -p /work_beegfs/sunam237/Day3/5_anvio-profiles/merged_profiles/PROFILE.db -c /work_beegfs/sunam237/Day3/5_anvio-profiles/contigs.db -C METABAT --driver metabat2 -T 12 --just-do-it --log-file log-metabat2

anvi-summarize -p /work_beegfs/sunam237/Day3/5_anvio-profiles/merged_profiles/PROFILE.db -c /work_beegfs/sunam237/Day3/5_anvio-profiles/contigs.db -o /work_beegfs/sunam237/Day3/5_anvio-profiles/SUMMARY_METABAT -C METABAT

#this prints the required resources into your logfile
jobinfo
```
Binning with CONCORD was performed with the following bash script.
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=concoct
#SBATCH --output=/work_beegfs/sunam237/Day3/5_anvio-profiles/scripts/concoct.out
#SBATCH --error=/work_beegfs/sunam237/Day3/5_anvio-profiles/scripts/concoct.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


module load miniconda3/4.7.12.1
source activate
conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

#concoct

anvi-cluster-contigs -p /work_beegfs/sunam237/Day3/5_anvio-profiles/merged_profiles/PROFILE.db -c /work_beegfs/sunam237/Day3/5_anvio-profiles/contigs.db -C CONCOCT --driver concoct -T 12 --just-do-it

anvi-summarize -p /work_beegfs/sunam237/Day3/5_anvio-profiles/merged_profiles/PROFILE.db -c /work_beegfs/sunam237/Day3/5_anvio-profiles/contigs.db -o /work_beegfs/sunam237/Day3/5_anvio-profiles/SUMMARY_CONCOCT -C CONCOCT

#this prints the required resources into your logfile
jobinfo
```
Binning consolidation by DASTool was performed with the following bash script.
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=dastool
#SBATCH --output=/work_beegfs/sunam237/Day3/5_anvio-profiles/scripts/dastool.out
#SBATCH --error=/work_beegfs/sunam237/Day3/5_anvio-profiles/scripts/dastool.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


module load miniconda3/4.7.12.1
source activate
conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

#dastool

anvi-cluster-contigs -p work_beegfs/sunam237/Day3/5_anvio-profiles/merged_profiles/PROFILE.db -c //work_beegfs/sunam237/Day3/5_anvio-profiles/contigs.db -C consolidated_bins --driver dastool -T 20 --search-engine diamond -S METABAT,CONCOCT --log-file log_consolidation_of_bins --just-do-it

anvi-summarize -p work_beegfs/sunam237/Day3/5_anvio-profiles/merged_profiles/PROFILE.db -c /work_beegfs/sunam237/Day3/5_anvio-profiles/contigs.db -o /work_beegfs/sunam237/Day3/5_anvio-profiles/SUMMARY_consolidated_bins -C consolidated_bins

#this prints the required resources into your logfile
jobinfo
```
Python errors, so the output .html files were copied and validated.

Questions:

Number of ARCHAEA bins you got from MetaBAT2?
- 3
  
Number of ARCHAEA bins you got from CONCOCT?
- 2

Number of ARCHAEA bins you got after consolidating the bins?
- 2

### MAGs Quality Estimation

See collections on login node (using sunam226):
```
anvi-estimate-genome-completeness -p /home/sunam226/Day3/5_anvio-profiles/merged_profiles/PROFILE.db -c /home/sunam226/Day3/5_anvio_profiles/contigs.db --list-collections
```
Continuation tomorrow...