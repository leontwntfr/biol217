# Day 8

## RNAseq

### Repitition of project publication script from yesterday

> see publication_data.md ([LINK](https://github.com/AammarTufail/Bioinformatics_Master_Module2023/blob/main/Day-7/publication_data.md))

Copy command for reads (in new directory in /Day8/pub_data_1):
```
cp ../../Day7/fastq/*.fastq ./READemption_analysis/input/reads
```

In the given bash script, the names of the read files were changed accordingly (now with load python module and .fa in deseq command).
```
#!/bin/bash
#SBATCH --job-name=pub_data
#SBATCH --output=pub_data.out
#SBATCH --error=pub_data.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH	--qos=long
#SBATCH --time=0-05:00:00
#SBATCH --partition=all
#SBATCH --reservation=biol217

source ~/.bashrc

module load miniconda3/4.7.12.1
module load python/3.7.4
conda activate /home/sunam226/.conda/envs/reademption
################################### ---CALCULATIONS---
#aligning:
reademption align --fastq -f READemption_analysis --poly_a_clipping --min_read_length 12 --segemehl_accuracy 95  

# coverage:
reademption coverage -p 4 --project_path READemption_analysis 

#gene wise quanty:
reademption gene_quanti -p 4 --features CDS,tRNA,rRNA --project_path READemption_analysis 

#differential gene expression:
reademption deseq -l sRNA_R1.fa,sRNA_R2.fa,wt_R1.fa,wt_R2.fa -c sRNA,sRNA,wt,wt \
	-r 1,2,1,2 --libs_by_species Methanosarcina=sRNA_R1,sRNA_R2,wt_R1,wt_R2 --project_path READemption_analysis

############################## ---PLOTS---
reademption viz_align --project_path READemption_analysis
reademption viz_gene_quanti --project_path READemption_analysis
reademption viz_deseq --project_path READemption_analysis
conda deactivate
jobinfo
```

Did not work again (deseq2 permission denied).

Installation of reademption environment on own sunam: [LINK](https://reademption.readthedocs.io/en/latest/installation.html)



## R

Guide for what plot to choose: [LINK](https://apandre.files.wordpress.com/2011/02/chartchooserincolor.jpg)

```
getwd() # to get working directory

# set working directory
setwd()

# create folders ----
dir.create('data')
dir.create('data/raw_data')




x = 2+2 # better not use =
x <- 2+2 # use this instead

x <- y <- z <- 5
class(x)

x <- 'Hello World'
class(x)

x <- TRUE
class(x)

x <- 1+2i
class(x)

x <- charToRaw('Hello World') # stored differently (without any kind of compression, easy to adjust)
class(x)

x <- data.frame(matrix(1:6, nrow = 2, ncol = 3))
class(x)



data()

data('iris')

View(iris)

class(iris$Sepal.Length)


plot(iris)
boxplot(data = iris, iris$Petal.Length~iris$Species) # ~ for every species individually



#### Plotting ----

# install and activate packages

#install.packages('ggplot2')
library(ggplot2)

install.packages(c('readxl', 'plotly'), dependencies = TRUE) # install multiple packages

install.packages('tidyverse')




## ggplot2

ggplot(iris, mapping = aes(Species, Sepal.Length, fill = Species)) + geom_violin()

# -> Andrew Abela chart guide

ggplot(iris, mapping = aes(Petal.Length, Sepal.Length, color = Species)) + geom_point()

ggplot(iris, mapping = aes(Petal.Length, Sepal.Length, shape = Species)) + geom_point()

ggplot(iris, mapping = aes(Petal.Length, Sepal.Length, size = Petal.Width)) + geom_point()

# save high quality plot

plot1 <- ggplot(iris, mapping = aes(Petal.Length, Sepal.Length, color = Species)) + geom_point()

ggsave('test_plot_comp.tiff', plot = plot1, height = 6, width = 8, unit = 'in', dpi = 300, compression = 'lzw')


# svg by hand possible (export)

#install.packages('svglite')

#ggsave('test_plot.svg', plot = plot1, height = 6, width = 8, unit = 'in')





#### ---- afternoon

library(tidyr)

sp_table <- spread(data = iris, key = iris$Species, value = c(iris$Sepal.Length, iris$Sepal.Width, iris$Petal.Length, iris$Petal.Width))

# https://uc-r.github.io/tidyr
# how to use this -> long to wide QUESTION FOR PROOCOL


hist(iris$Sepal.Length)

hist(iris$Sepal.Width) # -> parametric data (bell curve)

hist(iris$Petal.Length)

hist(iris$Petal.Width)

boxplot(iris$Sepal.Width) # parametric data if median is near center (mean)

# paired t-test between the same thing but different time points/conditions
# unpaired: comparing something different


# choose data set and make plots ----

data('chickwts')

ggplot(chickwts, aes(feed, weight)) + geom_boxplot()

# other data

data('midwest')




### continue

ggplot(iris, aes(Species, Petal.Length)) + geom_boxplot() +
  facet_wrap(~Species)



#setRepositories()
```

## RNAseq - continuation

Copy output from sunam226 to new /output directory.
```
cp /home/sunam226/Day7/example_rna_seq/READemption_analysis/output/* ./
```
Permission denied

### Integrated Genome Browser

Introduction
