# EMS VCF Processing Pipeline   
This repository contains the Jupyter notebook **_Ems_vcf_pipeline.ipynb**, a reproducible workflow for raw fastq samples to vcfs.   
The notebook is designed for reproducibility on the UC Davis Farm HPC cluster.  
  
## Overview  

The pipeline performs the following steps:  
1. **select reference fasta and read fastq files directory**  
2. **run read qc, trimming, alignment (bwa mem) and vcf calling (bcftools) **    
3. **Generate summary statistics**  
The notebook is intended for the bond lab and friends to speed up systematics workflows.  

## Requirements  
This project assumes the following:  
1. You are running this in Ondemand as a Jupyter notebook on the farm hpc at UCDavis  
2. You have selected the kernel "Python [conda env:Bondlab_phylo_env]"  This has all the needed software  

### Software  
- TO re-iterate the above - start with ondemand on the UCD Farm hpc, and select the kernel "Python [conda env:Bondlab_phylo_env]"  

## Repository Structure  
To run this it will assume you have put your reference genome in:  
/_Ems_vcf_pipeline/data/references   (if it is not already there)  
And that you will put your fastq reads in a new directory in  
/_Ems_vcf_pipeline/data/read_directories  (if they are not already there)  
results will be in:  
/_Ems_vcf_pipeline/data    
full instructions are in the jupyter notebook.  
