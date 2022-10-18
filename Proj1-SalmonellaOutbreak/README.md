# Computational Biology - Salmonella Outbreak

## Table of contents
* [General info](#general-info)
* [Requirements](#technologies)
* [Setup](#setup)
* [Results](#results)


## General info
The tool takes as input two different FASTA files (the Wild strain and the Variantstrain) and give in output a list of SNPs, which are a germline substitutions of a single nucleotide at a specific position in the genome.

	
## Requirements
Project is written in Python and works with the following version:
*  python_version >= '3.6'
All the other requirements can be installed through the terminal writing:
```
$ pip install -r requirements.txt
```
Main libraries:
* levenshtein==0.16.0
* matplotlib==3.5.0
* scipy==1.7.3
* termcolor==1.1.0
* biopython==1.79
* numpy==1.21.4

## Setup
Before running the programme, reading the  [finalreport.pdf](finalreport.pdf) is highly recommended.

To run this project, after installing the requirements, open the terminal and launch it as:
```
$ python main0.py -f1 "PATH TO salmonella-enterica.reads.fna" -f2 
"PATH TO salmonella-enterica-variant.reads.fna"
```
where:
* **salmonella-enterica.reads.fna** is the Wild Type file
* **salmonella-enterica-variant.reads.fna** is the Variant Type file

The files used can be downloaded at: https://cloud-ljk.imag.fr/index.php/s/HkxDLozHRcqBcqz

## Results
Results are presented as:
![myimage-alt-tag](https://github.com/HanningYang/ComputationalBiology/blob/main/Proj1-SalmonellaOutbreak/result.png)

[result.png](result.png)
```
SNP
    R: TTACATGCCAATACAATGTAGGCTGCTCTACACCTAGCTTCTGGGCGAGTTTACGGGTTGTTAAACCTTCGATTCCGACCTCATTAAGCAGCTCTAATGCG
    V: TTACATGCCAATACAATGTAGGCTGCTCTACACCTAGCTTCTGGGCGAGGGGACGGGTTGTTAAACCTTCGATTCCGACCTCATTAAGCAGCTCTAATGCG
   	 Distance: 3
------------------------------------------------------
SNP
    R: CGCATTAGAGCTGCTTAATGAGGTCGGAATCGAAGGTTTAACAACCCGTAAACTCGCCCAGAAGCTAGGTGTAGAGCAGCCTACATTGTATTGGCATGTAA
    V: CGCATTAGAGCTGCTTAATGAGGTCGGAATCGAAGGTTTAACAACCCGTCCCCTCGCCCAGAAGCTAGGTGTAGAGCAGCCTACATTGTATTGGCATGTAA
   	 Distance: 3
```



