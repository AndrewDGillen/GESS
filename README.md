# GESS - the Gene Expression Similarity Score

## Defining GESS

Gene expression patterns can be highly informative of a gene's biology - after all, temporal and spatial regulation is a key part of controlling gene function. However, comparing expression patterns is not trivial, and requires the user to define how "similar" two expression patterns are.

The Gene Expression Similarity Score (lovingly referred to as GESS) is a broadly-applicable solution to this problem. For any given dataset with multiple contexts (ie tissues, timepoints, cell typess...) available, each gene can be ranked by context from highest expression to lowest, thus defining that gene’s expression pattern. We considered the major defining features of this expression pattern to be the ordering of these contexts across the totality of the annotation, and the difference in gene enrichment within like contexts across datasets.

Thus, given two datasets Q (Query) and T (Target) with at least some comparable contexts, the similarity in expression between any two genes can be defined across datasets as:

![Screenshot from 2024-05-29 14-59-43](https://github.com/AndrewDGillen/GESS/assets/88687148/5a8264b6-2390-4f9f-b44e-f40e2863374c)

For a given gene in a set of like contexts, the expression difference between Q and T can be defined as:

![Screenshot from 2024-05-29 15-00-51](https://github.com/AndrewDGillen/GESS/assets/88687148/5f44b2ad-3857-4a23-9f8d-e4d7da3d9575)

Where EQ  and ET represent log2 expression of the gene in Q and T datasets respectively. Expression weight is applied as a modifier designed to emphasise expression differences where datasets differ in fold-change directionality within a context relative to baseline. 

Expression weight is defined as:

![Screenshot from 2024-05-29 15-01-58](https://github.com/AndrewDGillen/GESS/assets/88687148/b6629caf-635c-47d9-a896-1f5d1afab474)

For a given context, the position difference between Q and T can be defined as:

![Screenshot from 2024-05-29 15-02-36](https://github.com/AndrewDGillen/GESS/assets/88687148/b09f15d3-4001-46da-9feb-b199b7c4e8bb)

Where PQ  and PT represent the position of a given context within a ranked list of all contexts in datasets Q and T respectively. 

Position Uniqueness is applied as a modifier to alter the relative scale of position differences based on how unique the expression noted in a specific context is across both datasets  within the whole annotation level. Position Uniqueness is defined as:

![image](https://github.com/AndrewDGillen/GESS/assets/88687148/6faca0d1-125c-4e09-9c10-593f7c490da6)

Where UQ  and UT represent Context Uniqueness (U) of a given context in datasets Q and T respectively. 

Context Uniqueness is defined as:

![Screenshot from 2024-05-29 15-04-26](https://github.com/AndrewDGillen/GESS/assets/88687148/1550953c-3490-4574-8a6b-966b704a2d95)

## Using GESS

GESS, as you'll appreciate, is somewhat complicated to manually calculate. To facilitate general usage of GESS, we provide the programs in the current repository

### Installation

-> **For command line version** - Simply download this repository
-> **For python module** Downloadable from PyPI (https://pypi.org/project/GESS/) using the command pip install GESS

### Requirements
  Python3 >= 3.10
  numpy
  pandas
  seaborn
  fastcluster
  matplotlib
  scanpy
  
### Running GESS 
The following examples are for running GESS from the command line. For a python walkthrough, check "GESS_Tutorial.ipynb" for a Jupyter notebook tutorial

With a single-cell dataset of interest "Sample1.h5ad", which contains information of interest across "cell type" and "tissue" levels, GESS can be calculated using Expression data between two genes with the following command:

  python GESSfinder.py -qg <Gene1> -tg <Gene2> -qd <Sample1.h5ad> -anno <"cell type" "tissue"> -mode <"expression">

A second dataset can be defined using the -td argument.

Single cell data can be analysed using either "expression" (average expression in each defined annotation level) or "prevalence" (number of cells of each annotation expressing a gene) modes.

Bulk RNA sequencing data can also be analysed based on an Enrichment matrix (Expression of gene in each sample is normalised against a calibrator sample). This also requires manual definition of each annotation level in a seperate annotation.csv (see example in Example Data folder). Bulk RNA seq GESS can then be run with the command:

  python GESSfinder.py -qg <Gene1> -tg <Gene2> -qd <EnrichmentMatrix> -a <Annotation.csv> -qs <Query Selection> -ts <Target Selection> -mode <"bulk"> -anno <"Tissue" "Function">

This repository also contains code for automatically calculating GESS across multiple genes of interest, and hierarchically clustering results based on their GESS value. This function can be called using GESSmatricise.py, defining datasets as above. Gene lists can be defined in line-separated .txt files (see example in Example Data folder). For example, for single-cell data:

  python GESSmatricise.py -iq <GeneList1.txt> -qd <Sample1.h5ad> -anno <"cell type" "tissue"> -mode <"expression">
