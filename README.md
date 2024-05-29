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

->
