
Author : Ilango Guy
Maintainer : Ilango Guy
Draw : Adeline Cezard
Licence : CEPR
Affiliation: Research Center for Respiratory Diseases ; Team 4  (France)

-----------------------------------------------
# Analyse youR daTa Yourself (ARTY) V2
-----------------------------------------------


Analyze youR daTa Yourself V2 is now a shiny application. 


------------------------------------------------
# Prerequisite
------------------------------------------------ 
DESeq2
```
BiocManager::install("DESeq2")
```
tidyverse
```
install.packages("tidyverse")
```

Shiny
```
install.packages("shiny")
```
plotly
```
install.packages("plotly")
```


------------------------------------------------
# How to use:

------------------------------------------------ 

Run the script from R or Rstudio


Features : 

**Design** 
-----------

Create a design matrix from your data. 
For example  : if you have this matrix 

| id | A_1 | A_2 | A_3 | B_1 | B_2| B_3 |
| --- | --- | --- | --- | --- | --- | --- |
| gene1 | 10 | 2 | 3 | 8 | 100 | 6 |
| gene2 | 9 | 20 | 1 | 0 | 41 | 8 |



You will have to put "A,B" as groups and 3 as number of replicate to generate your design matrix

**DE** 
-----------

Create a differential expression matrix from your count matrix and the design matrix ( you can generate it from design tab) 

**Volcano** 
-----------

Generate a Volcano plot from DE matrix

**Heatmap**
-----------

Generate a plot of heatmap from the count matrix

**PCA** 
-----------

Add the count matrix then the design matrix and choose the annotation from design matrix
Choose if you want a 2d or 3d pca
