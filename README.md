# Trans-eQTL mapping in gene sets identifies network effects of genetic variants

All trans-eQTLs that are associated with gene co-expression networks and biological pathways can be found [here](http://www.networks-liulab.org/transPCO). 

Our preprint is online [here](https://www.biorxiv.org/content/10.1101/2022.11.11.516189v1). 



## Introduction

We developed a pipeline called trans-PCO to detect trans effects of gene networks (gene module). It combines careful read and gene filters with a principal component (PC)-based multivariate association test. 

Trans-PCO allows the use of many types of gene groups or sets. For example, genes with correlated expression levels in a co-expression gene network, or genes in the same pathway, or protein-protein interaction network.



![](docs/pipeline.png){width=60%}

Figure: Three main steps in trans-PCO pipeline. The first step of trans-PCO pre-processes RNA-seq data to reduce false positive trans-eQTL associations due to read alignment errors. The second step involves grouping genes into gene sets, such as co-expression modules or biological pathways. The last step tests for trans-eQTLs of each gene set by a PC-based multivariate association test.

<img src="docs/pipeline.png" width="5000">
