# CasPhy
Li Lin<br>

## Introduction
CasPhy is an ultrafast phylogenetic reconstruction algorithm designed for lineage-tracing data (CRISPR-Cas9 induced indels and base-editor induced point mutations) and conventional multiple sequence alignments (MSAs) of DNA.

## Requirements

-   Linux (developed on ubuntu 24.04) on an x86-64 system that supports AVX2 (most Intel and AMD CPUs)
-   CUDA Toolkit \>=12.9
-   R \>=4.3.1
-   R packages `ape`, `bigmemory`, and `jsonlite`

## Installation

Download CasPhy and compile the CUDA scripts. The `-arch` flag must be set according to the compute capability of your GPU.

```         
git clone https://github.com/LiLin-biosoft/CasPhy.git
cd CasPhy
for i in $(ls cu/|sed 's/\.cu//'); do nvcc cu/${i}.cu -o ${i} -arch=sm_120;done
```

In the R console, install CasPhy and its dependencies:

```         
install.packages("ape") ## developed on version 5.8.1
install.packages("bigmemory") ## developed on version 4.6.4
install.packages("jsonlite") ## developed on version 2.0.0
install.packages("CasPhy_1.5.7.tar.gz")
```
Installation will take several minutes.

To buid phylogenies from MSA using VeryFastTree, VeryFastTree (<https://github.com/citiususc/veryfasttree>; developed on version 4.0.5) is required.

To buid phylogenies from MSA using RAxML-ng, RAxML-ng (<https://github.com/amkozlov/raxml-ng/>; developed on version 1.2.2) and EPA-ng (<https://github.com/pierrebarbera/epa-ng>; developed on version 0.3.8) is required.

An HTML tutorial is available at `vignettes.html`. All code takes just several minutes. All Data files used in tutorial can be found in `examples/`.
