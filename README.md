# CasPhy
Li Lin<br>

## Introduction
CasPhy is an ultrafast phylogenetic reconstruction algorithm designed for lineage-tracing data and conventional multiple sequence alignment (MSAs).

## Requirements

-   Linux on an x86-64 system that supports AVX2 (most Intel and AMD CPUs)
-   CUDA Toolkit \>=12.9
-   R \>=4.0.0
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
install.packages("ape")
install.packages("bigmemory")
install.packages("jsonlite")
install.packages("CasPhy_1.5.7.tar.gz")
```

To buid phylogenies from MSA using VeryFastTree, VeryFastTree (<https://github.com/citiususc/veryfasttree>) is required.

To buid phylogenies from MSA using RAxML-ng, RAxML-ng (<https://github.com/amkozlov/raxml-ng/>) and EPA-ng (<https://github.com/pierrebarbera/epa-ng>) is required.

An HTML tutorial is available at `vignettes.html`. Data files used in tutorial can be found in `examples/`.
