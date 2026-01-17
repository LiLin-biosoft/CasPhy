# CasPhy

Li Lin<br>

CasPhy is a phylogeny reconstruction algorithm for lineage-tracing data. It can be installed from Github:<br>

```         
install.packages('devtools')  
devtools::install_github("LiLin-biosoft/CasPhy", build_vignettes = T)
browseVignettes('CasPhy')
```

To use GPU acceleration, CUDA script is required:<br>

```         
git clone https://github.com/LiLin-biosoft/CasPhy-cuda.git
cd CasPhy-cuda
mkdir build
cd build
cmake ..
make
```
or

```         
git clone https://github.com/LiLin-biosoft/CasPhy-cuda.git
cd CasPhy-cuda
mkdir build
cd build
cmake -DCUDA_ARCH=89 .. ## set compute capacity to 89 for RTX 40 series
make
```
To perform phylogenetic reconstruction, VeryFastTree is required(<https://github.com/citiususc/veryfasttree>).