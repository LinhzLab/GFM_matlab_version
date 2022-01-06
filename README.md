# GFM
Generalized factor model for ultra-high dimensional variables with mixed types --MATLAB version

Please see our new paper for model details:

[Wei Liu, Huazhen Lin, Shurong Zheng & Jin Liu (2021) . Generalized factor model for ultra-high dimensional mixed data. Journal of the American Statistics Association (Online).](https://www.tandfonline.com/doi/abs/10.1080/01621459.2021.1999818?journalCode=uasa20)

 
## Advantage of computation
We develop a two-step iterative procedure so that each update can be carried out parallelly across all variables and samples. The fast computation version is provided for ultra-high data, see gfm_example.m file.

## R version of GFM

We also develop the corresponding R verion of GFM that is available at [here](https://github.com/LinhzLab/GFM-1). However, we suggest users to use the MATLAB version since it is more computationally efficient than R version in handling large-scale high-dimensional data.
