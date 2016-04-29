# LDA Variational EM Algorithm 

This **R** package implements the variational expectation-maximization (VEM) algorithm for the latent Dirichlet allocation (LDA) model. This VEM algorithm is for the LDA full Bayesian model. Except for this additional feature, this implementation is quite similar to the original implementation of the LDA VEM algorithm by Dr. David Blei. 

For package documentation run 

``` help("ldavem") ```

in an R console. All major functions and datasets are documented and linked to the package index. Raw data files for each dataset are available in the **data-raw** folder. To load raw data see ``` demo/load_raw_data.R  ```.    

To see all demo R scripts available in this package, run 

``` demo(package="ldavem") ```

in an R console. Some scripts can be executed via running  

``` demo(file-name, package="ldavem") ```

in an R console. The rest of them may require commandline arguments for execution. Please see the documentation provided in each script before execution.    

Authors
----------------------------
* [Clint P. George](http://www.cise.ufl.edu/~cgeorge) (Please contact for questions and comments)


Dependencies
----------------------------

This package uses the following R packages, which are already included in this R package.   
* **Rcpp**
* **RcppArmadillo** based on the **Armadillo** C++ package 
* **lattice**
* Hyperparameter optimization is adapted from Blei (2004)'s [implementation](https://github.com/blei-lab/lda-c/blob/master/lda-alpha.c)  

Installation Guide 
------------------

* Download the package source from [Git Download Link](https://github.com/clintpgeorge/ldavem/archive/master.zip)
* Unzip the dowloaded file and rename the folder **ldavem-master** to **ldavem** 
* To install **ldavem** run ```R CMD INSTALL ldavem``` on the commandline 
* To uninstall **ldavem** run ```R CMD REMOVE ldavem``` on the commandline 

References
----------

1. Blei, D. M., Ng, A. Y. and Jordan, M. I. (2003). Latent Dirichlet 
allocation. Journal of Machine Learning Research 3 993-1022.
2. Blei, D. M. (2004). [C implementation of variational EM for latent Dirichlet allocation (LDA)](https://github.com/blei-lab/lda-c) 

Acknowledgements
----------------

Clint is supported by the NIH Grant #7 R21 GM101719-03 and the University of Florida Informatics Institute.

I would like to thank Dr. George Michaildis and Wei Xia for the valuable discussions and critics that helped the development.


