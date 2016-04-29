#' #############################################################################
#' Runs the LDA Variational Inference Algorithm on synthetic dataset
#' 
#' See gen_synth_data_multi_h.R to see how to generate one such dataset with 
#' a different configuration   
#'
#' Versions: 
#'    April 27, 2015 - Initial version 
#'    
#' Example: 
#'    Rscript run_lda_vi_synth.R
#' #############################################################################

rm(list = ls());
library(ldavem)

setwd('~') # sets the working directory 
prefix <- "multi-h-J2-K2-V20"
data("multi-h-J2-K2-V20")

base.alpha     <- .2
base.eta       <- .8

vi.max.iter    <- 100 # the maximum number of Gibbs iterations
em.max.iter    <- 10
vi.conv.thresh <- 1e-6 
em.conv.thresh <- 1e-4 
estimate.alpha <- 1
estimate.eta   <- 1
verbose        <- 2 
SEED           <- 2008 


fn.prefix      <- paste(
  "vi-", prefix, "-K", K, "-h(",
  base.eta, "-", base.alpha, ")-viter", vi.max.iter, "-eiter", em.max.iter, 
  "-", format(Sys.time(), "%Y%b%d%H%M%S"), sep = ""
)
fn.prefix      <- gsub("\\.", "d", fn.prefix)
rdata.file     <- paste(fn.prefix, ".RData", sep = "")[1]



# Variational Inference  -------------------------------------------------------

cat("Variational Inference...\n\n")
ptm            <- proc.time()
set.seed(SEED)
model          <-
  lda_vem( 
    K, V, ds$docs, base.alpha, base.eta, vi.max.iter, em.max.iter, 
    vi.conv.thresh, em.conv.thresh, estimate.alpha, estimate.eta, verbose 
  )
gs_ptm         <- proc.time() - ptm
cat("\nExecution time = ", gs_ptm[3], "\n\n")


# save.image(rdata.file)
# cat("\nThe R Session is saved to:", rdata.file, "\n")
# 


