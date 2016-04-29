#' #############################################################################
#' Runs the LDA Variational Inference Algorithm on Real Dataset 
#' 
#' See help(wt16)  
#'
#' Versions: 
#'    April 27, 2015 - Initial version 
#'    
#' Example: 
#'    Rscript run_lda_vem_synth.R
#' #############################################################################

rm(list = ls());
library(ldavem)

data(wt16)

base.alpha     <- .2
base.eta       <- .8

vi.max.iter    <- 20 # the maximum number of Gibbs iterations
em.max.iter    <- 100
vi.conv.thresh <- 1e-6 
em.conv.thresh <- 1e-4 
estimate.alpha <- 1
estimate.eta   <- 0
verbose        <- 2 
SEED           <- 1983 
K              <- 2 
V              <- length(vocab)  

fn.prefix      <- paste(
  "vi-", ds.name, "-K", K, "-h(",
  base.eta, "-", base.alpha, ")-viter", vi.max.iter, "-eiter", em.max.iter, 
  "-", format(Sys.time(), "%Y%b%d%H%M%S"), sep = ""
)
fn.prefix      <- gsub("\\.", "d", fn.prefix)
rdata.file     <- paste(fn.prefix, ".RData", sep = "")[1]



# Variational Inference  -------------------------------------------------------
set.seed(1983)
model <-
  lda_vem(
    K, V, docs, base.alpha, base.eta, vi.max.iter, em.max.iter, 
    vi.conv.thresh, em.conv.thresh, estimate.alpha, estimate.eta, verbose 
  )


# save.image(rdata.file)
# cat("\nThe R Session is saved to:", rdata.file, "\n")
# 


