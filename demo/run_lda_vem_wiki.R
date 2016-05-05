#' #############################################################################
#' Runs the LDA Variational Inference Algorithm on Wikipedia Datasets  
#' 
#' See help(wt16)  
#'
#' Versions: 
#'    May 05, 2016   - Tested on various Wikipedia datasets   
#'    April 27, 2015 - Initial version 
#'    
#' Example: 
#'    Rscript run_lda_vem_synth.R
#' #############################################################################

rm(list = ls());
library(ldavem)

data(bop) # change to appropriate dataset 

base.alpha     <- .1
base.eta       <- .1

vi.max.iter    <- 1 # the maximum number of Gibbs iterations
em.max.iter    <- 1
vi.conv.thresh <- 1e-6 
em.conv.thresh <- 1e-4 
estimate.alpha <- 0
estimate.eta   <- 0
verbose        <- 2 
SEED           <- 1983 
K              <- length(class.labels)  
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


# C-6: cats 
# opt-alpha: 0.2096168387 opt-eta: 1.146523621
# 
# C-7: felines 
# em_iter #54 vi-lb: -77080.55634 opt-alpha: 0.0562323916 opt-eta: 0.7102011106
# 
# C-8: Canis 
# em_iter #34 vi-lb: -54153.80766 opt-alpha: 0.1625960584 opt-eta: 0.8411708841
# 
# C-9: Birds of Prey 
# 

