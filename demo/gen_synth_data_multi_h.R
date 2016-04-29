rm(list = ls());
library(ldavem)
setwd("D:/data/lda-hp-data/synth-multi-h")

J                  <- 2
K                  <- 2
V                  <- 20
alpha.vec          <- c(.2, .8)
eta.h              <- .8
doc.size           <- 80
collection.size    <- c(100, 100) # number of documents in each collection
num.docs           <- sum(collection.size)
ds.name            <- paste("multi-h-J", J, "-K", K, "-V", V, "-D", num.docs, sep = "")
set.seed(1983)
ds                 <- gen_synth_corpus_multi_alpha(K, V, J, collection.size, doc.size, alpha.vec, eta.h)


# 
# alpha.h            <- .2
# num.docs           <- 40 # number of documents in each collection
# set.seed(1983)
# ds2                 <- gen_synth_corpus2(K, V, num.docs, doc.size, alpha.h, eta.h)




rda.file           <- paste(ds.name, ".rda", sep = "")
save(ds.name, ds, collection.size, doc.size, K, V, J, alpha.vec, eta.h, file = rda.file)
