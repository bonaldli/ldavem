#' @title LDA Variational EM Algorithm
#'
#' @description
#' This R package implements the variational expectation-maximization (VEM) 
#' algorithm for the latent Dirichlet allocation (LDA) model. This VEM algorithm 
#' is for the LDA full Bayesian model---Except for this, this implementation is 
#' quite similar to the original implementation of the LDA VEM algorithm by 
#' David Blei.
#' 
#' 
#' @references 
#' Latent Dirichlet Allocation. Blei, Ng, Jordan (2003) 
#' 
#' 
#' @docType package
#' 
#' @aliases
#' ldavem
#' package-ldavem
#' 
#' @useDynLib ldavem 
#' 
#' @name ldavem
#' 
#' @importFrom Rcpp evalCpp
NULL