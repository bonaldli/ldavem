# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>

using namespace Rcpp ;
using namespace std ;
using namespace arma ;


///////////////////////////////////////////////////////////////////////////////
// Helper functions
///////////////////////////////////////////////////////////////////////////////

/**
 * Extend division reminder to vectors
 *
 * @param   a       Dividend 
 * @param   n       Divisor
 */
template<typename T>
T mod(T a, int n)
{
  return a - floor(a/n)*n;
}   

/**
* Samples an integer from [0, K) uniformly at random
*
* Arguments:
* 		K - the upper interval
* Returns:
* 		the sampled integer
*/
unsigned int sample_uniform_int(unsigned int K){
  return (unsigned int) (runif(1)(0) * (double)K); // To speedup
}

//' A speedy sampling from a multimomial distribution
//'
//' @param theta a multinomial probability vector (K x 1 vector)
//'
//' @return returns a class index from [0, K)
//'
//' @note
//' Author: Clint P. George
//'
//' Created on: February 11, 2016
//'
//' @family utils
//'
//' @export
// [[Rcpp::export]]
unsigned int sample_multinomial (arma::vec theta) {
  
  unsigned int t = 0;
  double total_prob = accu(theta);
  double u = runif(1)(0) * total_prob;
  double cumulative_prob = theta(0);
  
  while(u > cumulative_prob){
    t++;
    cumulative_prob += theta(t);
  }
  
  return t;
  
}


//' Samples from a Dirichlet distribution given a hyperparameter
//'
//' @param num_elements the dimention of the Dirichlet distribution
//' @param alpha the hyperparameter vector (a column vector)
//'
//' @return returns a Dirichlet sample (a column vector)
//'
//' @note
//' Author: Clint P. George
//'
//' Created on: 2014
//'
//' @family utils
//'
//' @export
// [[Rcpp::export]]
arma::vec sample_dirichlet(unsigned int num_elements, arma::vec alpha){
  
  arma::vec dirichlet_sample = arma::zeros<arma::vec>(num_elements);
  
  for ( register unsigned int i = 0; i < num_elements; i++ )
    dirichlet_sample(i) = rgamma(1, alpha(i), 1.0)(0); // R::rgamma(1, alpha(i));
  
  dirichlet_sample /= accu(dirichlet_sample);
  
  return dirichlet_sample;
  
}

/**
* Samples from a Dirichlet distribution given a hyperparameter
*
* Aruguments:
* 		num_elements - the dimention of the Dirichlet distribution
* 		alpha - the hyperparameter vector (a column vector)
* Returns:
* 		the Dirichlet sample (a column vector)
*/
arma::rowvec sample_dirichlet_row_vec (unsigned int num_elements, arma::rowvec alpha){
  
  arma::rowvec dirichlet_sample = arma::zeros<arma::rowvec>(num_elements);
  
  for ( register unsigned int i = 0; i < num_elements; i++ )
    dirichlet_sample(i) = rgamma(1, alpha(i), 1.0)(0); // R::rgamma(1, alpha(i));
  
  dirichlet_sample /= accu(dirichlet_sample);
  
  return dirichlet_sample;
  
}


arma::vec log_gamma_vec(arma::vec x_vec){
  
  arma::vec lgamma_vec = arma::zeros<arma::vec>(x_vec.n_elem);
  
  for (unsigned int i = 0; i < x_vec.n_elem; i++)
    lgamma_vec(i) = lgamma(x_vec(i));
  
  return lgamma_vec;
  
}

arma::vec digamma_vec(arma::vec x_vec){
  // digamma(wrap()) will do, with comparable performance 
  arma::vec ret = arma::zeros<arma::vec>(x_vec.n_elem);
  for (unsigned int i = 0; i < x_vec.n_elem; i++)
    ret(i) = Rf_digamma(x_vec(i));
  return ret;
}

arma::rowvec digamma_rowvec(arma::rowvec x_vec){
  arma::rowvec ret = arma::zeros<arma::rowvec>(x_vec.n_elem);
  for (unsigned int i = 0; i < x_vec.n_elem; i++)
    ret(i) = Rf_digamma(x_vec(i));
  return ret;
}



arma::vec gamma_col_vec(arma::vec x_vec){
  // It took 2hrs in the April 19, 2014 morning to make this function work. 
  // The main issue was with accessing the R gamma function from the
  // RcppArmadillo namespace. See
  // http://dirk.eddelbuettel.com/code/rcpp/html/Rmath_8h_source.html
  // gamma(as<NumericVector>(wrap(x_vec))) is another option, but it seems to be
  // slow. See
  // http://stackoverflow.com/questions/14253069/convert-rcpparmadillo-vector-to-rcpp-vector
  
  arma::vec gamma_vec = arma::zeros<arma::vec>(x_vec.n_elem);
  for (unsigned int i = 0; i < x_vec.n_elem; i++)
    gamma_vec(i) = Rf_gammafn(x_vec(i));
  return gamma_vec;
}



