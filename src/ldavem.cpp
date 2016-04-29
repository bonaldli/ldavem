# include "ldavem.h"
# include "opt-hp.h"
# include <assert.h>


//' LDA: Variational EM   
//'
//' This implements the variational inference algorithm for the LDA (full 
//' Bayesian) model. This includes optimization routines for both \eqn{\alpha} 
//' and \eqn{\eta} hyperparameters. 
//' 
//' References: 
//'   * Latent Dirichlet Allocation. D. Blei, A. Ng, M.I. Jordan (2003)    
//'
//' @param num_topics Number of topics in the corpus
//' @param vocab_size  Vocabulary size
//' @param docs_tf A list of corpus documents read from the Blei corpus using 
//'   \code{\link{read_docs}} (term indices starts with 0)
//' @param alpha_h Hyperparameter for \eqn{\theta} sampling
//' @param eta_h Smoothing parameter for the \eqn{\beta} matrix
//' @param vi_max_iter Maximum number of iterations for variational inference
//' @param em_max_iter Maximum number of iterations for variational EM
//' @param vi_conv_threshold Convergence threshold for the document variational
//'   inference loop
//' @param em_conv_threshold Convergence threshold for the variational EM loop
//' @param estimate_alpha If true, run hyperparameter \eqn{\alpha} optimization
//' @param estimate_eta If true, run hyperparameter \eqn{\eta} optimization
//' @param verbose from {0, 1, 2}
//'
//' @return TBA 
//'
//' @export
//'
//' @family Variational Inference  
//' 
//' @note Created on April 26, 2016 
//'
// [[Rcpp::export]]
List lda_vem(
    unsigned int num_topics,
    unsigned int vocab_size,
    List docs_tf,
    double alpha_h,
    double eta_h,
    unsigned int vi_max_iter,
    unsigned int em_max_iter,
    double vi_conv_thresh, 
    double em_conv_thresh, 
    bool estimate_alpha, 
    bool estimate_eta,  
    int verbose
) {
  
  unsigned int num_docs = docs_tf.size(); // number of documents
  unsigned int d, em_iter, vi_iter, c, word_id, word_count, k;
  unsigned int num_words = 0; // number of words in the corpus
  unsigned int num_uwords; 
  double alpha_ss = 0; 
  double eta_ss = 0;  
  double alpha_ss_d; 
  double doc_vi_lb; 
  double doc_vi_cr; 
  double doc_vi_lb_old;
  double em_conv_ratio; 
  double em_lb_old; 
  double em_lb_current;

  vector < vector < unsigned int > > doc_uword_ids; // unique word ids 
  vector < vector < unsigned int > > doc_uword_counts; // unique word counts 
  vector < double > vi_lb; // variational lower-bound 
  
  vec phi_ui; 
  vec dig_digsum; 
  vec vi_gamma_ss; 
  vec vi_gamma_new; 
  uvec doc_word_counts = zeros<uvec>(num_docs); // doc lengths
  mat vi_mu = zeros<mat>(num_topics, vocab_size); // K x V matrix
  mat vi_gamma = zeros<mat>(num_topics, num_docs); // K x D matrix

  
  cout << endl << endl;
  if (verbose > 0){
    cout << "lda-vem (c++): Initializes variables and count statistics....";
  }

  //////////////////////////////////////////////////////////////////////////////
  // Reads corpus statistics  
  //////////////////////////////////////////////////////////////////////////////
  for (d = 0; d < num_docs; d++){
    umat document = as<umat>(docs_tf(d));
    vector < unsigned int > uword_ids;
    vector < unsigned int > uword_counts;

    for (c = 0; c < document.n_cols; c++){
      word_id = document(0,c);
      word_count = document(1,c);
      num_words += word_count; // increments number of words in the corpus
      doc_word_counts(d) += word_count; // increments doc word counts
      uword_ids.push_back(word_id); // saves unique word ids 
      uword_counts.push_back(word_count); // saves unique word counts 
    }

    doc_uword_ids.push_back(uword_ids);
    doc_uword_counts.push_back(uword_counts); 
  }

  if (verbose > 0){ cout << "DONE" << endl; }

  if (verbose > 1){
    cout << "lda-vem (c++): Number of docs: " << num_docs << endl;
    cout << "lda-vem (c++): Number of total words: " << num_words << endl;
    cout << "lda-vem (c++): Number of topics: " << num_topics << endl;
    cout << "lda-vem (c++): Vocabulary size: " << vocab_size << endl;
    cout << "lda-vem (c++): alpha_h: " << alpha_h << endl;
    cout << "lda-vem (c++): eta_h: " << eta_h << endl;
  }
  
  if (verbose > 0){
    cout << "lda-vem (c++): Variational-EM..." << endl << endl;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  // Initializes variational-EM parameters 
  //////////////////////////////////////////////////////////////////////////////
  
  // vi_mu - the variational Dirichlet parameter for \beta 
  // random initiliazation, adapted from 
  //     void random_initialize_ss(lda_suffstats* ss, lda_model* model) 
  // given in the LDA-C package  
  vi_mu.randu(); 
  vi_mu += (1.0 / (double) vocab_size); 
  vi_mu.each_col() /= sum(vi_mu, 1); // normalizes each row to sum to 1 
  
  em_iter = 0;
  em_conv_ratio = 1.0; 
  em_lb_old = 0; 
  cout.precision(10);

  while (((em_conv_ratio < 0) || (em_conv_ratio > em_conv_thresh) || (em_iter <= 2)) && (em_iter <= em_max_iter)) {
    
    em_iter++;
    
    ////////////////////////////////////////////////////////////////////////////
    // E Step 
    ////////////////////////////////////////////////////////////////////////////
    
    if (verbose > 1) { 
      cout << "lda-vem (c++): em_iter #" << em_iter; 
      cout << " alpha: " << alpha_h << " eta: " << eta_h << endl << endl; 
    } else if (verbose > 0) { 
      cout << endl; 
    }
    
    // resets variational Dirichlets in each EM iteration 
    // vi_gamma - the variational Dirichlet parameter for \theta 
    // vi_mu - the variational Dirichlet parameter for \beta 
    
    vi_gamma.fill(alpha_h); 
    for (d = 0; d < num_docs; d++){
      vi_gamma.col(d) += ((double) doc_word_counts(d) / (double) num_topics); 
    }
    mat vi_mu_new = zeros<mat>(num_topics, vocab_size); // K x V matrix
    vi_mu_new.fill(eta_h); 
    
    
    em_lb_current = 0.0; // corpus variational lowerbound 
    alpha_ss = 0; // alpha sufficient statistics 
    
    for (d = 0; d < num_docs; d++) { // for each document
      
      vector < unsigned int > uword_ids = doc_uword_ids[d];
      vector < unsigned int > uword_counts = doc_uword_counts[d];
      num_uwords = uword_ids.size(); 
      
      // Variational Inference 
      
      mat phi_d = zeros<mat>(num_topics, num_uwords);
      alpha_ss_d = 0; 
      doc_vi_cr = 1.0; 
      doc_vi_lb = 0; 
      doc_vi_lb_old = 0; 
      vi_iter = 0; 
      
      while ((doc_vi_cr > vi_conv_thresh) && (vi_iter < vi_max_iter)) { // for each VI iteration 
        
        vi_iter++;  
        if (verbose > 2) { 
          cout << "lda-vem (c++): doc #" << (d + 1) << " vi_iter # " << vi_iter; 
        } else if (verbose > 1) { 
          cout << "."; 
        }
        
        // vi updates 
        
        vi_gamma_new = ones<vec>(num_topics) * alpha_h; 
        vi_gamma_ss = exp(digamma_vec(vi_gamma.col(d)) - Rf_digamma(sum(vi_gamma.col(d)))); 
        
        for (c = 0; c < num_uwords; c++){ // for each unique word 
          word_id = uword_ids[c];
          word_count = uword_counts[c]; 
          
          phi_ui = vi_mu.col(word_id) % vi_gamma_ss;   
          phi_ui /= sum(phi_ui); // normalize to sum to 1. 
          phi_d.col(c) = phi_ui; 
          
          vi_gamma_new += word_count * phi_ui;   
        } // for each unique word 
        
        vi_gamma.col(d) = vi_gamma_new; 
        
        // computes log document d's variational lower-bound 
        // Reference: Blei, Ng, Jordan (2003), Appendix A.3, Equation 15.  
        
        dig_digsum = digamma_vec(vi_gamma_new) - Rf_digamma(sum(vi_gamma_new)); // K x 1 vector 
        alpha_ss_d = sum(dig_digsum); 
        doc_vi_lb = lgamma(alpha_h * num_topics) - num_topics * lgamma(alpha_h) + (alpha_h - 1.0) * alpha_ss_d; // 15.a 
        doc_vi_lb += (sum(log_gamma_vec(vi_gamma_new)) - lgamma(sum(vi_gamma_new)) - sum((vi_gamma_new - 1.0) % dig_digsum)); // 15.d 
        for (c = 0; c < num_uwords; c++) { // for each unique word 
          word_id = uword_ids[c];
          word_count = uword_counts[c];
          phi_ui = phi_d.col(c); 
          doc_vi_lb += (sum(phi_ui % dig_digsum) * word_count); // 15.b 
          doc_vi_lb += (sum(phi_ui % log(vi_mu.col(word_id))) * word_count); // 15.c
          doc_vi_lb += (sum(phi_ui % log(phi_ui)) * word_count); // 15.e
        }
        
        // checks for convergence
        
        assert(!isnan(doc_vi_lb));
        if (verbose > 2){ cout << " vi-lb: " << doc_vi_lb << " alpha-ss: " << alpha_ss_d; }
        if (vi_iter > 1){
          doc_vi_cr = (doc_vi_lb_old - doc_vi_lb) / doc_vi_lb_old;
          if (verbose > 2){ cout << " conv-ratio: " << doc_vi_cr; }
        } 
        doc_vi_lb_old = doc_vi_lb; 
        
        
        if (verbose > 2){ cout << endl; }
      
      
      } // for each VI iteration
      

      em_lb_current += doc_vi_lb; 
      alpha_ss += alpha_ss_d; // computes \alpha sufficient statistics 
      
      for (c = 0; c < num_uwords; c++) {
        vi_mu_new.col(uword_ids[c]) += uword_counts[c] * phi_d.col(c);  
      }
      
    } // for each document
    
    if (verbose > 1){ cout << endl; }
      

    ////////////////////////////////////////////////////////////////////////////
    // M Step 
    ////////////////////////////////////////////////////////////////////////////
    
    // updates variational Dirichlet parameter \mu, for topics (i.e. \beta_k's) 
    
    vi_mu = vi_mu_new; 
    vi_mu.each_col() /= sum(vi_mu, 1); // normalizes each row to sum to 1 
    
    // computes \eta sufficient statistics 
    
    eta_ss = 0;  
    for (k = 0; k < num_topics; k++) { // for each topic
      eta_ss += sum(digamma_rowvec(vi_mu.row(k)) - Rf_digamma(sum(vi_mu.row(k)))); 
    }
    
    // hyperparameter optimization for \alpha and \eta   
        
    if(estimate_alpha){
      alpha_h = opt_hp(100, alpha_ss, num_docs, num_topics);
      if (verbose > 1){ cout << endl; }
    }
  
    if(estimate_eta){
      eta_h = opt_hp(100, eta_ss, num_topics, vocab_size);
      if (verbose > 1){ cout << endl; }
    }
    

    ////////////////////////////////////////////////////////////////////////////
    // EM: check for convergence
    ////////////////////////////////////////////////////////////////////////////
    
    if (verbose > 1){  
      cout << "lda-vem (c++): em_iter #" << em_iter; 
      cout << " vi-lb: " << em_lb_current; 
      cout << " opt-alpha: " << alpha_h; 
      cout << " opt-eta: " << eta_h; 
    }
    if (em_iter > 1){
      em_conv_ratio = (em_lb_old - em_lb_current) / em_lb_old;
      if (verbose > 1){ cout << " em-conv-ratio: " << em_conv_ratio; }
      if (em_conv_ratio < 0.0) { vi_max_iter *= 2; }
    }
    if (verbose > 1){  cout << endl << endl; }
    
    em_lb_old = em_lb_current;
    vi_lb.push_back(em_lb_current); 

  
  } // for each EM iteration
  

  
  if (verbose > 0){
    cout << "lda-vem (c++): Completed Variational-EM." << endl;
  }
  
  return List::create(
    Named("vi_mu") = wrap(vi_mu), 
    Named("vi_gamma") = wrap(vi_gamma), 
    Named("alpha_h") = wrap(alpha_h), 
    Named("eta_h") = wrap(eta_h)
  );
  
}


