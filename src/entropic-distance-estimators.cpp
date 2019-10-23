#include <RcppArmadillo.h>
#include <iostream>

using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List et_distance_objective(arma::vec lambda_exact
                               , arma::vec theta_extended
                               , arma::mat return_matrix
                               , double mu_penalty){
  
  // Stack return matrix twice side by side in order to deal with extended parameter vector
  arma::mat return_matrix_stacked = arma::join_rows(return_matrix, -return_matrix);
  
  // Calculate value of penalty
  double penalty = mu_penalty * arma::accu(theta_extended);
  
  // Calculate returns on exact SDF portfolio
  arma::mat lambda_return = return_matrix * lambda_exact;
  
  // Caclulate expected utility of exact SDF portfolio
  double objective = - arma::as_scalar(arma::mean(arma::exp(lambda_return)));
  
  // Add (subtract, actually) contribution of expected utility of approximate SDF portfolio times exact return
  objective = objective - arma::as_scalar(arma::mean((arma::exp(return_matrix_stacked * theta_extended) - arma::ones(lambda_return.size())) % lambda_return));
  
  // Add penalty
  objective = objective + penalty;
  
  // Switch sign (this is the dual problem and normally you would maximize, but a solver minimizes)
  objective = -objective;
  
  //// Prepare gradient calculation
  // Caclulate exact SDF from portfolio
  // No minus here!
  arma::mat sdf_lambda = arma::exp(lambda_return);
  // Price assets with exact SDF
  arma::mat returns_prices_exact = return_matrix.each_col() % sdf_lambda;
  
  // Calculate model SDF
  arma::mat sdf_theta = arma::exp(return_matrix_stacked * theta_extended) - arma::ones(lambda_return.size());
  arma::mat returns_prices_model = return_matrix.each_col() % sdf_theta;
  
  // Average priced returns
  arma::vec gradient = arma::mean(returns_prices_exact, 0).t() + arma::mean(returns_prices_model, 0).t();
  
  // Return
  return Rcpp::List::create(Rcpp::Named("objective") = objective
                              , Rcpp::Named("gradient") = gradient
                              , Rcpp::Named("other") = 2.0);
}


// Currently unnecessary
arma::vec et_distance_lambda_gradient(arma::vec lambda_exact
                               , arma::vec theta_extended
                               , arma::mat return_matrix
                               , double mu_penalty){

  // Stack return matrix twice side by side in order to deal with extended parameter vector
  arma::mat return_matrix_stacked = arma::join_rows(return_matrix, -return_matrix);
  
  // Calculate returns on exact SDF portfolio
  arma::mat lambda_return = return_matrix * lambda_exact;
  
  // Caclulate exact SDF from portfolio
  // No minus here!
  arma::mat sdf_lambda = arma::exp(lambda_return);
  // Price assets with exact SDF
  arma::mat returns_prices_exact = return_matrix.each_col() % sdf_lambda;
  
  // Calculate model SDF
  arma::mat sdf_theta = arma::exp(return_matrix_stacked * theta_extended) - arma::ones(lambda_return.size());
  arma::mat returns_prices_model = return_matrix.each_col() % sdf_theta;
  
  // Average priced returns
  arma::vec res = arma::mean(returns_prices_exact, 0).t() + arma::mean(returns_prices_model, 0).t();
    
  return res;
} 
  

// Currently unnecessary
arma::mat et_distance_lambda_hessian(arma::vec lambda_exact
                              , arma::vec theta_extended
                              , arma::mat return_matrix
                              , double mu_penalty){
  
  // Calculate returns on exact SDF portfolio
  arma::mat lambda_return = return_matrix * lambda_exact;
  
  // Caclulate exact SDF from portfolio, to power 1/2
  // No minus here!
  arma::mat sdf_lambda = arma::exp(0.5 * lambda_return);
  
  return_matrix.each_col() %= sdf_lambda;
  
  arma::mat res = return_matrix.t() * return_matrix;
  
  res = res * 1.0/lambda_return.n_rows;
  
  return res;
}

// [[Rcpp::export]]
arma::mat et_distance_theta_gradient(arma::vec lambda_opt
                                       , arma::vec theta_extended
                                       , arma::mat return_matrix
                                       , double mu_penalty){
  
  // lambda_opt comes from solving the inner problem in R
  
  // Stack return matrix twice side by side in order to deal with extended parameter vector
  arma::mat return_matrix_stacked = arma::join_rows(return_matrix, -return_matrix);
  
  // Calculate returns on exact SDF portfolio
  arma::mat lambda_return = return_matrix * lambda_opt;
  
  // Calculate model SDF
  arma::mat sdf_theta = arma::exp(return_matrix_stacked * theta_extended);
  
  return_matrix_stacked.each_col() %= sdf_theta % lambda_return;
  
  arma::vec res = arma::mean(return_matrix_stacked, 0).t();
    
  // Invert sign, because you will be minimising and you inverted the sign in the inner function
  res = -res + mu_penalty * arma::ones(res.size());
    
  return res;
}

// Currently unnecessary
arma::mat et_distance_theta_hessian(arma::vec lambda_opt
                                       , arma::vec theta_extended
                                       , arma::mat return_matrix
                                       , double mu_penalty){
  
  // lambda_opt comes from solving the inner problem in R
  
  // Stack return matrix twice side by side in order to deal with extended parameter vector
  arma::mat return_matrix_stacked = arma::join_rows(return_matrix, -return_matrix);
  
  // Calculate returns on exact SDF portfolio
  arma::mat lambda_return = return_matrix * lambda_opt;
  
  // Calculate model SDF 
  arma::mat sdf_theta = arma::exp(return_matrix_stacked * theta_extended);
  
  arma::mat return_matrix_stacked_extended = return_matrix_stacked;
  return_matrix_stacked_extended.each_col() %= sdf_theta % lambda_return;
  
  // initialise cube for holding r * r^T * e^(theta^T r) * lambda^T r
  arma::cube variance_cube_extended(theta_extended.n_rows, theta_extended.n_rows, return_matrix.n_rows);
  
  // create temp matrix for holding intermediate results
  arma::mat temp_res(theta_extended.n_rows, theta_extended.n_rows);
  cout << temp_res.size() << "\n";
  // iterate over number of obs
  for(unsigned int obs_iter = 0; obs_iter < return_matrix.n_rows; obs_iter++){
    temp_res = return_matrix_stacked.row(obs_iter).t() * return_matrix_stacked_extended.row(obs_iter);
    // Assign
    variance_cube_extended.slice(obs_iter) = temp_res;
  }
  
  arma::mat res = arma::mean(variance_cube_extended, 2);
  
  return res;
}