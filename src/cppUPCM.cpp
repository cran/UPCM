
// [[Rcpp::depends(RcppArmadillo)]]
// #define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <Rcpp.h>

#include <vector>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(openmp)]]

// square root of a matrix
mat matsqrt2(mat A){
  vec eigval;
  mat eigvec;
  
  eig_sym(eigval, eigvec, A);
  
  colvec d = (eigval+abs(eigval))*0.5;
  colvec d2 = sqrt(d);
  mat B = eigvec*diagmat(d2)*trans(eigvec);
  
  return B;
}

// // response function for acat
// vec responseFun(vec eta){
//   int q = eta.n_rows;
//   // eta(find(eta>10)) = ones(size(find(eta>10)))*10;
//   // eta(find(eta<-10)) = -ones(size(find(eta<-10)))*10;
//   mat eta_help1 = ones(q+1)*trans(join_cols(zeros(1),eta));
//   eta_help1 = eta_help1 % trimatl(ones(eta_help1.n_rows,eta_help1.n_cols));
//   vec pi_vec = ones(q)/as_scalar(sum(prod(exp(eta_help1),1)));
//   for(int k=1; k<q ;k++){
//     pi_vec(k) =pi_vec(k-1)*exp(eta(k-1));
//   }
//   vec pi_vec2 = (pi_vec-0.5)*0.99999+0.5;
//   if(prod(pi_vec2)!=prod(pi_vec2)){
//     eta(find(eta>10)) = ones(size(find(eta>10)))*10;
//     eta(find(eta<-10)) = -ones(size(find(eta<-10)))*10;
//     eta_help1 = ones(q+1)*trans(join_cols(zeros(1),eta));
//     eta_help1 = eta_help1 % trimatl(ones(eta_help1.n_rows,eta_help1.n_cols));
//     pi_vec = ones(q)/as_scalar(sum(prod(exp(eta_help1),1)));
//     for(int k=1; k<q ;k++){
//       pi_vec(k) =pi_vec(k-1)*exp(eta(k-1));
//     }
//     pi_vec2 = (pi_vec-0.5)*0.99999+0.5;
//     Rcout<<trans(pi_vec2)<<endl;
//     Rcout<<trans(eta)<<endl;
//     Rcout<<eta_help1<<endl;
//   }
//   // pi_vec = pi_vec2/sum(pi_vec2)*sum(pi_vec);
//   return pi_vec2;
// }


// response function for acat
vec responseFun(vec eta){
  int q = eta.n_rows;
  // eta(find(eta>10)) = ones(size(find(eta>10)))*10;
  // eta(find(eta<-10)) = -ones(size(find(eta<-10)))*10;
  mat eta_help1 = ones(q+1)*trans(join_cols(zeros(1),eta));
  eta_help1 = eta_help1 % trimatl(ones(eta_help1.n_rows,eta_help1.n_cols));
  vec pi_vec = ones(q)/as_scalar(sum(prod(exp(eta_help1),1)));
  for(int k=1; k<q ;k++){
    pi_vec(k) =pi_vec(k-1)*exp(eta(k-1));
  }
  vec pi_vec2 = (pi_vec-0.5)*0.99999+0.5;
  // pi_vec = pi_vec2/sum(pi_vec2)*sum(pi_vec);
  return pi_vec2;
}

// response function for acat
vec responseFun2(vec eta){
  int q = eta.n_rows;
  // eta(find(eta>10)) = ones(size(find(eta>10)))*10;
  // eta(find(eta<-10)) = -ones(size(find(eta<-10)))*10;
  mat eta_help1 = ones(q+1)*trans(join_cols(zeros(1),eta));
  eta_help1 = eta_help1 % trimatl(ones(eta_help1.n_rows,eta_help1.n_cols));
  vec pi_vec = ones(q)/as_scalar(sum(prod(exp(eta_help1),1)));
  for(int k=1; k<q ;k++){
    pi_vec(k) =pi_vec(k-1)*exp(eta(k-1));
  }
  return pi_vec;
}

// create inverse sigma for acat
mat createSigmaInv(vec mu){
  mat Sigma = diagmat(mu) - mu * trans(mu);
  mat SigmaInv;
  try{
    SigmaInv = inv(Sigma);
  }
  catch(...){
    Rcout<<"Sigma"<<endl;
    Rcout<<mu<<endl;
    Rcout<<Sigma<<endl;
    SigmaInv = pinv(Sigma);
  }
  return SigmaInv;
}

// create derivative matrix for acat
mat createD(vec mu){
  int q = mu.n_rows;
  mat mu_k = join_cols(mu,1-sum(mu,1));
  mat D2 = zeros(q,q) - diagmat(1/mu);
  
  if(q==2){
    D2(1,0) = 1/mu(1);
  }else{
    D2(span(1,q-1),span(0,q-2)) = D2(span(1,q-1),span(0,q-2)) + diagmat(1/mu(span(1,q-1)));
  }

  D2(span::all,q-1) = -ones(q,1)/(1-as_scalar(sum(mu)));
  D2(q-1,q-1) = as_scalar(-(1-sum(mu(span(0,q-2))))/((1-sum(mu))*mu(q-1)));
  
  mat D;
  try{
    D = inv(D2);
  }
  catch(...){
    Rcout<<"D"<<endl;
    Rcout<<mu<<endl;
    Rcout<<D2<<endl;
    D = pinv(D2);
  }
  return D;
}

 // [[Rcpp::export]]
  double loglikUPCM(arma::vec alpha,
           arma::vec Y,
           int Q,
           int q,
           int n,
           int I,
           int pall,
           arma::mat GHweights,
           arma::vec GHnodes,
           int pX,
           arma::mat X,
           int cores,
           double lambda) { 
      
      // initialize loglikelihood       
      vec f = zeros(n);   

    double P2 = accu(alpha%alpha)*lambda;
    
      // initialize design for one person, without random effects
      mat Z = -join_rows(diagmat(ones(q*I)),zeros(q*I,2*pX+3));
      mat etai = Z*alpha;
      // current linear predictor without random effects, per item and category
      etai = reshape(etai,q,I);

      // create sigma matrix from current parameters
      double co_var = alpha(pall-2)*sqrt(alpha(pall-1))*sqrt(alpha(pall-3));
      mat sigma = zeros(2,2);        
      sigma(0,0) = alpha(pall-3);
      sigma(1,0) = co_var;
      sigma(0,1) = co_var;
      sigma(1,1) = alpha(pall-1);
      // square root of current sigma matrix
      
      // mat sigma12 = trans(chol(sigma));
      mat sigma12;
      try{
        sigma12 = trans(chol(sigma));
      }
      catch(...)
      {
        sigma = sigma + diagmat(ones(2)*0.0001);
        sigma12 = trans(chol(sigma));
      }

      vec betaX = alpha(span(q*I,q*I+pX-1));
      vec betaU = alpha(span(q*I+pX,q*I+2*pX-1));

  // initialize number of different threads
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
  
  mat Xi; vec yi; mat prods_i; int j; int jj; int k; int i; 
  vec rnd_act; mat etaij; mat mui; vec eta_k; vec mu_k;
  vec yi_k; 
  // loop through all persons
#pragma omp parallel for private(i, j, jj, k, yi, yi_k, Xi, prods_i, rnd_act, etaij, mui, eta_k, mu_k) shared(f)
      for(i=0; i<n ;i++){

        Xi = X(i,span::all);
        
        // get response of person i
        yi = Y(span(i*q*I,i*q*I+q*I-1));
        prods_i = ones(Q,Q);
          // loop through all knots, for both random effects
          for(j=0; j<Q ;j++){ 
            for(jj=0; jj<Q ;jj++){ 
              // initialize current random effects or knots
              rnd_act = zeros(2);
              rnd_act(0) = GHnodes(j);
              rnd_act(1) = GHnodes(jj);
              
              rnd_act = sigma12*rnd_act+join_cols(Xi*betaX,Xi*betaU);
              

              // loop through all items  
              for(k=0; k<I ;k++){  
                // response of person i for current item k
                yi_k = yi(span(k*q,q+k*q-1));
                yi_k = join_cols(yi_k,1-sum(yi_k,0));
                
                // get eta and mu of person i for current item k
                eta_k = (etai(span::all,k) + ones(q)*(rnd_act(0)))%(ones(q)*(exp(rnd_act(1))));
                mu_k = responseFun(eta_k);
                mu_k = join_cols(mu_k,1-sum(mu_k,0));
                // if(prod(mu_k % yi_k - (yi_k-1))<0){
                //   Rcout<<trans(mu_k)<<endl;
                // }
                // create prob for current item, update prob matrix for respective knots
                prods_i(j,jj) = prods_i(j,jj)*prod(mu_k % yi_k - (yi_k-1));
              }
            }
          }
          f(i)  = -log(accu(prods_i%GHweights));
          
         // accumulate all likelihood contributions, weights by respective weights and probs of knots
      }
  return sum(f)+P2;
}

// [[Rcpp::export]]
double loglikUPCM2(arma::vec alpha,
                  arma::mat Y,
                  int Q,
                  int q,
                  int n,
                  int I,
                  int pall,
                  arma::mat GHweights,
                  arma::vec GHnodes,
                  int pX,
                  arma::mat X,
                  int cores,
                  double lambda) { 
  
  
  // initialize loglikelihood       
  vec f = zeros(n);   
  
  double P2 = accu(alpha%alpha)*lambda;
  
  // initialize design for one person, without random effects
  mat Z = -join_rows(diagmat(ones(q*I)),zeros(q*I,2*pX+3));
  mat etai = Z*alpha;
  // current linear predictor without random effects, per item and category
  etai = reshape(etai,q,I);
  
  // create sigma matrix from current parameters
  double co_var = alpha(pall-2)*sqrt(alpha(pall-1))*sqrt(alpha(pall-3));
  mat sigma = zeros(2,2);        
  sigma(0,0) = alpha(pall-3);
  sigma(1,0) = co_var;
  sigma(0,1) = co_var;
  sigma(1,1) = alpha(pall-1);
  // square root of current sigma matrix
  
  mat sigma12;
  try{
    sigma12 = trans(chol(sigma));
  }
  catch(...)
  {
    sigma = sigma + diagmat(ones(2)*0.0001);
    sigma12 = trans(chol(sigma));
  }
  
  
  vec betaX = alpha(span(q*I,q*I+pX-1));
  vec betaU = alpha(span(q*I+pX,q*I+2*pX-1));
  
  // initialize number of different threads
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
  
  mat Xi; mat yi; mat prods_i; int j; int jj; int k; int i; 
  vec rnd_act; mat etaij; mat mui; vec eta_k; vec mu_k;
  // loop through all persons
#pragma omp parallel for private(i, j, jj, k, yi, Xi, prods_i, rnd_act, etaij, mui, eta_k, mu_k) shared(f)
  for(i=0; i<n ;i++){
    Xi = X(i,span::all);
    
    // get response of person i
    yi = Y(span::all,span(i*I,i*I+I-1));

    prods_i = ones(Q,Q);
    // loop through all knots, for both random effects
    for(j=0; j<Q ;j++){ 
      for(jj=0; jj<Q ;jj++){ 
        // initialize current random effects or knots
        rnd_act = zeros(2);
        rnd_act(0) = GHnodes(j);
        rnd_act(1) = GHnodes(jj);

        rnd_act = sigma12*rnd_act+join_cols(Xi*betaX,Xi*betaU);
        
        
        etaij = (etai + ones(q,I)*rnd_act(0))%(ones(q,I)*exp(rnd_act(1)));
         mui = ones(q+1,I);
        // loop through all items  
        for(k=0; k<I ;k++){
          
          // get eta and mu of person i for current item k
          eta_k = etaij(span::all,k);
          mu_k = responseFun2(eta_k);
          mui(span::all,k) = join_cols(mu_k,1-sum(mu_k,0));
        }
          // create prob for current item, update prob matrix for respective knots
          prods_i(j,jj) = prods_i(j,jj)*prod(prod(mui % yi - (yi-1)));
        
      }
    }
    f(i)  = - log(accu(prods_i%GHweights));
    // accumulate all likelihood contributions, weights by respective weights and probs of knots
  }
  
  
  return sum(f)+P2;
}

// // [[Rcpp::export]]
// double loglikUPCM3(arma::vec alpha,
//                    arma::mat Y,
//                    int Q,
//                    int q,
//                    int n,
//                    int I,
//                    arma::mat GHweights,
//                    arma::vec GHnodes,
//                    int scaled,
//                    arma::mat X,
//                    int cores,
//                    arma::vec betaX,
//                    arma::vec betaU,
//                    arma::vec delta) { 
//   
//   
//   // initialize loglikelihood       
//   vec f = zeros(n);   
// 
//   // initialize design for one person, without random effects
//   mat Z = -diagmat(ones(q*I));
//   mat etai = Z*delta;
//   // current linear predictor without random effects, per item and category
//   etai = reshape(etai,q,I);
// 
//   // create sigma matrix from current parameters
//   double co_var = alpha(1)*sqrt(alpha(0))*sqrt(alpha(2));
//   mat sigma = zeros(2,2);        
//   sigma(0,0) = alpha(0);
//   sigma(1,0) = co_var;
//   sigma(0,1) = co_var;
//   sigma(1,1) = alpha(2);
//   // square root of current sigma matrix
//   mat sigma12 = trans(chol(sigma));
// 
//   // initialize number of different threads
// #ifdef _OPENMP
//   if(cores > 1)
//     omp_set_num_threads(cores);
// #endif
//   
//   mat Xi; mat yi; mat prods_i; int j; int jj; int k; int i; 
//   vec rnd_act; mat etaij; mat mui; vec eta_k; vec mu_k;
//   // loop through all persons
// #pragma omp parallel for private(i, j, jj, k, yi, Xi, prods_i, rnd_act, etaij, mui, eta_k, mu_k) shared(f)
//   for(i=0; i<n ;i++){
//     Xi = X(i,span::all);
//     
//     // get response of person i
//     vec yi = Y(span(i*q*I,i*q*I+q*I-1));
//     
//     prods_i = ones(Q,Q);
//     // loop through all knots, for both random effects
//     for(j=0; j<Q ;j++){ 
//       for(jj=0; jj<Q ;jj++){ 
//         // initialize current random effects or knots
//         rnd_act = zeros(2);
//         rnd_act(0) = GHnodes(j);
//         rnd_act(1) = GHnodes(jj);
//         
//         rnd_act = sigma12*rnd_act+join_cols(Xi*betaX,Xi*betaU);
//         
//         
//         etaij = (etai + ones(q,I)*rnd_act(0))%(ones(q,I)*exp(rnd_act(1)));
//         mui = ones(q+1,I);
//         // loop through all items  
//         for(k=0; k<I ;k++){
//           
//           // get eta and mu of person i for current item k
//           eta_k = etaij(span::all,k);
//           mu_k = responseFun(eta_k);
//           mui(span::all,k) = join_cols(mu_k,1-sum(mu_k,0));
//         }
//         // create prob for current item, update prob matrix for respective knots
//         prods_i(j,jj) = prods_i(j,jj)*prod(prod(mui % yi - (yi-1)));
//         
//       }
//     }
//     f(i)  = - log(accu(prods_i%GHweights));
//     // accumulate all likelihood contributions, weights by respective weights and probs of knots
//   }
//   
//   
//   return sum(f);
// }

// [[Rcpp::export]]
double loglikUPCM4(arma::vec alpha,
                   arma::vec Y,
                   int Q,
                   int q,
                   int n,
                   int I,
                   arma::mat GHweights,
                   arma::vec GHnodes,
                   arma::mat X,
                   int cores,
                   arma::vec betaX,
                   arma::vec betaU,
                   arma::vec delta,
                   double lambda) {
  
  
  // initialize loglikelihood       
  vec f = zeros(n);   

  vec alpha2 = join_cols(join_cols(join_cols(delta,betaX),betaU),alpha);
  
  double P2 = accu(alpha2%alpha2)*lambda;
  // initialize design for one person, without random effects
  mat Z = -diagmat(ones(q*I));
  mat etai = Z*delta;
  // current linear predictor without random effects, per item and category
  etai = reshape(etai,q,I);
  
    // create sigma matrix from current parameters
    double co_var = alpha(1)*sqrt(alpha(0))*sqrt(alpha(2));
    mat sigma = zeros(2,2);
    sigma(0,0) = alpha(0);
    sigma(1,0) = co_var;
    sigma(0,1) = co_var;
    sigma(1,1) = alpha(2);
    // square root of current sigma matrix

    // mat sigma12 = trans(chol(sigma));
    mat sigma12;
    try{
      sigma12 = trans(chol(sigma));
    }
    catch(...)
    {
      sigma = sigma + diagmat(ones(2)*0.0001);
      sigma12 = trans(chol(sigma));
    }
  
  // initialize number of different threads
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
  
  
  mat Xi; vec yi; mat prods_i; int j; int jj; int k; int i; 
  vec rnd_act; mat etaij; mat mui; vec eta_k; vec mu_k;
  vec yi_k; 
  // loop through all persons
  #pragma omp parallel for private(i, j, jj, k, yi, yi_k, Xi, prods_i, rnd_act, etaij, mui, eta_k, mu_k) shared(f)
  for(i=0; i<n ;i++){

Xi = X(i,span::all);
    
    // get response of person i
   yi = Y(span(i*q*I,i*q*I+q*I-1));
   prods_i = ones(Q,Q);
    // loop through all knots, for both random effects
    for(j=0; j<Q ;j++){ 
      for(jj=0; jj<Q ;jj++){ 
        // initialize current random effects or knots
        rnd_act = zeros(2);
        rnd_act(0) = GHnodes(j);
        rnd_act(1) = GHnodes(jj);
        rnd_act = sigma12*rnd_act+join_cols(Xi*betaX,Xi*betaU);
        
        
        // loop through all items  
        for(k=0; k<I ;k++){  
          // response of person i for current item k
          yi_k = yi(span(k*q,q+k*q-1));
          yi_k = join_cols(yi_k,1-sum(yi_k,0));
          
          // get eta and mu of person i for current item k
          eta_k = (etai(span::all,k) + ones(q)*(rnd_act(0)))%(ones(q)*(exp(rnd_act(1))));
          mu_k = responseFun(eta_k);
          mu_k = join_cols(mu_k,1-sum(mu_k,0));
          
          // create prob for current item, update prob matrix for respective knots

          prods_i(j,jj) = prods_i(j,jj)*prod(mu_k % yi_k - (yi_k-1));
        }
      }
    }
    f(i)  = -log(accu(prods_i%GHweights));
    // accumulate all likelihood contributions, weights by respective weights and probs of knots
  }
  
  return sum(f)+P2;
}


// [[Rcpp::export]]
arma::vec scoreUPCM(arma::vec alpha,
                    arma::vec Y,
                    int Q,
                    int q,
                    int n,
                    int I,
                    int pall,
                    arma::mat GHweights,
                    arma::vec GHnodes,
                    int pX,
                    arma::mat X,
                    int cores,
                    arma::vec sigma,
                    double lambda) {


  vec P2 = 2*alpha*lambda;
  // initialize cij matrix for all persons and all knots,
  // will be needed for normalization per person afterwards
  // in contrast, for implementation of PCM this was calculated BEFORE looping through persons
  // should be faster this way
  mat cij_mat = zeros(n,Q*Q);

  // initialize matrix containing score contributions per person and per parameter
  // has to be kept because one has to use cij_mat before one can sum over persons
  mat help_mat = zeros(pall,n);

  // initialize design for one person, without random effects
  mat Z = -join_rows(diagmat(ones(q*I)),zeros(q*I,2*pX));
  mat etai = Z*alpha;
  // current linear predictor without random effects, per item and category
  etai = reshape(etai,q,I);

  // create sigma matrix from current parameters
  double co_var = sigma(1)*sqrt(sigma(0))*sqrt(sigma(2));
  mat sigmaMat = zeros(2,2);
  sigmaMat(0,0) = sigma(0);
  sigmaMat(1,0) = co_var;
  sigmaMat(0,1) = co_var;
  sigmaMat(1,1) = sigma(2);

  // square root of current sigma matrix
  // mat sigma12 = trans(chol(sigma));
  mat sigma12;
  try{
    sigma12 = trans(chol(sigmaMat));
  }
  catch(...)
  {
    sigmaMat = sigmaMat + diagmat(ones(2)*0.0001);
    sigma12 = trans(chol(sigmaMat));
  }


  vec betaX = alpha(span(q*I,q*I+pX-1));
  vec betaU = alpha(span(q*I+pX,q*I+2*pX-1));

  // initialize number of different threads
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif


  mat Xi; vec yi; int j; int jj; int k; int i; 
  vec rnd_act; vec eta_k; vec mu_k;
  vec yi_k; int pos_knot; mat Xihelp; mat Xihelp2; mat Zij; mat D_i; mat SigmaInv_i; vec mu_i;
  double cij_help;
  // loop through all persons
#pragma omp parallel for private(i, Xihelp, Xihelp2, Zij, D_i, SigmaInv_i, mu_i, pos_knot, j, jj, k, yi, yi_k, Xi, rnd_act, eta_k, mu_k, cij_help) shared(cij_mat, help_mat)
  for(i=0; i<n ;i++){

    // initialize a running vector over all Q*Q knots
    pos_knot = 0;

    // get response of person i
    yi = Y(span(i*q*I,i*q*I+q*I-1));

    Xi = X(i,span::all);
      Xihelp = Xi;
      if(q>1){
        for(int f=1; f<q ;f++){
          Xihelp = join_cols(Xihelp,Xi);
        }
      }
      
    // loop through all knots, for both random effects
    for(j=0; j<Q ;j++){
      for(jj=0; jj<Q ;jj++){
        Zij = Z;
        rnd_act = zeros(2);
        rnd_act(0) = GHnodes(j);
        rnd_act(1) = GHnodes(jj);
        
        rnd_act = sigma12*rnd_act+join_cols(Xi*betaX,Xi*betaU);
        
        // initialize derivatives, inverse sigma and mu for all items and current knots
        D_i = zeros(q*I,q*I);
        SigmaInv_i = zeros(q*I,q*I);
        mu_i = zeros(q*I);

        // initialize current cij value with weight and prob of current knot-combination
        cij_help = GHweights(j,jj);

        // loop through all items
        for(k=0; k<I ;k++){
          eta_k = (etai(span::all,k) + ones(q)*rnd_act(0));
          Xihelp2 = Xihelp%(eta_k*trans(ones(pX)));
          eta_k = eta_k%(ones(q)*exp(rnd_act(1)));
          
          // response of person i for current item k
          yi_k = yi(span(k*q,q+k*q-1));
          yi_k = join_cols(yi_k,1-sum(yi_k,0));

          // create current eta and mu
          mu_k = responseFun(eta_k);
          mu_i(span(k*q,q+k*q-1)) = mu_k;

          // update the respective part of design matrix

            Zij(span(k*q,q+k*q-1),span(pall-2*pX,pall-1)) = join_rows(Xihelp,Xihelp2);


          // update derivative and inverse sigma for current item
          D_i(span(k*q,q+k*q-1),span(k*q,q+k*q-1)) = createD(mu_k);
          SigmaInv_i(span(k*q,q+k*q-1),span(k*q,q+k*q-1)) = createSigmaInv(mu_k);

          // update current cij by probability of current response
          mu_k = join_cols(mu_k,1-sum(mu_k,0));

          cij_help = cij_help * prod(mu_k % yi_k - (yi_k-1));
        }
        Zij = Zij * exp(rnd_act(1));

        // after looping through all items, update cij_mat
        cij_mat(i,pos_knot) = cij_help;

          // Rcout<<cij_help<<endl;
          // vec x = trans(Zij)*D_i*SigmaInv_i*(yi-mu_i);
          // double xx =x(0);
          // if(xx!=xx){
          //   Rcout<<trans(trans(Zij)*D_i*SigmaInv_i*(yi-mu_i))<<endl;
            // Rcout<<Zij<<endl;
            // Rcout<<trans(trans(Zij)*D_i)<<endl;
            // Rcout<<SigmaInv_i<<endl;
            // Rcout<<yi-mu_i<<endl;
          // }
          // Rcout<< trans((trans(Zij)*D_i*SigmaInv_i*(yi-mu_i)))<<endl;


        // update all contribution of current knot-combination and person i to column of person i
        help_mat(span::all,i) = help_mat(span::all,i) + (trans(Zij)*D_i*SigmaInv_i*(yi-mu_i))*cij_help;

        pos_knot = pos_knot + 1;
      }}
  }

  // normalization value per person
  vec cij_norm = sum(cij_mat,1);

  // normalize row of each parameter by cij_norm
  for(int e=0; e<pall ;e++){
    help_mat(e,span::all) = help_mat(e,span::all)%trans(1/cij_norm);
  }
// Rcout<<cij_norm<<endl;
  // sum score contributions over all persons per covariate
  vec s = -sum(help_mat,1)+P2;

  return s;
}



// [[Rcpp::export]]
double loglikUGPCM(arma::vec alpha,
                  arma::vec Y,
                  int Q,
                  int q,
                  int n,
                  int I,
                  int pall,
                  arma::mat GHweights,
                  arma::vec GHnodes,
                  int pX,
                  arma::mat X,
                  int cores,
                  double lambda) { 
  
  // initialize loglikelihood       
  vec f = zeros(n);   
  
  double P2 = accu(alpha%alpha)*lambda;
  
  // initialize design for one person, without random effects
  mat Z = -join_rows(diagmat(ones(q*I)),zeros(q*I,2*pX+2+I));
  mat etai = Z*alpha;
  // current linear predictor without random effects, per item and category
  etai = reshape(etai,q,I);
  
  // create sigma matrix from current parameters
  double co_var = alpha(pall-1-I)*sqrt(alpha(pall-I))*sqrt(alpha(pall-2-I));
  mat sigma = zeros(2,2);        
  sigma(0,0) = alpha(pall-2-I);
  sigma(1,0) = co_var;
  sigma(0,1) = co_var;
  sigma(1,1) = alpha(pall-I);
  // square root of current sigma matrix

  // mat sigma12 = trans(chol(sigma));
  mat sigma12;
  try{
    sigma12 = trans(chol(sigma));
  }
  catch(...)
  {
    sigma = sigma + diagmat(ones(2)*0.0001);
    sigma12 = trans(chol(sigma));
  }
  
  vec betaX = alpha(span(q*I,q*I+pX-1));
  vec betaU = alpha(span(q*I+pX,q*I+2*pX-1));
  
  vec slopes = join_cols(ones<vec>(1), alpha(span(pall-I+1, pall-1)));
  
  // initialize number of different threads
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
  
  mat Xi; vec yi; mat prods_i; int j; int jj; int k; int i; 
  vec rnd_act; mat etaij; mat mui; vec eta_k; vec mu_k;
  vec yi_k; 
  // loop through all persons
#pragma omp parallel for private(i, j, jj, k, yi, yi_k, Xi, prods_i, rnd_act, etaij, mui, eta_k, mu_k) shared(f)
  for(i=0; i<n ;i++){
    
    Xi = X(i,span::all);
    
    // get response of person i
    yi = Y(span(i*q*I,i*q*I+q*I-1));
    prods_i = ones(Q,Q);
    // loop through all knots, for both random effects
    for(j=0; j<Q ;j++){ 
      for(jj=0; jj<Q ;jj++){ 
        // initialize current random effects or knots
        rnd_act = zeros(2);
        rnd_act(0) = GHnodes(j);
        rnd_act(1) = GHnodes(jj);
        
        rnd_act = sigma12*rnd_act+join_cols(Xi*betaX,Xi*betaU);
        
        
        // loop through all items  
        for(k=0; k<I ;k++){  
          // response of person i for current item k
          yi_k = yi(span(k*q,q+k*q-1));
          yi_k = join_cols(yi_k,1-sum(yi_k,0));
          
          // get eta and mu of person i for current item k
          eta_k = (etai(span::all,k) + ones(q)*(rnd_act(0)))%(slopes(k)*ones(q)*(exp(rnd_act(1))));
          mu_k = responseFun(eta_k);
          mu_k = join_cols(mu_k,1-sum(mu_k,0));
          // if(prod(mu_k % yi_k - (yi_k-1))<0){
          //   Rcout<<trans(mu_k)<<endl;
          // }
          // create prob for current item, update prob matrix for respective knots
          prods_i(j,jj) = prods_i(j,jj)*prod(mu_k % yi_k - (yi_k-1));
        }
      }
    }
    f(i)  = -log(accu(prods_i%GHweights));
    
    // accumulate all likelihood contributions, weights by respective weights and probs of knots
  }
  return sum(f)+P2;
}


// [[Rcpp::export]]
double loglikUGPCM2(arma::vec alpha,
                   arma::mat Y,
                   int Q,
                   int q,
                   int n,
                   int I,
                   int pall,
                   arma::mat GHweights,
                   arma::vec GHnodes,
                   int pX,
                   arma::mat X,
                   int cores,
                   double lambda) { 
  
  
  // initialize loglikelihood       
  vec f = zeros(n);   
  
  double P2 = accu(alpha%alpha)*lambda;
  
  // initialize design for one person, without random effects
  mat Z = -join_rows(diagmat(ones(q*I)),zeros(q*I,2*pX+2+I));
  mat etai = Z*alpha;
  // current linear predictor without random effects, per item and category
  etai = reshape(etai,q,I);
  
  // create sigma matrix from current parameters
  double co_var = alpha(pall-1-I)*sqrt(alpha(pall-I))*sqrt(alpha(pall-2-I));
  mat sigma = zeros(2,2);        
  sigma(0,0) = alpha(pall-2-I);
  sigma(1,0) = co_var;
  sigma(0,1) = co_var;
  sigma(1,1) = alpha(pall-I);
  // square root of current sigma matrix

  // mat sigma12 = trans(chol(sigma));
  mat sigma12;
  try{
    sigma12 = trans(chol(sigma));
  }
  catch(...)
  {
    sigma = sigma + diagmat(ones(2)*0.0001);
    sigma12 = trans(chol(sigma));
  }
  
  vec betaX = alpha(span(q*I,q*I+pX-1));
  vec betaU = alpha(span(q*I+pX,q*I+2*pX-1));
  
  vec slopes = join_cols(ones<vec>(1), alpha(span(pall-I+1, pall-1)));
  
  // initialize number of different threads
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
  
  mat Xi; mat yi; mat prods_i; int j; int jj; int k; int i; 
  vec rnd_act; mat etaij; mat mui; vec eta_k; vec mu_k;
  // loop through all persons
#pragma omp parallel for private(i, j, jj, k, yi, Xi, prods_i, rnd_act, etaij, mui, eta_k, mu_k) shared(f)
  for(i=0; i<n ;i++){
    Xi = X(i,span::all);
    
    // get response of person i
    yi = Y(span::all,span(i*I,i*I+I-1));
    
    prods_i = ones(Q,Q);
    // loop through all knots, for both random effects
    for(j=0; j<Q ;j++){ 
      for(jj=0; jj<Q ;jj++){ 
        // initialize current random effects or knots
        rnd_act = zeros(2);
        rnd_act(0) = GHnodes(j);
        rnd_act(1) = GHnodes(jj);
        
        rnd_act = sigma12*rnd_act+join_cols(Xi*betaX,Xi*betaU);
        
        
        etaij = (etai + ones(q,I)*rnd_act(0))%((ones(q)*trans(slopes))*exp(rnd_act(1)));
        mui = ones(q+1,I);
        // loop through all items  
        for(k=0; k<I ;k++){
          
          // get eta and mu of person i for current item k
          eta_k = etaij(span::all,k);
          mu_k = responseFun2(eta_k);
          mui(span::all,k) = join_cols(mu_k,1-sum(mu_k,0));
        }
        // create prob for current item, update prob matrix for respective knots
        prods_i(j,jj) = prods_i(j,jj)*prod(prod(mui % yi - (yi-1)));
        
      }
    }
    f(i)  = - log(accu(prods_i%GHweights));
    // accumulate all likelihood contributions, weights by respective weights and probs of knots
  }
  
  
  return sum(f)+P2;
}

// [[Rcpp::export]]
arma::vec scoreUGPCM(arma::vec alpha,
                    arma::vec Y,
                    int Q,
                    int q,
                    int n,
                    int I,
                    int pall,
                    arma::mat GHweights,
                    arma::vec GHnodes,
                    int pX,
                    arma::mat X,
                    int cores,
                    arma::vec sigma,
                    double lambda) {
  
  
  vec P2 = 2*alpha*lambda;
  // initialize cij matrix for all persons and all knots,
  // will be needed for normalization per person afterwards
  // in contrast, for implementation of PCM this was calculated BEFORE looping through persons
  // should be faster this way
  mat cij_mat = zeros(n,Q*Q);
  
  // initialize matrix containing score contributions per person and per parameter
  // has to be kept because one has to use cij_mat before one can sum over persons
  mat help_mat = zeros(pall,n);
  
  // initialize design for one person, without random effects
  mat Z = -join_rows(diagmat(ones(q*I)),zeros(q*I,2*pX+I-1));
  mat etai = Z*alpha;
  // current linear predictor without random effects, per item and category
  etai = reshape(etai,q,I);
  
  // create sigma matrix from current parameters
  double co_var = sigma(1)*sqrt(sigma(0))*sqrt(sigma(2));
  mat sigmaMat = zeros(2,2);
  sigmaMat(0,0) = sigma(0);
  sigmaMat(1,0) = co_var;
  sigmaMat(0,1) = co_var;
  sigmaMat(1,1) = sigma(2);
  
  // square root of current sigma matrix
  // mat sigma12 = trans(chol(sigma));
  mat sigma12;
  try{
    sigma12 = trans(chol(sigmaMat));
  }
  catch(...)
  {
    sigmaMat = sigmaMat + diagmat(ones(2)*0.0001);
    sigma12 = trans(chol(sigmaMat));
  }
  
  
  vec betaX = alpha(span(q*I,q*I+pX-1));
  vec betaU = alpha(span(q*I+pX,q*I+2*pX-1));
  
  vec slopes = join_cols(ones<vec>(1), alpha(span(pall-I+1, pall-1)));
  
  // initialize number of different threads
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
  
  
  mat Xi; vec yi; int j; int jj; int k; int i; 
  vec rnd_act; vec eta_k; vec mu_k;
  vec yi_k; int pos_knot; mat Xihelp; mat Xihelp2; mat Xihelp3; mat Zij; mat D_i; mat SigmaInv_i; vec mu_i;
  double cij_help;
  // loop through all persons
#pragma omp parallel for private(i, Xihelp, Xihelp2, Xihelp3, Zij, D_i, SigmaInv_i, mu_i, pos_knot, j, jj, k, yi, yi_k, Xi, rnd_act, eta_k, mu_k, cij_help) shared(cij_mat, help_mat)
  for(i=0; i<n ;i++){
    
    // initialize a running vector over all Q*Q knots
    pos_knot = 0;
    
    // get response of person i
    yi = Y(span(i*q*I,i*q*I+q*I-1));
    
    Xi = X(i,span::all);
    Xihelp = Xi;
    if(q>1){
      for(int f=1; f<q ;f++){
        Xihelp = join_cols(Xihelp,Xi);
      }
    }
    
    // loop through all knots, for both random effects
    for(j=0; j<Q ;j++){
      for(jj=0; jj<Q ;jj++){
        Zij = Z;
        rnd_act = zeros(2);
        rnd_act(0) = GHnodes(j);
        rnd_act(1) = GHnodes(jj);
        
        rnd_act = sigma12*rnd_act+join_cols(Xi*betaX,Xi*betaU);
        
        // initialize derivatives, inverse sigma and mu for all items and current knots
        D_i = zeros(q*I,q*I);
        SigmaInv_i = zeros(q*I,q*I);
        mu_i = zeros(q*I);
        
        // initialize current cij value with weight and prob of current knot-combination
        cij_help = GHweights(j,jj);
        
        // loop through all items
        for(k=0; k<I ;k++){
          eta_k = (etai(span::all,k) + ones(q)*rnd_act(0));
          Xihelp2 = (Xihelp%(eta_k*trans(ones(pX))));
          Xihelp3 = zeros(q,I-1);
          if(k>0){
            Xihelp3(span::all,k-1) = eta_k*(1/slopes(k));  
          }
          
          eta_k = eta_k%(slopes(k)*ones(q)*exp(rnd_act(1)));
          
          
          // response of person i for current item k
          yi_k = yi(span(k*q,q+k*q-1));
          yi_k = join_cols(yi_k,1-sum(yi_k,0));
          
          // create current eta and mu
          mu_k = responseFun(eta_k);
          mu_i(span(k*q,q+k*q-1)) = mu_k;
          
          // update the respective part of design matrix
          Zij(span(k*q,q+k*q-1),span(pall-2*pX-I+1,pall-1)) = join_rows(join_rows(Xihelp,Xihelp2),Xihelp3);
          Zij(span(k*q,q+k*q-1),span::all) = Zij(span(k*q,q+k*q-1),span::all)*slopes(k);
          
          
          // update derivative and inverse sigma for current item
          D_i(span(k*q,q+k*q-1),span(k*q,q+k*q-1)) = createD(mu_k);
          SigmaInv_i(span(k*q,q+k*q-1),span(k*q,q+k*q-1)) = createSigmaInv(mu_k);
          
          // update current cij by probability of current response
          mu_k = join_cols(mu_k,1-sum(mu_k,0));
          
          cij_help = cij_help * prod(mu_k % yi_k - (yi_k-1));
        }
        // double f1 = accu(D_i);
        // if(f1!=f1){
        //   Rcout<<"D_i"<<endl;
        //   // Rcout<<D_i<<endl;
        // }
        // double f2 = accu(SigmaInv_i);
        // if(f2!=f2){
        //   Rcout<<"Sigma_i"<<endl;
        //   // Rcout<<SigmaInv_i<<endl;
        // }
        
        Zij = Zij * exp(rnd_act(1));
        
        // after looping through all items, update cij_mat
        cij_mat(i,pos_knot) = cij_help;
        
        
        // update all contribution of current knot-combination and person i to column of person i
        help_mat(span::all,i) = help_mat(span::all,i) + (trans(Zij)*D_i*SigmaInv_i*(yi-mu_i))*cij_help;
        
        pos_knot = pos_knot + 1;
      }}
  }
  
  // normalization value per person
  vec cij_norm = sum(cij_mat,1);
  
  // normalize row of each parameter by cij_norm
  for(int e=0; e<pall ;e++){
    help_mat(e,span::all) = help_mat(e,span::all)%trans(1/cij_norm);
  }
  // sum score contributions over all persons per covariate
  vec s = -sum(help_mat,1)+P2;
  
  return s;
}

// [[Rcpp::export]]
double loglikUGPCM4(arma::vec alpha,
                   arma::vec Y,
                   int Q,
                   int q,
                   int n,
                   int I,
                   arma::mat GHweights,
                   arma::vec GHnodes,
                   arma::mat X,
                   int cores,
                   arma::vec betaX,
                   arma::vec betaU,
                   arma::vec delta,
                   arma::vec slopes,
                   double lambda) {
  
  
  // initialize loglikelihood       
  vec f = zeros(n);   
  
  vec alpha2 = join_cols(join_cols(join_cols(join_cols(delta,betaX),betaU),alpha),slopes);
  
  double P2 = accu(alpha2%alpha2)*lambda;
  // initialize design for one person, without random effects
  mat Z = -diagmat(ones(q*I));
  mat etai = Z*delta;
  // current linear predictor without random effects, per item and category
  etai = reshape(etai,q,I);
  
  // create sigma matrix from current parameters
  double co_var = alpha(1)*sqrt(alpha(0))*sqrt(alpha(2));
  mat sigma = zeros(2,2);
  sigma(0,0) = alpha(0);
  sigma(1,0) = co_var;
  sigma(0,1) = co_var;
  sigma(1,1) = alpha(2);
  // square root of current sigma matrix

  // mat sigma12 = trans(chol(sigma));
  mat sigma12;
  try{
    sigma12 = trans(chol(sigma));
  }
  catch(...)
  {
    sigma = sigma + diagmat(ones(2)*0.0001);
    sigma12 = trans(chol(sigma));
  }
  
  // initialize number of different threads
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
  
  
  mat Xi; vec yi; mat prods_i; int j; int jj; int k; int i; 
  vec rnd_act; mat etaij; mat mui; vec eta_k; vec mu_k;
  vec yi_k; 
  // loop through all persons
#pragma omp parallel for private(i, j, jj, k, yi, yi_k, Xi, prods_i, rnd_act, etaij, mui, eta_k, mu_k) shared(f)
  for(i=0; i<n ;i++){
    
    Xi = X(i,span::all);
    
    // get response of person i
    yi = Y(span(i*q*I,i*q*I+q*I-1));
    prods_i = ones(Q,Q);
    // loop through all knots, for both random effects
    for(j=0; j<Q ;j++){ 
      for(jj=0; jj<Q ;jj++){ 
        // initialize current random effects or knots
        rnd_act = zeros(2);
        rnd_act(0) = GHnodes(j);
        rnd_act(1) = GHnodes(jj);
        rnd_act = sigma12*rnd_act+join_cols(Xi*betaX,Xi*betaU);
        
        
        // loop through all items  
        for(k=0; k<I ;k++){  
          // response of person i for current item k
          yi_k = yi(span(k*q,q+k*q-1));
          yi_k = join_cols(yi_k,1-sum(yi_k,0));
          
          // get eta and mu of person i for current item k
          eta_k = (etai(span::all,k) + ones(q)*(rnd_act(0)))%(slopes(k)*ones(q)*(exp(rnd_act(1))));
          mu_k = responseFun(eta_k);
          mu_k = join_cols(mu_k,1-sum(mu_k,0));
          
          // create prob for current item, update prob matrix for respective knots
          
          prods_i(j,jj) = prods_i(j,jj)*prod(mu_k % yi_k - (yi_k-1));
        }
      }
    }
    f(i)  = -log(accu(prods_i%GHweights));
    // accumulate all likelihood contributions, weights by respective weights and probs of knots
  }
  
  return sum(f)+P2;
}





// [[Rcpp::export]]
double loglikUPCMnoX(arma::vec alpha,
                  arma::vec Y,
                  int Q,
                  int q,
                  int n,
                  int I,
                  int pall,
                  arma::mat GHweights,
                  arma::vec GHnodes,
                  int pX,
                  arma::mat X,
                  int cores,
                  double lambda) { 
  
  // initialize loglikelihood       
  vec f = zeros(n);   
  
  double P2 = accu(alpha%alpha)*lambda;
  
  // initialize design for one person, without random effects
  mat Z = -join_rows(diagmat(ones(q*I)),zeros(q*I,2*pX+3));
  mat etai = Z*alpha;
  // current linear predictor without random effects, per item and category
  etai = reshape(etai,q,I);
  
  // create sigma matrix from current parameters
  double co_var = alpha(pall-2)*sqrt(alpha(pall-1))*sqrt(alpha(pall-3));
  mat sigma = zeros(2,2);        
  sigma(0,0) = alpha(pall-3);
  sigma(1,0) = co_var;
  sigma(0,1) = co_var;
  sigma(1,1) = alpha(pall-1);
  // square root of current sigma matrix
  
  // mat sigma12 = trans(chol(sigma));
  mat sigma12;
  try{
    sigma12 = trans(chol(sigma));
  }
  catch(...)
  {
    sigma = sigma + diagmat(ones(2)*0.0001);
    sigma12 = trans(chol(sigma));
  }
  
  
  // initialize number of different threads
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
  
   vec yi; mat prods_i; int j; int jj; int k; int i; 
  vec rnd_act; mat etaij; mat mui; vec eta_k; vec mu_k;
  vec yi_k; 
  // loop through all persons
#pragma omp parallel for private(i, j, jj, k, yi, yi_k, prods_i, rnd_act, etaij, mui, eta_k, mu_k) shared(f)
  for(i=0; i<n ;i++){
    
    
    // get response of person i
    yi = Y(span(i*q*I,i*q*I+q*I-1));
    prods_i = ones(Q,Q);
    // loop through all knots, for both random effects
    for(j=0; j<Q ;j++){ 
      for(jj=0; jj<Q ;jj++){ 
        // initialize current random effects or knots
        rnd_act = zeros(2);
        rnd_act(0) = GHnodes(j);
        rnd_act(1) = GHnodes(jj);
        
        rnd_act = sigma12*rnd_act;
        
        
        // loop through all items  
        for(k=0; k<I ;k++){  
          // response of person i for current item k
          yi_k = yi(span(k*q,q+k*q-1));
          yi_k = join_cols(yi_k,1-sum(yi_k,0));
          
          // get eta and mu of person i for current item k
          eta_k = (etai(span::all,k) + ones(q)*(rnd_act(0)))%(ones(q)*(exp(rnd_act(1))));
          mu_k = responseFun2(eta_k);
          mu_k = join_cols(mu_k,1-sum(mu_k,0));
// if(prod(mu_k)<0){
//   Rcout<<i<<endl;
//   Rcout<<mu_k<<endl;
//   // Rcout<<yi_k<<endl;
//   // Rcout<<eta_k<<endl;
//   // Rcout<<rnd_act<<endl;
//   // Rcout<<etai(span::all,k)<<endl;
// }
          // create prob for current item, update prob matrix for respective knots
          prods_i(j,jj) = prods_i(j,jj)*prod(mu_k % yi_k - (yi_k-1));
        }
      }
    }
    f(i)  = -log(accu(prods_i%GHweights));
    // if(isnan(f(i))){
    //   Rcout<<prods_i<<endl;
    //   }
    // accumulate all likelihood contributions, weights by respective weights and probs of knots
  }
  // Rcout<<trans(f)<<endl;
  return sum(f)+P2;
}


// [[Rcpp::export]]
arma::vec scoreUPCMnoX(arma::vec alpha,
                    arma::vec Y,
                    int Q,
                    int q,
                    int n,
                    int I,
                    int pall,
                    arma::mat GHweights,
                    arma::vec GHnodes,
                    int cores,
                    arma::vec sigma,
                    double lambda) {
  

  vec P2 = 2*alpha*lambda;
  // initialize cij matrix for all persons and all knots,
  // will be needed for normalization per person afterwards
  // in contrast, for implementation of PCM this was calculated BEFORE looping through persons
  // should be faster this way
  mat cij_mat = zeros(n,Q*Q);
  
  // initialize matrix containing score contributions per person and per parameter
  // has to be kept because one has to use cij_mat before one can sum over persons
  mat help_mat = zeros(pall,n);
  
  // initialize design for one person, without random effects
  mat Z = -diagmat(ones(q*I));
  mat etai = Z*alpha;
  // current linear predictor without random effects, per item and category
  etai = reshape(etai,q,I);
  
  // create sigma matrix from current parameters
  double co_var = sigma(1)*sqrt(sigma(0))*sqrt(sigma(2));
  mat sigmaMat = zeros(2,2);
  sigmaMat(0,0) = sigma(0);
  sigmaMat(1,0) = co_var;
  sigmaMat(0,1) = co_var;
  sigmaMat(1,1) = sigma(2);
  
  // square root of current sigma matrix
  // mat sigma12 = trans(chol(sigma));
  mat sigma12;
  try{
    sigma12 = trans(chol(sigmaMat));
  }
  catch(...)
  {
    sigmaMat = sigmaMat + diagmat(ones(2)*0.0001);
    sigma12 = trans(chol(sigmaMat));
  }
  
  // initialize number of different threads
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
  
  
   vec yi; int j; int jj; int k; int i; 
  vec rnd_act; vec eta_k; vec mu_k;
  vec yi_k; int pos_knot; mat Zij; mat D_i; mat SigmaInv_i; vec mu_i;
  double cij_help;
  // loop through all persons
#pragma omp parallel for private(i,  Zij, D_i, SigmaInv_i, mu_i, pos_knot, j, jj, k, yi, yi_k,  rnd_act, eta_k, mu_k, cij_help) shared(cij_mat, help_mat)
  for(i=0; i<n ;i++){
    
    // initialize a running vector over all Q*Q knots
    pos_knot = 0;
    
    // get response of person i
    yi = Y(span(i*q*I,i*q*I+q*I-1));
    
    // loop through all knots, for both random effects
    for(j=0; j<Q ;j++){
      for(jj=0; jj<Q ;jj++){
        Zij = Z;
        rnd_act = zeros(2);
        rnd_act(0) = GHnodes(j);
        rnd_act(1) = GHnodes(jj);
        
        rnd_act = sigma12*rnd_act;
        
        // initialize derivatives, inverse sigma and mu for all items and current knots
        D_i = zeros(q*I,q*I);
        SigmaInv_i = zeros(q*I,q*I);
        mu_i = zeros(q*I);
        
        // initialize current cij value with weight and prob of current knot-combination
        cij_help = GHweights(j,jj);
        
        // loop through all items
        for(k=0; k<I ;k++){
          eta_k = (etai(span::all,k) + ones(q)*rnd_act(0));

          eta_k = eta_k%(ones(q)*exp(rnd_act(1)));
          
          // response of person i for current item k
          yi_k = yi(span(k*q,q+k*q-1));
          yi_k = join_cols(yi_k,1-sum(yi_k,0));
          
          // create current eta and mu
          mu_k = responseFun(eta_k);
          mu_i(span(k*q,q+k*q-1)) = mu_k;
          
          
          // update derivative and inverse sigma for current item
          D_i(span(k*q,q+k*q-1),span(k*q,q+k*q-1)) = createD(mu_k);
          SigmaInv_i(span(k*q,q+k*q-1),span(k*q,q+k*q-1)) = createSigmaInv(mu_k);
          
          // update current cij by probability of current response
          mu_k = join_cols(mu_k,1-sum(mu_k,0));
          
          cij_help = cij_help * prod(mu_k % yi_k - (yi_k-1));
        }
        Zij = Zij * exp(rnd_act(1));
        
        // after looping through all items, update cij_mat
        cij_mat(i,pos_knot) = cij_help;
        
        // update all contribution of current knot-combination and person i to column of person i
        help_mat(span::all,i) = help_mat(span::all,i) + (trans(Zij)*D_i*SigmaInv_i*(yi-mu_i))*cij_help;
        
        pos_knot = pos_knot + 1;
      }}
  }
  
  // normalization value per person
  vec cij_norm = sum(cij_mat,1);
  
  // normalize row of each parameter by cij_norm
  for(int e=0; e<pall ;e++){
    help_mat(e,span::all) = help_mat(e,span::all)%trans(1/cij_norm);
  }
  // Rcout<<cij_norm<<endl;
  // sum score contributions over all persons per covariate
  vec s = -sum(help_mat,1)+P2;
  
  return s;
}


// [[Rcpp::export]]
double loglikUPCM4noX(arma::vec alpha,
                   arma::vec Y,
                   int Q,
                   int q,
                   int n,
                   int I,
                   arma::mat GHweights,
                   arma::vec GHnodes,
                   int cores,
                   arma::vec delta,
                   double lambda) {
  
  
  // initialize loglikelihood       
  vec f = zeros(n);   
  
  vec alpha2 = join_cols(delta,alpha);
  
  double P2 = accu(alpha2%alpha2)*lambda;
  // initialize design for one person, without random effects
  mat Z = -diagmat(ones(q*I));
  mat etai = Z*delta;
  // current linear predictor without random effects, per item and category
  etai = reshape(etai,q,I);
  
  // create sigma matrix from current parameters
  double co_var = alpha(1)*sqrt(alpha(0))*sqrt(alpha(2));
  mat sigma = zeros(2,2);
  sigma(0,0) = alpha(0);
  sigma(1,0) = co_var;
  sigma(0,1) = co_var;
  sigma(1,1) = alpha(2);
  // square root of current sigma matrix
  
  // mat sigma12 = trans(chol(sigma));
  mat sigma12;
  try{
    sigma12 = trans(chol(sigma));
  }
  catch(...)
  {
    sigma = sigma + diagmat(ones(2)*0.0001);
    sigma12 = trans(chol(sigma));
  }
  
  // initialize number of different threads
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
  
  
 vec yi; mat prods_i; int j; int jj; int k; int i; 
  vec rnd_act; mat etaij; mat mui; vec eta_k; vec mu_k;
  vec yi_k; 
  // loop through all persons
#pragma omp parallel for private(i, j, jj, k, yi, yi_k,  prods_i, rnd_act, etaij, mui, eta_k, mu_k) shared(f)
  for(i=0; i<n ;i++){
    
    // get response of person i
    yi = Y(span(i*q*I,i*q*I+q*I-1));
    prods_i = ones(Q,Q);
    // loop through all knots, for both random effects
    for(j=0; j<Q ;j++){ 
      for(jj=0; jj<Q ;jj++){ 
        // initialize current random effects or knots
        rnd_act = zeros(2);
        rnd_act(0) = GHnodes(j);
        rnd_act(1) = GHnodes(jj);
        rnd_act = sigma12*rnd_act;
        
        
        // loop through all items  
        for(k=0; k<I ;k++){  
          // response of person i for current item k
          yi_k = yi(span(k*q,q+k*q-1));
          yi_k = join_cols(yi_k,1-sum(yi_k,0));
          
          // get eta and mu of person i for current item k
          eta_k = (etai(span::all,k) + ones(q)*(rnd_act(0)))%(ones(q)*(exp(rnd_act(1))));
          mu_k = responseFun2(eta_k);
          mu_k = join_cols(mu_k,1-sum(mu_k,0));
          
          // create prob for current item, update prob matrix for respective knots
          
          prods_i(j,jj) = prods_i(j,jj)*prod(mu_k % yi_k - (yi_k-1));
        }
      }
    }
    f(i)  = -log(accu(prods_i%GHweights));
    // accumulate all likelihood contributions, weights by respective weights and probs of knots
  }

  return sum(f)+P2;
}


// [[Rcpp::export]]
double loglikUPCM2noX(arma::vec alpha,
                   arma::mat Y,
                   int Q,
                   int q,
                   int n,
                   int I,
                   int pall,
                   arma::mat GHweights,
                   arma::vec GHnodes,
                   int pX,
                   arma::mat X,
                   int cores,
                   double lambda) { 
  
  
  // initialize loglikelihood       
  vec f = zeros(n);   
  
  double P2 = accu(alpha%alpha)*lambda;
  
  // initialize design for one person, without random effects
  mat Z = -join_rows(diagmat(ones(q*I)),zeros(q*I,2*pX+3));
  mat etai = Z*alpha;
  // current linear predictor without random effects, per item and category
  etai = reshape(etai,q,I);
  
  // create sigma matrix from current parameters
  double co_var = alpha(pall-2)*sqrt(alpha(pall-1))*sqrt(alpha(pall-3));
  mat sigma = zeros(2,2);        
  sigma(0,0) = alpha(pall-3);
  sigma(1,0) = co_var;
  sigma(0,1) = co_var;
  sigma(1,1) = alpha(pall-1);
  // square root of current sigma matrix
  
  mat sigma12;
  try{
    sigma12 = trans(chol(sigma));
  }
  catch(...)
  {
    sigma = sigma + diagmat(ones(2)*0.0001);
    sigma12 = trans(chol(sigma));
  }
  
  
  // initialize number of different threads
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
  
  mat yi; mat prods_i; int j; int jj; int k; int i; 
  vec rnd_act; mat etaij; mat mui; vec eta_k; vec mu_k;
  // loop through all persons
#pragma omp parallel for private(i, j, jj, k, yi, prods_i, rnd_act, etaij, mui, eta_k, mu_k) shared(f)
  for(i=0; i<n ;i++){
    
    // get response of person i
    yi = Y(span::all,span(i*I,i*I+I-1));
    
    prods_i = ones(Q,Q);
    // loop through all knots, for both random effects
    for(j=0; j<Q ;j++){ 
      for(jj=0; jj<Q ;jj++){ 
        // initialize current random effects or knots
        rnd_act = zeros(2);
        rnd_act(0) = GHnodes(j);
        rnd_act(1) = GHnodes(jj);
        
        rnd_act = sigma12*rnd_act;
        
        
        etaij = (etai + ones(q,I)*rnd_act(0))%(ones(q,I)*exp(rnd_act(1)));
        mui = ones(q+1,I);
        // loop through all items  
        for(k=0; k<I ;k++){
          
          // get eta and mu of person i for current item k
          eta_k = etaij(span::all,k);
          mu_k = responseFun2(eta_k);
          mui(span::all,k) = join_cols(mu_k,1-sum(mu_k,0));
        }
        // create prob for current item, update prob matrix for respective knots
        prods_i(j,jj) = prods_i(j,jj)*prod(prod(mui % yi - (yi-1)));
        
      }
    }
    f(i)  = - log(accu(prods_i%GHweights));
    // accumulate all likelihood contributions, weights by respective weights and probs of knots
  }
  
  
  return sum(f)+P2;
}


// [[Rcpp::export]]
double loglikUGPCMnoX(arma::vec alpha,
                   arma::vec Y,
                   int Q,
                   int q,
                   int n,
                   int I,
                   int pall,
                   arma::mat GHweights,
                   arma::vec GHnodes,
                   int pX,
                   arma::mat X,
                   int cores,
                   double lambda) { 
  
  // initialize loglikelihood       
  vec f = zeros(n);   
  
  double P2 = accu(alpha%alpha)*lambda;
  
  // initialize design for one person, without random effects
  mat Z = -join_rows(diagmat(ones(q*I)),zeros(q*I,2*pX+2+I));
  mat etai = Z*alpha;
  // current linear predictor without random effects, per item and category
  etai = reshape(etai,q,I);
  
  // create sigma matrix from current parameters
  double co_var = alpha(pall-1-I)*sqrt(alpha(pall-I))*sqrt(alpha(pall-2-I));
  mat sigma = zeros(2,2);        
  sigma(0,0) = alpha(pall-2-I);
  sigma(1,0) = co_var;
  sigma(0,1) = co_var;
  sigma(1,1) = alpha(pall-I);
  // square root of current sigma matrix
  
  // mat sigma12 = trans(chol(sigma));
  mat sigma12;
  try{
    sigma12 = trans(chol(sigma));
  }
  catch(...)
  {
    sigma = sigma + diagmat(ones(2)*0.0001);
    sigma12 = trans(chol(sigma));
  }
  
  vec slopes = join_cols(ones<vec>(1), alpha(span(pall-I+1, pall-1)));
  
  // initialize number of different threads
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
  
  vec yi; mat prods_i; int j; int jj; int k; int i; 
  vec rnd_act; mat etaij; mat mui; vec eta_k; vec mu_k;
  vec yi_k; 
  // loop through all persons
#pragma omp parallel for private(i, j, jj, k, yi, yi_k, prods_i, rnd_act, etaij, mui, eta_k, mu_k) shared(f)
  for(i=0; i<n ;i++){
    
    
    // get response of person i
    yi = Y(span(i*q*I,i*q*I+q*I-1));
    prods_i = ones(Q,Q);
    // loop through all knots, for both random effects
    for(j=0; j<Q ;j++){ 
      for(jj=0; jj<Q ;jj++){ 
        // initialize current random effects or knots
        rnd_act = zeros(2);
        rnd_act(0) = GHnodes(j);
        rnd_act(1) = GHnodes(jj);
        
        rnd_act = sigma12*rnd_act;
        
        
        // loop through all items  
        for(k=0; k<I ;k++){  
          // response of person i for current item k
          yi_k = yi(span(k*q,q+k*q-1));
          yi_k = join_cols(yi_k,1-sum(yi_k,0));
          
          // get eta and mu of person i for current item k
          eta_k = (etai(span::all,k) + ones(q)*(rnd_act(0)))%(slopes(k)*ones(q)*(exp(rnd_act(1))));
          mu_k = responseFun2(eta_k);
          mu_k = join_cols(mu_k,1-sum(mu_k,0));

          // create prob for current item, update prob matrix for respective knots
          prods_i(j,jj) = prods_i(j,jj)*prod(mu_k % yi_k - (yi_k-1));
        }
      }
    }
    f(i)  = -log(accu(prods_i%GHweights));
    
    // accumulate all likelihood contributions, weights by respective weights and probs of knots
  }
  return sum(f)+P2;
}

// [[Rcpp::export]]
arma::vec scoreUGPCMnoX(arma::vec alpha,
                     arma::vec Y,
                     int Q,
                     int q,
                     int n,
                     int I,
                     int pall,
                     arma::mat GHweights,
                     arma::vec GHnodes,
                     int cores,
                     arma::vec sigma,
                     double lambda) {
  
  
  vec P2 = 2*alpha*lambda;
  // initialize cij matrix for all persons and all knots,
  // will be needed for normalization per person afterwards
  // in contrast, for implementation of PCM this was calculated BEFORE looping through persons
  // should be faster this way
  mat cij_mat = zeros(n,Q*Q);
  
  // initialize matrix containing score contributions per person and per parameter
  // has to be kept because one has to use cij_mat before one can sum over persons
  mat help_mat = zeros(pall,n);
  
  // initialize design for one person, without random effects
  mat Z = -join_rows(diagmat(ones(q*I)),zeros(q*I,I-1));
  mat etai = Z*alpha;
  // current linear predictor without random effects, per item and category
  etai = reshape(etai,q,I);
  
  // create sigma matrix from current parameters
  double co_var = sigma(1)*sqrt(sigma(0))*sqrt(sigma(2));
  mat sigmaMat = zeros(2,2);
  sigmaMat(0,0) = sigma(0);
  sigmaMat(1,0) = co_var;
  sigmaMat(0,1) = co_var;
  sigmaMat(1,1) = sigma(2);
  
  // square root of current sigma matrix
  // mat sigma12 = trans(chol(sigma));
  mat sigma12;
  try{
    sigma12 = trans(chol(sigmaMat));
  }
  catch(...)
  {
    sigmaMat = sigmaMat + diagmat(ones(2)*0.0001);
    sigma12 = trans(chol(sigmaMat));
  }
  
  
  vec slopes = join_cols(ones<vec>(1), alpha(span(pall-I+1, pall-1)));
  
  // initialize number of different threads
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
  
vec yi; int j; int jj; int k; int i; 
  vec rnd_act; vec eta_k; vec mu_k;
  vec yi_k; int pos_knot; mat Zij; mat D_i; mat SigmaInv_i; vec mu_i;
  double cij_help;
  // loop through all persons
#pragma omp parallel for private(i, Zij, D_i, SigmaInv_i, mu_i, pos_knot, j, jj, k, yi, yi_k, rnd_act, eta_k, mu_k, cij_help) shared(cij_mat, help_mat)
  for(i=0; i<n ;i++){
    
    // initialize a running vector over all Q*Q knots
    pos_knot = 0;
    
    // get response of person i
    yi = Y(span(i*q*I,i*q*I+q*I-1));
    
    
    // loop through all knots, for both random effects
    for(j=0; j<Q ;j++){
      for(jj=0; jj<Q ;jj++){
        Zij = Z;
        rnd_act = zeros(2);
        rnd_act(0) = GHnodes(j);
        rnd_act(1) = GHnodes(jj);
        
        rnd_act = sigma12*rnd_act;
        
        // initialize derivatives, inverse sigma and mu for all items and current knots
        D_i = zeros(q*I,q*I);
        SigmaInv_i = zeros(q*I,q*I);
        mu_i = zeros(q*I);
        
        // initialize current cij value with weight and prob of current knot-combination
        cij_help = GHweights(j,jj);
        
        // loop through all items
        for(k=0; k<I ;k++){
          eta_k = (etai(span::all,k) + ones(q)*rnd_act(0));
          
          eta_k = eta_k%(slopes(k)*ones(q)*exp(rnd_act(1)));
          
          
          // response of person i for current item k
          yi_k = yi(span(k*q,q+k*q-1));
          yi_k = join_cols(yi_k,1-sum(yi_k,0));
          
          // create current eta and mu
          mu_k = responseFun(eta_k);
          mu_i(span(k*q,q+k*q-1)) = mu_k;
          
          // update the respective part of design matrix
          Zij(span(k*q,q+k*q-1),span::all) = Zij(span(k*q,q+k*q-1),span::all)*slopes(k);
          
          
          // update derivative and inverse sigma for current item
          D_i(span(k*q,q+k*q-1),span(k*q,q+k*q-1)) = createD(mu_k);
          SigmaInv_i(span(k*q,q+k*q-1),span(k*q,q+k*q-1)) = createSigmaInv(mu_k);
          
          // update current cij by probability of current response
          mu_k = join_cols(mu_k,1-sum(mu_k,0));
          
          cij_help = cij_help * prod(mu_k % yi_k - (yi_k-1));
        }
        
        Zij = Zij * exp(rnd_act(1));
        
        // after looping through all items, update cij_mat
        cij_mat(i,pos_knot) = cij_help;
        
        
        // update all contribution of current knot-combination and person i to column of person i
        help_mat(span::all,i) = help_mat(span::all,i) + (trans(Zij)*D_i*SigmaInv_i*(yi-mu_i))*cij_help;
        
        pos_knot = pos_knot + 1;
      }}
  }
  
  // normalization value per person
  vec cij_norm = sum(cij_mat,1);
  
  // normalize row of each parameter by cij_norm
  for(int e=0; e<pall ;e++){
    help_mat(e,span::all) = help_mat(e,span::all)%trans(1/cij_norm);
  }
  // sum score contributions over all persons per covariate
  vec s = -sum(help_mat,1)+P2;
  
  return s;
}




// [[Rcpp::export]]
double loglikUGPCM2noX(arma::vec alpha,
                    arma::mat Y,
                    int Q,
                    int q,
                    int n,
                    int I,
                    int pall,
                    arma::mat GHweights,
                    arma::vec GHnodes,
                    int pX,
                    arma::mat X,
                    int cores,
                    double lambda) { 
  
  
  // initialize loglikelihood       
  vec f = zeros(n);   
  
  double P2 = accu(alpha%alpha)*lambda;
  
  // initialize design for one person, without random effects
  mat Z = -join_rows(diagmat(ones(q*I)),zeros(q*I,2*pX+2+I));
  mat etai = Z*alpha;
  // current linear predictor without random effects, per item and category
  etai = reshape(etai,q,I);
  
  // create sigma matrix from current parameters
  double co_var = alpha(pall-1-I)*sqrt(alpha(pall-I))*sqrt(alpha(pall-2-I));
  mat sigma = zeros(2,2);        
  sigma(0,0) = alpha(pall-2-I);
  sigma(1,0) = co_var;
  sigma(0,1) = co_var;
  sigma(1,1) = alpha(pall-I);
  // square root of current sigma matrix
  
  // mat sigma12 = trans(chol(sigma));
  mat sigma12;
  try{
    sigma12 = trans(chol(sigma));
  }
  catch(...)
  {
    sigma = sigma + diagmat(ones(2)*0.0001);
    sigma12 = trans(chol(sigma));
  }
  
  
  vec slopes = join_cols(ones<vec>(1), alpha(span(pall-I+1, pall-1)));
  
  // initialize number of different threads
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
  
 mat yi; mat prods_i; int j; int jj; int k; int i; 
  vec rnd_act; mat etaij; mat mui; vec eta_k; vec mu_k;
  // loop through all persons
#pragma omp parallel for private(i, j, jj, k, yi, prods_i, rnd_act, etaij, mui, eta_k, mu_k) shared(f)
  for(i=0; i<n ;i++){

    
    // get response of person i
    yi = Y(span::all,span(i*I,i*I+I-1));
    
    prods_i = ones(Q,Q);
    // loop through all knots, for both random effects
    for(j=0; j<Q ;j++){ 
      for(jj=0; jj<Q ;jj++){ 
        // initialize current random effects or knots
        rnd_act = zeros(2);
        rnd_act(0) = GHnodes(j);
        rnd_act(1) = GHnodes(jj);
        
        rnd_act = sigma12*rnd_act;
        
        
        etaij = (etai + ones(q,I)*rnd_act(0))%((ones(q)*trans(slopes))*exp(rnd_act(1)));
        mui = ones(q+1,I);
        // loop through all items  
        for(k=0; k<I ;k++){
          
          // get eta and mu of person i for current item k
          eta_k = etaij(span::all,k);
          mu_k = responseFun2(eta_k);
          mui(span::all,k) = join_cols(mu_k,1-sum(mu_k,0));
        }
        // create prob for current item, update prob matrix for respective knots
        prods_i(j,jj) = prods_i(j,jj)*prod(prod(mui % yi - (yi-1)));
        
      }
    }
    f(i)  = - log(accu(prods_i%GHweights));
    // accumulate all likelihood contributions, weights by respective weights and probs of knots
  }
  
  
  return sum(f)+P2;
}



// [[Rcpp::export]]
double loglikUGPCM4noX(arma::vec alpha,
                    arma::vec Y,
                    int Q,
                    int q,
                    int n,
                    int I,
                    arma::mat GHweights,
                    arma::vec GHnodes,
                    int cores,
                    arma::vec delta,
                    arma::vec slopes,
                    double lambda) {
  
  
  // initialize loglikelihood       
  vec f = zeros(n);   
  
  vec alpha2 = join_cols(join_cols(delta,alpha),slopes);
  
  double P2 = accu(alpha2%alpha2)*lambda;
  // initialize design for one person, without random effects
  mat Z = -diagmat(ones(q*I));
  mat etai = Z*delta;
  // current linear predictor without random effects, per item and category
  etai = reshape(etai,q,I);
  
  // create sigma matrix from current parameters
  double co_var = alpha(1)*sqrt(alpha(0))*sqrt(alpha(2));
  mat sigma = zeros(2,2);
  sigma(0,0) = alpha(0);
  sigma(1,0) = co_var;
  sigma(0,1) = co_var;
  sigma(1,1) = alpha(2);
  // square root of current sigma matrix
  
  // mat sigma12 = trans(chol(sigma));
  mat sigma12;
  try{
    sigma12 = trans(chol(sigma));
  }
  catch(...)
  {
    sigma = sigma + diagmat(ones(2)*0.0001);
    sigma12 = trans(chol(sigma));
  }
  
  // initialize number of different threads
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
  
  
 vec yi; mat prods_i; int j; int jj; int k; int i; 
  vec rnd_act; mat etaij; mat mui; vec eta_k; vec mu_k;
  vec yi_k; 
  // loop through all persons
#pragma omp parallel for private(i, j, jj, k, yi, yi_k, prods_i, rnd_act, etaij, mui, eta_k, mu_k) shared(f)
  for(i=0; i<n ;i++){
    

    
    // get response of person i
    yi = Y(span(i*q*I,i*q*I+q*I-1));
    prods_i = ones(Q,Q);
    // loop through all knots, for both random effects
    for(j=0; j<Q ;j++){ 
      for(jj=0; jj<Q ;jj++){ 
        // initialize current random effects or knots
        rnd_act = zeros(2);
        rnd_act(0) = GHnodes(j);
        rnd_act(1) = GHnodes(jj);
        rnd_act = sigma12*rnd_act;
        
        
        // loop through all items  
        for(k=0; k<I ;k++){  
          // response of person i for current item k
          yi_k = yi(span(k*q,q+k*q-1));
          yi_k = join_cols(yi_k,1-sum(yi_k,0));
          
          // get eta and mu of person i for current item k
          eta_k = (etai(span::all,k) + ones(q)*(rnd_act(0)))%(slopes(k)*ones(q)*(exp(rnd_act(1))));
          mu_k = responseFun(eta_k);
          mu_k = join_cols(mu_k,1-sum(mu_k,0));
          
          // create prob for current item, update prob matrix for respective knots
          
          prods_i(j,jj) = prods_i(j,jj)*prod(mu_k % yi_k - (yi_k-1));
        }
      }
    }
    f(i)  = -log(accu(prods_i%GHweights));
    // accumulate all likelihood contributions, weights by respective weights and probs of knots
  }
  
  return sum(f)+P2;
}

