// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
double posterior_logit(const arma::vec &params, const Rcpp::List &xx){
  
  
  uvec ind_exp = xx["ind_exp"];
  int p = xx["p"];
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));

  
  mat Xzero = xx["Xzero"];
  mat Xrest = xx["Xrest"];
  mat Xlrest = xx["Xlrest"];
//  mat Xzerof = ss["Xzerof"];
//  mat Xrestf = ss["Xrestf"];
//  mat Xlrestf = ss["Xlrestf"];
  mat S = xx["S"];
  
  vec em1 = exp(-Xzero*bt);
  vec em2 = exp(-Xrest*bt);
  vec em3 = exp(-Xlrest*bt);
  double ret = sum(log(1/(1+em1))) + sum(log(1/(1+em2) - 1/(1+em3)))- as_scalar(0.5*params.t()*S*params);
  
  return(ret);
}


// [[Rcpp::export]]
double ll_logit(const arma::vec &params, const Rcpp::List &xx){
  
  
  uvec ind_exp = xx["ind_exp"];
  int p = xx["p"];
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  
  
  mat Xzero = xx["Xzero"];
  mat Xrest = xx["Xrest"];
  mat Xlrest = xx["Xlrest"];
  //  mat Xzerof = ss["Xzerof"];
  //  mat Xrestf = ss["Xrestf"];
  //  mat Xlrestf = ss["Xlrestf"];
  mat S = xx["S"];
  
  vec em1 = exp(-Xzero*bt);
  vec em2 = exp(-Xrest*bt);
  vec em3 = exp(-Xlrest*bt);
  double ret = sum(log(1/(1+em1))) + sum(log(1/(1+em2) - 1/(1+em3)));
  
  return(ret);
}

// 
// [[Rcpp::export]]
double posterior_logitf(const arma::vec &params, const Rcpp::List &xx){


  mat Xzero = xx["Xzero"];
  mat Xrest = xx["Xrest"];
  mat Xlrest = xx["Xlrest"];

  mat S = xx["S"];
  uvec ind_exp = xx["ind_exp"];
  double b1 = xx["b1"];
  int p = xx["p"];

  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  bt.resize(p);
  bt(p-1) = b1;
  vec em1 = exp(-Xzero*bt);
  vec em2 = exp(-Xrest*bt);
  vec em3 = exp(-Xlrest*bt);

    double ret = sum(log(1/(1+em1))) + sum(log(1/(1+em2) - 1/(1+em3)))- as_scalar(0.5*params.t()*S*params);

  return(ret);
}

// [[Rcpp::export]]
vec dlogis(const arma::vec &x){
  
  vec ret = exp(x)/square(1+exp(x));
  ret.replace(datum::nan, 0);   
  return(ret);
}

// [[Rcpp::export]]
vec plogis(const arma::vec &x){
  
  vec ret = 1/(1+exp(-x));
  
  return(ret);
}

// [[Rcpp::export]]
vec gradf_logit(const arma::vec &params, const Rcpp::List  &xx){
  
  mat Xzero = xx["Xzero"];
  mat Xrest = xx["Xrest"];
  mat Xlrest = xx["Xlrest"];
  //  mat Xzerof = ss["Xzerof"];
  //  mat Xrestf = ss["Xrestf"];
  //  mat Xlrestf = ss["Xlrestf"];
  mat S = xx["S"];
  uvec ind_exp = xx["ind_exp"];
  int p = xx["p"];

  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  
  int qb = bt.n_elem;
  vec c = vec(qb);
  c.ones();
  c(ind_exp) = bt(ind_exp);
  
  // mat C = diagmat(c);
  // Rcout << "The value is " << C;
  vec v1 = Xzero*bt;
  vec v2 = Xrest*bt;
  vec v3 = Xlrest*bt;
  vec ppd = plogis(v2)-plogis(v3);
  
  vec ret = ((Xzero.each_row() % c.t()).t())*(dlogis(v1)/plogis(v1))+
    ((Xrest.each_row() % c.t()).t())*((dlogis(v2))/
      (ppd)) -
        ((Xlrest.each_row() % c.t()).t())*((dlogis(v3))/
          (ppd)) -
            S*params;
  
  
  return(ret);
}
// 
// 
// 
// [[Rcpp::export]]
vec gradf_logitf(const arma::vec &params, const Rcpp::List &xx){

  mat Xzero = xx["Xzero"];
  mat Xrest = xx["Xrest"];
  mat Xlrest = xx["Xlrest"];
   mat Xzerof = xx["Xzerof"];
   mat Xrestf = xx["Xrestf"];
   mat Xlrestf = xx["Xlrestf"];
   mat S = xx["S"];
   uvec ind_exp = xx["ind_exp"];
   double b1 = xx["b1"];
   int p = xx["p"];

  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));

  int qb = bt.n_elem;
  vec c = vec(qb);
  c.ones();
  c(ind_exp) = bt(ind_exp);
  bt.resize(p);
  bt(p-1) = b1;
  // mat C = diagmat(c);

  vec v1 = Xzero*bt;
  vec v2 = Xrest*bt;
  vec v3 = Xlrest*bt;
  vec ppd = plogis(v2)-plogis(v3);

  vec ret = ((Xzerof.each_row() % c.t()).t())*(dlogis(v1)/plogis(v1))+
    ((Xrestf.each_row() % c.t()).t())*((dlogis(v2))/
      (ppd)) -
        ((Xlrestf.each_row() % c.t()).t())*((dlogis(v3))/
          (ppd)) -
            S*params;


  return(ret);
}


// [[Rcpp::export]]
double posterior_probit(const arma::vec &params, const Rcpp::List &xx){
  uvec ind_exp = xx["ind_exp"];
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));

  mat Xzero = xx["Xzero"];
  mat Xrest = xx["Xrest"];
  mat Xlrest = xx["Xlrest"];
  int p = xx["p"];
  
  //  mat Xzerof = ss["Xzerof"];
  //  mat Xrestf = ss["Xrestf"];
  //  mat Xlrestf = ss["Xlrestf"];
  mat S = xx["S"];

  vec m1 = normcdf(Xzero*bt);
  vec m2 = normcdf(Xrest*bt);
  vec m3 = normcdf(Xlrest*bt);
  // Rcout << "The value is " << as_scalar(0.5*params.t()*S*params);
  double ret = sum(log(m1)) + sum(log(m2 -m3))- as_scalar(0.5*params.t()*S*params);

  return(ret);
}

// [[Rcpp::export]]
double ll_probit(const arma::vec &params, const Rcpp::List &xx){
  uvec ind_exp = xx["ind_exp"];
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  
  mat Xzero = xx["Xzero"];
  mat Xrest = xx["Xrest"];
  mat Xlrest = xx["Xlrest"];
  int p = xx["p"];
  
  //  mat Xzerof = ss["Xzerof"];
  //  mat Xrestf = ss["Xrestf"];
  //  mat Xlrestf = ss["Xlrestf"];
  mat S = xx["S"];
  
  vec m1 = normcdf(Xzero*bt);
  vec m2 = normcdf(Xrest*bt);
  vec m3 = normcdf(Xlrest*bt);
  // Rcout << "The value is " << as_scalar(0.5*params.t()*S*params);
  double ret = sum(log(m1)) + sum(log(m2 -m3));
  
  return(ret);
}

// [[Rcpp::export]]
vec gradf_probit(const arma::vec &params, const Rcpp::List &xx){
  
  mat Xzero = xx["Xzero"];
  mat Xrest = xx["Xrest"];
  mat Xlrest = xx["Xlrest"];
  //  mat Xzerof = ss["Xzerof"];
  //  mat Xrestf = ss["Xrestf"];
  //  mat Xlrestf = ss["Xlrestf"];
  mat S = xx["S"];
  uvec ind_exp = xx["ind_exp"];
  int p = xx["p"];
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  
  int qb = bt.n_elem;
  vec c = vec(qb);
  c.ones();
  c(ind_exp) = bt(ind_exp);
  
  vec v1 = Xzero*bt;
  vec v2 = Xrest*bt;
  vec v3 = Xlrest*bt;
  vec ppd = normcdf(v2)-normcdf(v3);
  
  vec ret = ((Xzero.each_row() % c.t()).t())*(normpdf(v1)/normcdf(v1))+
    ((Xrest.each_row() % c.t()).t())*((normpdf(v2))/
      (ppd)) -
        ((Xlrest.each_row() % c.t()).t())*((normpdf(v3))/
          (ppd)) -
            S*params;
  
  
  return(ret);
}

// [[Rcpp::export]]
double posterior_probitf(const arma::vec &params, const Rcpp::List &xx){
  uvec ind_exp = xx["ind_exp"];
  double b1 = xx["b1"];
  int p = xx["p"];
  
  
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  bt.resize(p);
  bt(p-1) = b1;
  
  mat Xzero = xx["Xzero"];
  mat Xrest = xx["Xrest"];
  mat Xlrest = xx["Xlrest"];

  mat S = xx["S"];
  
  vec m1 = normcdf(Xzero*bt);
  vec m2 = normcdf(Xrest*bt);
  vec m3 = normcdf(Xlrest*bt);
  // Rcout << "The value is " << as_scalar(0.5*params.t()*S*params);
  double ret = sum(log(m1)) + sum(log(m2 -m3))- as_scalar(0.5*params.t()*S*params);
  
  return(ret);
}


// [[Rcpp::export]]
vec gradf_probitf(const arma::vec &params, const Rcpp::List &xx){
  
  mat Xzero = xx["Xzero"];
  mat Xrest = xx["Xrest"];
  mat Xlrest = xx["Xlrest"];
  mat Xzerof = xx["Xzerof"];
  mat Xrestf = xx["Xrestf"];
  mat Xlrestf = xx["Xlrestf"];
  mat S = xx["S"];
  uvec ind_exp = xx["ind_exp"];
  double b1 = xx["b1"];
  int p = xx["p"];
  
  
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  
  vec c = vec(p-1);
  c.ones();
  c(ind_exp) = bt(ind_exp);
  bt.resize(p);
  bt(p-1) = b1;
  // mat C = diagmat(c);
  
  vec v1 = Xzero*bt;
  vec v2 = Xrest*bt;
  vec v3 = Xlrest*bt;
  vec ppd = normcdf(v2)-normcdf(v3);
  
  vec ret = ((Xzerof.each_row() % c.t()).t())*(normpdf(v1)/normcdf(v1))+
    ((Xrestf.each_row() % c.t()).t())*((normpdf(v2))/
      (ppd)) -
        ((Xlrestf.each_row() % c.t()).t())*((normpdf(v3))/
          (ppd)) -
            S*params;
  
  
  return(ret);
}




// [[Rcpp::export]]
vec dcll(const arma::vec &x){
  
  // vec ret = x-exp(x);
  vec ret = exp(x-exp(x));
 // vec ret = exp(x);
  return(ret);
}

// [[Rcpp::export]]
vec pcll(const arma::vec &x){
  
  vec ret = 1-exp(-exp(x));
  
  return(ret);
}



// [[Rcpp::export]]
double posterior_cll(const arma::vec &params, const Rcpp::List &xx){
  uvec ind_exp = xx["ind_exp"];
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));

  
  mat Xzero = xx["Xzero"];
  mat Xrest = xx["Xrest"];
  mat Xlrest = xx["Xlrest"];
  //  mat Xzerof = ss["Xzerof"];
  //  mat Xrestf = ss["Xrestf"];
  //  mat Xlrestf = ss["Xlrestf"];
  mat S = xx["S"];
  
  vec m1 = pcll(Xzero*bt);
  vec m2 = pcll(Xrest*bt);
  vec m3 = pcll(Xlrest*bt);
  

  double ret = sum(log(m1)) + sum(log(m2 -m3))- as_scalar(0.5*params.t()*S*params);
  
  return(ret);
}


// [[Rcpp::export]]
double ll_cll(const arma::vec &params, const Rcpp::List &xx){
  uvec ind_exp = xx["ind_exp"];
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  
  
  mat Xzero = xx["Xzero"];
  mat Xrest = xx["Xrest"];
  mat Xlrest = xx["Xlrest"];
  //  mat Xzerof = ss["Xzerof"];
  //  mat Xrestf = ss["Xrestf"];
  //  mat Xlrestf = ss["Xlrestf"];
  mat S = xx["S"];
  
  vec m1 = pcll(Xzero*bt);
  vec m2 = pcll(Xrest*bt);
  vec m3 = pcll(Xlrest*bt);
  
  
  double ret = sum(log(m1)) + sum(log(m2 -m3));
  
  return(ret);
}

// [[Rcpp::export]]
vec gradf_cll(const arma::vec &params, const Rcpp::List &xx){
  
  mat Xzero = xx["Xzero"];
  mat Xrest = xx["Xrest"];
  mat Xlrest = xx["Xlrest"];

  mat S = xx["S"];
  uvec ind_exp = xx["ind_exp"];
  int p = xx["p"];
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  
  int qb = bt.n_elem;
  vec c = vec(qb);
  c.ones();
  c(ind_exp) = bt(ind_exp);
  
  // Rcout << "The value is " << bt;
  vec v1 = Xzero*bt;
  vec v2 = Xrest*bt;
  vec v3 = Xlrest*bt;
  vec ppd = pcll(v2)-pcll(v3);
  
  vec ret = ((Xzero.each_row() % c.t()).t())*(dcll(v1)/pcll(v1))+
    ((Xrest.each_row() % c.t()).t())*((dcll(v2))/
      (ppd)) -
        ((Xlrest.each_row() % c.t()).t())*((dcll(v3))/
          (ppd)) -
            S*params;
  
  
  return(ret);
}

// [[Rcpp::export]]
double posterior_cllf(const arma::vec &params, const Rcpp::List &xx){
  uvec ind_exp = xx["ind_exp"];
  double b1 = xx["b1"];
  int p = xx["p"];
  
  
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  bt.resize(p);
  bt(p-1) = b1;
  
  mat Xzero = xx["Xzero"];
  mat Xrest = xx["Xrest"];
  mat Xlrest = xx["Xlrest"];
  
  mat S = xx["S"];
  
  vec m1 = pcll(Xzero*bt);
  vec m2 = pcll(Xrest*bt);
  vec m3 = pcll(Xlrest*bt);
  // Rcout << "The value is " << as_scalar(0.5*params.t()*S*params);
  double ret = sum(log(m1)) + sum(log(m2 -m3))- as_scalar(0.5*params.t()*S*params);
  
  return(ret);
}


// [[Rcpp::export]]
vec gradf_cllf(const arma::vec &params, const Rcpp::List &xx){
  
  mat Xzero = xx["Xzero"];
  mat Xrest = xx["Xrest"];
  mat Xlrest = xx["Xlrest"];
  mat Xzerof = xx["Xzerof"];
  mat Xrestf = xx["Xrestf"];
  mat Xlrestf = xx["Xlrestf"];
  mat S = xx["S"];
  uvec ind_exp = xx["ind_exp"];
  double b1 = xx["b1"];
  int p = xx["p"];
  
  
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  
  int qb = bt.n_elem;
  vec c = vec(qb);
  c.ones();
  c(ind_exp) = bt(ind_exp);
  bt.resize(p);
  bt(p-1) = b1;
  // mat C = diagmat(c);
  
  vec v1 = Xzero*bt;
  vec v2 = Xrest*bt;
  vec v3 = Xlrest*bt;
  
  vec ppd = pcll(v2)-pcll(v3);
  
  
  vec ret = ((Xzerof.each_row() % c.t()).t())*(dcll(v1)/pcll(v1))+
    ((Xrestf.each_row() % c.t()).t())*((dcll(v2))/
      (ppd)) -
        ((Xlrestf.each_row() % c.t()).t())*((dcll(v3))/
          (ppd)) -
            S*params;
  
  
  return(ret);
}



// 
// 
// // [[Rcpp::export]]
// vec dcll(const arma::vec &x){
//   
//   // vec ret = x-exp(x);
//   vec ret = exp(x-exp(x));
//   // vec ret = exp(x);
//   return(ret);
// }

// // [[Rcpp::export]]
// vec cdflogit(const arma::vec &x, const arma::vec &y){
//   int n = y.n_elem;
//   vec ret = vec(n);
//   for (int i=0; i<n; i++) {
//     if(y[i]==-1){
//       ret[i] = 0;
//     } else if(y[i]==2){
//       ret[i] = 1;
//     } else {
//       ret[i] =1-exp(-exp(x[i]));
//     }
//     }
//   return(ret);
// }
// 
// 
// // [[Rcpp::export]]
// vec pdflogit(const arma::vec &x, const arma::vec &y){
//   int n = y.n_elem;
//   vec ret = vec(n);
//   
//   for (int i=0; i<n; i++) {
//     if(y[i]==-1) {
//     ret[i] = 0;
//     } else if(y[i]==2){
//       ret[i] = 0;
//   } else {
//     double ex= exp(x[i]);
//     ret[i] = ex/pow(ex+1,2);
//   }
//   }
//   return(ret);
// }

// // [[Rcpp::export]]
// double posterior_logitd(const arma::vec &params, const Rcpp::List &xx){
//   
//   
//   uvec ind_exp = xx["ind_exp"];
//   int p = xx["p"];
//   
//   arma::vec bt= params;
//   bt.elem(ind_exp) = exp(params.elem(ind_exp));
//   
//   vec y = xx["y"];
//   mat X = xx["X"];
//   mat Xl = xx["Xl"];
// 
//   mat S = xx["S"];
//   
//   // vec em1 = exp(-Xzero*bt);
//   // vec em2 = exp(-Xrest*bt);
//   // vec em3 = exp(-Xlrest*bt);
//   // double ret =  sum(log(1/(1+em2) - 1/(1+em3)))- as_scalar(0.5*params.t()*S*params);
//   double ret = sum(log(cdflogit(X*bt, y) - cdflogit(Xl*bt, y)));
//   return(ret);
// }
// 
// // [[Rcpp::export]]
// vec gradf_logitd(const arma::vec &params, const Rcpp::List  &xx){
//   
//   vec y = xx["y"];
//   mat X = xx["X"];
//   mat Xl = xx["Xl"];
//   int n = y.n_elem;
//   mat S = xx["S"];
//   
//   uvec ind_exp = xx["ind_exp"];
//   int p = xx["p"];
//   
//   arma::vec bt= params;
//   bt.elem(ind_exp) = exp(params.elem(ind_exp));
//   
//   int qb = bt.n_elem;
//   vec c = vec(qb);
//   c.ones();
//   c(ind_exp) = bt(ind_exp);
//   
//   // mat C = diagmat(c);
//   // Rcout << "The value is " << C;
//   vec v2 = X*bt;
//   vec v3 = Xl*bt;
//   vec ppd = cdflogit(v2, y)-cdflogit(v3, y-1);
//   
//   
//   vec pdf2 = pdflogit(v2, y);
//   vec pdf3 = pdflogit(v3, y-1);
//   
//  
//   vec ret = ((X.each_row() % c.t()).t())*((pdf2)/
//     (ppd)) -
//       ((Xl.each_row() % c.t()).t())*((pdf3)/
//         (ppd)) -
//           S*params;
//     
// 
//   
//   
//   return(ret);
// }



// [[Rcpp::export]]
double posterior_logitd(const arma::vec &params, const Rcpp::List &xx){
  
  
  mat S = xx["S"];
  mat Xzero = xx["Xzero"];
  mat Xlrest = xx["Xlrest"];
  // mat Xlone = xx["Xlone"];
  // mat Xone = xx["Xone"];
  mat Xrest = xx["Xrest"];
  
  uvec ind_exp = xx["ind_exp"];
  // int p = xx["p"];
  int p = xx["p"];
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  // bt.resize(p);
  // bt(2) = 1000;
  // bt.resize(p);
  bt.insert_rows( 2, 1 );
  bt(2) = 1000;  
    
    vec em1 = Xzero*bt;
    vec em2 = Xrest*bt;
    vec em3 = Xlrest*bt;
    // vec em4 = Xlone*bt;
    
    double ret = sum(log(plogis(em1))) + sum(log(plogis(em2) - plogis(em3))) - as_scalar(0.5*params.t()*S*params);
    

    return(ret);

}


// [[Rcpp::export]]
vec gradf_logitd(const arma::vec &params, const Rcpp::List  &xx){
  
  mat Xzero = xx["Xzero"];
  mat Xrest = xx["Xrest"];
  mat Xlrest = xx["Xlrest"];
  mat Xzerof = xx["Xzerof"];
  mat Xrestf = xx["Xrestf"];
  mat Xlrestf = xx["Xlrestf"];
  int p = xx["p"];
  uvec ind_exp = xx["ind_exp"];
  mat S = xx["S"];
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  
  // int qb = bt.n_elem;
  int qb = bt.n_elem;
  vec c = vec(qb);
  c.ones();
  c(ind_exp) = bt(ind_exp);
  // bt.resize(p);
  // Rcout << "The value is " << c;
  
  bt.insert_rows( 2, 1 );
  bt(2) = 1000;
  // mat C = diagmat(c);
  // Rcout << "The value is " << bt;
  
  vec v1 = Xzero*bt;
  vec v2 = Xrest*bt;
  vec v3 = Xlrest*bt;
  
  vec ppd = plogis(v2)-plogis(v3);
  // Rcout << "The value is " << v2;
  
                 
  vec ret = ((Xzerof.each_row() % c.t()).t())*(dlogis(v1)/plogis(v1))+
    ((Xrestf.each_row() % c.t()).t())*((dlogis(v2))/
      (ppd)) -
        ((Xlrestf.each_row() % c.t()).t())*((dlogis(v3))/
          (ppd)) - S*params;

  return(ret);
  
}



// [[Rcpp::export]]
double posterior_logitd2(const arma::vec &params, const Rcpp::List &xx){
  
  
  mat S = xx["S"];
  mat Xzero = xx["Xzero"];
  mat Xlrest = xx["Xlrest"];
  // mat Xlone = xx["Xlone"];
  // mat Xone = xx["Xone"];
  mat Xrest = xx["Xrest"];
  
  uvec ind_exp = xx["ind_exp"];
  // int p = xx["p"];
  int p = xx["p"];
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  // bt.resize(p);
  // bt(2) = 1000;
  // bt.resize(p);
  bt.insert_rows( 2, 1 );
  bt(2) = 1000;
  bt.insert_rows( 5, 1 );
  bt(5) = 1000;
  // bt.insert_rows( 8, 1 );
  // bt(8) = 1000;
  // bt.insert_rows( 11, 1 );
  // bt(11) = 1000;
  // Rcout << "The value is " << bt;
  
  vec em1 = Xzero*bt;
  vec em2 = Xrest*bt;
  vec em3 = Xlrest*bt;
  // vec em4 = Xlone*bt;
  // Rcout << "The value is " << em2;
  
  double ret = sum(log(plogis(em1))) + sum(log(plogis(em2) - plogis(em3)));  //- as_scalar(0.5*params.t()*S*params);
  
  
  return(ret);
  
}


// [[Rcpp::export]]
vec gradf_logitd2(const arma::vec &params, const Rcpp::List  &xx){
  
  mat Xzero = xx["Xzero"];
  mat Xrest = xx["Xrest"];
  mat Xlrest = xx["Xlrest"];
  mat Xzerof = xx["Xzerof"];
  mat Xrestf = xx["Xrestf"];
  mat Xlrestf = xx["Xlrestf"];
  int p = xx["p"];
  uvec ind_exp = xx["ind_exp"];
  mat S = xx["S"];
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  
  // int qb = bt.n_elem;
  int qb = bt.n_elem;
  vec c = vec(qb);
  c.ones();
  c(ind_exp) = bt(ind_exp);
  // bt.resize(p);
  // Rcout << "The value is " << c;
  
  bt.insert_rows( 2, 1 );
  bt(2) = 100000;
  bt.insert_rows( 5, 1 );
  bt(5) = 1000;
  // bt.insert_rows( 8, 1 );
  // bt(8) = 1000;
  // bt.insert_rows( 11, 1 );
  // bt(11) = 1000;
  // mat C = diagmat(c);
  // Rcout << "The value is " << bt;
  
  vec v1 = Xzero*bt;
  vec v2 = Xrest*bt;
  vec v3 = Xlrest*bt;
  
  vec ppd = plogis(v2)-plogis(v3);
  // Rcout << "The value is " <<plogis(v2);
  
  
  vec ret = ((Xzerof.each_row() % c.t()).t())*(dlogis(v1)/plogis(v1))+
    ((Xrestf.each_row() % c.t()).t())*((dlogis(v2))/
      (ppd)) -
        ((Xlrestf.each_row() % c.t()).t())*((dlogis(v3))/
          (ppd)); 
    // - S*params;
  
  return(ret);
  
}



// [[Rcpp::export]]
double posterior_logitd3(const arma::vec &params, const Rcpp::List &xx){
  
  
  mat S = xx["S"];
  mat X = xx["Xnew"];
  mat Xl = xx["Xlnew"];

  uvec zero_ind = xx["zero_ind"];
  uvec max_ind = xx["max_ind"];
  uvec ind_exp = xx["ind_exp"];
  // int p = xx["p"];
  int p = xx["p"];
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));

  vec Fh = plogis(X*bt);
  vec Fhl = plogis(Xl*bt);
  
  Fhl(zero_ind).fill(0);
  Fh(max_ind).fill(1);
  
  double ret = sum(log(Fh - Fhl))- as_scalar(0.5*params.t()*S*params);
  
  
  return(ret);
  
}


// [[Rcpp::export]]
double ll_logitd3(const arma::vec &params, const Rcpp::List &xx){
  
  
  mat S = xx["S"];
  mat X = xx["Xnew"];
  mat Xl = xx["Xlnew"];
  
  uvec zero_ind = xx["zero_ind"];
  uvec max_ind = xx["max_ind"];
  uvec ind_exp = xx["ind_exp"];
  // int p = xx["p"];
  int p = xx["p"];
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  
  vec Fh = plogis(X*bt);
  vec Fhl = plogis(Xl*bt);
  
  Fhl(zero_ind).fill(0);
  Fh(max_ind).fill(1);
  
  double ret = sum(log(Fh - Fhl));
  
  
  return(ret);
  
}

// [[Rcpp::export]]
vec gradf_logitd3(const arma::vec &params, const Rcpp::List  &xx){
  
  mat X = xx["Xnew"];
  mat Xl = xx["Xlnew"];

  int p = xx["p"];
  uvec ind_exp = xx["ind_exp"];
  mat S = xx["S"];
  uvec zero_ind = xx["zero_ind"];
  uvec max_ind = xx["max_ind"];
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  
  // int qb = bt.n_elem;
  int qb = bt.n_elem;
  vec c = vec(qb);
  c.ones();
  c(ind_exp) = bt(ind_exp);

  vec v1 = X*bt;
  vec v2 = Xl*bt;

  vec Fh = plogis(v1);
  vec Fhl = plogis(v2);
  
//need to correct for conditions
  Fhl(zero_ind).fill(0);
  Fh(max_ind).fill(1);

 vec ppd = Fh - Fhl;
  
  vec ret = (((X.each_row() % c.t()).t())*(dlogis(v1)/ppd )-
    ((Xl.each_row() % c.t()).t())*(dlogis(v2)/ppd ))-
    S*params;
  

  
  return(ret);
  
}