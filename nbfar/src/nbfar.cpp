// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <stdlib.h>
#include <stdexcept>
#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
//
// Task: converggence = |dev - dev_{old}|/(|dev| + 0.1) < epsilon  = 1e-8, initial value of thhe null model;
// likelihood =  2|Saturated model - 2fit moodel|
//
//
// nbrrr_likelihood, update_mu_alpha, grad_mu_nb, get_sv
//

// Model diagnostic to verify deviance calculation
// i = 30
// fit = MASS::glm.nb(Y[,i, drop = FALSE]~  offset(ofset[,i, drop = FALSE]) )
//   th = fit$theta
//   mu = exp(coef(fit))
//   fit$null.deviance
//   y = Y[,i, drop = FALSE]
// loglik <- function(th, mu, y) sum((lgamma(th + y) - lgamma(th) +
//   th * log(th) + y * log(mu + (y == 0)) -
//   (th + y) * log(th + mu)))
//   2*loglik(th,y, y) - 2*loglik(th,mu, y)
//   fit$null.deviance - (2*loglik(th,y, y) - 2*loglik(th,mu, y))
// [[Rcpp::export]]
double nb_dev(arma::mat Y, arma::mat MU, arma::vec Phi,
              arma::mat naind) {
  arma::mat T1 = Y.each_row() + Phi.t();
  T1 = Y%log(Y+(Y==0))-Y%log(MU) - T1%(log(T1) - log(MU.each_row() + Phi.t()));
  T1.elem( find_nonfinite(T1) ).zeros();
  return(2*accu(T1%naind)/accu(naind));
  // return(2*accu(T1.elem(find(naind ==1)) ));
}



// tem = 1:q
//   for (i in 1:ncol(Y)) {
//   tem[i] = norm(crossprod(X, X*(Y[,i] + 1)),"2")/2
// }
//  sanity chheck  all(get_sc(X,Y) == X*Y[,1])
// [[Rcpp::export]]
double get_sc(arma::mat X, arma::mat Y){
  int iter, q = Y.n_cols; Y = Y + 1;
  Y.elem(find_nonfinite(Y)).zeros();
  arma::vec tem(Y.n_cols);
  for(iter = 0; iter < q; iter++){
    tem(iter) = norm(X.t()*(X.each_col()%Y.col(iter)),2);
  }
  return(tem.max()/2.0);
}



// eta = cbind(1,X)%*%C0
//   Phi <- rep(0.3,q)
//   b = grad_eta_nb(Y, eta,Phi)
//   eta1 = exp(eta);
// tem = t(eta1) + Phi;
// tem1 = t(eta1) * Phi;
// tem1 = tem1/tem
//   tem = Phi/tem
//   a = -Y*t(tem) + t(tem1)
//   all(a == b)

// help to compute: Grad XTgrad_eta_nb; ZTgrad_eta_nb;
// [[Rcpp::export]]
arma::mat grad_eta_nb(arma::mat Y, arma::mat eta, arma::vec Phi){
  eta = exp(eta);
  arma::mat tem = eta.each_row() + Phi.t();
  arma::mat tem1 = eta.each_row()%Phi.t();
  tem1 = tem1/tem;
  tem = Phi.t()/tem.each_row();
  return(-Y%tem + tem1);
}


// [[Rcpp::export]]
arma::mat grad_mu_nb(const arma::mat &Y, const arma::mat &mu,
                     const arma::vec &Phi){
  // arma::mat grad_mu_nb(arma::mat Y, arma::mat mu,
  //                      arma::vec Phi){
  arma::mat tem  = mu, tem1 = mu;
  tem.each_row() += Phi.t();
  tem1.each_row() %= Phi.t();
  tem1 /= tem;
  tem = Phi.t()/tem.each_row();
  // tem = -Y%tem  + tem1;
  tem %= -Y;  tem += tem1;
  tem.elem(find_nonfinite(tem) ).zeros();
  return(tem);
}

// [[Rcpp::export]]
arma::mat grad_mu_nb_uv(const arma::mat &Y, const arma::mat &mu,
                     const arma::vec &Phi, arma::mat &d2l){
  // arma::mat grad_mu_nb(arma::mat Y, arma::mat mu,
  //                      arma::vec Phi){
  arma::mat tem  = mu, tem1 = mu;
  tem.each_row() += Phi.t();
  tem1.each_row() %= Phi.t();
  tem1 /= tem;
  d2l = (Y+1)%(tem1/tem);
  tem = Phi.t()/tem.each_row();
  // tem = -Y%tem  + tem1;
  tem %= -Y;  tem += tem1;
  tem.elem(find_nonfinite(tem) ).zeros();
  return(tem);
}




// // [[Rcpp::export]]
// arma::vec xxx(arma::vec Phi){
//   Phi.subvec(0,3) = linspace(2,5,4);
//   return(Phi);
// }



// mu = exp(X0 %*% C0); phi = rep(20, q)
// phi2 = update_mu_phi(Y, mu, phi)
// phi2 = update_mu_phi(Y, mu, phi2)


// Newton-Raphson update for the dispersion parameter
// use need to pass mean parameter mu and current update
// [[Rcpp::export]]
arma::vec update_mu_phi(arma::mat Y, arma::mat mu, arma::vec Phi){
  // compute T1
  arma::vec tem,a,b;
  arma::mat T1_g = zeros(size(Y));
  arma::mat T1_h = zeros(size(Y));
  int i,j,k,q = Y.n_cols;
  for(i=0; i < q; i++){
    k = max(Y.col(i)); tem = ones(k+1);
    tem.subvec(1,k) = linspace(0,k-1,k) + Phi(i);
    a = 1/tem; b = a/tem; a(0) = 0; b(0) = 0;
    a = cumsum(a); b = -1*cumsum(b);
    for( j=0; j < (int) Y.n_rows; j++){
      T1_g(j,i) = a(Y(j,i));
      T1_h(j,i) = b(Y(j,i));
    }
  }
  // cout << sum(T1_h < 0,0)  << endl;

  // compute T2
  arma::mat n1 = mu - Y;
  arma::mat d1 = mu.each_row() + Phi.t();
  arma::mat T2 = n1/d1;

  // compute T3
  arma::mat T3_g = -1*log(d1);
  T3_g.each_row() += log(Phi.t());
  arma::mat T3_h = (mu.each_row()/Phi.t())/d1;

  // arma::mat grad_phi = T1_g + T2 + T3_g;
  // arma::mat hess_phi = T1_h + -T2/d1 + T3_h;
  // cout << "a " << endl;
  tem = conv_to<arma::vec>::from(sum(T1_g + T2 + T3_g,0)/sum(T1_h -T2/d1 + T3_h,0));
  // cout << sum(T1_g + T2 + T3_g,0) << endl;
  // cout << sum(T1_h -T2/d1 + T3_h,0) << endl;
  a = Phi - tem;
  // if(a.has_nan()) {
  //   cout << a << "  "<< Phi << "  "<< accu(Y) << "  " <<accu(mu)<< std::endl;
  //   throw std::runtime_error("error aditya");
  // }
  return(a);
}


// phi2 = update_mu_alpha(Y, mu, phi)
//   phi2 = update_mu_alpha(Y, mu, phi2)


// Newton-Raphson update for the dispersion parameter
// after reparameterization of the variable
// use need to pass mean parameter mu and current update
// [[Rcpp::export]]
arma::vec update_mu_alpha(const arma::mat &Y, const arma::mat &mu,
                          const arma::vec &Phi, const arma::mat &naind){
  // arma::vec update_mu_alpha( arma::mat Y,  arma::mat mu,
  //                            arma::vec Phi,  arma::mat naind){
  // compute T1
  arma::vec alpha = 1/Phi;
  arma::vec tem,a,b;
  arma::mat T1_g = zeros(size(Y));
  arma::mat T1_h = zeros(size(Y));
  int i,j,k,q = Y.n_cols;
  for(i=0; i < q; i++){
    k = max(Y.col(i)); tem = ones(k+1);
    tem.subvec(1,k) = linspace(0,k-1,k)/(linspace(0,k-1,k)*alpha(i) + 1);
    a = tem; b = pow(tem,2); a(0) = 0; b(0) = 0;
    a = cumsum(a); b = -1*cumsum(b);

    for( j=0; j < (int) Y.n_rows; j++){
      T1_g(j,i) = a(Y(j,i));
      T1_h(j,i) = b(Y(j,i));
    }
  }
  // cout << a  << endl;
  // cout << b  << endl;

  // compute T2
  arma::mat n1 = Y.each_row()%alpha.t();
  arma::mat d1 = mu.each_row()%alpha.t();
  arma::mat d2 = d1 + 1;
  arma::mat n2 = mu.each_row()/alpha.t();
  arma::mat T2_g = n2%(n1+1)/d2;

  n2 = d2.each_row()%alpha.t();
  n2 = pow(n2,2);
  arma::mat T2_h = ((d1%(n1 + 3) + 2)%mu)/n2;


  // compute T3
  arma::mat t1 = log(d2);
  arma::mat T3_g = t1.each_row()/pow(alpha.t(),2);
  arma::mat T3_h = 2*(t1.each_row()/pow(alpha.t(),3));

  // arma::mat grad_phi = T1_g + T2 + T3_g;
  // arma::mat hess_phi = T1_h + -T2/d1 + T3_h;
  // cout << "a " << endl;
  n1 = T1_g - T2_g + T3_g;
  d1 = T1_h +T2_h - T3_h;
  tem = conv_to<arma::vec>::from(sum(n1%naind,0)/sum(d1%naind,0));
  // cout << sum(T1_g + T2 + T3_g,0) << endl;
  // cout << sum(T1_h -T2/d1 + T3_h,0) << endl;
  // a = 1/(alpha - tem);
  // if(a.has_nan()) {
  //   cout << a << "  "<< Phi << "  "<< alpha << "  " <<accu(mu)<< std::endl;
  //   throw std::runtime_error("error aditya");
  // }
  a = 1/(alpha - tem);
  a.elem(find(a<0)).fill(1e-6);
  return(a);
}


// y = Y[,i, drop = FALSE]
// loglik <- function(th, mu, y) sum((lgamma(th + y) - lgamma(th) + th * log(th) + y *
//   log(mu + (y == 0)) - (th + y) * log(th + mu)))
//
//   fit$null.deviance
//   loglik(th,mu, Y[,i, drop = FALSE])
//
//
//   nbrrr_likelihood(y,mu*matrix(1,nrow = nrow(Y),1), log(mu)*matrix(1,nrow = nrow(Y),1),
//                    matrix(th,1,1), (y !=0)+0)
// [[Rcpp::export]]
arma::vec nbrrr_likelihood(const arma::mat &Y, const arma::mat &MU,
                           const arma::mat &ETA,
                           const arma::vec &Phi, const arma::mat &naind) {
  // arma::vec nbrrr_likelihood( arma::mat Y,  arma::mat MU,
  //                             arma::mat ETA,
  //                             arma::vec Phi,  arma::mat naind) {
  // compute T1
  arma::vec tem,a, out(3);
  arma::mat T1 = zeros(size(Y));
  int i,j,k, q = Y.n_cols;
  for(i=0; i < q; i++){
    k = max(Y.col(i)); tem = ones(k+1);
    tem.subvec(1,k) = log(linspace(0,k-1,k) + Phi(i));
    a = tem; a(0) = 0;
    a = cumsum(a);
    for( j=0; j < (int) Y.n_rows; j++)
      T1(j,i) = a(Y(j,i));
  }

  // compute T3
  arma::mat T3 = (Y.each_row()+Phi.t())%(log(1+ MU.each_row()/Phi.t())); // last
  arma::mat b = T1 + Y%(ETA.each_row() - log(Phi.t())) - T3; // likelohood
  b.elem( find_nonfinite(b) ).zeros();
  b= b%naind;
  out(2) = accu(naind);
  out(0) = accu(b)/out(2);
  out(1) = accu(pow(b - out(0),2)%naind)/out(2);
  return(out);
}


// // [[Rcpp::export]]
// arma::vec nbrrr_likelihoodx(arma::mat Y, arma::mat MU, arma::mat ETA,
//                             arma::vec Phi, arma::mat naind) {
//   // compute T1
//   arma::vec tem,a, out(3);
//   arma::mat T1 = zeros(size(Y)),tem_mat;
//   // int i,j,k;
//   tem_mat = Y.each_row()+Phi.t();
//   T1 = lgamma(tem_mat); // - Y_l;
//   T1.each_row() -= lgamma(Phi.t());
//   // for(i=0; i < Y.n_cols; i++){
//   //   k = max(Y.col(i)); tem = ones(k+1);
//   //   tem.subvec(1,k) = log(linspace(0,k-1,k) + Phi(i));
//   //   a = tem; a(0) = 0;
//   //   a = cumsum(a);
//   //   for( j=0; j < Y.n_rows; j++)
//   //     T1(j,i) = a(Y(j,i));
//   // }
//
//   // compute T3
//   arma::mat T3 = tem_mat%(log(1+ MU.each_row()/Phi.t())); // last
//   arma::mat b = T1 + Y%(ETA.each_row() - log(Phi.t())) - T3; // likelohood
//   b.elem( find_nonfinite(b) ).zeros();
//   b= b%naind;
//   out(2) = accu(naind);
//   out(0) = accu(b)/out(2);
//   out(1) = accu(pow(b - out(0),2)%naind)/out(2);
//   return(out);
// }



// [[Rcpp::export]]
arma::uvec mySdiff(arma::uvec x, arma::uvec y){
  // cout<<zeros<mat>(4,5);
  for (int j = 0; j < (int) y.n_elem; j++)
    x = x.elem(find(x != y(j)));
  return(x);
}





// X;Y;Z;O;r; Zini, PhiIni; ndev[ for convergence critteria]
// [[Rcpp::export]]
Rcpp::List nbrrr_cpp(arma::mat Y, arma::mat X0, int rnk, arma::vec cindex,
                     arma::mat ofset,  arma::mat  Zini, arma::vec PhiIni,
                     arma::mat  Cini,
                     Rcpp::List control, int msind, arma::mat naind){
  bool converged=false;
  int pt = X0.n_cols, q = Y.n_cols,  maxit = control["initmaxit"];
  int cObj = control["objI"];
  Rcpp::List out;
  // define epsilon to see convergence of the data
  double epsilon = control["initepsilon"];  //epsilon = epsilon*ndev;


  arma::uvec cIndex =  arma::conv_to< uvec >::from(cindex-1);
  arma::uvec cIndexC = mySdiff(linspace<uvec>(0,pt-1,pt), cIndex);
  arma::uvec t4=find(naind==0);


  arma::mat MU(size(Y)), ETA(size(Y)),X2 = X0.cols(cIndexC), grad_mu;
  arma::mat Res,Ct,X1 = X0.cols(cIndex); // X1 is the control variable
  double sb = get_sc(X1, Y), sc = get_sc(X2, Y);
  X2 = X2/sc; X1 = X1/sb;


  arma::mat C = zeros<mat>(pt,q);
  C.rows(cIndex) = Zini; C.rows(cIndexC)  = Cini;
  arma::vec Phi = PhiIni, Phi2 = 1/PhiIni;
  rnk = rnk-1; // adjust input rank to extract

  wall_clock timer;
  // defining required variable used in the loop
  int j,iter;
  double elp;
  arma::mat C_temp1=C.rows(cIndexC),C_temp2=C.rows(cIndex);
  arma::mat C_temp,Ut,Vt;
  arma::vec diffobj, obj,dt, objtem(3);
  diffobj.zeros(maxit);obj.zeros(maxit+1);



  // MU = familyLinkinv3(ofset + X3*C,t1,t2,t3);
  ETA = ofset + X0*C;
  MU = exp(ETA);
  if(msind == 1) MU.elem(t4).zeros();

  // nbrrr_likelihood: compute likelihood function
  if(cObj!=0){
    obj(0) = nb_dev(Y, MU, Phi, naind);
  } else {
    objtem = nbrrr_likelihood(Y, MU, ETA, Phi,naind);
    obj(0) = objtem(0);
  }

  //
  // cout << obj(0) << ' ' << accu(MU) << ' ' << accu(Phi) << ' ' << accu(naind) << std::endl;

  timer.tic();
  for(iter = 1; iter < maxit; iter++){
    C_temp = C;

    // C-step:
    grad_mu = grad_mu_nb(Y, MU, Phi);
    if(msind == 1) grad_mu.elem(t4).zeros();
    Ct = C.rows(cIndexC) -   X2.t()*grad_mu;
    // cout << " 2 " << std::endl;
    // cout << accu(Ct) << ' '<< accu(MU) << " "<< accu(Phi) << " "<< accu(grad_mu) <<   std::endl;
    svd(Ut,dt,Vt,Ct);
    // cout << " 1 " << std::endl;
    for(j = 0; j <= rnk; j++)
      Ut.col(j) = Ut.col(j)*dt(j);
    C.rows(cIndexC) = Ut.cols(0,rnk)*Vt.cols(0,rnk).t();

    // Z-step:
    ETA  = ofset + X0*C;
    MU = exp(ETA);
    if(msind == 1) MU.elem(t4).zeros();
    grad_mu = grad_mu_nb(Y, MU, Phi);
    if(msind == 1) grad_mu.elem(t4).zeros();
    C.rows(cIndex) = C.rows(cIndex) -   X1.t()*grad_mu;


    // Phi-step:
    ETA  = ofset + X0*C;
    MU = exp(ETA);
    if(msind == 1) MU.elem(t4).zeros();
    Phi = update_mu_alpha(Y, MU, Phi,naind);
    // Phi = update_mu_phi(Y, MU, Phi);

    if(cObj!=0){
      obj(iter) = nb_dev(Y, MU, Phi, naind);
    } else {
      objtem = nbrrr_likelihood(Y, MU, ETA, Phi,naind);
      obj(iter) = objtem(0);
    }
    diffobj(iter) = abs((obj(iter) - obj(iter- 1)))/(abs(obj(iter)) + 0.1);

    if (diffobj(iter) < epsilon ) {
      converged = true;
      break;
    }
  }

  elp = timer.toc();

  out["C"] = C;
  out["PHI"] = Phi;
  out["sc"] = sc;
  out["sb"] = sb;
  out["eta"] = ETA;
  out["mu"] = exp(ETA);
  out["objval"] = obj;
  out["diffobj"] = diffobj;
  out["converged"] = converged;
  out["ExecTimekpath"] = elp;
  out["maxit"] = iter;
  out["converge"] = 0;
  if(converged) out["converge"] = 1;
  return(out);
}



// [[Rcpp::export]]
double get_sv1(arma::cube xyx, arma::vec ue, int q){
  int iter; arma::vec tem(q);
  for(iter = 0; iter < q; iter++)
    tem(iter) = as_scalar(ue.t()*xyx.slice(iter)*ue);
  return(tem.max()/2.0);
}

// [[Rcpp::export]]
double get_sv2(const arma::mat &xyx, const arma::mat &Y, int q){
  int iter; arma::vec tem(q);
  for(iter = 0; iter < q; iter++)
    tem(iter) = as_scalar(xyx.t()*((Y.col(iter)+1)%xyx.each_col()));
  return(tem.max()/2.0);
}



// [[Rcpp::export]]
double get_sv(const arma::cube &xyx, const arma::vec &ue, int q, arma::uvec tem_uvec){
  // double get_sv( arma::cube xyx,  arma::vec ue, int q, arma::uvec tem_uvec){
  int iter; arma::vec tem(q); arma::mat X2X2;
  // arma::uvec tem_uvec = find(ue);
  for(iter = 0; iter < q; iter++){
    X2X2 = xyx.slice(iter);
    tem(iter) = as_scalar((ue(tem_uvec).t()*(X2X2.submat(tem_uvec,tem_uvec))*ue(tem_uvec)));
  }
  return(tem.max()/2.0);
}
// [[Rcpp::export]]
double softThres(double x, double lambda) {
  return((x > lambda) ? x - lambda :
           (x < -lambda) ? x + lambda : 0.);
}

// [[Rcpp::export]]
arma::vec softT(arma::vec x, arma::vec lambda) {
  arma::vec y; y.zeros(x.n_elem);
  for(int i=0; i < (int) x.n_elem; i++)
    y(i) = softThres(x(i),lambda(i));
  return(y);
}

// [[Rcpp::export]]
int nzcount(arma::vec x) {
  arma::vec y = nonzeros(x) ;
  return y.n_elem;
}






// ctrl <- nbfar_control( objI = 0, initmaxit = 5000, initepsilon = 1e-10)
//   naind <- (!is.na(Y)) + 0 # matrix(1,n,q)
//   msind <- any(naind == 0) + 0
// zerosoltest <- nbzerosol_cpp(Y, X0[,1,drop = FALSE], offset, ctrl, msind, naind)
//   zerosoltestx <- nbZeroSol(Y, X0, c_index = 1, offset, naind)
//   zerosoltest$Z - zerosoltestx$Z
//   zerosoltest$PHI -  zerosoltestx$PHI
// X;Y;Z;O;r; Zini, PhiIni; ndev[ for convergence critteria]
// [[Rcpp::export]]
Rcpp::List nbzerosol_cpp(arma::mat Y, arma::mat X0,
                         arma::mat ofset,
                         Rcpp::List control, int msind, arma::mat naind){
  bool converged=false;
  int pt = X0.n_cols, q = Y.n_cols,  maxit = control["initmaxit"];
  int cObj = control["objI"];
  Rcpp::List out;
  // define epsilon to see convergence of the data
  double epsilon = control["initepsilon"];

  arma::uvec t4=find(naind==0);


  arma::mat MU(size(Y)), ETA(size(Y)), grad_mu;
  arma::mat Res,Ct,X1 = X0; // X1 is the control variable
  double sb = get_sc(X1, Y);X1 = X1/sb;

  arma::mat C = zeros<mat>(pt,q);
  arma::vec Phi = ones<vec>(q);

  wall_clock timer;
  // defining required variable used in the loop
  int iter;
  double elp;
  arma::vec diffobj, obj, objtem(3);
  diffobj.zeros(maxit);obj.zeros(maxit+1);


  ETA = ofset + X0*C;
  MU = exp(ETA);
  if(msind == 1) MU.elem(t4).zeros();


  // nbrrr_likelihood:
  if(cObj!=0){
    obj(0) = nb_dev(Y, MU, Phi, naind);
  } else {
    objtem = nbrrr_likelihood(Y, MU, ETA, Phi,naind);
    obj(0) = objtem(0);
  }

  timer.tic();
  for(iter = 1; iter < maxit; iter++){
    // Z-step:
    if(msind == 1) MU.elem(t4).zeros();
    grad_mu = grad_mu_nb(Y, MU, Phi);
    if(msind == 1) grad_mu.elem(t4).zeros();
    C = C -   X1.t()*grad_mu;

    // Phi-step:
    ETA  = ofset + X0*C;
    MU = exp(ETA);
    if(msind == 1) MU.elem(t4).zeros();
    Phi = update_mu_alpha(Y, MU, Phi,naind);

    if(cObj!=0){
      obj(iter) = nb_dev(Y, MU, Phi, naind);
    } else {
      objtem = nbrrr_likelihood(Y, MU, ETA, Phi,naind);
      obj(iter) = objtem(0);
    }
    diffobj(iter) = abs((obj(iter) - obj(iter- 1)))/(abs(obj(iter)) + 0.1);

    if (diffobj(iter) < epsilon ) {
      converged = true;
      break;
    }
  }
  elp = timer.toc();

  out["Z"] = C;
  out["PHI"] = Phi;
  // out["objval"] = obj;
  // out["diffobj"] = diffobj;
  out["converged"] = converged;
  out["ExecTimekpath"] = elp;
  out["maxit"] = iter;
  out["converge"] = 0;
  if(converged) out["converge"] = 1;
  return(out);
}





// [[Rcpp::export]]
Rcpp::List nbfar_cpp(arma::mat Y, arma::mat Xm,int nlam, arma::vec cindex,
                     arma::mat ofset, Rcpp::List initw,
                     double Dini, arma::mat  Zini, arma::vec PhiIni,
                     arma::mat  Uini, arma::vec Vini,
                     double lmax, Rcpp::List control,
                     int msind, arma::mat naind,
                     int maxit, double epsilon){

  arma::mat X0 = Xm, sd0 = stddev(Xm);
  sd0.elem( find(sd0 == 0) ).ones(); sd0.ones();
  sd0 = 1/sd0;
  X0.each_row() %= sd0;
  arma::vec sd1  = conv_to< arma::vec >::from(sd0);
  arma::vec tem_v;


  // Extract parameter
  // bool converged=false;
  int pt = X0.n_cols, q = Y.n_cols, n = Y.n_rows,ii=0;
  int p = pt - cindex.n_elem, eea = 0;
  int cObj = control["objI"];
  double alp = control["elnetAlpha"];
  double spu = control["spU"], spv = control["spV"];
  double wd =  initw["wd"], gamma0 = control["gamma0"];
  arma::vec wu = initw["wu"], wv = initw["wv"];
  arma::uvec cIndex =  arma::conv_to< uvec >::from(cindex-1);
  arma::uvec cIndexC = mySdiff(linspace<uvec>(0,pt-1,pt), cIndex);
  arma::mat X2 = X0.cols(cIndexC),X1 = X0.cols(cIndex), X3 =X0;
  double sb = get_sc(X1, Y); X1 =  X1/sb;
  arma::mat xx = X2.t()*X2, X2X2, Y_l = lgamma(Y+1); X2X2 = xx/n;
  arma::cube xyx = zeros(p,p,q);
  for(ii=0; ii < q; ii++) xyx.slice(ii) = X2.t()*(X2.each_col()%(Y.col(ii)+1));
  if(eea==1) X2X2 = diagmat(ones<vec>(X2.n_cols));
  // generate sequence of lambda
  double tem = control["lamMaxFac"];
  double lmx = tem*lmax; tem  = control["lamMinFac"];
  double lmn = lmax*tem;
  arma::uvec t4=find(naind==0);
  arma::vec lamSeq =  exp(linspace<vec>(log(lmx),  log(lmn), nlam));


  // ------------ define  path vaariable
  arma::mat uklam = zeros(p,nlam+1),vklam = zeros(q,nlam+1);
  arma::mat philam = zeros(q,nlam+1), BIClam = zeros(4,nlam+1);
  arma::mat objval = zeros(maxit+1,nlam+1);
  arma::vec dklam = zeros(nlam+1),lselectSeq = zeros(nlam+1);
  arma::vec indlam = zeros(nlam+1),execTime = zeros(nlam+1);
  arma::cube zpath = zeros(pt-p,q,nlam+1),Ckpath= zeros(pt,q,nlam+1);
  arma::cube MUkpath = zeros(n,q,nlam+1), ETAkpath = zeros(n,q,nlam+1);

  //  -------- Inittialization, weight for penalty
  //  Update u^0, v^0, d^0,z^0, phi^0 estimate
  Uini = Uini%(1/sd0(cIndexC)); wu = pow(abs(Uini),-1.0*gamma0);
  Zini.each_col() %= (1/sd0(cIndex));
  // Solution for largest lambda [ue = ve = de = 0]
  Rcpp::List zerosolx = nbzerosol_cpp(Y, X0.cols(cIndex),ofset,control,msind,naind);
  arma::mat C = zeros(pt,q), Ct, Z00 = zerosolx["Z"];
  // Z00.each_col() %= (1/sd0(cIndex));
  C.rows(cIndex) = Z00;
  // C.each_col() %= (1/sd0.t());
  arma::vec Phi = zerosolx["PHI"], PHI00 = zerosolx["PHI"];
  arma::vec ue = zeros(p),ve = zeros(q);
  double de = 0;
  C.rows(cIndexC) = de*(ue*ve.t());
  arma::mat ETA = ofset + X0*C;
  arma::mat MU = exp(ETA);
  if(msind == 1) MU.elem(t4).zeros();
  // tem_v = nbrrr_likelihoodx(Y, MU, ETA, Phi, naind, Y_l);
  tem_v = nbrrr_likelihood(Y, MU, ETA, Phi, naind);
  double obj0 = tem_v(0);

  // Model fit intital
  int ik = 0;
  double df = (pt-p)*q, qn = accu(naind);
  double SSE = nb_dev(Y, MU, Phi, naind)/qn;
  BIClam(0,ik) = SSE + (df*log((double)qn))/(qn); //BIC
  BIClam(1,ik) = SSE + 2*df*log((double) pt*q)/(qn); //BICP
  BIClam(2,ik) = SSE + log(log( (double) qn))*df*log((double) pt*q)/(qn); //GIC
  BIClam(3,ik) = SSE + (2/qn)*(df); //AIC
  Ct = C; Ct = Ct.each_col()%sd0.t();
  Ckpath.slice(ik) = Ct; MUkpath.slice(ik) = MU;
  ETAkpath.slice(ik) = ETA;
  zpath.slice(ik) = Ct.rows(cIndex);
  philam.col(ik) = Phi;



  // defining auxilary variable for the loop
  int iter;
  double svk=0,suk=0,m1=1,lam,elp,fac,sv,su,suf = norm(xx,2);
  // double de_temp;
  // arma::mat C_temp,MU_temp;
  // arma::vec Phi_temp, ue_temp, ve_temp;
  arma::vec diffobj, obj, obj2, plfacv, plfacu,PhiI, relerror = zeros(nlam);
  arma::vec vest, uest, maxitc = zeros(nlam),convval = zeros(nlam), xtyv, xtyu;
  arma::mat facW = alp*wd*(wu*wv.t()),facL, d2l = zeros(n,q);
  arma::mat X2Y = X2.t()*Y,X1Y = X1.t()*Y,zb,xuv;
  arma::mat grad_mu, ofset_exp = exp(ofset);
  arma::uvec tem_uvec;
  arma::vec time_prof = zeros<arma::vec>(5);

  // sv = (Y.max()+1)*n/2;
  // su = sqrt(Y.max()+1)*norm(xx,2)/2;
  // converggence: dev/obj
  wall_clock timer;
  for(ii=0; ii < nlam; ii++){
    lam = lamSeq(ii);
    facL = lam*facW;
    fac=(1-alp)*lam;
    diffobj.zeros(maxit);obj.zeros(maxit+1);obj2.zeros(5*(maxit+1));

    // set initialization for the optimization;
    if(norm(C.rows(cIndexC),"fro") == 0){
      ue = Uini; ve = Vini; de = Dini;
      C.rows(cIndex) = Zini;
      Phi = PhiIni;
    }
    C.rows(cIndexC) = de*(ue*ve.t());
    zb =ofset + X0.cols(cIndex)*(C.rows(cIndex));
    xuv = X0.cols(cIndexC)*(C.rows(cIndexC));
    // MU = ofset_exp%exp(X0*C);
    MU = exp(zb + xuv);
    // MU = ofset_exp%exp(X0*C);
    obj(0) = obj0;

    // timer.tic();
    for(iter = 1; iter < maxit; iter++){
      // C_temp = C;
      // ue_temp = ue; ve_temp = ve; de_temp = de;
      // MU_temp = MU;
      // Phi_temp = Phi;
      // if(msind == 1) MU.elem(t4).zeros(); /////////

      timer.tic();
      // update  ve
      tem_uvec = find(ue);
      // sv = get_sv(xyx,ue,q, tem_uvec);
      // grad_mu = grad_mu_nb(Y, MU, Phi);
      grad_mu = grad_mu_nb_uv(Y, MU, Phi, d2l);
      if(msind == 1) grad_mu.elem(t4).zeros();
      // cout << " C " << accu(grad_mu) <<  " "<< accu(MU)<<  " "<< accu(C)<<  " "<< accu(ue)<<  " "<< accu(ve) <<  " "<< accu(de) << std::endl;
      xuv =X2.cols(tem_uvec)*ue(tem_uvec);
      // sv = get_sv2(xuv,Y,q);
      sv = n*d2l.max();
      xtyv = de*ve - (grad_mu.t()*xuv)/sv;
      plfacv =  (facL.t()*abs(ue))/sv;
      vest = softT(xtyv,plfacv)/(1+2*fac*accu(square(ue))/sv);
      svk = norm( vest,2);
      if(svk==0){
        de = 0;ue.zeros(p);ve.zeros(q);
        C.rows(cIndexC) = 0*C.rows(cIndexC);
        // C.rows(cIndex) = Z00;
        // Phi = PHI00;
        break;
      } else {
        de = svk;
        ve = vest/svk;
      }
      time_prof(0)  += timer.toc();


      // update ue
      tem_uvec = find(ve); ETA  = zb;
      ETA.cols(tem_uvec) = ETA.cols(tem_uvec) +  xuv*vest(tem_uvec).t();
      MU = exp(ETA);
      grad_mu = grad_mu_nb_uv(Y, MU, Phi, d2l);

      timer.tic();
      su = suf*sqrt((d2l.cols(tem_uvec)*square(ve(tem_uvec))).max());
      // su = norm(xx + X2.t()*(X2.each_col()%(Y.cols(tem_uvec)*square(ve(tem_uvec)))),2)/2;
      // su = norm(xx + X2.t()*(X2.each_col()%(Y.cols(tem_uvec)*square(ve(tem_uvec)))),2)/2;

      // grad_mu = grad_mu_nb(Y, MU, Phi);
      if(msind == 1) grad_mu.elem(t4).zeros();
      xtyu = de*ue - (X2.t()*(grad_mu.cols(tem_uvec)*ve(tem_uvec)))/su;
      plfacu = (facL*abs(ve))/su;
      uest = softT(xtyu,plfacu)/(1+2*fac/su);
      suk = norm(uest,2);
      if(suk==0){
        de = 0;ue.zeros(p);ve.zeros(q);
        C.rows(cIndexC) = 0*C.rows(cIndexC);
        // C.rows(cIndex) = Z00;
        // Phi = PHI00;
        xuv = X2.cols(tem_uvec)*uest(tem_uvec);
        break;
      } else {
        tem_uvec = find(uest);
        xuv = X2.cols(tem_uvec)*uest(tem_uvec);
        de = norm(xuv,2)/sqrt(n);
        // de = as_scalar(sqrt((uest(tem_uvec).t()*X2X2.submat(tem_uvec,tem_uvec)*uest(tem_uvec))));
        ue = uest/de;
      }
      C.rows(cIndexC) = (ue*de)*ve.t();
      time_prof(1)  += timer.toc();


      timer.tic();
      // Z-step:
      // xuv = X2.cols(tem_uvec)*uest(tem_uvec);
      ETA  = zb;
      tem_uvec = find(ve); xuv = xuv*ve(tem_uvec).t();
      ETA.cols(tem_uvec) = ETA.cols(tem_uvec) +  xuv;
      MU = exp(ETA);
      // if(msind == 1) MU.elem(t4).zeros(); /////// 2
      grad_mu = grad_mu_nb(Y, MU, Phi);
      if(msind == 1) grad_mu.elem(t4).zeros();
      C.rows(cIndex) = C.rows(cIndex) -   X1.t()*grad_mu;
      zb =ofset + X0.cols(cIndex)*C.rows(cIndex);
      time_prof(2)  += timer.toc();

      timer.tic();
      // Phi-step:
      ETA  = zb;  ETA.cols(tem_uvec) = ETA.cols(tem_uvec) +  xuv;
      MU = exp(ETA);
      // if(msind == 1) MU.elem(t4).zeros(); ///////
      Phi = update_mu_alpha(Y, MU, Phi,naind);
      time_prof(3)  += timer.toc();

      timer.tic();
      if(cObj!=0){
        obj(iter) = nb_dev(Y, MU, Phi, naind);
      } else {
        tem_v = nbrrr_likelihood(Y, MU, ETA, Phi, naind);
        obj(iter) = -1*tem_v(0) + alp*as_scalar(de*lam*wd*(wu.t()*abs(ue))*(wv.t()*abs(ve))) +
          (1-alp)*lam*(de*de)*as_scalar((ue.t()*ue)*(ve.t()*ve));
      }
      time_prof(4)  += timer.toc();
      diffobj(iter) = abs((obj(iter) - obj(iter- 1)))/(abs(obj(iter)) + 0.1);
      relerror(ii) = diffobj(iter);
      if (abs(diffobj(iter)) < epsilon ) {
        // converged = true;
        convval(ii) = 1;
        break;
      }
    }
    maxitc(ii) = iter;
    timer.tic();
    elp = timer.toc();

    indlam(ii) = ii;
    uklam.col(ik) = ue%sd1(cIndexC);
    vklam.col(ik) = ve;
    philam.col(ik) = Phi;
    dklam(ik) = de;
    C.rows(cIndexC) = de*(ue*ve.t());
    if(eea==1){
      m1 = norm(uklam.col(ik),2);
      uklam.col(ik) = uklam.col(ik)/m1;
      dklam(ik) = m1*dklam(ik);
    }
    lselectSeq(ik) = lam;
    objval.col(ik) = obj;//
    execTime(ik) = elp;


    Ct = C;
    Ct.each_col() %= sd0.t();
    zpath.slice(ik) = Ct.rows(cIndex);
    Ckpath.slice(ik) = Ct;

    ETA = ofset + X0*C; ///////
    MU = exp(ETA); ///////
    SSE = nb_dev(Y, MU, Phi, naind)/qn;
    df = nzcount(ue) + nzcount(ve)  +(pt-p)*q -1;

    BIClam(0,ik) = SSE + (df*log((double)qn))/(qn); //BIC
    BIClam(1,ik) = SSE + 2*df*log((double) pt*q)/(qn); //BICP
    BIClam(2,ik) = SSE + log(log( (double) qn))*df*log((double) pt*q)/(qn); //GIC
    BIClam(3,ik) = SSE + (2/qn)*(df); //AIC
    MUkpath.slice(ik) = MU;
    ETAkpath.slice(ik) = ETA;
    ik = ik+1;

    if( (iter>1) && ((nzcount(ue) > (p*spu)) || (nzcount(ve) > (q*spv)))  ) {break;}
  }
  ik = ik-1;

  Rcpp::List out;
  out["ukpath"] = uklam.cols(0,ik);
  out["vkpath"] = vklam.cols(0,ik);
  out["dkpath"] = arma::conv_to<arma::vec>::from(dklam.head(ik+1)) ;
  out["phipath"] = philam.cols(0,ik);
  out["zpath"] = arma::conv_to<arma::cube>::from(zpath.head_slices(ik+1));
  out["mukpath"] = arma::conv_to<arma::cube>::from(MUkpath.head_slices(ik+1));
  out["etapath"] = arma::conv_to<arma::cube>::from(ETAkpath.head_slices(ik+1));
  out["ICKpath"] = BIClam.cols(0,ik);
  out["nkpath"] = ik+1;
  out["objkval"] = objval.cols(0,ik);
  out["Ckpath"] = arma::conv_to<arma::cube>::from(Ckpath.head_slices(ik+1));
  out["lamKpath"] = arma::conv_to<arma::vec>::from(lselectSeq.head(ik+1));
  out["ExecTimekpath"] = arma::conv_to<arma::vec>::from(execTime.head(ik+1));
  out["lamseq"] = lamSeq;
  out["maxit"] = maxitc;
  out["converge"] = convval;
  out["converge_error"] = relerror;
  out["time_prof"] = time_prof;
  return(out);
}









// [[Rcpp::export]]
Rcpp::List cv_nbfar_cpp(arma::mat Y, arma::mat Xm,int nlam, arma::vec cindex,
                        arma::mat ofset, Rcpp::List initw,
                        arma::mat  Zini, arma::vec PhiIni,
                        Rcpp::List xx,
                        double lmax, Rcpp::List control,
                        int maxit, double epsilon,
                        int nfold){
  arma::mat Ytr(size(Y)), Yte(size(Y)), natr(size(Y)), nate(size(Y));
  arma::uvec indna = find_finite(Y),teIndex, trIndex;
  arma::vec ID = linspace(0, indna.n_elem-1, indna.n_elem);
  ID = ID - nfold*floor(ID/nfold); ID = ID(randperm(indna.n_elem));
  int misind,insel;
  Rcpp::List fitF;
  double Dini = xx["D"],lam;
  arma::mat  Uini = xx["U"]; arma::vec Vini  = xx["V"];
  arma::vec phiest, tttval(3), tec(nfold);
  arma::mat mutest, etatest;
  arma::mat dev(nfold, nlam), sdcal(nfold, nlam);
  dev.fill(datum::nan); sdcal.fill(datum::nan);

  for(int ifold=0; ifold < nfold; ifold++){
    // cout << "Fold " << ifold << std::endl;
    teIndex = indna(find(ID == ifold));
    trIndex = indna(find(ID != ifold));
    Ytr.zeros(); Yte.zeros();natr.zeros(); nate.zeros();
    Ytr.elem(trIndex) = Y.elem(trIndex);
    Yte.elem(teIndex) = Y.elem(teIndex);
    natr.elem(trIndex) = ones<vec>(trIndex.n_elem);
    nate.elem(teIndex) = ones<vec>(teIndex.n_elem);
    misind = 1;
    fitF =  nbfar_cpp(Ytr, Xm,  nlam, cindex,
                      ofset, initw , Dini, Zini, PhiIni,
                      Uini, Vini, lmax , control, misind,
                      natr,
                      maxit, epsilon);
    int nfit = fitF["nkpath"]; arma::vec  lamkpath = fitF["lamKpath"];
    arma::vec lamseq = fitF["lamseq"]; arma::cube mukpath = fitF["mukpath"];
    arma::cube etapath = fitF["etapath"]; arma::mat phipath = fitF["phipath"];
    for (int im = 0; im < nfit; im++) {
      lam  = lamkpath(im);
      insel =  conv_to< int >::from(find(lamseq == lam));
      mutest = mukpath.slice(im); etatest = etapath.slice(im);
      mutest.elem(trIndex) =  ones<vec>(trIndex.n_elem);
      etatest.elem(trIndex) =  zeros<vec>(trIndex.n_elem);
      phiest = phipath.col(im);
      tttval = nbrrr_likelihood(Yte, mutest, etatest, phiest, nate);
      dev(ifold, insel) = tttval(0);
      sdcal(ifold,insel) = tttval(1);
      tec(ifold) = tttval(2);
    }
  }
  Rcpp::List out;
  out["lamseq"] = fitF["lamseq"];
  out["dev"] = dev;
  out["sdcal"] = sdcal;
  out["tec"] = tec;
  // out["lamseq"] = arma::conv_to<arma::vec>::from(lamseq);
  return out;
}










struct ParCv : public RcppParallel::Worker
{
  const arma::vec ID; const arma::uvec indna;
  arma::vec lamseq; arma::mat dev; arma::mat sdcal;
  arma::vec tec;
  const arma::mat Y; const arma::mat Xm;
  int nlam; const arma::vec  cindex;
  const arma::mat ofset; const Rcpp::List initw;
  double Dini; const arma::mat Zini; const arma::vec PhiIni;
  const arma::mat Uini; const arma::vec Vini;
  double lmax; const Rcpp::List control;
  int maxit; double epsilon;

  // // destination matrix
  // arma::mat dev, sdcal;
  // arma::vec tec, lamseq;
  //
  // // source matrix
  // const arma::uvec indna;
  // double Dini, lmax,epsilon;
  // const arma::mat  Uini, Y, Xm,ofset, Zini;
  // const arma::vec Vini, cindex, PhiIni, ID;
  // const Rcpp::List control, zerosol, initw;
  // int maxit,nlam;

  // initialize with source and destination
  ParCv(const arma::vec ID, const arma::uvec indna,
        arma::vec lamseq, arma::mat dev,arma::mat sdcal,
        arma::vec tec,
        const arma::mat Y, const arma::mat Xm,
        int nlam, const arma::vec  cindex,
        const arma::mat ofset, const Rcpp::List initw,
        double Dini, const arma::mat Zini, const arma::vec PhiIni,
        const arma::mat Uini, const arma::vec Vini,
        double lmax , const Rcpp::List control,
        int maxit, double epsilon)
    : ID(ID),indna(indna),lamseq(lamseq),dev(dev),sdcal(sdcal),tec(tec),
      Y(Y), Xm(Xm),  nlam(nlam), cindex(cindex),
      ofset(ofset), initw(initw) , Dini(Dini), Zini(Zini),
      PhiIni(PhiIni), Uini(Uini), Vini(Vini), lmax(lmax) , control(control),
      maxit(maxit), epsilon(epsilon) {}


  // take the square root of the range of elements requested
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t k = begin; k < end; k++) {
      int ifold = k;
      arma::uvec teIndex, trIndex;
      arma::mat Ytr(size(Y)), Yte(size(Y)), natr(size(Y)), nate(size(Y));
      int misind,insel;
      Rcpp::List fitF;
      double lam;
      arma::vec phiest, tttval(3);
      arma::mat mutest, etatest;

      teIndex = indna(find(ID == ifold));
      trIndex = indna(find(ID != ifold));
      Ytr.zeros(); Yte.zeros();natr.zeros(); nate.zeros();
      Ytr.elem(trIndex) = Y.elem(trIndex);
      Yte.elem(teIndex) = Y.elem(teIndex);
      natr.elem(trIndex) = ones<vec>(trIndex.n_elem);
      nate.elem(teIndex) = ones<vec>(teIndex.n_elem);
      misind = 1;
      fitF =  nbfar_cpp(Ytr, Xm,  nlam, cindex,
                        ofset, initw , Dini, Zini, PhiIni,
                        Uini, Vini, lmax , control, misind,
                        natr,
                        maxit, epsilon);
      int nfit = fitF["nkpath"]; arma::vec  lamkpath = fitF["lamKpath"];
      arma::vec lamseqx = fitF["lamseq"]; arma::cube mukpath = fitF["mukpath"];
      arma::cube etapath = fitF["etapath"]; arma::mat phipath = fitF["phipath"];
      for (int im = 0; im < nfit; im++) {
        lam  = lamkpath(im);
        insel =  conv_to< int >::from(find(lamseqx == lam));
        mutest = mukpath.slice(im); etatest = etapath.slice(im);
        mutest.elem(trIndex) =  ones<vec>(trIndex.n_elem);
        etatest.elem(trIndex) =  zeros<vec>(trIndex.n_elem);
        phiest = phipath.col(im);
        tttval = nbrrr_likelihood(Yte, mutest, etatest, phiest, nate);
        dev(ifold, insel) = tttval(0);
        sdcal(ifold,insel) = tttval(1);
        tec(ifold) = tttval(2);
      }
      lamseq = lamseqx;
    }
  }
};



// [[Rcpp::export]]
Rcpp::List cv_nbfar_par(arma::mat Y, arma::mat Xm,int nlam, arma::vec cindex,
                        arma::mat ofset, Rcpp::List initw,
                        arma::mat  Zini, arma::vec PhiIni,
                        Rcpp::List xx,
                        double lmax, Rcpp::List control,
                        int maxit, double epsilon,
                        int nfold){
  arma::uvec indna = find_finite(Y);
  arma::vec ID = linspace(0, indna.n_elem-1, indna.n_elem);
  ID = ID - nfold*floor(ID/nfold); ID = ID(randperm(indna.n_elem));
  // input
  double Dini = xx["D"];
  arma::mat  Uini = xx["U"]; arma::vec Vini  = xx["V"];
  // Output
  arma::mat dev(nfold, nlam), sdcal(nfold, nlam);
  dev.fill(datum::nan); sdcal.fill(datum::nan);
  arma::vec tec(nfold), lamseq(nlam);
  // Constructer of the class
  ParCv cv_nbfar(ID,indna,lamseq,dev,sdcal,tec, Y, Xm,  nlam, cindex,
                 ofset, initw , Dini, Zini, PhiIni, Uini, Vini, lmax , control,
                 maxit, epsilon);
  parallelFor(0, nfold, cv_nbfar);
  Rcpp::List out;
  out["lamseq"] = cv_nbfar.lamseq;
  out["dev"] = cv_nbfar.dev;
  out["sdcal"] = cv_nbfar.sdcal;
  out["tec"] = cv_nbfar.tec;
  return out;
}




