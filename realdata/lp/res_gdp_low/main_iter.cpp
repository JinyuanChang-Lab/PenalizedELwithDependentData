#include <RcppArmadillo.h>
#include <math.h>
#include <iomanip>
#include <float.h>
#include <iostream>
#include <algorithm>
#include <random>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec log_star(arma::vec X, double thresh){
  int m = X.n_elem;
  arma::vec logX(m);
  for(int i=0; i<m; i++){
    if (X(i) < thresh) {
      logX(i) = log(thresh) - 1.5 + 2.0 * X(i)/thresh - 0.5 * pow((X(i)/thresh), 2);
    }
    else logX(i) = log(X(i));
  }
  return logX;
}

// [[Rcpp::export]]
arma::vec der_log_star(arma::vec X, double thresh){
  int m = X.n_elem;
  arma::vec derlogX(m);
  for(int i=0; i<m; i++){
    if (X(i) < thresh) {
      derlogX(i) = 2.0/thresh - X(i)/pow(thresh,2);
    }
    else derlogX(i) = 1.0/X(i);
  }
  return derlogX;
}

// [[Rcpp::export]]
arma::mat scad_penalty(arma::mat beta, double lam){
  double a = 3.7;
  int s = beta.n_rows;
  int t = beta.n_cols; 
  arma::mat scadval(s,t);
  for(int i=0; i<s; i++){
    for(int j=0; j<t; j++){
      if (abs(beta(i,j) ) <= lam) scadval(i,j) = lam * abs(beta(i,j));
      else if ((abs(beta(i,j) ) > lam) & (abs(beta(i,j) ) <= a * lam)) 
        scadval(i,j) = (2*a*lam*abs(beta(i,j)) - pow(abs(beta(i,j)),2) - pow(lam,2)) / (2*(a-1));
      else if (abs(beta(i,j) ) > a * lam) scadval(i,j) = pow(lam,2)*(a + 1)/2;
    }
  }
  return scadval;
}

// [[Rcpp::export]]
arma::mat scad_grad(arma::mat beta, double lam){
  double a = 3.7;
  int s = beta.n_rows;
  int t = beta.n_cols; 
  arma::mat scadgrad(s,t);
  for(int i=0; i<s; i++){
    for(int j=0; j<t; j++){
      if (abs(beta(i,j) ) <= lam) scadgrad(i,j) = lam * (abs(beta(i,j)) - beta(i,j) > 0? -1:1);
      else if ((abs(beta(i,j) ) > lam) & (abs(beta(i,j) ) <= a * lam)) 
        scadgrad(i,j) = max(a*lam - abs(beta(i,j)), 0.0) / (a - 1)*(abs(beta(i,j)) - beta(i,j) > 0? -1:1);
    }
  }
  return scadgrad;
}

// [[Rcpp::export]]
Rcpp::List scad_quad(arma::mat beta, double lam){
  int s = beta.n_rows;
  int t = beta.n_cols;
  
  double thresh = lam / 4;
  arma::mat init_beta = beta;
  for(int i=0; i<s; i++){
    for(int j=0; j<t; j++){
      if(abs(init_beta(i,j)) < 1e-4) init_beta(i,j) = thresh;
      else if(abs(init_beta(i,j)) < thresh) init_beta(i,j) = sign(init_beta(i,j)) * thresh;
    }
  }
  
  arma::mat temp_1 = scad_penalty(init_beta, lam);
  arma::mat temp_2 = scad_grad(abs(init_beta), lam);
  arma::mat temp_3 = (beta % beta - init_beta % init_beta) / 2;
  
  arma::mat scad_val_mat = temp_1 + temp_2 % temp_3 / abs(init_beta);
  double scad_val = accu(scad_val_mat);
  
  arma::mat scad_der = temp_2 % beta / abs(init_beta);
  
  return List::create(Named("scad_val") = scad_val, Named("scad_der") = scad_der);
}

// [[Rcpp::export]]
double ee_lambda(mat lambda, mat g_ee, mat theta, double nu, double tau){
  int n = g_ee.n_rows;
  
  vec lambda_g = g_ee * lambda + 1; 
  vec log_g = log_star(lambda_g, 1/n);
  double obj_1 = sum(log_g);
  
  double scad_lambda = scad_quad(lambda, nu)["scad_val"]; 
  double obj_2 = n * scad_lambda;
  
  double scad_theta = scad_quad(theta, tau)["scad_val"];
  double obj_3 = n * scad_theta; 
  
  double obj = obj_1 - obj_2 + obj_3;
  
  return obj; 
}

// [[Rcpp::export]]
double ee_fun(mat lambda, mat g_ee){
  int n = g_ee.n_rows;
  
  vec lambda_g = g_ee * lambda + 1; 
  vec log_g = log_star(lambda_g, 1/n);
  double obj_1 = sum(log_g);
  
  return obj_1;
}

// [[Rcpp::export]]
Rcpp::List grad_theta(mat theta, uvec supp_theta, mat lambda, mat g_ee, cube gradG, double tau){
  int n = gradG.n_rows;
  double alpha = 0.001;
  int p  = theta.n_rows;
  
  vec lambda_g = 1.0 + g_ee * lambda;
 
  uvec prox_set(p, fill::zeros);
  int count = 0;
  
  double len_theta = supp_theta.n_elem;
  vec grad_vect(p, fill::zeros);
  for(int i=0; i<len_theta; i++){
    int index_theta = supp_theta(i);
    vec lambda_partial_g = gradG.slice(index_theta) * lambda;
    // double scad_der = scad_quad(theta.row(index_theta), tau)["scad_der"];
    double scad_der = scad_grad(theta.row(index_theta), tau)(0,0);
    double temp = sum(der_log_star(lambda_g, 1.0/n) % lambda_partial_g);
    
    if(abs(theta.row(index_theta)(0) - alpha * temp) <= n * tau * alpha){
      prox_set(count) = index_theta;
      count = count + 1;
    }
    
    grad_vect(index_theta) = temp + n * scad_der;
  }
  
  uvec proxset = prox_set.head(count);
  
  //return grad_vect;
  return List::create(Named("grad_vect") = grad_vect, Named("prox_set") = proxset);
}

arma::vec logstar_arma(arma::vec X, double m, double thresh){
  arma::vec logX(m);
  for(int i=0; i<m; i++){
    if (X(i) < thresh) {
      logX(i) = log(thresh) - 1.5 + 2.0 * X(i)/thresh - 0.5 * pow((X(i)/thresh), 2);
    }
    else logX(i) = log(X(i));
  }
  return logX;
}

arma::vec der2logstar_arma(arma::vec X, double m, double thresh){
  arma::vec der2logX(m);
  for(int i=0; i<m; i++){
    if (X(i) < thresh) {
      der2logX(i) = -1.0/pow(thresh,2);
    }
    else der2logX(i) = -1.0/pow(X(i),2);
  }
  return der2logX;
}

arma::vec derlogstar_arma(arma::vec X, double m, double thresh){
  arma::vec derlogX(m);
  for(int i=0; i<m; i++){
    if (X(i) < thresh) { 
      derlogX(i) = 2.0/thresh - X(i)/pow(thresh,2);
    }
    else derlogX(i) = 1.0/X(i);
  }
  return derlogX;
}

arma::vec Hessian_x(arma::vec dxu, arma::mat Gee, arma::mat Geet,
                    arma::vec d0, arma::vec d1, arma::vec d2){
  // double n = Gee.n_rows;
  int r = (dxu.n_elem)/2;
  arma::vec x1 = dxu.head(r);
  arma::vec x2 = dxu.tail(r);
  arma::vec Hes_x(dxu.n_elem);
  Hes_x.head(r) = Geet*(d0%((Gee*x1)))+d1%x1+d2%x2;
  Hes_x.tail(r) = d2%x1+d1%x2;
  return Hes_x;
}

arma::vec Pinv_x(arma::vec dxu,arma::vec p1,arma::vec p2,
                 arma::vec p3){
  int r = (dxu.n_elem)/2;
  arma::vec x1 = dxu.head(r);
  arma::vec x2 = dxu.tail(r);
  arma::vec Pr(dxu.n_elem);
  Pr.head(r) = p1%x1-p2%x2;
  Pr.tail(r) = -p2%x1+p3%x2;
  return Pr;
}

Rcpp::List pcg_(arma::vec dxu, arma::vec b, double tol, int maxit,
                arma::vec d0, arma::vec d1, arma::vec d2,
                arma::vec p1, arma::vec p2, arma::vec p3,
                arma::mat A, arma::mat At) {
  int n = dxu.n_elem;
  
  int iter,imin=0, flag=1;
  double normr,normr_act, normrmin, rho=1.0;
  double tolb, n2b = arma::norm(b,2);
  
  double rho1,beta,pq,alpha,relres;
  arma::vec x=dxu;
  
  arma::vec xmin=x, r, z, p, q;
  
  if(n2b==0){
    iter=0;
    flag=0;
    Rcout<<"n2b=0"<<endl;
    return List::create(Named("dlambda") = x, Named("iter") = iter, Named("flag")=flag);
  }
  
  tolb = tol*n2b;
  
  r = b - Hessian_x(x,A,At,d0,d1,d2);
  normr = arma::norm(r,2);
  normr_act = normr;
  if(normr<=tolb){
    flag=0;
    relres = normr/n2b;
    iter = 0;
    return List::create(Named("dlambda") = x, Named("iter") = iter, Named("flag")=flag);
  }
  
  
  normrmin = normr;
  int stag=0;
  int moresteps = 0;
  int maxstagsteps = 3;
  int maxmsteps = min(min(floor(n/50),5.0),double(n-maxit));
  int ii;
  bool verbose=false;
  if(verbose){
    Rcpp::Rcout<<n<<endl;
    Rcpp::Rcout<<n2b<<endl;
    Rcpp::Rcout<<x<<endl;
    Rcpp::Rcout<<r<<endl;
    Rcpp::Rcout<<tolb<<endl;
    Rcpp::Rcout<<maxmsteps<<endl;
    Rcpp::Rcout<<"true"<<endl;
  }
  
  for(int ii=1;ii<=maxit;ii++){
    z = Pinv_x(r, p1,p2,p3);
    if(!z.is_finite()){
      flag=2;
      break;
    }
    rho1 = rho;
    rho = arma::dot(r,z);
    if ((rho == 0) || std::isinf(rho))
    {flag = 4;
      break;}
    
    if(ii==1)
      p=z;
    else{
      beta = rho/rho1;
      if ((beta == 0) || std::isinf(beta))
      {flag = 4;
        break;}
      p = z+beta*p;
    }
    
    q = Hessian_x(p,A,At,d0,d1,d2);
    pq = arma::dot(p,q);
    
    
    
    if(pq<=0 || std::isinf(pq)){
      flag =4;
      break;
    }
    else
      alpha = rho/pq;
    
    if(std::isinf(alpha)){
      flag=4;
      break;
    }
    
    if((arma::norm(p,2)*abs(alpha))< ((2e-16)*arma::norm(x,2)))stag++;
    else stag = 0;
    
    x = x+alpha*p;
    r = r-alpha*q;
    normr = arma::norm(r,2);
    normr_act = normr;
    //if(cout)Rcout<<normr<<"  "<<tolb<<"  "<<stag<<" "<<maxstagsteps<<" "<<moresteps<<endl;
    
    if(normr<=tolb || stag >=maxstagsteps || moresteps){
      r = b-Hessian_x(x,A,At,d0,d1,d2);
      normr_act = arma::norm(r,2);
      if(normr_act<=tolb){
        flag = 0;
        iter = ii;
        break;
      }
      else{
        if(stag>=maxstagsteps &&moresteps==0)stag=0;
        moresteps++;
        if(moresteps>=maxmsteps){
          flag = 3;
          iter=ii;
          break;
        }
        //  Rcpp::Rcout<<"pcg:tooSmallTolerance";
        
      }
    }
    if (normr_act<normrmin){
      normrmin = normr_act;
      xmin = x;
      imin=ii;
    }
    if(stag>=maxstagsteps){flag=3;break;}
    
  }
  
  if(flag==0)relres = normr_act/n2b;
  else{
    arma::vec r_comp = b-Hessian_x(xmin,A,At,d0,d1,d2);
    if(arma::norm(r_comp,2)<=normr_act){
      x = xmin;
      iter = imin;
      relres = arma::norm(r_comp,2)/n2b;
    }
    else{
      iter = ii;
      relres = normr_act/n2b;
    }
  }
  return List::create(Named("dlambda") = x, Named("iter") = iter, Named("flag")=flag);
}

// [[Rcpp::export]]
arma::mat optim_lambda(arma::mat Gee, double nu,
                      double tar_gap = 1e-9, double eta = 1e-3, double pcgmaxi=5000,
                      bool verbose = false, string preconmat="identity"){
  
  // IPM PARAMETERS
  double MU = 2;         // updating parameter of t
  int MAX_NT_ITER= 400;     // maximum IPM (Newton) iteration
  
  // LINE SEARCH PARAMETERS
  double ALPHA           = 0.01;     // minimum fraction of decrease in the objective
  double BETA            = 0.5;      // stepsize decrease factor
  int MAX_LS_ITER     = 100;      // maximum backtracking line search iteration
  string status;
  
  
  
  int pitr = 0;  //pcg iter
  int pflag = 0 ; //pcg convergence
  
  int r = Gee.n_cols;
  double n = Gee.n_rows;
  
  arma::mat H,Geet = Gee.t();
  
  double t = min(max(1.0,1.0/nu),2*n/1e-3);
  double reltol = tar_gap;
  double thresh = 1.0/n;
  
  double pobj=1e100, dobj=-pobj, step = 1000;
  arma::vec lambda(r, arma::fill::ones);
  lambda = 0.1 * lambda;
  arma::vec u(r, arma::fill::ones);
  arma::vec f(r*2);
  f.head(r) = lambda-u;
  f.tail(r)= -lambda-u;
  int ntiter; //number of iteration
  int lsiter; //
  arma::vec d0;
  arma::vec z= Gee*lambda + 1.0;
  arma::vec dlambda_u(2*r, arma::fill::zeros);
  arma::vec diagxtx = 2.0*MU*arma::vec(r,arma::fill::ones);
  d0 = (-1.0/n)* der2logstar_arma(z, n,thresh);
  if (preconmat=="diag")
    arma::vec diagxtx = 2*arma::diagvec(Geet*diagmat(d0)*Gee);
  
  if(verbose==true){
    Rcpp::Rcout<<"Solving a proplem of size" <<n<<r<<", with nu="<<nu<<endl;
    Rcpp::Rcout<<"-------------------------------------------"<<endl;
    Rcpp::Rcout<<"iter  step len  pcgiters    gap          primobj       dualobj"<<endl;
    
  }
  
  //main loop
  double gap, normg, phi,newphi,gdx,pcgtol;
  arma::vec newlambda,newf(2*r),newz,newu,dlambda,du;
  
  arma::vec q1;
  arma::vec q2;
  arma::vec d1;
  arma::vec d2;
  arma::vec gradphi;
  arma::vec prb,prs;
  Rcpp::List result;
  arma::mat lower;
  arma::wall_clock timer;
  // double time=0.0,time1=0.0;
  // 
  int bb=1;
  
  for (ntiter=0; ntiter<=MAX_NT_ITER; ntiter++){
    
    //arma::vec z= Gee*lambda + 1.0;
    z= Gee*lambda + 1.0;
    
    double s = min(nu * n / arma::norm(Geet * derlogstar_arma(z, n,thresh),"inf"), 1.0);
    arma::vec vv = s/n*derlogstar_arma(z,n,thresh);
    //Rcout<<lambda<<endl<<u<<endl;
    
    //double priobj = -1.0 * (1.0 / n) * logstar(lambda_g, n).sum() + nu * lambda.array().abs().sum();
    //double duobj = (n * vv).array().log().sum() / n + 1.0 - vv.sum();
    
    
    //double tmppobj = pobj;
    pobj = -1.0/n * arma::sum(logstar_arma(z, n, thresh))  + nu * arma::norm(lambda,1);
    dobj = max((arma::sum(log(n *vv)) / n + 1.0 - arma::sum(vv)), dobj);
    
    
    gap = pobj-dobj;
    
    if(verbose==true){
      Rcpp::Rcout<<" "<<ntiter<<"  ";
      Rcpp::Rcout<<step<<"     ";
      Rcpp::Rcout<<pitr<<"    ";
      Rcpp::Rcout<<setiosflags(ios::scientific)<<gap<<"  ";
      Rcpp::Rcout<<setiosflags(ios::scientific)<<pobj<<"  ";
      Rcpp::Rcout<<setiosflags(ios::scientific)<<dobj<<" ";
      Rcpp::Rcout<<" "<<MU<<endl;
    }
    if (gap/(-dobj)<reltol){
      status = "Solved";
      if(verbose==true)
        Rcpp::Rcout<<"Absolute tolerance reached"<<endl;
      // Rcout<<"number of seconds: "<<time<<endl;
      return lambda;
      //return Rcpp::List::create(Named("lambda")=lambda, Named("status")=status);
    }
    
    // if(ntiter==10){
    //   status = "Solved";
    //   //Rcpp::Rcout<<"Absolute tolerance reached"<<endl;
    //   return Rcpp::List::create(Named("lambda")=lambda, Named("status")=status);
    // }
    
    if(step >= 0.5)
      t = max(8.0*min(2.0*n/gap,t),t);
    else if (step < 1e-5) {
      t = 1.5*t;
    }
    //Rcout<<t<<endl;
    
    q1 = 1.0/(u+lambda);
    q2 = 1.0/(u-lambda);
    
    // Rcout<<t<<endl;
    timer.tic();
    
    d0 = (-1.0/n)* der2logstar_arma(z, n,thresh);
    d1 = (q1%q1 + q2%q2)/t;
    d2 = (q1%q1 - q2%q2)/t;
    
    //Rcout<<d1<<endl<<d2;
    
    gradphi =  arma::join_vert((-1.0/n)*Geet*derlogstar_arma(z, n,thresh)-(q1-q2)/t,
                               nu*arma::vec(r,arma::fill::ones)-(q1+q2)/t);
    
    prb = diagxtx+d1;
    prs = prb%d1-(d2%d2);
    
    normg = arma::norm(gradphi,2);
    pcgtol = min(1e-1,eta*gap/min(1.0,normg));
    if (ntiter !=0 && pitr ==0)
      pcgtol = pcgtol*0.1;
    
    arma::vec tmp = dlambda_u;
    result = pcg_(dlambda_u,-gradphi,pcgtol,pcgmaxi,d0,d1,d2,d1/prs, d2/prs, prb/prs,Gee,Geet);
    
    arma::vec dlambda_u = result["dlambda"];
    if(dlambda_u.has_nan())Rcout<<"555"<<endl;
    //if(dlambda_u.has_nan()||(d1/prs).has_nan())
    //Rcout<<dlambda_u.has_nan()<<" "<<(d1/prs).has_nan()<<endl;
    //Rcout<<dlambda_u.has_nan()<<endl;
    pflag = result["flag"];
    pitr = result["iter"];
    //    if(pflag==2||pflag==4){
    //      gradphi = -1.0*gradphi;
    //      arma::mat H = Hessian(Gee, Geet, d0, d1, d2);
    //      //lower = arma::chol(H,"lower");
    //      //dlambda_u = solve(trimatl(lower),gradphi);
    //      //dlambda_u = solve(trimatu(lower.t()),dlambda_u);
    //      
    //      //if(bb==1)Rcout<<H<<endl;
    //      //Rcout<<H.is_sympd()<<endl;
    //      //dlambda_u = solve(H,gradphi,arma::solve_opts::likely_sympd);
    //      dlambda_u = arma::inv_sympd(H)*gradphi;
    //      //dlambda_u = arma::pinv(H)*gradphi;
    //      //dlambda_u = dlambda_u;
    //    }
    // if(dlambda_u.has_nan()){
    //   Rcout<<pflag<<endl;
    //   Rcout<<H.is_symmetric()<<endl;
    //   MU++;
    //   Rcout<<"111"<<endl;
    //   goto init;
    // }
    // while(dlambda_u.has_nan()){
    //   diagxtx = diagxtx+1.0;
    //   prb = diagxtx+d1;
    //   prs = prb%d1-(d2%d2);
    //   
    //   result = pcg_(tmp,-gradphi,pcgtol,pcgmaxi,d0,d1,d2,d1/prs, d2/prs, prb/prs,Gee,Geet);
    //   arma::vec dlambda_u = result["dlambda"];
    //   Rcout<<dlambda_u.has_nan()<<"  "<<gradphi.has_nan()<<endl;
    //   pflag = result["flag"];
    //   pitr = result["iter"];
    //   
    // }
    // Rcout<<"while down"<<endl;
    
    // //Rcout<<dlambda_u<<endl;
    // arma::vec tmp = dlambda_u;
    
    
    
    
    // arma::mat H = Hessian(Gee, Geet, d0, d1, d2);
    // dlambda_u = arma::inv_sympd(H)*gradphi;
    // dlambda_u = -1.0*dlambda_u;
    //Rcout<<norm((dlambda_u-tmp),"inf")<<endl;
    // time = time+timer.toc();
    
    
    
    //Rcout<<H<<endl<<gradphi;
    //Rcout<<t<<endl;
    
    //Rcout<<gradphi<<endl<<dlambda_u<<endl;
    //Rcout<<z<<endl;
    
    //Rcout<<Hessian(Gee, Geet, d0, d1, d2)<<endl;
    
    
    if(pflag==1)pitr=pcgmaxi;
    
    dlambda = dlambda_u.head(r);
    du = dlambda_u.tail(r);
    
    
    phi = -1.0/n * arma::sum(logstar_arma(z, n,thresh)) + nu * arma::sum(u) - arma::sum(log(-f))/t;
    step = 1.00;
    gdx = arma::dot(gradphi,dlambda_u);
    for (lsiter = 1;lsiter<=MAX_LS_ITER;lsiter++){
      newlambda = lambda+step*dlambda; newu=u+step*du;
      newf.head(r) = newlambda-newu;
      //Rcout<<1<<endl;
      newf.tail(r) = -1.0*newlambda-newu;
      //newf = arma::join_vert(newlambda-newu,-newlambda-newu);
      if (arma::max(newf)<0){
        newz = Gee*newlambda + 1.0;
        
        newphi = -1.0/n * arma::sum(logstar_arma(newz, n,thresh))+nu*arma::sum(newu)-arma::sum(log(-newf))/t;
        if((newphi-phi)<=(ALPHA*step*gdx))
          break;
      }
      step=BETA*step;
    }
    
    if(lsiter == MAX_LS_ITER)break;
    
    lambda = newlambda; u = newu; f=newf;
  }
  
  
  
  if(lsiter == MAX_LS_ITER){
    Rcpp::Rcout<<"MAX_LS_ITER exceeded in BLS";
    status = "Failed";}
  else if(ntiter == MAX_NT_ITER)
    status = "Failed";
  
  arma::mat lambda_matrix = arma::mat(lambda.memptr(), lambda.n_elem, 1, false, true);
  // return Rcpp::List::create(Named("lambda")=lambda, Named("status")=status);
  return lambda_matrix;
}

arma::mat bias(arma::mat theta, arma::mat lambda, mat XX, arma::uvec supp_theta, arma::uvec supp_lambda,
               arma::mat g_ee, Rcpp::Function grad_g){
  int n = XX.n_rows;
  int r = lambda.n_rows;
  
  cube gradG = as<arma::cube>(grad_g(theta, XX));
  
  mat Temp = mean(gradG, 0); // dim = 1?
  mat G = Temp.submat(supp_lambda, supp_theta);
  
  mat Temp2 = g_ee.t() * g_ee / n;
  mat V = Temp2.submat(supp_lambda, supp_lambda);
  
  mat Down = repmat(g_ee * lambda, 1, r) + 1;
  mat Eta = mean(g_ee / Down, 0).t();
  
  if(det(V) == 0){V = V + 0.001 * eye(size(V));}
  
  mat M = G.t() * inv(V) * G;
  
  if(det(M) == 0){M = M + 0.001 * eye(size(M));}
  
  mat bias = inv(M) * G.t() * inv(V) * Eta.rows(supp_lambda);
  
  return(bias);
}

// [[Rcpp::export]]
Rcpp::List main_iter(Rcpp::Function auxi_fun, mat XX, mat init_theta, Rcpp::Function grad_g,
                     double tau, double nu, double eps_tol = 0.005, int iter_num = 1000){
  mat lambda, theta, g_ee, min_theta;
  
  int n = XX.n_rows;
  int p = init_theta.n_rows;
  g_ee = as<arma::mat>(auxi_fun(init_theta, XX));
  int r = g_ee.n_cols;
  
  uvec supp_theta = linspace<uvec>(0,p-1,p);
  
  int iternum = iter_num + 1;
  mat theta_iter(iternum, p);
  theta_iter.row(0) = init_theta.col(0).t();
  
  vec obj_vec(iter_num);
  
  lambda = optim_lambda(g_ee, nu);
  theta = init_theta;
  
  cube gradG;
  double beta1 = 0.9, beta2 = 0.999, epsilon = 1e-6, alpha = 0.001;
  vec m(p, fill::zeros);
  vec v(p, fill::zeros);
  vec m_hat, v_hat;
  
  for(int count_num = 0; count_num < iter_num; count_num++){
    gradG = as<arma::cube>(grad_g(theta, XX));
    Rcpp::List res = grad_theta(theta, supp_theta, lambda, g_ee, gradG, tau);
    vec driv = res["grad_vect"];
    uvec prox_set = res["prox_set"];
    
    theta.rows(prox_set).fill(0);
    supp_theta = find(vectorise(theta) != 0);
    uvec zero_theta = find(vectorise(theta) == 0);
    
    m = beta1 * m + (1.0 - beta1) * driv;
    v = beta2 * v + (1.0 - beta2) * driv % driv;
    m_hat = m / (1.0 - pow(beta1, count_num + 1));
    v_hat = v / (1.0 - pow(beta2, count_num + 1));
    
    theta.col(0) = theta.col(0) - (alpha / (sqrt(v_hat) + epsilon)) % m_hat;
    theta.rows(zero_theta).fill(0);
    
    g_ee = as<arma::mat>(auxi_fun(theta, XX));
    lambda = optim_lambda(g_ee, nu);
    
    theta_iter.row(count_num + 1) = theta.col(0).t();
    obj_vec(count_num) = ee_lambda(lambda, g_ee, theta, tau, nu);
    
    if(max(abs(theta_iter.row(count_num + 1) - theta_iter.row(count_num))) < 1e-6 || count_num == iter_num - 1){
      //if(count_num == iter_num - 1) cout << "Max iterion reached!" << endl; 
      break;
      }

    //Rcout<<count_num<<endl;
  }
  
  uvec supp_lambda = find(abs(vectorise(lambda)) > 0.001);
  mat theta_debiased;
  
  if(supp_lambda.n_elem == 0 || supp_theta.n_elem == 0){
    theta_debiased = theta;
  }
  else{
    theta_debiased = theta;
    mat bias_mat = bias(theta, lambda, XX, supp_theta, supp_lambda, g_ee, grad_g);
    theta_debiased.rows(supp_theta) = theta.rows(supp_theta) - bias_mat;
  }
  
  return List::create(Named("theta") = theta, Named("theta.iter") = theta_iter, 
                      Named("lambda") = lambda, Named("g.ee") = g_ee,
                      Named("obj_vec") = obj_vec, Named("supp_theta") = supp_theta,
                      Named("theta.debiased") = theta_debiased,
                      Named("supp_lambda") = supp_lambda, Named("supp_theta") = supp_theta); 
}

