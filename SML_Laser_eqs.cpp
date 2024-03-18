
double fermi(double confine, double EF, double T){
  return ( 1.0/(1.0 + exp((confine-EF)*sm::e0/(sm::kB*T))) );
}

double get_rho_eq(double n, double nc, double confine, double temp){
  
  //Joyce-Dixon approximation coefficients
  const double A1 = 0.3535533905932737622004221810; //1.0/sqrt(8.0);
  const double A2 = -4.95009e-3;
//   const double A3 = 1.48386e-4;
//   const double A4 = -4.42563e-6;
  
  //Pade approximation coefficients - build upon the first two Joyce-Dixon coefficients
  const double K1 = 4.7;
  const double K2 = sqrt(2.0 * fabs(A2) / K1);
  
  const double r = n/nc;
//   const double v = log(r) + A1*r + (K1*log(1 + K2*r) - K1*K2*r);
//   const double EFapprox = v * kB * TEMP / e0;
  const double EFapprox = (sm::kB * temp / sm::e0 ) * ( log(r) + A1*r + (K1*log(1 + K2*r) - K1*K2*r) );

  return fermi(confine, EFapprox, temp);
}
  

auto SMLL = [](vars *x, vars_vec_wdX *Xhist, vars *d, parameters *p){
  std::complex<double>i = std::complex<double>(0.0,1.0); // imaginary i
  
//   const double Rrel = p->R * ( get_rho_eq(x->N, p->NC, p->confine, p->temp) - x->rho);
  const double Rrel = p->R0 * (x->N * x->N / (p->n_tr * p->n_tr) ) * ( get_rho_eq(x->N, p->NC, p->confine, p->temp) - x->rho);
  
  
//   d->E = 0.5 * (p->g*(2.0*x->rho - 1.0)*(1.0-i*p->alpha) - p->kappa) * x->E + p->noise + p->fakeSpE; //no rotating frame adjustment
  d->E = 0.5 * ((p->g*(2.0*x->rho - 1.0) - p->kappa)*(1.0-i*p->alpha)) * x->E + p->noise + p->fakeSpE; //adjust rotating frame by i*kappa to counter the effect of alpha -> optical frequency \approx 0.0 in the no feedback steady state
  
  d->rho = -p->gamma_rho * x->rho + Rrel - p->g* (2.0*x->rho - 1.0) * norm(x->E);
  
  d->N = -p->gamma_N * x->N + p->J/p->h_bulk - 2.0 * (p->n_SML/p->h_bulk) * Rrel; //factor 2 to account for spin degeneracy of the SML states
  
  //feedback
  d->E += + 0.5*p->kappa * p->K * exp(-i*p->C) * Xhist->at_C_cubicHermite(p->tau, IND::E);
	  	  
};
