// Andriy Zatserklyaniy <zatserkl@fnal.gov> Jan 7, 2015. See http://mathworld.wolfram.com/Erf.html

#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <iomanip>

using std::cout;     using std::endl;

Double_t erfa_sum(Double_t x)
{
   //
   // Returns sum from the asymptotic expansion of erf(x) for x >> 1
   // 
   // erf(x) ~ 1 - exp(-x*x)/sqrt(pi) * sum
   //
   // sum = 1/x - (1/2)/x^3 + (3/4)/x^5 - (15/8)/x^7 + (105/16)/x^9 - ...
   //

   // asymptotic expansion of the erf(x) for x >> 1
   //
   // precision (NB: seems to be better than std for x > 5)
   //
   // x     erfa-erf    iterations  terminated
   // ----------------------------------------
   // 1     +5.35e-02   3           diverged
   // 2     -6.64e-05   6           diverged
   // 3     +2.02e-09   11          diverged
   // 3.5   -2.58e-12   14          diverged
   // 4     -1.22e-15   18          diverged
   // 4.5    0.00e+00   22          diverged
   // 5      0.00e+00   27          diverged
   // 6      0.00e+00   31          abs(term) < eps
   // 10     0.00e+00   12          abs(term) < eps
   // 20     0.00e+00   8           abs(term) < eps
   // 40     0.00e+00   6           abs(term) < eps
   // 60     0.00e+00   6           abs(term) < eps

   // bool debug = true;
   bool debug = false;

   Double_t erf = 0;
   if (debug) erf = TMath::Erf(x);
   if (debug) cout<< "TMath::Erf(" << x << ") = " << std::setprecision(18) << erf <<endl;

   Double_t eps = 1e-16;

   if (TMath::Abs(x) < eps) return 0.5;

   Double_t sign_function = 1.;

   if (x < 0) {
      x *= -1.;
      sign_function = -1.;
   }

   const Int_t NDIM = 100;
   Double_t term[NDIM];
   Double_t nominator[NDIM];
   Double_t denominator[NDIM];

   Int_t n = 0;
   nominator[n] = 1;
   denominator[n] = 1;
   term[n] = 1./x;

   Double_t sign_term = -1.;
   Double_t x1power = 1./x;
   Double_t x2 = x*x;

   Double_t erfa_local = 0;
   if (debug) {
      erfa_local = sign_function * (1. - TMath::Exp(-x*x)/TMath::Sqrt(TMath::Pi()) * term[0]);
      cout<< n << "\tnominator = " << nominator[n] << " denominator = " << denominator[n] << " x1power = " << x1power << " term = " << term[n] << "   erfa = " << erfa_local << " erfa-erf = " << erfa_local - erf <<endl;
   }

   while (TMath::Abs(term[n]) > eps)
   {
      if (n == NDIM-1) break;
      ++n;

      nominator[n] = nominator[n-1] * (2*n - 1);
      nominator[n] = nominator[n-1] * (2*n - 1);
      denominator[n] = denominator[n-1] * 2;
      x1power /= x2;
      term[n] = sign_term*x1power*nominator[n]/denominator[n];
      sign_term *= -1.;

      if (debug) {
         Double_t sum_local = 0;
         for (int i=n; i>=0; --i) sum_local += term[i];
         erfa_local = sign_function * (1. - TMath::Exp(-x*x)/TMath::Sqrt(TMath::Pi()) * sum_local);
         cout<< n << "\tnominator = " << nominator[n] << " denominator = " << denominator[n] << " x1power = " << x1power << " term = " << term[n] << "   sum/TMath::Pi() = " << sum_local/TMath::Pi() << "   erfa = " << erfa_local << " erfa-erf = " << erfa_local - erf <<endl;
      }

      if (TMath::Abs(term[n]) > TMath::Abs(term[n-1])) {
         if (debug) cout<< "start to diverge: terminate iterations" << endl;
         --n;
         break;
      }
   }
   if (debug) cout<< "n = " << n <<endl;

   Double_t sum = 0;
   for (int i=n; i>=0; --i) sum += term[i];

   //-- Double_t erfa = sign_function * (1. - TMath::Exp(-x*x)/TMath::Sqrt(TMath::Pi()) * sum);

   return sum / TMath::Sqrt(TMath::Pi());
}

/*
Double_t x = 3; erfa_sum_test(&x,0) - TMath::Erf(x)
*/
Double_t erfa_sum_test(Double_t* xx, Double_t*)
{
   //
   // Calculates sum in the asimptotic expansion of erf(x) for x >> 1
   // 
   // erf(x) ~ 1 - exp(-x*x)/sqrt(pi) * sum
   //
   // sum = 1/x - (1/2)/x^3 + (3/4)/x^5 - (15/8)/x^7 + (105/16)/x^9 - ...
   //

   Double_t x = *xx;

   Double_t sign_function = 1.;

   if (x < 0) {
      x *= -1.;
      sign_function = -1.;
   }

   Double_t sum_pi = erfa_sum(x);
   return sign_function * (1. - TMath::Exp(-x*x)*sum_pi);
}

Double_t fxexp(Double_t* xx, Double_t* par) {
   Double_t A = par[0];
   Double_t x0 = par[1];
   Double_t tau = par[2];
   Double_t sigma = TMath::Abs(par[3]);

   Double_t x = *xx - x0;

   bool debug = true;
   // bool debug = false;

   Double_t fun = 0;
   
   Double_t eps = 1e-12;

   if (tau < eps) return fun;

   // unsmeared value
   fun = 0;
   if (x > eps) {
      Double_t frac = x/tau;
      Double_t arg = -frac;
      Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

      // Double_t integral = 2*tau*tau*tau;     // integral_0_inf x*x*exp(-x/tau)
      fun += (A/tau) * frac * exp;
   }

   if (sigma < tau*1e-3) return fun;     // return unsmeared value

   Double_t mu = (x/sigma - sigma/tau)/TMath::Sqrt(2);
   Double_t factor = 2.*sigma*sigma;

   Double_t exp_x2 = 0;
   if (x*x/2/sigma/sigma < 500) exp_x2 = TMath::Exp(-x*x/2/sigma/sigma);
   else return fun;                    // return unsmeared value

   Double_t erf_exp = 0;
   if (mu < -5.) {
      erf_exp = erfa_sum(mu) * exp_x2;
   }
   else {
      erf_exp = mu >= 0? 1 + TMath::Erf(mu): TMath::Erfc(-mu);
      erf_exp *= TMath::Exp(mu*mu - x*x/2/sigma/sigma);
   }

   Double_t term1 = TMath::Sqrt(TMath::Pi()) * mu * erf_exp;
   Double_t term2 = exp_x2;

   fun = A/TMath::Sqrt(2*TMath::Pi())/sigma * factor * 1/tau/tau * (term1 + term2)/2.;

   if (debug) cout<< "x = " << x << " mu = " << mu << " factor = " << factor << " term_erf = " << term_erf << " term1 = " << term1 << " term2 = " << term2 << " fun = " << fun <<endl;
   return fun;
}

TF1* xexp(Double_t xmin, Double_t xmax
           , Double_t A=1
           , Double_t x0=0
           , Double_t tau=2
           , Double_t sigma=5
           )
{
   Double_t par[10];
   Int_t npar = 0;
   par[npar++] = A;
   par[npar++] = x0;
   par[npar++] = tau;
   par[npar++] = sigma;

   TF1* xexp = new TF1("fxexp", fxexp, xmin, xmax, npar);
   xexp->SetParameters(par);
   xexp->SetNpx(1000);

   npar = 0;
   xexp->SetParName(npar++, "A");
   xexp->SetParName(npar++, "x0");
   xexp->SetParName(npar++, "#tau");
   xexp->SetParName(npar++, "#sigma");

   new TCanvas;
   xexp->Draw();

   return xexp;
}

Double_t fx2exp(Double_t* xx, Double_t* par) {
   Double_t A = par[0];
   Double_t x0 = par[1];
   Double_t tau = par[2];
   Double_t sigma = TMath::Abs(par[3]);

   Double_t x = *xx - x0;

   bool debug = true;   // goes negative (mu < -60) for: f2 = x2exp(-5,10, 1,0,0.01,1)
   // bool debug = false;

   Double_t fun = 0;
   
   Double_t eps = 1e-12;

   if (tau < eps) return fun;

   // unsmeared value
   fun = 0;
   if (x > eps) {
      Double_t frac = x/tau;
      Double_t arg = -frac;
      Double_t exp = TMath::Abs(arg) < 500.? TMath::Exp(arg): 0;

      // Double_t integral = 2*tau*tau*tau;     // integral_0_inf x*x*exp(-x/tau)
      fun += (A/tau) * frac*frac/2 * exp;
   }

   if (sigma < tau*1e-3) return fun;     // return unsmeared value

   //
   // smearing
   //

   Double_t mu = (x/sigma - sigma/tau)/TMath::Sqrt(2);
   Double_t factor = sigma*sigma*sigma*2*TMath::Sqrt(2);

   Double_t exp_x2 = 0;
   if (x*x/2/sigma/sigma < 500.) exp_x2 = TMath::Exp(-x*x/2/sigma/sigma);
   else return fun;                    // return unsmeared value

   Double_t erf_exp = 0;
   if (mu < -5.) {
      erf_exp = erfa_sum(mu) * exp_x2;
   }
   else {
      erf_exp = mu >= 0? 1 + TMath::Erf(mu): TMath::Erfc(-mu);
      erf_exp *= TMath::Exp(mu*mu - x*x/2/sigma/sigma);
   }

   Double_t term1 = TMath::Sqrt(TMath::Pi()) * (1. + 2.*mu*mu) * erf_exp;
   Double_t term2 = 2*mu * exp_x2;

   fun = A/TMath::Sqrt(2*TMath::Pi())/sigma * factor * 1/tau/tau/tau/2 * (term1 + term2)/4.;

   if (debug) cout<< "x = " << x << " mu = " << mu << " exp_x2 = " << exp_x2 << " erf_exp = " << erf_exp << " term1 = " << term1 << " term2 = " << term2 << " fun = " << fun <<endl;
   return fun;
}

TF1* x2exp(Double_t xmin, Double_t xmax
           , Double_t A=1
           , Double_t x0=0
           , Double_t tau=2
           , Double_t sigma=5
           )
{
   Double_t par[10];
   Int_t npar = 0;
   par[npar++] = A;
   par[npar++] = x0;
   par[npar++] = tau;
   par[npar++] = sigma;

   TF1* x2exp = new TF1("fx2exp", fx2exp, xmin, xmax, npar);
   x2exp->SetParameters(par);
   x2exp->SetNpx(1000);

   npar = 0;
   x2exp->SetParName(npar++, "A");
   x2exp->SetParName(npar++, "x0");
   x2exp->SetParName(npar++, "#tau");
   x2exp->SetParName(npar++, "#sigma");

   new TCanvas;
   x2exp->Draw();

   return x2exp;
}
