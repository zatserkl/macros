// Andriy Zatserklyaniy <zatserkl@fnal.gov> Jan 8, 2015

#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

using std::cout;     using std::endl;

// ------------------------------------------------
//
//    FunRC
//
// ------------------------------------------------

class FunRC {
   static Double_t erfa_sum(Double_t x)
   {
      //
      // Returns sum from the asymptotic expansion of erf(x) for x >> 1
      // 
      // erf(x) ~ 1 - exp(-x*x)/sqrt(pi) * sum
      //
      // sum = 1/x - (1/2)/x^3 + (3/4)/x^5 - (15/8)/x^7 + (105/16)/x^9 - ...
      //
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
//private:
public:
   static Double_t fRC0(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];
      Double_t sigma = par[npar++];

      Double_t x = *xx - x0;

      Double_t fRC0 = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t arg1 = -x/tau;
         Double_t exp1 = TMath::Abs(arg1) < 50? TMath::Exp(arg1): 0;

         Double_t arg2 = sigma*sigma/(2.*tau*tau);
         Double_t exp2 = TMath::Abs(arg2) < 50? TMath::Exp(arg2): 0;

         Double_t arg3 = sigma/tau/TMath::Sqrt(2.);
         if (TMath::Abs(sigma) > eps) arg3 -= x/sigma/TMath::Sqrt(2.);
         Double_t erfc = TMath::Erfc(arg3);
         fRC0 = A/(2.*tau) * exp1 * exp2 * erfc;
      }

      return fRC0;
   }
   static Double_t fRC0b(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];
      Double_t sigma = par[npar++];
      Double_t b = par[npar++];

      Double_t x = *xx - x0;

      Double_t fRC0b = b;
      if (TMath::Abs(tau) > eps) {
         Double_t arg1 = -x/tau;
         Double_t exp1 = TMath::Abs(arg1) < 50? TMath::Exp(arg1): 0;

         Double_t arg2 = sigma*sigma/(2.*tau*tau);
         Double_t exp2 = TMath::Abs(arg2) < 50? TMath::Exp(arg2): 0;

         Double_t arg3 = sigma/tau/TMath::Sqrt(2.);
         if (TMath::Abs(sigma) > eps) arg3 -= x/sigma/TMath::Sqrt(2.);
         Double_t erfc = TMath::Erfc(arg3);
         fRC0b = b + A/(2.*tau) * exp1 * exp2 * erfc;
      }

      return fRC0b;
   }
   static Double_t fRC1(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];

      Double_t x = *xx - x0;

      if (x < eps) return 0;

      Double_t fRC1 = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t arg = -x/tau;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         Double_t integral = tau*tau;     // integral_0_inf x*exp(-x/tau)
         fRC1 = A * (x/integral) * exp;
      }

      return fRC1;
   }
   static Double_t fRC1gb(Double_t* xx, Double_t* par)
   {
      Double_t A = par[0];
      Double_t x0 = par[1];
      Double_t tau = par[2];
      Double_t sigma = par[3];
      Double_t bkg = par[4];

      Double_t x = *xx - x0;

      Double_t fun = 0;

      Double_t eps = 1e-12;

      if (tau < eps) return fun;

      // unsmeared value
      fun = bkg;

      if (x > eps) {
         Double_t frac = x/tau;
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         // Double_t integral = 2*tau*tau*tau;     // integral_0_inf x*x*exp(-x/tau)
         fun += (A/tau) * frac * exp;
      }

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

      fun = bkg + A/TMath::Sqrt(2*TMath::Pi())/sigma * factor * 1/tau/tau * (term1 + term2)/2.;

      // cout<< "x = " << x << " mu = " << mu << " factor = " << factor << " term_erf = " << term_erf << " term1 = " << term1 << " term2 = " << term2 << " fun = " << fun <<endl;
      return fun;
   }
   static Double_t fRC1g(Double_t *xx, Double_t *par)
   {
      // needs a room for an extra parameter (bkg)

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];
      Double_t sigma = par[npar++];

      // set bkg to 0
      par[npar] = 0;
      return fRC1gb(xx,par);
   }
   static Double_t fRC1b(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];
      Double_t b = par[npar++];

      Double_t x = *xx - x0;

      if (x < eps) return b;

      Double_t fRC1b = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t frac = x/tau;
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         // Double_t integral = tau*tau;     // integral_0_inf x*exp(-x/tau)
         fRC1b = b + (A/tau) * frac * exp;
      }

      return fRC1b;
   }
   static Double_t fRC1bd(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];
      Double_t b = par[npar++];
      Double_t d = par[npar++];

      Double_t x = *xx - x0;

      Double_t fRC1bd = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t frac = x/tau;
         Double_t fracd = (x+d)/(tau+d);
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         // Double_t integral = tau*tau;     // integral_0_inf x*exp(-x/tau)
         fRC1bd = b + (A/tau) * fracd * exp;

         if (fRC1bd < 0) return 0;
      }

      return fRC1bd;
   }
   static Double_t fRC2(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];

      Double_t x = *xx - x0;

      if (x < eps) return 0;

      Double_t fRC2 = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t frac = x/tau;
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         // Double_t integral = 2*tau*tau*tau;     // integral_0_inf x*x*exp(-x/tau)
         fRC2 = (A/tau) * frac*frac/2 * exp;
      }

      return fRC2;
   }
   static Double_t dfRC2(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];

      Double_t x = *xx - x0;

      if (x < eps) return 0;

      Double_t fun = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t frac = x/tau;
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         // Double_t integral = 2*tau*tau*tau;     // integral_0_inf x*x*exp(-x/tau)
         fun = (2./x - 1./tau) * (A/tau) * frac*frac/2 * exp;
      }

      return fun;
   }
   static Double_t fRC2b(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];
      Double_t b = par[npar++];

      Double_t x = *xx - x0;

      if (x < eps) return b;

      Double_t fRC2b = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t frac = x/tau;
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         // Double_t integral = 2*tau*tau*tau;     // integral_0_inf x*x*exp(-x/tau)
         fRC2b = b + (A/tau) * frac*frac/2 * exp;
      }

      return fRC2b;
   }
   static Double_t fRC2gb(Double_t* xx, Double_t* par)
   {
      Double_t A = par[0];
      Double_t x0 = par[1];
      Double_t tau = par[2];
      Double_t sigma = par[3];
      Double_t bkg = par[4];

      Double_t x = *xx - x0;

      Double_t eps = 1e-12;

      // bool debug = true;   // goes negative (mu < -60) for: f2 = x2exp(-5,10, 1,0,0.01,1)
      bool debug = false;

      // unsmeared value
      Double_t fun = bkg;
      if (x > eps) {
         Double_t frac = x/tau;
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 500.? TMath::Exp(arg): 0;
         // Double_t integral = 2*tau*tau*tau;     // integral_0_inf x*x*exp(-x/tau)
         fun += (A/tau) * frac*frac/2 * exp;
      }

      //-- if (sigma < tau*1e-3) return fun;     // return unsmeared value

      //
      // smearing
      //

      Double_t mu = (x/sigma - sigma/tau)/TMath::Sqrt(2);
      Double_t factor = sigma*sigma*sigma*2*TMath::Sqrt(2);

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

      Double_t term1 = TMath::Sqrt(TMath::Pi()) * (1. + 2.*mu*mu) * erf_exp;
      Double_t term2 = 2*mu * exp_x2;

      fun = bkg + A/TMath::Sqrt(2*TMath::Pi())/sigma * factor * 1/tau/tau/tau/2 * (term1 + term2)/4.;

      if (debug) cout<< "x = " << x << " mu = " << mu << " exp_x2 = " << exp_x2 << " erf_exp = " << erf_exp << " term1 = " << term1 << " term2 = " << term2 << " fun = " << fun <<endl;
      return fun;
   }
   static Double_t fRC3(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];

      Double_t x = *xx - x0;

      if (x < eps) return 0;

      Double_t fRC3 = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t frac = x/tau;
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         // Double_t integral = (2*3)*tau*tau*tau*tau;   // integral_0_inf x*x*x*exp(-x/tau)
         fRC3 = (A/tau) * frac*frac*frac/(2*3) * exp;
      }

      return fRC3;
   }
   static Double_t fRC3b(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];
      Double_t b = par[npar++];

      Double_t x = *xx - x0;

      if (x < eps) return b;

      Double_t fRC3b = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t frac = x/tau;
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         // Double_t integral = (2*3)*tau*tau*tau*tau;   // integral_0_inf x*x*x*exp(-x/tau)
         fRC3b = b + (A/tau) * frac*frac*frac/(2*3) * exp;
      }

      return fRC3b;
   }
   static Double_t fRC4(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];

      Double_t x = *xx - x0;

      if (x < eps) return 0;

      Double_t fRC4 = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t frac = x/tau;
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         // Double_t integral = (2*3*4)*tau*tau*tau*tau*tau;   // integral_0_inf x*x*x*x*exp(-x/tau)
         fRC4 = (A/frac) * frac*frac*frac*frac/(2*3*4) * exp;
      }

      return fRC4;
   }
   static Double_t fRC4b(Double_t *xx, Double_t *par)
   {
      const Double_t eps = 1e-12;

      Int_t npar = 0;
      Double_t A = par[npar++];
      Double_t x0 = par[npar++];
      Double_t tau = par[npar++];
      Double_t b = par[npar++];

      Double_t x = *xx - x0;

      if (x < eps) return b;

      Double_t fRC4b = 0;
      if (TMath::Abs(tau) > eps) {
         Double_t frac = x/tau;
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 50? TMath::Exp(arg): 0;

         // Double_t integral = (2*3*4)*tau*tau*tau*tau*tau;   // integral_0_inf x*x*x*x*exp(-x/tau)
         fRC4b = b + (A/frac) * frac*frac*frac*frac/(2*3*4) * exp;
      }

      return fRC4b;
   }
public:
   static TF1* RC0(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau, Double_t sigma
         , const char* name="fRC0"
         )
   {
      TF1* fps0 = new TF1(name, fRC0, xmin, xmax, 4);
      fps0->SetNpx(1024);
      Int_t npar = 0;
      fps0->SetParName(npar++, "A");
      fps0->SetParName(npar++, "x0");
      fps0->SetParName(npar++, "#tau");
      fps0->SetParName(npar++, "#sigma");
      npar = 0;
      fps0->SetParameter(npar++, A);
      fps0->SetParameter(npar++, x0);
      fps0->SetParameter(npar++, tau);
      fps0->SetParameter(npar++, sigma);
      return fps0;
   }
   static TF1* RC0b(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau, Double_t sigma
         , Double_t b=0
         , const char* name="fRC0b"
         )
   {
      TF1* fps0b = new TF1(name, fRC0b, xmin, xmax, 5);
      fps0b->SetNpx(1024);
      Int_t npar = 0;
      fps0b->SetParName(npar++, "A");
      fps0b->SetParName(npar++, "x0");
      fps0b->SetParName(npar++, "#tau");
      fps0b->SetParName(npar++, "#sigma");
      fps0b->SetParName(npar++, "b");
      npar = 0;
      fps0b->SetParameter(npar++, A);
      fps0b->SetParameter(npar++, x0);
      fps0b->SetParameter(npar++, tau);
      fps0b->SetParameter(npar++, sigma);
      fps0b->SetParameter(npar++, b);
      return fps0b;
   }
   static TF1* RC1bd(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , Double_t b=0
         , Double_t d=0
         , const char* name="fRC1bd"
         )
   {
      TF1* fps1bd = new TF1(name, fRC1bd, xmin, xmax, 5);
      fps1bd->SetNpx(1024);
      Int_t npar = 0;
      fps1bd->SetParName(npar++, "A");
      fps1bd->SetParName(npar++, "x0");
      fps1bd->SetParName(npar++, "#tau");
      fps1bd->SetParName(npar++, "b");
      fps1bd->SetParName(npar++, "d");
      npar = 0;
      fps1bd->SetParameter(npar++, A);
      fps1bd->SetParameter(npar++, x0);
      fps1bd->SetParameter(npar++, tau);
      fps1bd->SetParameter(npar++, b);
      fps1bd->SetParameter(npar++, d);
      return fps1bd;
   }
   static TF1* RC1(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , const char* name="fRC1"
         )
   {
      TF1* fps1 = new TF1(name, fRC1, xmin, xmax, 3);
      fps1->SetNpx(1024);
      Int_t npar = 0;
      fps1->SetParName(npar++, "A");
      fps1->SetParName(npar++, "x0");
      fps1->SetParName(npar++, "#tau");
      npar = 0;
      fps1->SetParameter(npar++, A);
      fps1->SetParameter(npar++, x0);
      fps1->SetParameter(npar++, tau);
      return fps1;
   }
   static TF1* RC1b(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , Double_t b=0
         , const char* name="fRC1b"
         )
   {
      TF1* fps1b = new TF1(name, fRC1b, xmin, xmax, 4);
      fps1b->SetNpx(1024);
      Int_t npar = 0;
      fps1b->SetParName(npar++, "A");
      fps1b->SetParName(npar++, "x0");
      fps1b->SetParName(npar++, "#tau");
      fps1b->SetParName(npar++, "b");
      npar = 0;
      fps1b->SetParameter(npar++, A);
      fps1b->SetParameter(npar++, x0);
      fps1b->SetParameter(npar++, tau);
      fps1b->SetParameter(npar++, b);
      return fps1b;
   }
   static TF1* RC1g(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , Double_t sigma
         , const char* name="fRC1g"
         )
   {
      TF1* fps1g = new TF1(name, fRC1g, xmin, xmax, 4);
      fps1g->SetNpx(1024);
      Int_t npar = 0;
      fps1g->SetParName(npar++, "A");
      fps1g->SetParName(npar++, "x0");
      fps1g->SetParName(npar++, "#tau");
      fps1g->SetParName(npar++, "#sigma");
      npar = 0;
      fps1g->SetParameter(npar++, A);
      fps1g->SetParameter(npar++, x0);
      fps1g->SetParameter(npar++, tau);
      fps1g->SetParameter(npar++, sigma);
      return fps1g;
   }
   static TF1* RC1gb(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , Double_t sigma
         , Double_t bkg = 0
         , const char* name="fRC1gb"
         )
   {
      TF1* fps1gb = new TF1(name, fRC1gb, xmin, xmax, 5);
      //--------------------------------------------------fps1gb->SetNpx(1024);
      Int_t npar = 0;
      fps1gb->SetParName(npar++, "A");
      fps1gb->SetParName(npar++, "x0");
      fps1gb->SetParName(npar++, "#tau");
      fps1gb->SetParName(npar++, "#sigma");
      fps1gb->SetParName(npar++, "bkg");
      npar = 0;
      fps1gb->SetParameter(npar++, A);
      fps1gb->SetParameter(npar++, x0);
      fps1gb->SetParameter(npar++, tau);
      fps1gb->SetParameter(npar++, sigma);
      fps1gb->SetParameter(npar++, bkg);
      return fps1gb;
   }
   static TF1* RC2(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , const char* name="fRC2"
         )
   {
      TF1* fps2 = new TF1(name, fRC2, xmin, xmax, 3);
      fps2->SetNpx(1024);
      Int_t npar = 0;
      fps2->SetParName(npar++, "A");
      fps2->SetParName(npar++, "x0");
      fps2->SetParName(npar++, "#tau");
      npar = 0;
      fps2->SetParameter(npar++, A);
      fps2->SetParameter(npar++, x0);
      fps2->SetParameter(npar++, tau);
      return fps2;
   }
   static TF1* RC2b(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , Double_t b=0
         , const char* name="fRC2b"
         )
   {
      TF1* fps2b = new TF1(name, fRC2b, xmin, xmax, 4);
      fps2b->SetNpx(1024);
      Int_t npar = 0;
      fps2b->SetParName(npar++, "A");
      fps2b->SetParName(npar++, "x0");
      fps2b->SetParName(npar++, "#tau");
      fps2b->SetParName(npar++, "b");
      npar = 0;
      fps2b->SetParameter(npar++, A);
      fps2b->SetParameter(npar++, x0);
      fps2b->SetParameter(npar++, tau);
      fps2b->SetParameter(npar++, b);
      return fps2b;
   }
   static TF1* RC2gb(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , Double_t sigma=0
         , Double_t b=0
         , const char* name="fRC2gb"
         )
   {
      TF1* fps2gb = new TF1(name, fRC2gb, xmin, xmax, 5);
      fps2gb->SetNpx(1024);
      Int_t npar = 0;
      fps2gb->SetParName(npar++, "A");
      fps2gb->SetParName(npar++, "x0");
      fps2gb->SetParName(npar++, "#tau");
      fps2gb->SetParName(npar++, "#sigma");
      fps2gb->SetParName(npar++, "b");
      npar = 0;
      fps2gb->SetParameter(npar++, A);
      fps2gb->SetParameter(npar++, x0);
      fps2gb->SetParameter(npar++, tau);
      fps2gb->SetParameter(npar++, sigma);
      fps2gb->SetParameter(npar++, b);
      return fps2gb;
   }
   static TF1* RC3(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , const char* name="fRC3"
         )
   {
      TF1* fps3 = new TF1(name, fRC3, xmin, xmax, 3);
      fps3->SetNpx(1024);
      Int_t npar = 0;
      fps3->SetParName(npar++, "A");
      fps3->SetParName(npar++, "x0");
      fps3->SetParName(npar++, "#tau");
      npar = 0;
      fps3->SetParameter(npar++, A);
      fps3->SetParameter(npar++, x0);
      fps3->SetParameter(npar++, tau);
      return fps3;
   }
   static TF1* RC3b(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , Double_t b=0
         , const char* name="fRC3b"
         )
   {
      TF1* fps3b = new TF1(name, fRC3b, xmin, xmax, 4);
      fps3b->SetNpx(1024);
      Int_t npar = 0;
      fps3b->SetParName(npar++, "A");
      fps3b->SetParName(npar++, "x0");
      fps3b->SetParName(npar++, "#tau");
      fps3b->SetParName(npar++, "b");
      npar = 0;
      fps3b->SetParameter(npar++, A);
      fps3b->SetParameter(npar++, x0);
      fps3b->SetParameter(npar++, tau);
      fps3b->SetParameter(npar++, b);
      return fps3b;
   }
   static TF1* RC4(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , const char* name="fRC4"
         )
   {
      TF1* fps4 = new TF1(name, fRC4, xmin, xmax, 3);
      fps4->SetNpx(1024);
      Int_t npar = 0;
      fps4->SetParName(npar++, "A");
      fps4->SetParName(npar++, "x0");
      fps4->SetParName(npar++, "#tau");
      npar = 0;
      fps4->SetParameter(npar++, A);
      fps4->SetParameter(npar++, x0);
      fps4->SetParameter(npar++, tau);
      return fps4;
   }
   static TF1* RC4b(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , Double_t b=0
         , const char* name="fRC4b"
         )
   {
      TF1* fps4b = new TF1(name, fRC4b, xmin, xmax, 4);
      fps4b->SetNpx(1024);
      Int_t npar = 0;
      fps4b->SetParName(npar++, "A");
      fps4b->SetParName(npar++, "x0");
      fps4b->SetParName(npar++, "#tau");
      fps4b->SetParName(npar++, "b");
      npar = 0;
      fps4b->SetParameter(npar++, A);
      fps4b->SetParameter(npar++, x0);
      fps4b->SetParameter(npar++, tau);
      fps4b->SetParameter(npar++, b);
      return fps4b;
   }
   //FunRC() {}
};
