#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TVirtualFFT.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TObject.h>

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
   static Double_t fRC1g(Double_t* xx, Double_t* par)
   {
      Double_t A = par[0];
      Double_t x0 = par[1];
      Double_t tau = par[2];
      Double_t sigma = TMath::Abs(par[3]);

      Double_t x = *xx - x0;

      bool debug = false;
      // bool debug = true;   // goes negative (mu < -60) for: f2 = x2exp(-5,10, 1,0,0.01,1)

      Double_t eps = 1e-12;

      Double_t fun = 0;
      if (tau < eps) return fun;

      // unsmeared value
      fun = 0;
      if (x > eps) {
         Double_t frac = x/tau;
         Double_t arg = -frac;
         Double_t exp = TMath::Abs(arg) < 500.? TMath::Exp(arg): 0;
         // Double_t integral = 2*tau*tau*tau;     // integral_0_inf x*x*exp(-x/tau)
         fun += (A/tau) * frac * exp;
      }

      if (sigma < tau*1e-3) return fun;     // return unsmeared value

      Double_t mu = (x/sigma - sigma/tau)/TMath::Sqrt(2);
      Double_t factor = 2.*sigma*sigma;

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

      Double_t term1 = TMath::Sqrt(TMath::Pi()) * mu * erf_exp;
      Double_t term2 = exp_x2;

      fun = A/TMath::Sqrt(2*TMath::Pi())/sigma * factor * 1/tau/tau * (term1 + term2)/2.;

      if (debug) cout<< "x = " << x << " mu = " << mu << " factor = " << factor << " erf_exp = " << erf_exp << " term1 = " << term1 << " term2 = " << term2 << " fun = " << fun <<endl;
      return fun;
   }
   static Double_t fRC1gb(Double_t* xx, Double_t* par)
   {
      // Double_t A = par[0];
      // Double_t x0 = par[1];
      // Double_t tau = par[2];
      // Double_t sigma = TMath::Abs(par[3]);
      Double_t bkg = par[4];
      return bkg + fRC1g(xx,par);
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
   static Double_t fRC1flat(Double_t* xx, Double_t* par) {
      Double_t A = par[0];
      Double_t x0 = par[1];
      Double_t tau = par[2];
      Double_t d = par[3];
      Double_t bkg = par[4];

      Double_t x = *xx - x0;

      Double_t eps = 1e-12;

      Double_t fun = bkg;

      if (x < eps) return fun;
      if (d < eps) return fun;
      if (tau < eps) return fun;

      Double_t limit = x > d? d: x;

      Double_t factor = A/d/tau/tau;

      fun += factor * TMath::Exp(-x/tau) * (tau*TMath::Exp(limit/tau)*(x - limit + tau) - tau*(tau+x));
      return fun;
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
   static Double_t fRC2g(Double_t* xx, Double_t* par)
   {
      Double_t A = par[0];
      Double_t x0 = par[1];
      Double_t tau = par[2];
      Double_t sigma = TMath::Abs(par[3]);

      Double_t x = *xx - x0;

      bool debug = false;
      // bool debug = true;   // goes negative (mu < -60) for: f2 = x2exp(-5,10, 1,0,0.01,1)

      Double_t eps = 1e-12;

      Double_t fun = 0;
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
   static Double_t fRC2gb(Double_t* xx, Double_t* par)
   {
      // Double_t A = par[0];
      // Double_t x0 = par[1];
      // Double_t tau = par[2];
      // Double_t sigma = TMath::Abs(par[3]);
      Double_t bkg = par[4];
      return bkg + fRC2g(xx,par);
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
   static TF1* RC1flat(Double_t xmin, Double_t xmax
         , Double_t A, Double_t x0, Double_t tau
         , Double_t d=0
         , Double_t bkg=0
         , const char* name="fRC1flat"
         )
   {
      TF1* fflat = new TF1(name, fRC1flat, xmin, xmax, 5);
      fflat->SetNpx(1024);
      Int_t npar = 0;
      fflat->SetParName(npar++, "A");
      fflat->SetParName(npar++, "x0");
      fflat->SetParName(npar++, "#tau");
      fflat->SetParName(npar++, "d");
      fflat->SetParName(npar++, "bkg");
      npar = 0;
      fflat->SetParameter(npar++, A);
      fflat->SetParameter(npar++, x0);
      fflat->SetParameter(npar++, tau);
      fflat->SetParameter(npar++, d);
      fflat->SetParameter(npar++, bkg);
      return fflat;
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
      fps1gb->SetNpx(1024);
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

//--------- Fast Fourier Transform -----------

Double_t fft_filter(Int_t N, Double_t* x, Double_t* y, Double_t fcut=0.2)  // fcut in GHz
{
   // NB: array y is in use for both input and output

   TVirtualFFT* fft = TVirtualFFT::FFT(1, &N, "R2C ES K");

   fft->SetPoints(y);
   fft->Transform();

   Double_t re[20000];
   Double_t im[20000];
   fft->GetPointsComplex(re, im);

   Double_t baseline = re[0]/N;
   // cout<< "re[0]/N = " << re[0]/N <<endl;

   // cut frequency band between the real world frequencies fcut1 and fcut2
   Double_t L = x[N-1] - x[0];
   Int_t scut = fcut*L;             // convert from real world frequency to frequency sample number

   for (int isample=0; isample<N/2+1; ++isample) {
      if (isample > scut) {
         re[isample] = 0;
         im[isample] = 0;
      }
   }

   // inverse transform
   TVirtualFFT* fft_inv = TVirtualFFT::FFT(1, &N, "C2R M K");  // inverse trasform
   fft_inv->SetPointsComplex(re, im);
   fft_inv->Transform();

   // get output of the inverse tranform and scale it to 1/N

   Double_t* output_inv = fft_inv->GetPointsReal();
   for (int isample=0; isample<N; ++isample) {
      y[isample] = output_inv[isample]/Double_t(N);
   }

   return baseline;     // re[0]/N
}

//--------------- moving average -------------------

void moving_average(Int_t naver, Int_t np, const Double_t* y, Double_t* yaver)   // Double_t version
{
   // assume naver odd
   if (naver % 2 == 0) {
      naver += 1;
      cout<< "\n--> set naver to odd number " << naver <<endl;
   }

   Double_t acc = 0;
   for (int i=0; i<naver; ++i) acc += y[i];

   yaver[0] = acc/naver;
   Int_t icurr = naver/2;
   for (int i=0; i<=icurr; ++i) yaver[i] = y[0];

   Int_t ifirst = 0;
   Int_t inext = naver;

   while (inext<np) {
      acc += y[inext++] - y[ifirst++];
      yaver[++icurr] = acc/naver;
   }

   for (int i=icurr+1; i<np; ++i) yaver[i] = yaver[icurr];
}

void moving_average(Int_t naver, Int_t np, const Float_t* y, Float_t* yaver)  // Float_t version
{
   // assume naver odd
   if (naver % 2 == 0) {
      naver += 1;
      cout<< "\n--> set naver to odd number " << naver <<endl;
   }

   Float_t acc = 0;
   for (int i=0; i<naver; ++i) acc += y[i];

   yaver[0] = acc/naver;
   Int_t icurr = naver/2;
   for (int i=0; i<=icurr; ++i) yaver[i] = y[0];

   Int_t ifirst = 0;
   Int_t inext = naver;

   while (inext<np) {
      acc += y[inext++] - y[ifirst++];
      yaver[++icurr] = acc/naver;
   }

   for (int i=icurr+1; i<np; ++i) yaver[i] = yaver[icurr];
}

//-------------------------- end of tools --------------------------

void wfmRun(const char* ifname0, const char* ifname1, const char* ifname2, const char* ifname3
         , Int_t& np0, Float_t* t0, Float_t* v0
         , Int_t& np1, Float_t* t1, Float_t* v1
         , Int_t& np2, Float_t* t2, Float_t* v2
         , Int_t& np3, Float_t* t3, Float_t* v3
         )
{
   std::ifstream ifile0(ifname0);
   if (!ifile0) {
      cout<< "Could not open file " << ifname0 <<endl;
      return;
   }
   std::ifstream ifile1(ifname1);
   if (!ifile1) {
      cout<< "Could not open file " << ifname1 <<endl;
      return;
   }
   std::ifstream ifile2(ifname2);
   if (!ifile2) {
      cout<< "Could not open file " << ifname2 <<endl;
      return;
   }
   std::ifstream ifile3(ifname3);
   if (!ifile3) {
      cout<< "Could not open file " << ifname3 <<endl;
      return;
   }

   bool debug = false;
   if (std::strcmp(ifname0, "TB_19Nov/Run9_600V_500V/C1_Ch3G10_CSA_Ch2G10_BB_Ch1G3_BB_Ch4Trigger00000.txt") == 0) debug = true;

   std::string line;
   char comma;

   for (int i=0; i<5; ++i) std::getline(ifile0,line);   // skip header
   np0 = 0;
   while (ifile0 >> t0[np0] >> comma >> v0[np0]) {
      //if (np0 < 10) cout<< "   t0[" << np0 << "] = " << t0[np0] << " v0[" << np0 << "] = " << v0[np0] <<endl;
      t0[np0] *= 1e9;
      v0[np0] *= 1e3;
      ++np0;
   }
   //cout<< "np0 = " << np0 << " line = " << line <<endl;

   for (int i=0; i<5; ++i) std::getline(ifile1,line);   // skip header
   np1 = 0;
   while (ifile1 >> t1[np1] >> comma >> v1[np1]) {
      //if (np1 < 10) cout<< "   t1[" << np1 << "] = " << t1[np1] << " v1[" << np1 << "] = " << v1[np1] <<endl;
      t1[np1] *= 1e9;
      v1[np1] *= 1e3;
      ++np1;
   }
   //cout<< "np1 = " << np1 << " line = " << line <<endl;

   for (int i=0; i<5; ++i) std::getline(ifile2,line);   // skip header
   np2 = 0;
   while (ifile2 >> t2[np2] >> comma >> v2[np2]) {
      //if (np2 < 10) cout<< "   t2[" << np2 << "] = " << t2[np2] << " v2[" << np2 << "] = " << v2[np2] <<endl;
      t2[np2] *= 1e9;
      v2[np2] *= 1e3;
      ++np2;
   }
   //cout<< "np2 = " << np2 << " line = " << line <<endl;

   for (int i=0; i<5; ++i) std::getline(ifile3,line);   // skip header
   np3 = 0;
   while (ifile3 >> t3[np3] >> comma >> v3[np3]) {
      //if (np3 < 10) cout<< "   t3[" << np3 << "] = " << t3[np3] << " v3[" << np3 << "] = " << v3[np3] <<endl;
      //if (debug) cout<< "   t3[" << np3 << "] = " << t3[np3] << " v3[" << np3 << "] = " << v3[np3] <<endl;
      t3[np3] *= 1e9;
      v3[np3] *= 1e3;
      ++np3;
   }
   //cout<< "np3 = " << np3 << " line = " << line <<endl;
}

// wfmRun("TB_19Nov/Run9_600V_500V","Ch3G10_CSA_Ch2G10_BB_Ch1G3_BB_Ch4Trigger",23772)
void wfmRun(const char* dir, const char* pattern, Int_t nlast)
{
   TFile* ofile = new TFile(Form("%s_%s.root",dir,pattern), "recreate");
   TTree* tree = new TTree("wfm", Form("%s.root",pattern));
   tree->SetMarkerColor(602);

   Float_t t0[100000], v0[100000];
   Float_t t1[100000], v1[100000];
   Float_t t2[100000], v2[100000];
   Float_t t3[100000], v3[100000];
   Int_t np0 = 0;
   Int_t np1 = 0;
   Int_t np2 = 0;
   Int_t np3 = 0;

   // read the first file to get the number of points
   std::string fname[4];
   for (int i=0; i<4; ++i) {
      fname[i] = Form("%s/C%d_%s%05d.txt",dir,i+1,pattern,0);
      cout<< "fname[" << i << "] = " << fname[i] <<endl;
   }

   wfmRun(fname[0].c_str(),fname[1].c_str(), fname[2].c_str(), fname[3].c_str()
            , np0, t0, v0
            , np1, t1, v1
            , np2, t2, v2
            , np3, t3, v3
            );

   if (np0 != np1) {
      cout<< "the number of points mismatch in files " << fname[0] << " and " << fname[1] <<endl;
      return;
   }
   cout<< "the number of points is " << np0 <<endl;

   tree->Branch("t1", &t0, Form("t1[%d]/F",np0));
   tree->Branch("v1", &v0, Form("v1[%d]/F",np0));
   tree->Branch("t2", &t1, Form("t2[%d]/F",np1));
   tree->Branch("v2", &v1, Form("v2[%d]/F",np1));
   tree->Branch("t3", &t2, Form("t3[%d]/F",np2));
   tree->Branch("v3", &v2, Form("v3[%d]/F",np2));
   tree->Branch("t4", &t3, Form("t4[%d]/F",np3));
   tree->Branch("v4", &v3, Form("v4[%d]/F",np3));

   tree->Fill();                                                  // fill the first event

   //cout<< "read the rest of the data" <<endl;
   //cout<< "before fname4 = " << fname4 <<endl;
   //ifile3 >> fname4;
   //cout<< "after fname4 = " << fname4 <<endl;

   for (int ievent=0; ievent<=nlast; ++ievent)
   {
      for (int ifile=0; ifile<4; ++ifile) {
         fname[ifile] = Form("%s/C%d_%s%05d.txt",dir,ifile+1,pattern,ievent);
      }

      wfmRun(fname[0].c_str(),fname[1].c_str(),fname[2].c_str(),fname[3].c_str(), np0,t0,v0, np1,t1,v1, np2,t2,v2, np3,t3,v3);
      tree->Fill();
      if (tree->GetEntries() % 100 == 0) cout<< "read " << tree->GetEntries() << " events" <<endl;
   }

   cout<< "Write " << tree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
   ofile->Write();
   // cout<< "Wrote " << tree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
}

namespace TreeRes {
   // trigger: channel 4
   Int_t npar[4];
   Double_t par0[4];
   Double_t par1[4];
   Double_t par2[4];
   Double_t par3[4];
   Double_t par4[4];
   Double_t par5[4];
   Double_t tstamp[4];
   Double_t pmax[4];
   Double_t bkg_raw[4];
   Double_t ebkg_raw[4];
   Double_t bkg[4];
   Double_t ebkg[4];
   Double_t xm[4];
   Double_t ym[4];
   Double_t x5mV[4];
   Double_t x10mV[4];
   Double_t tx[4];

   void clear() {
      for (int ich=0; ich<4; ++ich) {
         npar[ich] = 0;
         tstamp[ich] = 0;
         pmax[ich] = 0;
         bkg_raw[ich] = 0;
         ebkg_raw[ich] = 0;
         bkg[ich] = 0;
         ebkg[ich] = 0;
         xm[ich] = 0;
         ym[ich] = 0;
         x5mV[ich] = 0;
         x10mV[ich] = 0;
         tx[ich] = 0;
         par0[ich] = 0;
         par1[ich] = 0;
         par2[ich] = 0;
         par3[ich] = 0;
         par4[ich] = 0;
         par5[ich] = 0;
      }
   }

   void book(TTree* tree) {
      tree->Branch("npar",          &npar,         "npar[4]/I");
      tree->Branch("par0",          &par0,         "par0[4]/D");
      tree->Branch("par1",          &par1,         "par1[4]/D");
      tree->Branch("par2",          &par2,         "par2[4]/D");
      tree->Branch("par3",          &par3,         "par3[4]/D");
      tree->Branch("par4",          &par4,         "par4[4]/D");
      tree->Branch("par5",          &par5,         "par5[4]/D");
      tree->Branch("tstamp",          &tstamp,         "tstamp[4]/D");
      tree->Branch("pmax",          &pmax,         "pmax[4]/D");
      tree->Branch("bkg_raw",          &bkg_raw,         "bkg_raw[4]/D");
      tree->Branch("ebkg_raw",          &ebkg_raw,         "ebkg_raw[4]/D");
      tree->Branch("bkg",          &bkg,         "bkg[4]/D");
      tree->Branch("ebkg",          &ebkg,         "ebkg[4]/D");
      tree->Branch("xm",          &xm,         "xm[4]/D");
      tree->Branch("ym",          &ym,         "ym[4]/D");
      tree->Branch("x5mV",          &x5mV,         "x5mV[4]/D");
      tree->Branch("x10mV",          &x10mV,         "x10mV[4]/D");
      tree->Branch("tx",          &tx,         "tx[4]/D");
   }

   void connect(TTree* tree) {
      tree->SetBranchAddress("npar",          &npar);
      tree->SetBranchAddress("par0",          &par0);
      tree->SetBranchAddress("par1",          &par1);
      tree->SetBranchAddress("par2",          &par2);
      tree->SetBranchAddress("par3",          &par3);
      tree->SetBranchAddress("par4",          &par4);
      tree->SetBranchAddress("par5",          &par5);
      tree->SetBranchAddress("tstamp",          &tstamp);
      tree->SetBranchAddress("pmax",          &pmax);
      tree->SetBranchAddress("bkg_raw",          &bkg_raw);
      tree->SetBranchAddress("ebkg_raw",          &ebkg_raw);
      tree->SetBranchAddress("bkg",          &bkg);
      tree->SetBranchAddress("ebkg",          &ebkg);
      tree->SetBranchAddress("xm",          &xm);
      tree->SetBranchAddress("ym",          &ym);
      tree->SetBranchAddress("x5mV",          &x5mV);
      tree->SetBranchAddress("x10mV",          &x10mV);
      tree->SetBranchAddress("tx",          &tx);
   }
}  // namespace TreeRes

void trigger_event(Int_t event, Int_t chan=4, Double_t thres=100., Float_t sign=-1., Int_t N=0)
{
   TTree* wfm = (TTree*) gDirectory->Get("wfm");
   if (!wfm) {
      cout<< "Could not find tree 'wfm'" <<endl;
      return;
   }

   Float_t tFloat[20000];
   Float_t vFloat[20000];
   Double_t t[20000];
   Double_t v[20000];

   wfm->SetBranchAddress(Form("t%d",chan), &tFloat);
   wfm->SetBranchAddress(Form("v%d",chan), &vFloat);
   Int_t npoints = wfm->GetLeaf(Form("t%d",chan))->GetLen();

   if (N == 0) N = npoints;

   if (wfm->LoadTree(event) < 0) {
      cout<< "No such entry: " << event <<endl;
      return;
   }
   wfm->GetEntry(event);

   for (int i=0; i<npoints; ++i) {
      t[i] = tFloat[i];
      v[i] = sign*vFloat[i];
   }

   TGraph* g = new TGraph(npoints,t,v);
   g->SetNameTitle(Form("g_evt_%d_ch_%d",event,chan), Form("g_evt_%d_ch_%d",event,chan));
   g->SetMarkerStyle(7);
   g->SetMarkerColor(602);
   g->SetLineColor(602);

   new TCanvas;
   g->GetXaxis()->SetRangeUser(-2,1);
   g->Draw("ap");

   // select three points to fit
   Int_t thres_i = 0;
   while (v[thres_i] < thres && thres_i < npoints) ++thres_i;
   if (thres_i == npoints) return;            // no signal

   Int_t maximum_i = thres_i;
   for (int i=thres_i+1; i<npoints; ++i) {
      if (v[i] > v[maximum_i]) maximum_i = i;
      if (v[i] < v[maximum_i]) break;              // stop at the first maximum
   }

   Int_t halfmax_i = maximum_i - 1;
   while (v[halfmax_i] > 0.5*v[maximum_i] && halfmax_i > 0) --halfmax_i;

   // point halfmax_i is below the halfmax, the point halfmax+1 is above the maximum
   // add third point with max distance
   Double_t fit_x1 = 0;
   Double_t fit_x2 = 0;
   Double_t eps = 1e-7;
   if (v[halfmax_i+2] - v[halfmax_i+1] > v[halfmax_i] - v[halfmax_i-1]) {
      fit_x1 = t[halfmax_i] - eps;
      fit_x2 = t[halfmax_i+2] + eps;
   }
   else {
      fit_x1 = t[halfmax_i-1] - eps;
      fit_x2 = t[halfmax_i+1] + eps;
   }
   //cout<< "fit_x1 = " << fit_x1 << " fit_x2 = " << fit_x2 <<endl;

   g->Fit("pol1", "", "", fit_x1,fit_x2);
   Double_t delta = 0.1;
   g->GetFunction("pol1")->SetRange(fit_x1-delta, fit_x2+delta);
   Double_t p0 = g->GetFunction("pol1")->GetParameter(0);
   Double_t p1 = g->GetFunction("pol1")->GetParameter(1);
   Double_t tstamp = -p0/p1;
   cout<< "tstamp = " << tstamp <<endl;

   wfm->ResetBranchAddresses();
}

void trigger_fit(Int_t N, Double_t* t, Double_t* v, Double_t thres, Double_t* linepar, Double_t& tstamp, Int_t chan=0, Int_t event=-1, bool debug=false)
{
   linepar[0] = 0;
   linepar[1] = 0;
   tstamp = t[0];

   TGraph* g = new TGraph(N,t,v);
   g->SetNameTitle(Form("g_evt_%d_ch_%d",event,chan), Form("g_evt_%d_ch_%d",event,chan));
   g->SetMarkerStyle(7);
   g->SetMarkerColor(602);
   g->SetLineColor(602);

   if (debug) {
      new TCanvas;
      g->GetXaxis()->SetRangeUser(-2,1);
      g->Draw("ap");
   }

   // select three points to fit
   Int_t thres_i = 0;
   while (v[thres_i] < thres && thres_i < N) ++thres_i;
   if (thres_i == N) return;            // no signal

   Int_t maximum_i = thres_i;
   for (int i=thres_i+1; i<N; ++i) {
      if (v[i] > v[maximum_i]) maximum_i = i;
      if (v[i] < v[maximum_i]) break;              // stop at the first maximum
   }

   Int_t halfmax_i = maximum_i - 1;
   while (v[halfmax_i] > 0.5*v[maximum_i] && halfmax_i > 0) --halfmax_i;

   // point halfmax_i is below the halfmax, the point halfmax+1 is above the maximum
   // add third point with max distance
   Double_t fit_x1 = 0;
   Double_t fit_x2 = 0;
   Double_t eps = 1e-7;
   if (v[halfmax_i+2] - v[halfmax_i+1] > v[halfmax_i] - v[halfmax_i-1]) {
      fit_x1 = t[halfmax_i] - eps;
      fit_x2 = t[halfmax_i+2] + eps;
   }
   else {
      fit_x1 = t[halfmax_i-1] - eps;
      fit_x2 = t[halfmax_i+1] + eps;
   }
   //cout<< "fit_x1 = " << fit_x1 << " fit_x2 = " << fit_x2 <<endl;

   g->Fit("pol1", "Q", "goff", fit_x1,fit_x2);
   Double_t delta = 0.1;
   g->GetFunction("pol1")->SetRange(fit_x1-delta, fit_x2+delta);
   linepar[0] = g->GetFunction("pol1")->GetParameter(0);
   linepar[1] = g->GetFunction("pol1")->GetParameter(1);

   tstamp = TMath::Abs(linepar[1]) > eps? -linepar[0]/linepar[1]: t[0];
   if (debug) cout<< "tstamp = " << tstamp <<endl;
}

void trigger_run(const char* ifname, Double_t thres, Int_t event1=0, Int_t event2=-1, bool debug=false)
{
   TFile* ifile = new TFile(ifname);
   if (!ifile) {
      cout<< "Could not find file " << ifname <<endl;
      return;
   }
   TTree* wfm = (TTree*) ifile->Get("wfm");
   if (!wfm) {
      cout<< "Could not find tree 'wfm'" <<endl;
      return;
   }

   TFile* ofile = new TFile(Form("%s.time.root",ifname), "recreate");
   TTree* otree = new TTree("res", Form("Resulting tree for %s",ifname));
   otree->SetMarkerColor(46);
   otree->SetLineColor(46);

   TreeRes::book(otree);

   Float_t tFloat[20000];
   Float_t vFloat[20000];
   Double_t t[20000];
   Double_t v[20000];

   Int_t chan = 4;

   wfm->SetBranchAddress(Form("t%d",chan), &tFloat);
   wfm->SetBranchAddress(Form("v%d",chan), &vFloat);
   Int_t N = wfm->GetLeaf(Form("t%d",chan))->GetLen();

   if (event2 < event1) event2 = wfm->GetEntries() - 1;

   for (int ientry=event1; ientry<=event2; ++ientry)
   {
      if (wfm->LoadTree(ientry) < 0) break;
      wfm->GetEntry(ientry);
      if ((ientry-event1)%1000 == 0) cout<< "processing entry " << ientry <<endl;
      if (ientry-event1+1 > 10) debug = false;

      for (int i=0; i<N; ++i) {
         t[i] = tFloat[i];
         v[i] = -vFloat[i];
      }

      TreeRes::clear();
      TreeRes::npar[chan-1] = 2;
      Double_t par[10];
      trigger_fit(N, t, v, thres, par, TreeRes::tstamp[chan-1], chan, ientry, debug);
      TreeRes::par0[chan-1] = par[0];
      TreeRes::par1[chan-1] = par[1];
      otree->Fill();
   }

   cout<< "Writing " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
   otree->Write();

   wfm->ResetBranchAddresses();
}

//
//    channel 1
//

void channel1(Int_t event, Int_t chan=1, Double_t thres=10., Double_t sign=-1., Double_t cutoff_GHz=0.750, Int_t N=0)
{
   TTree* wfm = (TTree*) gDirectory->Get("wfm");
   if (!wfm) {
      cout<< "Could not find tree 'wfm'" <<endl;
      return;
   }

   Float_t tFloat[20000];
   Float_t vFloat[20000];
   Double_t t[20000];
   Double_t v[20000];

   wfm->SetBranchAddress(Form("t%d",chan), &tFloat);
   wfm->SetBranchAddress(Form("v%d",chan), &vFloat);
   Int_t npoints = wfm->GetLeaf(Form("t%d",chan))->GetLen();

   if (N == 0) N = npoints;

   if (wfm->LoadTree(event) < 0) {
      cout<< "No such entry: " << event <<endl;
      return;
   }
   wfm->GetEntry(event);

   for (int i=0; i<npoints; ++i) {
      t[i] = tFloat[i];
      v[i] = sign*vFloat[i];
   }

   TGraph* graw = new TGraph(npoints,t,v);
   graw->SetNameTitle(Form("graw_evt_%d_ch_%d",event,chan), Form("graw_evt_%d_ch_%d",event,chan));
   graw->SetMarkerStyle(7);
   graw->SetMarkerColor(602);
   graw->SetLineColor(602);

   new TCanvas;
   //g->GetXaxis()->SetRangeUser(-2,1);
   graw->Draw("ap");

   if (cutoff_GHz > 0) {
      Double_t re0 = fft_filter(N, t, v, cutoff_GHz);
      cout<< "re0 = " << re0 <<endl;
   }

   // for (int i=0; i<N; ++i) {
   //    v[i] -= re0;
   //    graw->SetPoint(i, graw->GetX()[i], graw->GetY()[i] - re0);
   // }

   TGraph* g = new TGraph(npoints,t,v);
   g->SetNameTitle(Form("g_evt_%d_ch_%d",event,chan), Form("g_evt_%d_ch_%d",event,chan));
   g->SetMarkerStyle(7);
   g->SetMarkerColor(46);
   g->SetLineColor(46);

   g->Draw("p");

   Int_t thres_i = 0;
   while (v[thres_i] < thres && thres_i < npoints) ++thres_i;
   if (thres_i == npoints) {
      cout<< "no signal above the threshold " << thres <<endl;
      return;            // no signal
   }

   /// // find the beginning of the pulse
   /// Int_t i5 = thres_i;
   /// Int_t local_minimum_i = thres_i;
   /// //while (v[i5] > 0.05*v[thres_i] && i5 > 0) --i5;
   /// for (i5=thres_i; i5>0; --i5) {
   ///    if (v[i5] < v[local_minimum_i]) local_minimum_i = i5;
   ///    if (v[i5] < 0.25*v[thres_i] && v[i5] > v[local_minimum_i]) break;
   ///    //-- if (v[i5] < 0) break;
   ///    if (v[i5] > 0.05*v[thres_i] && i5 > 0) continue;
   /// }

   /// // baseline
   /// Double_t baseline = 0;
   /// for (int i=0; i<i5; ++i) baseline += v[i];
   /// baseline /= i5;

   /// cout<< "thres_i = " << thres_i << " t[thres_i] = " << t[thres_i] << " i5 = " << i5 << " t[" << i5 << "] = " << t[i5] << " local_minimum_i = " << local_minimum_i << " t[local_minimum_i] = " << t[local_minimum_i] << " baseline = " << baseline <<endl;

   /// //for (int i=0; i<N; ++i) v[i] -= baseline;
   /// for (int i=0; i<N; ++i) {
   ///    graw->SetPoint(i, graw->GetX()[i], graw->GetY()[i] - baseline);
   ///    g->SetPoint(i, g->GetX()[i], g->GetY()[i] - baseline);
   /// }

   Int_t maximum_i = thres_i;
   for (int i=thres_i+1; i<npoints; ++i) {
      if (v[i] > v[maximum_i]) maximum_i = i;
      if (v[i] < 0.9*v[maximum_i]) break;              // stop at the first maximum
   }

   Int_t i88 = maximum_i;
   //-- while (v[i88] > 0.88*v[maximum_i] && i88 > 0) --i88;
   while (v[i88] > 0.75*v[maximum_i] && i88 > 0) --i88;

   Int_t i12 = i88;
   while (v[i12] > 0.12*v[maximum_i] && i12 > 0) --i12;

   Double_t fit_x1 = 0;
   Double_t fit_x2 = 0;
   Double_t eps = 1e-7;
   fit_x1 = t[i12] - eps;
   fit_x2 = t[i88] + eps;

   fit_x1 = -15.5;
   fit_x2 = -14;
   cout<< "fit_x1 = " << fit_x1 << " fit_x2 = " << fit_x2 <<endl;

   Double_t tau = 0.5*(t[maximum_i] - t[i12]);
   Double_t x0 = t[i12];
   Double_t A = v[maximum_i]*(tau/2.)*TMath::Exp(2);
   Double_t sigma = 0.25*tau;
   TF1* f = FunRC::RC2b(fit_x1,fit_x2, A,x0,tau);
   //TF1* f = FunRC::RC1gb(fit_x1,fit_x2, A,x0,tau,sigma);        // NB: RC1gb
   //-- TF1* f = FunRC::RC2gb(fit_x1,fit_x2, A,x0,tau,sigma);        // NB: RC2gb
   // Double_t d = tau;
   // TF1* f = FunRC::RC1flat(fit_x1,fit_x2, A,x0,tau,d);        // NB: RC1gb

   graw->Fit(f,"","",fit_x1,fit_x2);

   wfm->ResetBranchAddresses();
}

void channel1_fit(Int_t N, Double_t* t, Double_t* v, Double_t thres, Double_t cutoff_GHz, Int_t chan=1, Int_t event=-1, bool debug=false)
{
   if (debug) cout<< "channel1_fit event = " << event <<endl;

   TGraph* g = new TGraph(N,t,v);
   g->SetNameTitle(Form("g_evt_%d_ch_%d",event,chan), Form("g_evt_%d_ch_%d",event,chan));
   g->SetMarkerStyle(7);
   g->SetMarkerColor(602);
   g->SetLineColor(602);

   if (debug) {
      // cout<< "channel1_fit: plot g" <<endl;
      new TCanvas;
      // g->GetXaxis()->SetRangeUser(-2,1);
      g->Draw("ap");
   }

   if (cutoff_GHz > 0) {
      // apply low pass filter
      Double_t re0 = 0;
      if (cutoff_GHz > 0) re0 = fft_filter(N, t, v, cutoff_GHz);
      if (cutoff_GHz > 0 && debug) cout<< "re0 = " << re0 <<endl;
   }

   /// for (int i=0; i<N; ++i) {
   ///    g->SetPoint(i, g->GetX()[i], g->GetY()[i]-baseline);     // subtract from the raw plot as well
   ///    v[i] -= baseline;
   /// }

   Int_t thres_i = 0;
   while (v[thres_i] < thres && thres_i < N) ++thres_i;
   if (debug) cout<< "thres_i = " << thres_i << " t[thres_i] = " << t[thres_i] <<endl;
   if (thres_i == N) return;            // no signal

   Int_t maximum_i = thres_i;
   for (int i=thres_i+1; i<N; ++i) {
      if (v[i] > v[maximum_i]) maximum_i = i;
      if (v[i] < 0.9*v[maximum_i]) break;              // stop at the first maximum
   }
   if (debug) cout<< "maximum_i = " << maximum_i << " t[maximum_i] = " << t[maximum_i] <<endl;

   Int_t halfmax_i1 = maximum_i - 1;
   while (halfmax_i1 > 0 && v[halfmax_i1] > 0.5*v[maximum_i]) --halfmax_i1;

   Int_t halfmax_i2 = maximum_i + 1;
   while (halfmax_i2 > 0 && v[halfmax_i2] > 0.5*v[maximum_i]) ++halfmax_i2;

   if (debug) cout<< "t[halfmax_i1] = " << t[halfmax_i1] << " t[halfmax_i2] = " << halfmax_i2 <<endl;

   // signal region

   Int_t signal_i = maximum_i - 3*(maximum_i - halfmax_i1);
   if (signal_i < 0) signal_i = 0;
   if (debug) cout<< "t[signal_i] = " << t[signal_i] <<endl;

   Double_t baseline = 0;
   if (signal_i > 7) {
      // Double_t W = 0;      // weight sum
      // Double_t mean = 0;
      // for (int i=0; i<signal_i; ++i) {
      //    baseline += v[i];
      //    Double_t w = 1.;                    // equal weights
      //    mean = mean + w*(v[i] - mean)/(W+w);
      //    W = W + w;
      // }
      // baseline /= signal_i;
      // TreeRes::bkg[chan-1] = mean;
      // TreeRes::ebkg[chan-1] = mean/TMath::Sqrt(W);
      TreeRes::bkg_raw[chan-1] = TMath::Mean(signal_i, v);
      TreeRes::ebkg_raw[chan-1] = TMath::RMS(signal_i, v);
      TreeRes::bkg[chan-1] = TMath::Mean(signal_i, v);
      TreeRes::ebkg[chan-1] = TMath::RMS(signal_i, v);
      baseline = TreeRes::bkg[chan-1];
   }

   for (int i=0; i<N; ++i) {
      v[i] -= baseline;
      g->SetPoint(i, g->GetX()[i], g->GetY()[i]-baseline);
   }

   Double_t ey[20000];
   for (int i=0; i<N; ++i) ey[i] = TreeRes::ebkg[chan-1];

   //TGraph* gout = new TGraph(N,t,v);
   TGraphErrors* gout = new TGraphErrors(N,t,v,0,ey);
   gout->SetNameTitle(Form("gout_evt_%d_ch_%d",event,chan), Form("gout_evt_%d_ch_%d",event,chan));
   gout->SetMarkerStyle(7);
   gout->SetMarkerColor(46);
   gout->SetLineColor(46);

   if (debug) {
      // cout<< "channel2_fit: plot gout" <<endl;
      // new TCanvas;
      // gout->Draw("ap");
      gout->Draw("p");
   }

   // Int_t i88 = maximum_i;
   // while (v[i88] > 0.88*v[maximum_i] && i88 > 0) --i88;

   // Int_t i12 = i88;
   // while (v[i12] > 0.12*v[maximum_i] && i12 > 0) --i12;

   // if (debug) cout<< "i88 = " << i88 << " t[i88] = " << t[i88] << " i12 = " << i12 << " t[i12] = " << t[i12] <<endl;

   Int_t i80 = maximum_i;
   while (v[i80] > 0.75*v[maximum_i] && i80 > 0) --i80;

   Int_t i20 = i80;
   while (v[i20] > 0.25*v[maximum_i] && i20 > 0) --i20;

   if (debug) cout<< "i80 = " << i80 << " t[i80] = " << t[i80] << " i20 = " << i20 << " t[i20] = " << t[i20] <<endl;

   Double_t fit_x1 = 0;
   Double_t fit_x2 = 0;
   Double_t eps = 1e-7;
   // fit_x1 = t[i12] - eps;
   // fit_x2 = t[i88] + eps;
   fit_x1 = t[i20] - eps;
   fit_x2 = t[i80] + eps;
   fit_x1 -= 6.;  // ns
   if (debug) cout<< "fit_x1 = " << fit_x1 << " fit_x2 = " << fit_x2 <<endl;

   // Double_t tau = 0.5*(t[maximum_i] - t[i12]);
   // Double_t x0 = t[i12];
   // Double_t A = v[maximum_i]*(tau/2.)*TMath::Exp(2);
   Double_t tau = 0.5*(t[maximum_i] - t[i20]);
   Double_t x0 = t[i20];
   Double_t A = v[maximum_i]*(tau/2.)*TMath::Exp(2);

   TF1* f = FunRC::RC2(fit_x1,fit_x2, A,x0,tau);
   //-- Double_t sigma = 0.5;
   // TF1* f = FunRC::RC1gb(fit_x1,fit_x2, A,x0,tau,sigma);        // NB: RC1gb
   //-- TF1* f = FunRC::RC2gb(fit_x1,fit_x2, A,x0,tau,sigma);        // NB: RC2gb

   //-- Double_t d = tau;
   //-- TF1* f = FunRC::RC1flat(fit_x1,fit_x2, A,x0,tau,d);        // NB: RC1flat

   if (debug) gout->Fit(f,"","",fit_x1,fit_x2);
   else gout->Fit(f,"q","goff",fit_x1,fit_x2);

   Double_t delta = 0.1;
   TF1* fun = gout->GetFunction(f->GetName());
   if (!fun) return;

   fun->SetRange(fit_x1-delta, fit_x2+delta);
   TreeRes::npar[chan-1] = fun->GetNpar();
   if (fun->GetNpar() > 0) TreeRes::par0[chan-1] = fun->GetParameter(0);
   if (fun->GetNpar() > 1) TreeRes::par1[chan-1] = fun->GetParameter(1);
   if (fun->GetNpar() > 2) TreeRes::par2[chan-1] = fun->GetParameter(2);
   if (fun->GetNpar() > 3) TreeRes::par3[chan-1] = fun->GetParameter(3);
   if (fun->GetNpar() > 4) TreeRes::par4[chan-1] = fun->GetParameter(4);
   if (fun->GetNpar() > 5) TreeRes::par5[chan-1] = fun->GetParameter(5);

   // get baseline from the fit
   //-- Double_t bkg = fun->GetParameter("bkg");
   Double_t bkg = fun->GetParameter(4);

   TreeRes::pmax[chan-1] = v[maximum_i];
   TreeRes::ym[chan-1] = bkg + 0.5*(v[maximum_i] - bkg);
   TreeRes::xm[chan-1] = fun->GetX(TreeRes::ym[chan-1], t[i20], t[i80]);
   //TreeRes::tstamp = fun->GetX(0.5*v[maximum_i], t[i20], t[i80]);
   //-- Double_t derivative = (2./(TreeRes::xm[chan-1]-TreeRes::par[chan-1][1]) - 1/TreeRes::par[chan-1][2]) * TreeRes::ym[chan-1];
   Double_t derivative = fun->Derivative(TreeRes::xm[chan-1]);
   TreeRes::tstamp[chan-1] = derivative > eps? TreeRes::xm[chan-1] - (TreeRes::ym[chan-1] - bkg)/derivative: TreeRes::xm[chan-1];
   if (debug) cout<< "TreeRes::tstamp[chan-1] = " << TreeRes::tstamp[chan-1] << " derivative = " << derivative <<endl;

   TreeRes::tx[chan-1] = fit_x1;
   for (int i=i20; i<i80; ++i) {
      derivative = fun->Derivative(t[i]);
      Double_t tx = derivative > eps? t[i] - (fun->Eval(t[i]) - bkg)/derivative: t[i];
      if (tx > TreeRes::tx[chan-1]) TreeRes::tx[chan-1] = tx;
   }

   TreeRes::x5mV[chan-1] = fun->GetX(5+bkg, fit_x1, t[i80]);
   TreeRes::x10mV[chan-1] = fun->GetX(10+bkg, fit_x1, t[i80]);
}

void channel1_run(const char* ifname, Double_t thres=30, Double_t cutoff_GHz=0., Int_t event1=0, Int_t event2=-1, bool debug=false)
{
   TFile* ifile = new TFile(ifname);
   if (!ifile) {
      cout<< "Could not find file " << ifname <<endl;
      return;
   }
   TTree* wfm = (TTree*) ifile->Get("wfm");
   if (!wfm) {
      cout<< "Could not find tree 'wfm'" <<endl;
      return;
   }

   Int_t chan = 1;

   TFile* ofile = new TFile(Form("%s-chan%d.time.root",ifname,chan), "recreate");
   TTree* otree = new TTree("res", Form("Resulting tree for %s",ifname));
   otree->SetMarkerColor(46);
   otree->SetLineColor(46);

   TreeRes::book(otree);

   Float_t tFloat[20000];
   Float_t vFloat[20000];
   Double_t t[20000];
   Double_t v[20000];

   wfm->SetBranchAddress(Form("t%d",chan), &tFloat);
   wfm->SetBranchAddress(Form("v%d",chan), &vFloat);
   Int_t N = wfm->GetLeaf(Form("t%d",chan))->GetLen();

   if (event2 < event1) event2 = wfm->GetEntries() - 1;

   for (int ientry=event1; ientry<=event2; ++ientry)
   {
      if (wfm->LoadTree(ientry) < 0) break;
      wfm->GetEntry(ientry);

      if (false
          || (ientry-event1) < 10
          || (ientry-event1)%1000 == 0
          ) cout<< "--- processing entry " << ientry <<endl;
      if (ientry-event1+1 >= 10) debug = false;

      for (int i=0; i<N; ++i) {
         t[i] = tFloat[i];
         v[i] = -1.*vFloat[i];         // inverse
      }

      TreeRes::clear();
      channel1_fit(N, t, v, thres, cutoff_GHz, chan, ientry, debug);
      otree->Fill();
   }

   cout<< "Writing " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
   otree->Write();

   wfm->ResetBranchAddresses();
}

void channel2(Int_t event, Int_t chan=2, Double_t thres=20., Double_t sign=-1., Double_t cutoff_GHz=0.500, Int_t N=0)
{
   TTree* wfm = (TTree*) gDirectory->Get("wfm");
   if (!wfm) {
      cout<< "Could not find tree 'wfm'" <<endl;
      return;
   }

   Float_t tFloat[20000];
   Float_t vFloat[20000];
   Double_t t[20000];
   Double_t v[20000];

   wfm->SetBranchAddress(Form("t%d",chan), &tFloat);
   wfm->SetBranchAddress(Form("v%d",chan), &vFloat);
   Int_t npoints = wfm->GetLeaf(Form("t%d",chan))->GetLen();

   if (N == 0) N = npoints;

   if (wfm->LoadTree(event) < 0) {
      cout<< "No such entry: " << event <<endl;
      return;
   }
   wfm->GetEntry(event);

   for (int i=0; i<npoints; ++i) {
      t[i] = tFloat[i];
      v[i] = sign*vFloat[i];
   }

   TGraph* graw = new TGraph(npoints,t,v);
   graw->SetNameTitle(Form("graw_evt_%d_ch_%d",event,chan), Form("graw_evt_%d_ch_%d",event,chan));
   graw->SetMarkerStyle(7);
   graw->SetMarkerColor(602);
   graw->SetLineColor(602);

   new TCanvas;
   //g->GetXaxis()->SetRangeUser(-2,1);
   graw->Draw("ap");

   fft_filter(N, t, v, cutoff_GHz);

   TGraph* g = new TGraph(npoints,t,v);
   g->SetNameTitle(Form("g_evt_%d_ch_%d",event,chan), Form("g_evt_%d_ch_%d",event,chan));
   g->SetMarkerStyle(7);
   g->SetMarkerColor(46);
   g->SetLineColor(46);

   g->Draw("p");

   Int_t thres_i = 0;
   while (v[thres_i] < thres && thres_i < npoints) ++thres_i;
   if (thres_i == npoints) return;            // no signal

   Int_t maximum_i = thres_i;
   for (int i=thres_i+1; i<npoints; ++i) {
      if (v[i] > v[maximum_i]) maximum_i = i;
      if (v[i] < 0.9*v[maximum_i]) break;              // stop at the first maximum
   }

   Int_t i25 = maximum_i;
   while (v[i25] > 0.25*v[maximum_i] && i25 > 0) --i25;

   Int_t i12 = i25;
   while (v[i12] > 0.12*v[maximum_i] && i12 > 0) --i12;

   Double_t fit_x1 = 0;
   Double_t fit_x2 = 0;
   Double_t eps = 1e-7;
   fit_x1 = t[i12] - eps;
   fit_x2 = t[i25] + eps;
   cout<< "fit_x1 = " << fit_x1 << " fit_x2 = " << fit_x2 <<endl;

   Double_t tau = 0.5*(t[maximum_i] - t[i12]);
   Double_t x0 = t[i12];
   Double_t A = v[maximum_i]*(tau/2.)*TMath::Exp(2);
   Double_t sigma = 0.25*tau;
   TF1* f1g = FunRC::RC1g(fit_x1,fit_x2, A,x0,tau,sigma);

   g->Fit(f1g,"","",fit_x1,fit_x2);

   wfm->ResetBranchAddresses();
}

void channel2_fit(Int_t N, Double_t* t, Double_t* v, Double_t thres, Double_t cutoff_GHz, Int_t chan=0, Int_t event=-1, bool debug=false)
{
   if (debug) cout<< "channel2_fit event = " << event <<endl;

   TGraph* g = new TGraph(N,t,v);
   g->SetNameTitle(Form("g_evt_%d_ch_%d",event,chan), Form("g_evt_%d_ch_%d",event,chan));
   g->SetMarkerStyle(7);
   g->SetMarkerColor(602);
   g->SetLineColor(602);

   if (debug) {
      // cout<< "channel2_fit: plot g" <<endl;
      new TCanvas;
      // g->GetXaxis()->SetRangeUser(-2,1);
      g->Draw("ap");
   }

   // apply low pass filter
   Double_t re0 = 0;
   if (cutoff_GHz > 0) re0 = fft_filter(N, t, v, cutoff_GHz);
   if (cutoff_GHz > 0 && debug) cout<< "re0 = " << re0 <<endl;

   /// for (int i=0; i<N; ++i) {
   ///    g->SetPoint(i, g->GetX()[i], g->GetY()[i]-baseline);     // subtract from the raw plot as well
   ///    v[i] -= baseline;
   /// }

   Int_t thres_i = 0;
   while (v[thres_i] < thres && thres_i < N) ++thres_i;
   if (debug) cout<< "thres_i = " << thres_i << " t[thres_i] = " << t[thres_i] <<endl;
   if (thres_i == N) return;            // no signal

   Int_t maximum_i = thres_i;
   for (int i=thres_i+1; i<N; ++i) {
      if (v[i] > v[maximum_i]) maximum_i = i;
      if (v[i] < 0.9*v[maximum_i]) break;              // stop at the first maximum
   }
   if (debug) cout<< "maximum_i = " << maximum_i << " t[maximum_i] = " << t[maximum_i] <<endl;

   Int_t halfmax_i1 = maximum_i - 1;
   while (halfmax_i1 > 0 && v[halfmax_i1] > 0.5*v[maximum_i]) --halfmax_i1;

   Int_t halfmax_i2 = maximum_i + 1;
   while (halfmax_i2 > 0 && v[halfmax_i2] > 0.5*v[maximum_i]) ++halfmax_i2;

   if (debug) cout<< "t[halfmax_i1] = " << t[halfmax_i1] << " t[halfmax_i2] = " << halfmax_i2 <<endl;

   // signal region

   Int_t signal_i = maximum_i - 3*(maximum_i - halfmax_i1);
   if (signal_i < 0) signal_i = 0;
   if (debug) cout<< "t[signal_i] = " << t[signal_i] <<endl;

   Double_t baseline = 0;
   if (signal_i > 7) {
      // Double_t W = 0;      // weight sum
      // Double_t mean = 0;
      // for (int i=0; i<signal_i; ++i) {
      //    baseline += v[i];
      //    Double_t w = 1.;                    // equal weights
      //    mean = mean + w*(v[i] - mean)/(W+w);
      //    W = W + w;
      // }
      // baseline /= signal_i;
      // TreeRes::bkg[chan-1] = mean;
      // TreeRes::ebkg[chan-1] = mean/TMath::Sqrt(W);
      TreeRes::bkg_raw[chan-1] = TMath::Mean(signal_i, g->GetY());
      TreeRes::ebkg_raw[chan-1] = TMath::RMS(signal_i, g->GetY());
      TreeRes::bkg[chan-1] = TMath::Mean(signal_i, v);
      TreeRes::ebkg[chan-1] = TMath::RMS(signal_i, v);
      baseline = TreeRes::bkg[chan-1];
   }

   for (int i=0; i<N; ++i) {
      v[i] -= baseline;
      g->SetPoint(i, g->GetX()[i], g->GetY()[i]-baseline);
   }

   Double_t ey[20000];
   for (int i=0; i<N; ++i) ey[i] = TreeRes::ebkg[chan-1];

   //TGraph* gout = new TGraph(N,t,v);
   TGraphErrors* gout = new TGraphErrors(N,t,v,0,ey);
   gout->SetNameTitle(Form("gout_evt_%d_ch_%d",event,chan), Form("gout_evt_%d_ch_%d",event,chan));
   gout->SetMarkerStyle(7);
   gout->SetMarkerColor(46);
   gout->SetLineColor(46);

   if (debug) {
      // cout<< "channel2_fit: plot gout" <<endl;
      // new TCanvas;
      // gout->Draw("ap");
      gout->Draw("p");
   }

   // Int_t i88 = maximum_i;
   // while (v[i88] > 0.88*v[maximum_i] && i88 > 0) --i88;

   // Int_t i12 = i88;
   // while (v[i12] > 0.12*v[maximum_i] && i12 > 0) --i12;

   // if (debug) cout<< "i88 = " << i88 << " t[i88] = " << t[i88] << " i12 = " << i12 << " t[i12] = " << t[i12] <<endl;

   Int_t i80 = maximum_i;
   while (v[i80] > 0.75*v[maximum_i] && i80 > 0) --i80;

   Int_t i20 = i80;
   while (v[i20] > 0.25*v[maximum_i] && i20 > 0) --i20;

   if (debug) cout<< "i80 = " << i80 << " t[i80] = " << t[i80] << " i20 = " << i20 << " t[i20] = " << t[i20] <<endl;

   Double_t fit_x1 = 0;
   Double_t fit_x2 = 0;
   Double_t eps = 1e-7;
   // fit_x1 = t[i12] - eps;
   // fit_x2 = t[i88] + eps;
   fit_x1 = t[i20] - eps;
   fit_x2 = t[i80] + eps;
   //--------------------- fit_x1 -= 6.;  // ns
   fit_x1 -= 3.;  // ns
   if (debug) cout<< "fit_x1 = " << fit_x1 << " fit_x2 = " << fit_x2 <<endl;

   // Double_t tau = 0.5*(t[maximum_i] - t[i12]);
   // Double_t x0 = t[i12];
   // Double_t A = v[maximum_i]*(tau/2.)*TMath::Exp(2);
   Double_t tau = 0.5*(t[maximum_i] - t[i20]);
   Double_t x0 = t[i20];
   Double_t A = v[maximum_i]*(tau/2.)*TMath::Exp(2);

   //-- TF1* f = FunRC::RC2(fit_x1,fit_x2, A,x0,tau);
   Double_t sigma = 0.5;
   // TF1* f = FunRC::RC1gb(fit_x1,fit_x2, A,x0,tau,sigma);        // NB: RC1gb
   TF1* f = FunRC::RC2gb(fit_x1,fit_x2, A,x0,tau,sigma);        // NB: RC2gb

   //-- Double_t d = tau;
   //-- TF1* f = FunRC::RC1flat(fit_x1,fit_x2, A,x0,tau,d);        // NB: RC1flat

   if (debug) gout->Fit(f,"","",fit_x1,fit_x2);
   else gout->Fit(f,"q","goff",fit_x1,fit_x2);

   Double_t delta = 0.1;
   TF1* fun = gout->GetFunction(f->GetName());
   if (!fun) return;

   fun->SetRange(fit_x1-delta, fit_x2+delta);
   TreeRes::npar[chan-1] = fun->GetNpar();
   if (fun->GetNpar() > 0) TreeRes::par0[chan-1] = fun->GetParameter(0);
   if (fun->GetNpar() > 1) TreeRes::par1[chan-1] = fun->GetParameter(1);
   if (fun->GetNpar() > 2) TreeRes::par2[chan-1] = fun->GetParameter(2);
   if (fun->GetNpar() > 3) TreeRes::par3[chan-1] = fun->GetParameter(3);
   if (fun->GetNpar() > 4) TreeRes::par4[chan-1] = fun->GetParameter(4);
   if (fun->GetNpar() > 5) TreeRes::par5[chan-1] = fun->GetParameter(5);

   // get baseline from the fit
   //-- Double_t bkg = fun->GetParameter("bkg");
   Double_t bkg = fun->GetParameter(4);

   TreeRes::pmax[chan-1] = v[maximum_i];
   TreeRes::ym[chan-1] = bkg + 0.5*(v[maximum_i] - bkg);
   TreeRes::xm[chan-1] = fun->GetX(TreeRes::ym[chan-1], t[i20], t[i80]);
   //TreeRes::tstamp = fun->GetX(0.5*v[maximum_i], t[i20], t[i80]);
   //-- Double_t derivative = (2./(TreeRes::xm[chan-1]-TreeRes::par[chan-1][1]) - 1/TreeRes::par[chan-1][2]) * TreeRes::ym[chan-1];
   Double_t derivative = fun->Derivative(TreeRes::xm[chan-1]);
   TreeRes::tstamp[chan-1] = TreeRes::xm[chan-1] - (TreeRes::ym[chan-1] - bkg)/derivative;
   if (debug) cout<< "TreeRes::tstamp[chan-1] = " << TreeRes::tstamp[chan-1] << " derivative = " << derivative <<endl;

   TreeRes::tx[chan-1] = fit_x1;
   for (int i=i20; i<i80; ++i) {
      derivative = fun->Derivative(t[i]);
      Double_t tx = t[i] - (fun->Eval(t[i]) - bkg)/derivative;
      if (tx > TreeRes::tx[chan-1]) TreeRes::tx[chan-1] = tx;
   }

   TreeRes::x5mV[chan-1] = fun->GetX(5+bkg, fit_x1, t[i80]);
   TreeRes::x10mV[chan-1] = fun->GetX(10+bkg, fit_x1, t[i80]);
}

void channel2_run(const char* ifname, Double_t thres=30, Double_t cutoff_GHz=0.500, Int_t event1=0, Int_t event2=-1, bool debug=false)
{
   TFile* ifile = new TFile(ifname);
   if (!ifile) {
      cout<< "Could not find file " << ifname <<endl;
      return;
   }
   TTree* wfm = (TTree*) ifile->Get("wfm");
   if (!wfm) {
      cout<< "Could not find tree 'wfm'" <<endl;
      return;
   }

   Int_t chan = 2;

   TFile* ofile = new TFile(Form("%s-chan%d.time.root",ifname,chan), "recreate");
   TTree* otree = new TTree("res", Form("Resulting tree for %s",ifname));
   otree->SetMarkerColor(46);
   otree->SetLineColor(46);

   TreeRes::book(otree);

   Float_t tFloat[20000];
   Float_t vFloat[20000];
   Double_t t[20000];
   Double_t v[20000];

   wfm->SetBranchAddress(Form("t%d",chan), &tFloat);
   wfm->SetBranchAddress(Form("v%d",chan), &vFloat);
   Int_t N = wfm->GetLeaf(Form("t%d",chan))->GetLen();

   if (event2 < event1) event2 = wfm->GetEntries() - 1;

   for (int ientry=event1; ientry<=event2; ++ientry)
   {
      if (wfm->LoadTree(ientry) < 0) break;
      wfm->GetEntry(ientry);

      if (false
          || (ientry-event1) < 10
          || (ientry-event1)%1000 == 0
          ) cout<< "--- processing entry " << ientry <<endl;
      if (ientry-event1+1 >= 10) debug = false;

      for (int i=0; i<N; ++i) {
         t[i] = tFloat[i];
         v[i] = -1. * vFloat[i];         // invert
      }

      TreeRes::clear();
      channel2_fit(N, t, v, thres, cutoff_GHz, chan, ientry, debug);
      otree->Fill();
   }

   cout<< "Writing " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
   otree->Write();

   wfm->ResetBranchAddresses();
}

void channel3(Int_t event, Int_t chan=3, Double_t thres=40., Float_t sign=1., Double_t cutoff_GHz=0.2, Int_t N=0)
{
   TTree* wfm = (TTree*) gDirectory->Get("wfm");
   if (!wfm) {
      cout<< "Could not find tree 'wfm'" <<endl;
      return;
   }

   Float_t tFloat[20000];
   Float_t vFloat[20000];
   Double_t t[20000];
   Double_t v[20000];

   wfm->SetBranchAddress(Form("t%d",chan), &tFloat);
   wfm->SetBranchAddress(Form("v%d",chan), &vFloat);
   Int_t npoints = wfm->GetLeaf(Form("t%d",chan))->GetLen();

   if (N == 0) N = npoints;

   if (wfm->LoadTree(event) < 0) {
      cout<< "No such entry: " << event <<endl;
      return;
   }
   wfm->GetEntry(event);

   for (int i=0; i<npoints; ++i) {
      t[i] = tFloat[i];
      v[i] = sign*vFloat[i];
   }

   TGraph* graw = new TGraph(npoints,t,v);
   graw->SetNameTitle(Form("graw_evt_%d_ch_%d",event,chan), Form("graw_evt_%d_ch_%d",event,chan));
   graw->SetMarkerStyle(7);
   graw->SetMarkerColor(602);
   graw->SetLineColor(602);

   new TCanvas;
   //g->GetXaxis()->SetRangeUser(-2,1);
   graw->Draw("ap");

   if (cutoff_GHz > 0) fft_filter(N, t, v, cutoff_GHz);

   TGraph* g = new TGraph(npoints,t,v);
   g->SetNameTitle(Form("g_evt_%d_ch_%d",event,chan), Form("g_evt_%d_ch_%d",event,chan));
   g->SetMarkerStyle(7);
   g->SetMarkerColor(46);
   g->SetLineColor(46);

   g->Draw("p");

   Int_t thres_i = 0;
   while (v[thres_i] < thres && thres_i < npoints) ++thres_i;
   if (thres_i == npoints) return;            // no signal

   Int_t maximum_i = thres_i;
   for (int i=thres_i+1; i<npoints; ++i) {
      if (v[i] > v[maximum_i]) maximum_i = i;
      if (v[i] < 0.9*v[maximum_i]) break;              // stop at the first maximum
   }

   Int_t i88 = maximum_i;
   //-- while (v[i88] > 0.88*v[maximum_i] && i88 > 0) --i88;
   while (v[i88] > 0.75*v[maximum_i] && i88 > 0) --i88;

   Int_t i12 = i88;
   while (v[i12] > 0.12*v[maximum_i] && i12 > 0) --i12;

   Double_t fit_x1 = 0;
   Double_t fit_x2 = 0;
   Double_t eps = 1e-7;
   fit_x1 = t[i12] - eps;
   fit_x2 = t[i88] + eps;
   cout<< "fit_x1 = " << fit_x1 << " fit_x2 = " << fit_x2 <<endl;

   Double_t tau = 0.5*(t[maximum_i] - t[i12]);
   Double_t x0 = t[i12];
   Double_t A = v[maximum_i]*(tau/2.)*TMath::Exp(2);

   // fit_x1 = -10.;
   fit_x1 = -15.;
   // fit_x2 = 6.;

   Double_t sigma = 0.25*tau;
   Double_t bkg = 0;

   cout<< "A = " << A << " x0 = " << x0 << " tau = " << tau << " sigma = " << sigma << " bkg = " << bkg <<endl;

   // TF1* f1g = FunRC::RC1gb(fit_x1,fit_x2, A,x0,tau,sigma,bkg);
   // g->Fit(f1g,"","",fit_x1,fit_x2);

   TF1* f2g = FunRC::RC2gb(fit_x1,fit_x2, A,x0,tau,sigma,bkg);
   g->Fit(f2g,"","",fit_x1,fit_x2);

   // TF1* f1g = FunRC::RC0(fit_x1,fit_x2, A,x0,tau,sigma);
   // g->Fit(f1g,"","",fit_x1,fit_x2);

   // TF1* f2 = FunRC::RC2(fit_x1,fit_x2, A,x0,tau);
   // g->Fit(f2,"","",fit_x1,fit_x2);

   // TF1* f1 = FunRC::RC1(fit_x1,fit_x2, A,x0,tau);
   // g->Fit(f1,"","",fit_x1,fit_x2);

   wfm->ResetBranchAddresses();
}

void channel3_fit(Int_t N, Double_t* t, Double_t* v, Double_t thres, Double_t cutoff_GHz, Int_t chan=3, Int_t event=-1, bool debug=false)
{
   if (debug) cout<< "channel3_fit event = " << event <<endl;

   TGraph* g = new TGraph(N,t,v);
   g->SetNameTitle(Form("g_evt_%d_ch_%d",event,chan), Form("g_evt_%d_ch_%d",event,chan));
   g->SetMarkerStyle(7);
   g->SetMarkerColor(602);
   g->SetLineColor(602);

   if (debug) {
      // cout<< "channel3_fit: plot g" <<endl;
      new TCanvas;
      // g->GetXaxis()->SetRangeUser(-2,1);
      g->Draw("ap");
   }

   // apply low pass filter
   Double_t re0 = 0;
   if (cutoff_GHz > 0) re0 = fft_filter(N, t, v, cutoff_GHz);
   if (cutoff_GHz > 0 && debug) cout<< "re0 = " << re0 <<endl;

   /// for (int i=0; i<N; ++i) {
   ///    g->SetPoint(i, g->GetX()[i], g->GetY()[i]-baseline);     // subtract from the raw plot as well
   ///    v[i] -= baseline;
   /// }

   Int_t thres_i = 0;
   while (v[thres_i] < thres && thres_i < N) ++thres_i;
   if (debug) cout<< "thres_i = " << thres_i << " t[thres_i] = " << t[thres_i] <<endl;
   if (thres_i == N) return;            // no signal

   Int_t maximum_i = thres_i;
   for (int i=thres_i+1; i<N; ++i) {
      if (v[i] > v[maximum_i]) maximum_i = i;
      if (v[i] < 0.9*v[maximum_i]) break;              // stop at the first maximum
   }
   if (debug) cout<< "maximum_i = " << maximum_i << " t[maximum_i] = " << t[maximum_i] <<endl;

   Int_t halfmax_i1 = maximum_i - 1;
   while (halfmax_i1 > 0 && v[halfmax_i1] > 0.5*v[maximum_i]) --halfmax_i1;

   Int_t halfmax_i2 = maximum_i + 1;
   while (halfmax_i2 > 0 && v[halfmax_i2] > 0.5*v[maximum_i]) ++halfmax_i2;

   if (debug) cout<< "t[halfmax_i1] = " << t[halfmax_i1] << " t[halfmax_i2] = " << halfmax_i2 <<endl;

   // signal region

   Int_t signal_i = maximum_i - 3*(maximum_i - halfmax_i1);
   if (signal_i < 0) signal_i = 0;
   if (debug) cout<< "t[signal_i] = " << t[signal_i] <<endl;

   Double_t baseline = 0;
   if (signal_i > 7) {
      // Double_t W = 0;      // weight sum
      // Double_t mean = 0;
      // for (int i=0; i<signal_i; ++i) {
      //    baseline += v[i];
      //    Double_t w = 1.;                    // equal weights
      //    mean = mean + w*(v[i] - mean)/(W+w);
      //    W = W + w;
      // }
      // baseline /= signal_i;
      // TreeRes::bkg[chan-1] = mean;
      // TreeRes::ebkg[chan-1] = mean/TMath::Sqrt(W);
      TreeRes::bkg_raw[chan-1] = TMath::Mean(signal_i, v);
      TreeRes::ebkg_raw[chan-1] = TMath::RMS(signal_i, v);
      TreeRes::bkg[chan-1] = TMath::Mean(signal_i, v);
      TreeRes::ebkg[chan-1] = TMath::RMS(signal_i, v);
      baseline = TreeRes::bkg[chan-1];
   }

   for (int i=0; i<N; ++i) {
      v[i] -= baseline;
      g->SetPoint(i, g->GetX()[i], g->GetY()[i]-baseline);
   }

   Double_t ey[20000];
   for (int i=0; i<N; ++i) ey[i] = TreeRes::ebkg[chan-1];

   //TGraph* gout = new TGraph(N,t,v);
   TGraphErrors* gout = new TGraphErrors(N,t,v,0,ey);
   gout->SetNameTitle(Form("gout_evt_%d_ch_%d",event,chan), Form("gout_evt_%d_ch_%d",event,chan));
   gout->SetMarkerStyle(7);
   gout->SetMarkerColor(46);
   gout->SetLineColor(46);

   if (debug) {
      // cout<< "channel3_fit: plot gout" <<endl;
      // new TCanvas;
      // gout->Draw("ap");
      gout->Draw("p");
   }

   // Int_t i88 = maximum_i;
   // while (v[i88] > 0.88*v[maximum_i] && i88 > 0) --i88;

   // Int_t i12 = i88;
   // while (v[i12] > 0.12*v[maximum_i] && i12 > 0) --i12;

   // if (debug) cout<< "i88 = " << i88 << " t[i88] = " << t[i88] << " i12 = " << i12 << " t[i12] = " << t[i12] <<endl;

   Int_t i80 = maximum_i;
   while (v[i80] > 0.75*v[maximum_i] && i80 > 0) --i80;

   Int_t i20 = i80;
   while (v[i20] > 0.25*v[maximum_i] && i20 > 0) --i20;

   if (debug) cout<< "i80 = " << i80 << " t[i80] = " << t[i80] << " i20 = " << i20 << " t[i20] = " << t[i20] <<endl;

   Double_t fit_x1 = 0;
   Double_t fit_x2 = 0;
   Double_t eps = 1e-7;
   // fit_x1 = t[i12] - eps;
   // fit_x2 = t[i88] + eps;
   fit_x1 = t[i20] - eps;
   fit_x2 = t[i80] + eps;
   fit_x1 -= 6.;  // ns
   if (debug) cout<< "fit_x1 = " << fit_x1 << " fit_x2 = " << fit_x2 <<endl;

   // Double_t tau = 0.5*(t[maximum_i] - t[i12]);
   // Double_t x0 = t[i12];
   // Double_t A = v[maximum_i]*(tau/2.)*TMath::Exp(2);
   Double_t tau = 0.5*(t[maximum_i] - t[i20]);
   Double_t x0 = t[i20];
   Double_t A = v[maximum_i]*(tau/2.)*TMath::Exp(2);

   //-- TF1* f = FunRC::RC2(fit_x1,fit_x2, A,x0,tau);
   Double_t sigma = 0.5;
   // TF1* f = FunRC::RC1gb(fit_x1,fit_x2, A,x0,tau,sigma);        // NB: RC1gb
   TF1* f = FunRC::RC2gb(fit_x1,fit_x2, A,x0,tau,sigma);        // NB: RC2gb

   //-- Double_t d = tau;
   //-- TF1* f = FunRC::RC1flat(fit_x1,fit_x2, A,x0,tau,d);        // NB: RC1flat

   if (debug) gout->Fit(f,"","",fit_x1,fit_x2);
   else gout->Fit(f,"q","goff",fit_x1,fit_x2);

   Double_t delta = 0.1;
   TF1* fun = gout->GetFunction(f->GetName());
   if (!fun) return;

   fun->SetRange(fit_x1-delta, fit_x2+delta);
   TreeRes::npar[chan-1] = fun->GetNpar();
   if (fun->GetNpar() > 0) TreeRes::par0[chan-1] = fun->GetParameter(0);
   if (fun->GetNpar() > 1) TreeRes::par1[chan-1] = fun->GetParameter(1);
   if (fun->GetNpar() > 2) TreeRes::par2[chan-1] = fun->GetParameter(2);
   if (fun->GetNpar() > 3) TreeRes::par3[chan-1] = fun->GetParameter(3);
   if (fun->GetNpar() > 4) TreeRes::par4[chan-1] = fun->GetParameter(4);
   if (fun->GetNpar() > 5) TreeRes::par5[chan-1] = fun->GetParameter(5);

   // get baseline from the fit
   //-- Double_t bkg = fun->GetParameter("bkg");
   Double_t bkg = fun->GetParameter(4);

   TreeRes::pmax[chan-1] = v[maximum_i];
   TreeRes::ym[chan-1] = bkg + 0.5*(v[maximum_i] - bkg);
   TreeRes::xm[chan-1] = fun->GetX(TreeRes::ym[chan-1], t[i20], t[i80]);
   //TreeRes::tstamp = fun->GetX(0.5*v[maximum_i], t[i20], t[i80]);
   //-- Double_t derivative = (2./(TreeRes::xm[chan-1]-TreeRes::par[chan-1][1]) - 1/TreeRes::par[chan-1][2]) * TreeRes::ym[chan-1];
   Double_t derivative = fun->Derivative(TreeRes::xm[chan-1]);
   TreeRes::tstamp[chan-1] = TreeRes::xm[chan-1] - (TreeRes::ym[chan-1] - bkg)/derivative;
   if (debug) cout<< "TreeRes::tstamp[chan-1] = " << TreeRes::tstamp[chan-1] << " derivative = " << derivative <<endl;

   TreeRes::tx[chan-1] = fit_x1;
   for (int i=i20; i<i80; ++i) {
      derivative = fun->Derivative(t[i]);
      Double_t tx = t[i] - (fun->Eval(t[i]) - bkg)/derivative;
      if (tx > TreeRes::tx[chan-1]) TreeRes::tx[chan-1] = tx;
   }

   TreeRes::x5mV[chan-1] = fun->GetX(5+bkg, fit_x1, t[i80]);
   TreeRes::x10mV[chan-1] = fun->GetX(10+bkg, fit_x1, t[i80]);
}

void channel3_run(const char* ifname, Double_t thres=30, Double_t cutoff_GHz=0.200, Int_t event1=0, Int_t event2=-1, bool debug=false)
{
   TFile* ifile = new TFile(ifname);
   if (!ifile) {
      cout<< "Could not find file " << ifname <<endl;
      return;
   }
   TTree* wfm = (TTree*) ifile->Get("wfm");
   if (!wfm) {
      cout<< "Could not find tree 'wfm'" <<endl;
      return;
   }

   Int_t chan = 3;

   TFile* ofile = new TFile(Form("%s-chan%d.time.root",ifname,chan), "recreate");
   TTree* otree = new TTree("res", Form("Resulting tree for %s",ifname));
   otree->SetMarkerColor(46);
   otree->SetLineColor(46);

   TreeRes::book(otree);

   Float_t tFloat[20000];
   Float_t vFloat[20000];
   Double_t t[20000];
   Double_t v[20000];

   wfm->SetBranchAddress(Form("t%d",chan), &tFloat);
   wfm->SetBranchAddress(Form("v%d",chan), &vFloat);
   Int_t N = wfm->GetLeaf(Form("t%d",chan))->GetLen();

   if (event2 < event1) event2 = wfm->GetEntries() - 1;

   for (int ientry=event1; ientry<=event2; ++ientry)
   {
      if (wfm->LoadTree(ientry) < 0) break;
      wfm->GetEntry(ientry);

      if (false
          || (ientry-event1) < 10
          || (ientry-event1)%1000 == 0
          ) cout<< "--- processing entry " << ientry <<endl;
      if (ientry-event1+1 >= 10) debug = false;

      for (int i=0; i<N; ++i) {
         t[i] = tFloat[i];
         v[i] = vFloat[i];
      }

      TreeRes::clear();
      channel3_fit(N, t, v, thres, cutoff_GHz, chan, ientry, debug);
      otree->Fill();
   }

   cout<< "Writing " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
   otree->Write();

   wfm->ResetBranchAddresses();
}
