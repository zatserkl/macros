// Andriy Zatserklyaniy <zatserkl@fnal.gov> Jan 7, 2015. See http://mathworld.wolfram.com/Erf.html

#include <TROOT.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>

#include <iostream>
#include <iomanip>

using std::cout;     using std::endl;

/*
Double_t x = 3; erfa(&x,0) - TMath::Erf(x)
*/
Double_t erfa(Double_t* xx, Double_t*)
{
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

   bool debug = true;
   // bool debug = false;

   Double_t x = *xx;

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

   if (debug) cout<< "sum/TMath::Sqrt(TMath::Pi()) = " << sum/TMath::Sqrt(TMath::Pi()) <<endl;

   Double_t erfa = sign_function * (1. - TMath::Exp(-x*x)/TMath::Sqrt(TMath::Pi()) * sum);
   return erfa;
}

Double_t erfa_sum(Double_t x)
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

   bool debug = true;
   // bool debug = false;

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

