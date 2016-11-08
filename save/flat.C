#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>

#include <iostream>

using std::cout;     using std::endl;

Double_t fflat(Double_t* xx, Double_t* par) {
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

TF1* flat(Double_t xmin, Double_t xmax
          , Double_t A
          , Double_t x0
          , Double_t tau
          , Double_t d
          , Double_t bkg=0
          )
{
   Double_t par[10];
   Int_t npar = 0;

   par[npar++] = A;
   par[npar++] = x0;
   par[npar++] = tau;
   par[npar++] = d;
   par[npar++] = bkg;

   TF1* flat = new TF1("flat", fflat, xmin, xmax, npar);
   flat->SetNpx(1000);
   flat->SetParameters(par);

   npar = 0;
   flat->SetParName(npar++, "A");
   flat->SetParName(npar++, "x0");
   flat->SetParName(npar++, "tau");
   flat->SetParName(npar++, "d");
   flat->SetParName(npar++, "bkg");

   new TCanvas;
   flat->Draw();

   return flat;
}
