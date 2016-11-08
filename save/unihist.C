#include "DRS4Event.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TH1.h>
#include <TMath.h>
#include <TCanvas.h>

#include <iostream>
#include <string>

using std::cout;     using std::endl;

void uniform(Int_t N, const Float_t* x, const Float_t* y, Float_t* xuni, Float_t* yuni)
{
   // interpolate to the tempered (униформ) scale

   Float_t step = (x[N-1] - x[0]) / (N-1);

   xuni[0] = x[0];
   yuni[0] = y[0];

   Int_t idrs = 1;
   Int_t iuni = 1;
   while (iuni < N) {
      xuni[iuni] = xuni[iuni-1] + step;
      while (x[idrs] < xuni[iuni] && idrs < N-1) idrs++;

      Float_t dx_drs = x[idrs] - x[idrs-1];
      Float_t dy_drs = y[idrs] - y[idrs-1];
      Float_t d = xuni[iuni] - x[idrs-1];
      yuni[iuni] = y[idrs-1] + d*dy_drs/dx_drs;
      ++iuni;
   }
}

TH1F* unihist(Int_t N, const Float_t* x, const Float_t* y)
{
   // interpolate to the tempered (униформ) scale

   Float_t step = (x[N-1] - x[0]) / (N-1);
   TH1F* h = new TH1F("h_unihist","unihist", N, x[0]-0.5*step, x[N-1]+0.5*step);
   h->SetDirectory(0);

   Int_t idrs = 1;
   Int_t iuni = 1;
   while (iuni <= N) {
      while (x[idrs] < h->GetBinCenter(iuni) && idrs <= N-1) idrs++;

      Float_t dx_drs = x[idrs] - x[idrs-1];
      Float_t dy_drs = y[idrs] - y[idrs-1];
      Float_t d = h->GetBinCenter(iuni) - x[idrs-1];
      h->SetBinContent(iuni, y[idrs-1] + d*dy_drs/dx_drs);
      //cout<< "iuni = " << iuni << " h->GetBinCenter(iuni) = " << h->GetBinCenter(iuni) << " h->GetBinContent(iuni) = " << h->GetBinContent(iuni) <<endl;
      ++iuni;
   }
   return h;
}

TH1F* unihist(Int_t event, Int_t chan, TTree* tree=0)   // channel: 0..11
{
   DRS4Event drs4Event(tree);

   if (drs4Event.tree->LoadTree(event) < 0) {
      cout<< "Event out of range" <<endl;
      return 0;
   }
   drs4Event.tree->GetEntry(event);

   Float_t* x = drs4Event.chanT(chan);
   Float_t* y = drs4Event.chanV(chan);
   for (int i=0; i<1024; ++i) y[i] *= -1.;

   // original graph

   TGraph* gr = new TGraph(1024, x, y);
   gr->SetNameTitle(Form("gr_evt_%d_ch_%d",event,chan+1), Form("event %d channel %d",event,chan+1));
   gr->SetMarkerStyle(7);
   gr->SetMarkerColor(46);
   gr->SetLineColor(46);

   new TCanvas;
   gr->Draw("alp");

   // show DRS4's dx

   Float_t dx[1024];

   for (int i=0; i<1024-1; ++i) {
      dx[i] = x[i+1] - x[i];
   }
   dx[1024-1] = dx[1024-2];

   TGraph* gdx = new TGraph(1024, x, dx);
   gdx->SetNameTitle(Form("gdx_evt_%d_ch_%d",event,chan+1), Form("dx for event %d channel %d",event,chan+1));
   gdx->SetMarkerStyle(7);
   gdx->SetMarkerColor(8);
   gdx->SetLineColor(8);

   new TCanvas;
   gdx->Draw("alp");

   // interpolate to the tempered (униформ) scale

   TH1F* h = unihist(1024, x, y);
   new TCanvas;
   h->Draw();

   return h;
}
