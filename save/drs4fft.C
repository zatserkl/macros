#include "DRS4Event.h"

#include <TH1.h>
#include <TVirtualFFT.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph.h>

#include <iostream>

using std::cout;     using std::endl;

class DataFFT {
public:
   TVirtualFFT *fftK;            // own TVirtualFFT object, created with option K
   TVirtualFFT *fftK_back;       // backward trasform
   Int_t N;                      // should not be changed for next use of the class object
   TH1* hm;
   TH1* hp;
   TH1* hb;
   TH1* hRe;
   TH1* hIm;
   Double_t input[100000];
   Double_t output[100000];
   Double_t reX[100000];
   Double_t imX[100000];
private:
   DataFFT(): fftK(0), fftK_back(0), N(0), hm(0), hp(0), hb(0), hRe(0), hIm(0) {}
public:
   DataFFT(Int_t np): fftK(0), fftK_back(0), N(np), hm(0), hp(0), hb(0), hRe(0), hIm(0) {
      // Make our own TVirtualFFT object (using option "K")
      // Third parameter (option) consists of 3 parts:
      // -transform type:
      //  real input/complex output in our case
      // -transform flag: 
      //  the amount of time spent in planning
      //  the transform (see TVirtualFFT class description)
      // -to create a new TVirtualFFT object (option "K") or use the global (default)
      fftK = TVirtualFFT::FFT(1, &N, "R2C M K");
      assert(fftK);

      hm = new TH1F("hm", "Magnitude", N,0,N);
      hm->SetMarkerStyle(7);
      hm->SetDirectory(0);
      hp = new TH1F("hp", "Phase", N,0,N);
      hp->SetMarkerStyle(7);
      hp->SetDirectory(0);
      hb = new TH1F("hb", "The backward trasform", N,0,N);
      hb->SetMarkerStyle(7);
      hb->SetDirectory(0);
      hRe = new TH1F("hRe", "Real part", N,0,N);
      hRe->SetMarkerStyle(7);
      hRe->SetDirectory(0);
      hIm = new TH1F("hIm", "Imaginary part", N,0,N);
      hIm->SetMarkerStyle(7);
      hIm->SetDirectory(0);
   }
   ~DataFFT() {
      if (fftK) delete fftK;
      if (fftK_back) delete fftK_back;
   }
   void TransformData(Float_t* x) {
      for (int i=0; i<N; ++i) input[i] = x[i];  // convert to Double_t
      fftK->SetPoints(input);
      fftK->Transform();
      fftK->GetPointsComplex(reX,imX);
   }
   void TransformData(Double_t* x) {
      for (int i=0; i<N; ++i) input[i] = x[i];  // convert to Double_t
      fftK->SetPoints(input);
      fftK->Transform();
      fftK->GetPointsComplex(reX,imX);
   }
   void TransformBack(Double_t* re=0, Double_t* im=0) {
      if (!fftK_back) fftK_back = TVirtualFFT::FFT(1, &N, "C2R M K");
      assert(fftK_back);
      if (!re) re = reX;
      if (!im) im = imX;
      fftK_back->SetPointsComplex(re,im);
      fftK_back->Transform();
   }
   TH1* GetBack(Double_t* re=0, Double_t* im=0) {
      if (!fftK_back) fftK_back = TVirtualFFT::FFT(1, &N, "C2R M K");
      assert(fftK_back);
      if (!re) re = reX;
      if (!im) im = imX;
      fftK_back->SetPointsComplex(re,im);
      fftK_back->Transform();
      if (hb) TH1::TransformHisto(fftK_back,hb,"Re");
      else hb = TH1::TransformHisto(fftK_back,hb,"Re");
      assert(hb);
      hb->Scale(1./N);
      return hb;
   }
   void Mag(Double_t* mag) const {
      for (int i=0; i<N/2+1; ++i) {
         mag[i] = TMath::Sqrt(reX[i]*reX[i]+imX[i]*imX[i]) / N;
      }
   }
   TH1* GetMag() {
      if (hm) TH1::TransformHisto(fftK, hm, "MAG");
      else hm = TH1::TransformHisto(fftK, hm, "MAG");
      assert(hm);
      return hm;
   }
   TH1* GetPhase() {
      if (hp) TH1::TransformHisto(fftK, hp, "PH");               //-- substitution for hsin->FFT(hp,"PH")
      else hp = TH1::TransformHisto(fftK, hp, "PH");               //-- substitution for hsin->FFT(hp,"PH")
      assert(hp);
      return hp;
   }
   TH1* GetRe() {
      if (hRe) TH1::TransformHisto(fftK, hRe, "Re");
      else hRe = TH1::TransformHisto(fftK, hRe, "Re");
      assert(hRe);
      return hRe;
   }
   TH1* GetIm() {
      if (hIm) TH1::TransformHisto(fftK, hIm, "Im");
      else hIm = TH1::TransformHisto(fftK, hIm, "Im");
      assert(hIm);
      return hIm;
   }
};

void moving_average(Int_t naver, Int_t np, const Double_t* y, Double_t* yaver)
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

void magnitude(Int_t chan, Int_t N, Int_t offset=8, Int_t event1=0, Int_t event2=-1, bool debug=false)
{
   TTree* tree = 0;
   DRS4Event* drs4Event = new DRS4Event(tree);
   cout<< "tree = " << tree <<endl;
   if (!tree) return;

   cout<< "tree->GetEntries() = " << tree->GetEntries() <<endl;

   DataFFT* dataFFT = new DataFFT(N);

   if (event2 < event1) event2 = tree->GetEntries()-1;
   if (event2 - event1 + 1 > 100) {debug = false; cout<< "--> turn off debug" << endl;}

   Double_t mag[10000];
   for (int i=0; i<10000; ++i) mag[i] = 0;
   Int_t nevents = 0;

   for (int ientry=event1; ientry<=event2; ++ientry) {
      if (false
          || ientry - event1 < 10 
          || (ientry - event1) % 1000 == 0
      ) cout<< "processing entry " << ientry <<endl;
      drs4Event->GetEntry(ientry);

      Float_t* x = drs4Event->chanT(chan);
      Float_t* y = drs4Event->chanV(chan);
      Double_t xuni[10000];
      Double_t yuni[10000];
      drs4Event->uniform(1024,x,y, xuni,yuni);

      if (debug) {
         TGraph* guni = new TGraph(1024,xuni,yuni);
         guni->SetNameTitle(Form("guni_evt_%d_ch_%d_N_%d",ientry,chan,N), Form("Event %d ch %d N %d",ientry,chan,N));
         guni->SetMarkerStyle(7);
         guni->SetMarkerColor(2);
         guni->SetLineColor(2);
         new TCanvas;
         guni->Draw("apl");
      }

      Double_t yaver[10000];
      moving_average(9, 1024,yuni, yaver);

      Double_t input[10000];
      for (int i=0; i<N; ++i) {
         input[i] = -1.*yuni[offset+i];
         // input[i] = -1.*yaver[offset+i];
      }
      dataFFT->TransformData(input);

      Double_t mag_event[10000];
      dataFFT->Mag(mag_event);
      for (int i=0; i<=N/2; ++i) mag[i] += mag_event[i];
      ++nevents;
   }

   for (int i=0; i<=N/2; ++i) mag[i] /= nevents;

   Double_t samples[10000];
   for (int i=0; i<10000; ++i) samples[i] = i/Double_t(N);

   TGraph* gmag = new TGraph(N/2+1, samples, mag);
   gmag->SetNameTitle("gmag", Form("Magnitude average on %d events, N = %d offset = %d",nevents,N,offset));
   gmag->SetMarkerStyle(7);
   gmag->SetMarkerColor(2);
   gmag->SetLineColor(2);

   new TCanvas;
   gmag->Draw("apl");
}

void sigbkg_study()
{
   // sampling rate
   Double_t fsampling = 5e3;            // MS/s
   Int_t Nsampling = 1024;
   Double_t fres = fsampling/Nsampling;
   Double_t tres = 1./fres;            // 200 ps at 5 GS/s

   //-- Double_t noise1_T = 2.8;            // ns
   // Double_t noise1_T = 2.8e3;            // ns
   // Double_t noise1_T = 2.8e2;            // ns
   //-- Double_t noise1_T = 100;            // ns
   Double_t noise1_T = 40;            // ns
   Double_t noise1_f = 1/noise1_T;     // MHz
   Double_t noise1_omega = 2.*TMath::Pi()*noise1_f;

   // background
   TF1* fnoise1 = new TF1("fnoise1", Form("sin(%lg*x)",noise1_omega), 0,200);
   fnoise1->SetNpx(4000);

   //new TCanvas;
   //fnoise1->Draw();

   // signal
   TF1* fsig = new TF1("fsig", "[0]*TMath::Gaus(x,[1],[2])", 0,200);
   fsig->SetNpx(4000);
   Double_t sig_ampl = 10;
   Double_t sig_mean = 100;
   Double_t sig_sigma = 1;
   fsig->SetParameters(sig_ampl,sig_mean,sig_sigma);

   new TCanvas;
   fsig->Draw();

   // data arrays
   Double_t x[10000];
   Double_t y[10000];

   // Fill data arrays
   for (int i=0; i<Nsampling; ++i) {
      Double_t time = i*tres;
      x[i] = time;
      y[i] = 0;
      //-- y[i] += fsig->Eval(time);
      y[i] += fnoise1->Eval(time);
      // cout<< "x[" << i << "] = " << x[i] << "\t y[i] = " << y[i] <<endl;
   }

   TGraph* g = new TGraph(Nsampling,x,y);
   g->SetNameTitle("g","data");
   g->SetMarkerStyle(7);
   g->SetMarkerColor(2);
   g->SetLineColor(2);

   new TCanvas;
   g->Draw("apl");

   //
   // the length of the DFT
   //
   //-- Int_t N = 64;
   //-- Int_t N = 128;
   Int_t N = 256;
   // Int_t N = 512;
   // Int_t N = 1024;
   DataFFT* dataFFT = new DataFFT(N);

   Int_t offset = 0;

   Double_t input[10000];
   for (int i=0; i<N; ++i) {
      input[i] = y[offset+i];
   }
   dataFFT->TransformData(input);

   Double_t mag[10000];
   dataFFT->Mag(mag);

   Double_t samples[10000];
   for (int i=0; i<10000; ++i) samples[i] = fsampling*i/Double_t(N);

   TGraph* gmag = new TGraph(N/2+1, samples, mag);
   gmag->SetNameTitle("gmag", Form("Magnitude N = %d offset = %d",N,offset));
   gmag->SetMarkerStyle(7);
   gmag->SetMarkerColor(2);
   gmag->SetLineColor(2);

   new TCanvas;
   gmag->Draw("apl");

   // Re and Im parts

   //-- TGraph* gre = new TGraph(N/2+1, samples, dataFFT->reX);
   TGraph* gre = new TGraph(N, samples, dataFFT->reX);
   gre->SetNameTitle("gre","reX");
   gre->SetMarkerStyle(7);
   gre->SetMarkerColor(2);
   gre->SetLineColor(2);

   new TCanvas;
   gre->Draw("apl");

   //-- TGraph* gim = new TGraph(N/2+1, samples, dataFFT->imX);
   TGraph* gim = new TGraph(N, samples, dataFFT->imX);
   gim->SetNameTitle("gim","imX");
   gim->SetMarkerStyle(7);
   gim->SetMarkerColor(4);
   gim->SetLineColor(4);

   new TCanvas;
   gim->Draw("apl");

   // Re and Im histos

   new TCanvas;
   dataFFT->GetRe()->Draw();

   new TCanvas;
   dataFFT->GetIm()->Draw();

   // do backward transform by hand



   // // build magnitude and phase out of Re and Im
   //
   // Double_t mag[10000], phase[10000];
   // for (int i=2; i<N/2+1; ++i)               // loop over N/2+1 elements
   // {
   //    mag[i] = TMath::Sqrt(dataFFT->reX[i]*dataFFT->reX[i]+dataFFT->imX[i]*dataFFT->imX[i]);
   //    phase[i] = TMath::ATan2(dataFFT->imX[i], dataFFT->reX[i]);
   // }

   // // rebuild Re and Im out of magnitude and phase
   //
   // Double_t reMagPhase[10000], imMagPhase[10000];
   // for (int i=0; i<N/2+1; ++i)               // loop over N/2+1 elements
   // {
   //    reMagPhase[i] = mag[i]*TMath::Cos(phase[i]);
   //    imMagPhase[i] = mag[i]*TMath::Sin(phase[i]);
   // }

   // backward transform

   new TCanvas;
   dataFFT->GetBack()->Draw();
}

void sigbkg_functions()
{
   Double_t fsampling = 5e3;            // MS/s
   Int_t Nsampling = 1024;
   // Int_t Nsampling = 8192;
   Double_t fres = fsampling/Nsampling;
   Double_t tres = 1./fres;            // 200 ps at 5 GS/s

   Double_t noise1_T = 2.8;            // ns
   // Double_t noise1_T = 2.8e2;            // ns
   // Double_t noise1_T = 2.8e3;            // ns
   Double_t noise1_f = 1/noise1_T;     // MHz
   Double_t noise1_omega = 2.*TMath::Pi()*noise1_f;

   // background
   TF1* fnoise1 = new TF1("fnoise1", Form("sin(%lg*x)",noise1_omega), 0,200);
   fnoise1->SetNpx(4000);

   //new TCanvas;
   //fnoise1->Draw();

   // signal
   TF1* fsig = new TF1("fsig", "[0]*TMath::Gaus(x,[1],[2])", 0,200);
   fsig->SetNpx(4000);
   Double_t sig_ampl = 10;
   Double_t sig_mean = 100;
   Double_t sig_sigma = 1;
   fsig->SetParameters(sig_ampl,sig_mean,sig_sigma);

   new TCanvas;
   fsig->Draw();

   // data arrays
   Double_t x[10000];
   Double_t y[10000];

   // Fill data arrays
   for (int i=0; i<Nsampling; ++i) {
      Double_t time = i*tres;
      x[i] = time;
      y[i] = 0;
      y[i] += fsig->Eval(time);           // comment out to see noise only
      y[i] += fnoise1->Eval(time);
   }

   TGraph* g = new TGraph(Nsampling,x,y);
   g->SetNameTitle("g","data");
   g->SetMarkerStyle(7);
   g->SetMarkerColor(2);
   g->SetLineColor(2);

   new TCanvas;
   g->Draw("apl");

   // Int_t N = 128;
   Int_t N = 256;
   // Int_t N = 512;
   // Int_t N = 1024;
   DataFFT* dataFFT = new DataFFT(N);

   Int_t offset = 475;

   Double_t input[10000];
   for (int i=0; i<N; ++i) {
      input[i] = y[offset+i];
   }
   dataFFT->TransformData(input);

   Double_t mag[10000];
   dataFFT->Mag(mag);

   Double_t samples[10000];
   for (int i=0; i<10000; ++i) samples[i] = fsampling*i/Double_t(N);

   TGraph* gmag = new TGraph(N/2+1, samples, mag);
   gmag->SetNameTitle("gmag", Form("Magnitude N = %d offset = %d",N,offset));
   gmag->SetMarkerStyle(7);
   gmag->SetMarkerColor(2);
   gmag->SetLineColor(2);

   new TCanvas;
   gmag->Draw("apl");

   // Re and Im parts

   TGraph* gre = new TGraph(N/2+1, samples, dataFFT->reX);
   gre->SetNameTitle("gre","reX");
   gre->SetMarkerStyle(7);
   gre->SetMarkerColor(2);
   gre->SetLineColor(2);

   new TCanvas;
   gre->Draw("apl");

   TGraph* gim = new TGraph(N/2+1, samples, dataFFT->imX);
   gim->SetNameTitle("gim","imX");
   gim->SetMarkerStyle(7);
   gim->SetMarkerColor(4);
   gim->SetLineColor(4);

   new TCanvas;
   gim->Draw("apl");
}

void sigbkg()
{
   Double_t fsampling = 5e3;            // MS/s
   Int_t Nsampling = 1024;
   Double_t fres = fsampling/Nsampling;
   Double_t tres = 1./fres;            // 200 ps at 5 GS/s

   Double_t noise1_T = 2.8;            // ns
   Double_t noise1_f = 1/noise1_T;     // MHz
   Double_t noise1_omega = 2.*TMath::Pi()*noise1_f;

   // background
   TF1* fnoise1 = new TF1("fnoise1", Form("sin(%lg*x)",noise1_omega), 0,200);
   fnoise1->SetNpx(4000);

   //new TCanvas;
   //fnoise1->Draw();

   // signal
   TF1* fsig = new TF1("fsig", "[0]*TMath::Gaus(x,[1],[2])", 0,200);
   fsig->SetNpx(4000);
   Double_t sig_ampl = 10;
   Double_t sig_mean = 100;
   Double_t sig_sigma = 1;
   fsig->SetParameters(sig_ampl,sig_mean,sig_sigma);

   new TCanvas;
   fsig->Draw();

   // data arrays
   Double_t x[10000];
   Double_t y[10000];

   // Fill data arrays
   for (int i=0; i<Nsampling; ++i) {
      Double_t time = i*tres;
      x[i] = time;
      y[i] = 0;
      y[i] += fsig->Eval(time);
      y[i] += fnoise1->Eval(time);
   }

   TGraph* g = new TGraph(Nsampling,x,y);
   g->SetNameTitle("g","data");
   g->SetMarkerStyle(7);
   g->SetMarkerColor(2);
   g->SetLineColor(2);

   new TCanvas;
   g->Draw("apl");

   // Int_t N = 128;
   Int_t N = 256;
   // Int_t N = 512;
   // Int_t N = 1024;
   DataFFT* dataFFT = new DataFFT(N);

   Double_t input[10000];

   // background only

   Int_t offsetBkg = 220;
   // Int_t offsetBkg = 210;
   // Int_t offsetBkg = 200;
   //--------------- offsetBkg = 224;
   // Int_t offsetBkg = 232;                       // really bad because covers part of the signal
   // Int_t offsetBkg = 0;

   for (int i=0; i<N; ++i) {
      input[i] = y[offsetBkg+i];
   }
   dataFFT->TransformData(input);

   Double_t reXbkg[10000], imXbkg[10000];
   for (int i=0; i<N; ++i) {
      reXbkg[i] = dataFFT->reX[i];
      imXbkg[i] = dataFFT->imX[i];
   }

   Double_t magBkg[10000];
   //Double_t magBkg[10000], phaseBkg[10000];
   for (int i=0; i<N; ++i) {
      magBkg[i] = TMath::Sqrt(dataFFT->reX[i]*dataFFT->reX[i]+dataFFT->imX[i]*dataFFT->imX[i]);
      //phaseBkg[i] = TMath::ATan2(dataFFT->imX[i], dataFFT->reX[i]);
   }

   // background + signal

   Int_t offsetSig = 260;
   // Int_t offsetSig = 264;
   // Int_t offsetSig = 268;
   // Int_t offsetSig = 272;
   // Int_t offsetSig = 276;
   // Int_t offsetSig = 280;              // offset more than the signal width, but result almost the same

   // Int_t offsetSig = offsetBkg;

   for (int i=0; i<N; ++i) {
      input[i] = y[offsetSig+i];
   }
   dataFFT->TransformData(input);

   Double_t reXsig[10000], imXsig[10000];
   for (int i=0; i<N; ++i) {
      reXsig[i] = dataFFT->reX[i];
      imXsig[i] = dataFFT->imX[i];
   }

   // build magnitude and phase out of Re and Im

   Double_t magSig[10000], phaseSig[10000];
   for (int i=0; i<N; ++i) {
      magSig[i] = TMath::Sqrt(dataFFT->reX[i]*dataFFT->reX[i]+dataFFT->imX[i]*dataFFT->imX[i]);
      phaseSig[i] = TMath::ATan2(dataFFT->imX[i], dataFFT->reX[i]);
   }

   // subtract background

   Double_t reX[10000], imX[10000];
   for (int i=0; i<N; ++i) {
      reX[i] = reXsig[i] - reXbkg[i];
      imX[i] = imXsig[i] - imXbkg[i];
   }

   Double_t magSigBkg[10000];
   for (int i=0; i<N; ++i) {
      magSigBkg[i] = magSig[i] - magBkg[i];
   }
   Double_t reSigBkg[10000], imSigBkg[10000];
   for (int i=0; i<N; ++i) {
      // rebuild Re and Im out of magnitude and phase
      reSigBkg[i] = magSigBkg[i]*TMath::Cos(phaseSig[i]);
      imSigBkg[i] = magSigBkg[i]*TMath::Sin(phaseSig[i]);
   }

   // backward transorm

   new TCanvas;
   //-- dataFFT->GetBack(reX,imX)->Draw();
   dataFFT->GetBack(reSigBkg,imSigBkg)->Draw();

   /// Double_t samples[10000];
   /// for (int i=0; i<10000; ++i) samples[i] = fsampling*i/Double_t(N);

   /// TGraph* gmag = new TGraph(N/2+1, samples, mag);
   /// gmag->SetNameTitle("gmag", Form("Magnitude N = %d offset = %d",N,offset));
   /// gmag->SetMarkerStyle(7);
   /// gmag->SetMarkerColor(2);
   /// gmag->SetLineColor(2);

   /// new TCanvas;
   /// gmag->Draw("apl");

   /// // Re and Im parts

   /// TGraph* gre = new TGraph(N/2+1, samples, dataFFT->reX);
   /// gre->SetNameTitle("gre","reX");
   /// gre->SetMarkerStyle(7);
   /// gre->SetMarkerColor(2);
   /// gre->SetLineColor(2);

   /// new TCanvas;
   /// gre->Draw("apl");

   /// TGraph* gim = new TGraph(N/2+1, samples, dataFFT->imX);
   /// gim->SetNameTitle("gim","imX");
   /// gim->SetMarkerStyle(7);
   /// gim->SetMarkerColor(4);
   /// gim->SetLineColor(4);

   /// new TCanvas;
   /// gim->Draw("apl");
}

void sigbkg_re_im()
{
   Double_t fsampling = 5e3;            // MS/s
   Int_t Nsampling = 1024;
   Double_t fres = fsampling/Nsampling;
   Double_t tres = 1./fres;            // 200 ps at 5 GS/s

   Double_t noise1_T = 2.8;            // ns
   Double_t noise1_f = 1/noise1_T;     // MHz
   Double_t noise1_omega = 2.*TMath::Pi()*noise1_f;

   // background
   TF1* fnoise1 = new TF1("fnoise1", Form("sin(%lg*x)",noise1_omega), 0,200);
   fnoise1->SetNpx(4000);

   //new TCanvas;
   //fnoise1->Draw();

   // signal
   TF1* fsig = new TF1("fsig", "[0]*TMath::Gaus(x,[1],[2])", 0,200);
   fsig->SetNpx(4000);
   Double_t sig_ampl = 10;
   Double_t sig_mean = 100;
   Double_t sig_sigma = 1;
   fsig->SetParameters(sig_ampl,sig_mean,sig_sigma);

   new TCanvas;
   fsig->Draw();

   // data arrays
   Double_t x[10000];
   Double_t y[10000];

   // Fill data arrays
   for (int i=0; i<Nsampling; ++i) {
      Double_t time = i*tres;
      x[i] = time;
      y[i] = 0;
      y[i] += fsig->Eval(time);
      y[i] += fnoise1->Eval(time);
   }

   TGraph* g = new TGraph(Nsampling,x,y);
   g->SetNameTitle("g","data");
   g->SetMarkerStyle(7);
   g->SetMarkerColor(2);
   g->SetLineColor(2);

   new TCanvas;
   g->Draw("apl");

   // Int_t N = 128;
   Int_t N = 256;
   // Int_t N = 512;
   // Int_t N = 1024;
   DataFFT* dataFFT = new DataFFT(N);

   Int_t offset = 0;
   Double_t input[10000];

   // background only

   //-- offset = 220;
   offset = 210;

   for (int i=0; i<N; ++i) {
      input[i] = y[offset+i];
   }
   dataFFT->TransformData(input);

   Double_t reXbkg[10000], imXbkg[10000];
   for (int i=0; i<N; ++i) {
      reXbkg[i] = dataFFT->reX[i];
      imXbkg[i] = dataFFT->imX[i];
   }

   // background + signal

   offset = 260;

   for (int i=0; i<N; ++i) {
      input[i] = y[offset+i];
   }
   dataFFT->TransformData(input);

   Double_t reXsig[10000], imXsig[10000];
   for (int i=0; i<N; ++i) {
      reXsig[i] = dataFFT->reX[i];
      imXsig[i] = dataFFT->imX[i];
   }

   // ************* subtract background: does not work well, depends on the offset ***************

   Double_t reX[10000], imX[10000];
   for (int i=0; i<N; ++i) {
      reX[i] = reXsig[i] - reXbkg[i];
      imX[i] = imXsig[i] - imXbkg[i];
   }

   // backward transorm

   new TCanvas;
   dataFFT->GetBack(reX,imX)->Draw();

   /// Double_t samples[10000];
   /// for (int i=0; i<10000; ++i) samples[i] = fsampling*i/Double_t(N);

   /// TGraph* gmag = new TGraph(N/2+1, samples, mag);
   /// gmag->SetNameTitle("gmag", Form("Magnitude N = %d offset = %d",N,offset));
   /// gmag->SetMarkerStyle(7);
   /// gmag->SetMarkerColor(2);
   /// gmag->SetLineColor(2);

   /// new TCanvas;
   /// gmag->Draw("apl");

   /// // Re and Im parts

   /// TGraph* gre = new TGraph(N/2+1, samples, dataFFT->reX);
   /// gre->SetNameTitle("gre","reX");
   /// gre->SetMarkerStyle(7);
   /// gre->SetMarkerColor(2);
   /// gre->SetLineColor(2);

   /// new TCanvas;
   /// gre->Draw("apl");

   /// TGraph* gim = new TGraph(N/2+1, samples, dataFFT->imX);
   /// gim->SetNameTitle("gim","imX");
   /// gim->SetMarkerStyle(7);
   /// gim->SetMarkerColor(4);
   /// gim->SetLineColor(4);

   /// new TCanvas;
   /// gim->Draw("apl");
}

void analyse(Int_t event, Int_t chan, TTree* tree=0)
{
   DRS4Event drs4Event(tree);
   drs4Event.plot(event, chan);

   drs4Event.GetEntry(event);

   Float_t* x = drs4Event.chanT(chan);
   Float_t* y = drs4Event.chanV(chan);

   Double_t xuni[10000];
   Double_t yuni[10000];
   // Int_t N=1024;
   Int_t N=128;
   // Int_t N=64;

   drs4Event.uniform(1024,x,y, xuni,yuni);

   Double_t sample[10000];
   for (int i=0; i<1024; ++i) sample[i] = i;

   Int_t offset = 8+N;
   // Int_t offset = 8;
   Double_t input[10000];
   for (int i=0; i<N; ++i) {
      input[i] = -1.*yuni[offset+i];
   }

   DataFFT* dataFFT = new DataFFT(N);
   // dataFFT->TransformData(yuni);
   dataFFT->TransformData(input);

   //Use the following method to get just one point of the output
   Double_t re, im;
   dataFFT->fftK->GetPointComplex(0, re, im);
   cout<< "--> constant level is " << re/N <<endl;

   new TCanvas;
   dataFFT->GetMag()->Draw();

   new TCanvas;
   dataFFT->GetPhase()->Draw();

   new TCanvas;
   dataFFT->GetRe()->Draw();

   new TCanvas;
   dataFFT->GetIm()->Draw();

   new TCanvas;
   dataFFT->GetBack()->Draw("pl");

   TGraph* gre = new TGraph(N/2+1, sample, dataFFT->reX);
   gre->SetNameTitle("gre","reX");
   gre->SetMarkerStyle(7);
   gre->SetMarkerColor(2);
   gre->SetLineColor(2);

   new TCanvas;
   gre->Draw("apl");

   TGraph* gim = new TGraph(N/2+1, sample, dataFFT->imX);
   gim->SetNameTitle("gim","imX");
   gim->SetMarkerStyle(7);
   gim->SetMarkerColor(4);
   gim->SetLineColor(4);

   new TCanvas;
   gim->Draw("apl");

   Double_t mag[10000];
   for (int i=0; i<N/2+1; ++i) {
      mag[i] = TMath::Sqrt(dataFFT->reX[i]*dataFFT->reX[i]+dataFFT->imX[i]*dataFFT->imX[i]);
      //cout<< mag[i] << "\t " << dataFFT->GetMag()->GetBinContent(i+1) <<endl;
   }
   delete dataFFT;
}
