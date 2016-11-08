#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
using std::cout;        using std::endl;

// utils.C stuff
void cd(const char* filename="");                        // provides default parameter
TF1* ftemp(TCanvas* can=0);

class FitFunctions: public TNamed {
   // container for the fit functions
   // Eval proper function with correspondent range
public:
   std::vector<TF1> functions;
   FitFunctions(): TNamed() {this->SetName("objFitFunctions");}
   FitFunctions(const char* fname): TNamed() {
      this->SetName("objFitFunctions");
      TFile* file = TFile::Open(fname);
      if (!file) {
         cout<< "FitFunctions: Could not open file " << fname <<endl;
         return;
      }
   }
   void AddFunction(const TF1* f) {
      functions.push_back(*f);
   }
   Double_t Eval(Double_t x) {
      for (std::vector<TF1>::const_iterator it=functions.begin(); it<functions.end(); ++it) {
         if (x > it->GetXmin() && x <= it->GetXmax()) return it->Eval(x);
      }
      return 0;
   }

   ClassDef(FitFunctions, 2);
};

ClassImp(FitFunctions);

void wetfit(TGraph* gwet_sum, bool save=false)
{
   TFile* ofile = 0;
   FitFunctions* fitFunctions = new FitFunctions();

   if (save) ofile = TFile::Open("wetfit.root", "recreate");
   TF1* f;
   gwet_sum->Fit("pol2", "", "", 2000, 21000);
   f = ftemp();
   f->SetName("pol2_2000_21000");
   f->SetRange(2000, 21000);
   f->SetLineColor(2);
   if (save) {
      f->Write();
      fitFunctions->AddFunction(f);
   }
   gwet_sum->Fit("pol2", "+", "", 21000, 28650);
   f = ftemp();
   f->SetName("pol2_21000_28650");
   f->SetRange(21000, 28650);
   f->SetLineColor(4);
   if (save) {
      f->Write();
      fitFunctions->AddFunction(f);
   }
   gwet_sum->Fit("pol1", "+", "", 28650, 28710);
   f = ftemp();
   f->SetName("pol1_28650_28710");
   f->SetRange(28650, 28710);
   f->SetLineColor(1);
   if (save) {
      f->Write();
      fitFunctions->AddFunction(f);
   }
   gwet_sum->Fit("pol2", "+", "", 28710, 32500);
   f = ftemp();
   f->SetName("pol2_28710_32500");
   f->SetRange(28710, 32500);
   f->SetLineColor(8);
   if (save) {
      f->Write();
      fitFunctions->AddFunction(f);
   }
   gwet_sum->Fit("pol2", "+", "", 32000, 37500);
   f = ftemp();
   f->SetName("pol2_32500_37500");
   f->SetRange(32500, 37500);
   if (save) {
      f->Write();
      fitFunctions->AddFunction(f);
   }
   if (save) {
      cout<< "Wrote fit functions fitFunctions->GetName() = " << fitFunctions->GetName() << " into file " << ofile->GetName() <<endl;
      fitFunctions->Write(fitFunctions->GetName());
      ofile->Write();
   }
}

void scisum(Double_t thres=100, bool save=false)
{
   Double_t wet[100];
   Double_t sum[100];
   Int_t np = 0;

   const char* lstname = "t.lst";
   ifstream lst(lstname);
   if (!lst) {
      cout<< "Could not open file " << lstname <<endl;
      return;
   }

   TTree* T = 0;

   // TH1F* h_sum = new TH1F("h_sum", "h_sum", 10000,0,100000);
   TH1F* h_sum = 0;

   std::string filename;
   while (lst >> filename) {
      cd(filename.c_str());
      T = (TTree*) gDirectory->Get("T");
      Int_t iwet = 0;
      sscanf(filename.c_str(), "t-%d.root", &iwet);
      h_sum = new TH1F(Form("h_sum_%d",iwet), Form("h_sum for %s",filename.c_str()), 10000,0,100000);
      for (int ientry=0; ientry<T->GetEntries(); ++ientry)
      // for (int ientry=0; ientry<100000; ++ientry)
      {
         if (T->LoadTree(ientry) < 0) break;
         T->GetEntry(ientry);
         Double_t PhCh[5];
         // Double_t PhCh0 = T->GetLeaf("PhCh0")->GetValue();
         // Double_t PhCh1 = T->GetLeaf("PhCh1")->GetValue();
         // Double_t PhCh2 = T->GetLeaf("PhCh2")->GetValue();
         // Double_t PhCh3 = T->GetLeaf("PhCh3")->GetValue();
         // Double_t PhCh4 = T->GetLeaf("PhCh4")->GetValue();
         PhCh[0] = T->GetLeaf("PhCh0")->GetValue();
         PhCh[1] = T->GetLeaf("PhCh1")->GetValue();
         PhCh[2] = T->GetLeaf("PhCh2")->GetValue();
         PhCh[3] = T->GetLeaf("PhCh3")->GetValue();
         PhCh[4] = T->GetLeaf("PhCh4")->GetValue();
         Double_t NhitLyr0 = T->GetLeaf("NhitLyr0")->GetValue();
         Double_t NhitLyr5 = T->GetLeaf("NhitLyr5")->GetValue();
         //cout<< "NhitLyr0 = " << NhitLyr0 <<endl;
         if (true
               && NhitLyr0 == 1
               && NhitLyr5 == 1
         ) {
            Double_t sumchan = 0;
            // find the last channel
            //if (PhCh0 > thres) sumchan += PhCh0;
            //if (PhCh1 > thres) sumchan += PhCh1;
            //if (PhCh2 > thres) sumchan += PhCh2;
            //if (PhCh3 > thres) sumchan += PhCh3;
            //if (PhCh4 > thres) sumchan += PhCh4;
            for (int i=0; i<5; ++i) {
               if (PhCh[i] > thres) sumchan += PhCh[i];
               else break;
            }
            h_sum->Fill(sumchan);
         }
      }
      // break;
      wet[np] = iwet;
      sum[np] = h_sum->GetMean();
      ++np;

      new TCanvas;
      if (h_sum) h_sum->Draw();
   }

   TGraph* gsum_wet = new TGraph(np, wet, sum);
   gsum_wet->SetNameTitle("gsum_wet", "sum vs WET;WET, mm;sum, ADC counts");
   gsum_wet->SetMarkerStyle(20);
   gsum_wet->SetMarkerColor(4);
   new TCanvas;
   gsum_wet->Draw("ap");

   TGraph* gwet_sum = new TGraph(np, sum, wet);
   gwet_sum->SetNameTitle("gwet_sum", "WET vs sum;sum, ADC counts;WET, mm");
   gwet_sum->SetMarkerStyle(20);
   gwet_sum->SetMarkerColor(2);
   new TCanvas;
   gwet_sum->Draw("ap");

   // gwet_sum->Fit("pol2", "", "", 2000, 21000);
   // ftemp()->SetLineColor(2);
   // gwet_sum->Fit("pol2", "", "", 21000, 28650);
   // ftemp()->SetLineColor(4);
   // gwet_sum->Fit("pol1", "", "", 28650, 28710);
   // ftemp()->SetLineColor(1);
   // gwet_sum->Fit("pol3", "", "", 28710, 32500);
   // ftemp()->SetLineColor(8);
   // gwet_sum->Fit("pol2", "", "", 32500, 37500);
   wetfit(gwet_sum, save);
}
/*
// fitp(2000,21000,2)

****************************************
Minimizer is Linear
Chi2                      =      5.50923
NDf                       =           13
p0                        =      263.402   +/-   0.836681  
p1                        =  -0.00211025   +/-   0.000155485
p2                        = -1.55657e-07   +/-   6.45263e-09

// fitp(21000,28500,2,"+")

****************************************
Minimizer is Linear
Chi2                      =    0.0594618
NDf                       =            5
p0                        =      117.325   +/-   5.07463     
p1                        =   0.00762587   +/-   0.000415297 
p2                        = -2.88021e-07   +/-   8.43691e-09 

// fitp(28500,32500,3,"+")

****************************************
Minimizer is Linear
Chi2                      =      17.0376
NDf                       =            6
p0                        =      28624.3   +/-   17557.6     
p1                        =     -2.63398   +/-   1.74165     
p2                        =  8.09791e-05   +/-   5.75435e-05 
p3                        = -8.30702e-10   +/-   6.33241e-10

// fitp(32500,37500,2,"+")

****************************************
Minimizer is Linear
Chi2                      =    0.0113083
NDf                       =            3
p0                        =     -201.115   +/-   22.1756     
p1                        =    0.0205546   +/-   0.00125493  
p2                        = -4.06795e-07   +/-   1.77334e-08 

//---------------------------------------------------------------------
//------------------------- intermediate ------------------------------
//---------------------------------------------------------------------

// fitp(28650,28710,1,"+")

****************************************
Minimizer is Linear
Chi2                      =  5.58589e-20
NDf                       =            0
p0                        =      4174.88   +/-   1154.11     
p1                        =    -0.142262   +/-   0.0402377   

// fitp(28710,32100,3,"+")

****************************************
Minimizer is Linear
Chi2                      =      2.20836
NDf                       =            4
p0                        =     -1.94522   +/-   10493.5     
p1                        =     0.177744   +/-   1.03613     
p2                        = -1.10193e-05   +/-   3.40802e-05 
p3                        =  1.72027e-10   +/-   3.734e-10   
*/

void sciwet(const char* filename, Double_t thres=100)
{
   // open coefficients
   TFile* filefit = TFile::Open("wetfit.root");
   if (!filefit) {
      cout<< "Could not open FitFunctions file \"wetfit.root\"" <<endl;
      return;
   }
   FitFunctions* fitFunctions = (FitFunctions*) filefit->Get("objFitFunctions");
   if (!fitFunctions) {
      cout<< "Could not retieve objFitFunctions from \"wetfit.root\"" <<endl;
      return;
   }

   cd(filename);
   TTree* T = (TTree*) gDirectory->Get("T");
   Int_t iwet = 0;
   sscanf(filename, "t-%d.root", &iwet);
   TH1F* h_wet = new TH1F(Form("h_wet_%d",iwet), Form("h_wet for %s",filename), 500,-100,400);
   TH1F* h_sumchan = new TH1F(Form("h_sumchan_%d",iwet), Form("h_sumchan for %s",filename), 10000,0,100000);
   for (int ientry=0; ientry<T->GetEntries(); ++ientry)
   {
      if (T->LoadTree(ientry) < 0) break;
      T->GetEntry(ientry);
      Double_t PhCh[5];
      PhCh[0] = T->GetLeaf("PhCh0")->GetValue();
      PhCh[1] = T->GetLeaf("PhCh1")->GetValue();
      PhCh[2] = T->GetLeaf("PhCh2")->GetValue();
      PhCh[3] = T->GetLeaf("PhCh3")->GetValue();
      PhCh[4] = T->GetLeaf("PhCh4")->GetValue();
      Double_t NhitLyr0 = T->GetLeaf("NhitLyr0")->GetValue();
      Double_t NhitLyr5 = T->GetLeaf("NhitLyr5")->GetValue();
      if (true
            && NhitLyr0 == 1
            && NhitLyr5 == 1
         ) {
         Double_t sumchan = 0;
         // find the last channel
         for (int i=0; i<5; ++i) {
            if (PhCh[i] > thres) sumchan += PhCh[i];
            else break;
         }
         h_sumchan->Fill(sumchan);
         h_wet->Fill(fitFunctions->Eval(sumchan));
      }
   }

   new TCanvas;
   h_sumchan->Draw();

   new TCanvas;
   h_wet->Draw();
}

void wet_ch1_ch2(const char* filename)
{
   TFile* funfile = TFile::Open("fit_functions_ch1_ch2.root");
   TF1* f1 = (TF1*) funfile->Get("wet_chan1");
   TF1* f2 = (TF1*) funfile->Get("wet_chan2");
   // new TCanvas;
   // f1->Draw();
   // new TCanvas;
   // f2->Draw();

   Double_t xmin1 = 1500;
   Double_t xmax1 = 14300;
   Double_t xmin2 = 1500;
   Double_t xmax2 = 17000;

   cd(filename);
   TTree* T = (TTree*) gDirectory->Get("T");
   Int_t iwet = 0;
   sscanf(filename, "t-%d.root", &iwet);
   TH1F* h_wetrec1 = new TH1F(Form("h_wetrec1_%d",iwet), Form("h_wetrec1 for %s",filename), 500,-100,400);
   TH1F* h_wetrec2 = new TH1F(Form("h_wetrec2_%d",iwet), Form("h_wetrec2 for %s",filename), 500,-100,400);
   for (int ientry=0; ientry<T->GetEntries(); ++ientry)
   // for (int ientry=0; ientry<100; ++ientry)
   {
      if (T->LoadTree(ientry) < 0) break;
      T->GetEntry(ientry);
      Double_t PhCh[5];
      PhCh[0] = T->GetLeaf("PhCh0")->GetValue();
      PhCh[1] = T->GetLeaf("PhCh1")->GetValue();
      PhCh[2] = T->GetLeaf("PhCh2")->GetValue();
      PhCh[3] = T->GetLeaf("PhCh3")->GetValue();
      PhCh[4] = T->GetLeaf("PhCh4")->GetValue();
      Double_t NhitLyr0 = T->GetLeaf("NhitLyr0")->GetValue();
      Double_t NhitLyr5 = T->GetLeaf("NhitLyr5")->GetValue();
      if (true
            && NhitLyr0 == 1
            && NhitLyr5 == 1
         ) {
         if ((PhCh[1] > xmin1 && PhCh[2] < 100) && PhCh[1] < xmax1) {
            Double_t wetval1 = f1->GetX(PhCh[1]);
            //cout<< "wetval1 = " << wetval1 <<endl;
            h_wetrec1->Fill(wetval1);
            continue;
         }
         if ((PhCh[2] > xmin2 && PhCh[3] < 100) && PhCh[2] < xmax2) {
            Double_t wetval2 = f2->GetX(PhCh[2]);
            //cout<< "wetval2 = " << wetval2 <<endl;
            h_wetrec2->Fill(wetval2);
            continue;
         }
      }
   }

   new TCanvas;
   h_wetrec1->Draw();

   new TCanvas;
   h_wetrec2->Draw();
}

/*
cd("t-203.root");
T->Draw("PhCh0","PhCh0>100&&PhCh1<100 &&PhCh0>12000&&PhCh0<15000 &&NhitLyr0==1&&NhitLyr5==1");
187979
// fitgr(13600,14200)
 FCN=230.01 FROM MIGRAD    STATUS=CONVERGED      85 CALLS          86 TOTAL
                     EDM=5.06728e-07    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     2.86440e+03   1.03136e+01   6.00333e-02   6.02046e-05
   2  Mean         1.37974e+04   9.27989e-01   6.57909e-03   2.37051e-05
   3  Sigma        1.85378e+02   8.44872e-01   7.95371e-06   6.97873e-01


cd("t-211.root")
// T->Draw("PhCh0","PhCh0>100&&PhCh1<100 &&PhCh0>10000&&PhCh0<15000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)676725
// fitgr(12600,14000)
 FCN=1563.52 FROM MIGRAD    STATUS=CONVERGED      84 CALLS          85 TOTAL
                     EDM=1.281e-06    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     8.91100e+03   1.50978e+01   2.32883e-01  -5.83285e-05
   2  Mean         1.30880e+04   9.67612e-01   1.40592e-02   6.97619e-04
   3  Sigma        4.25136e+02   8.44889e-01   9.39881e-06   1.53347e+00

cd("t-216.root")
T->Draw("PhCh0","PhCh0>100&&PhCh1<100 &&PhCh0>10000&&PhCh0<15000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)641026
// fitgr(11700,13000)
 FCN=378.038 FROM MIGRAD    STATUS=CONVERGED      76 CALLS          77 TOTAL
                     EDM=2.56332e-08    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   2.7 per cent
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     7.83432e+03   1.49876e+01   5.62314e-02  -1.13593e-05
   2  Mean         1.22734e+04   9.84748e-01   2.14318e-03  -1.06938e-04
   3  Sigma        4.68725e+02   1.26511e+00  -3.98218e-06   1.07047e-01

cd("t-224.root")
// T->Draw("PhCh0","PhCh0>100&&PhCh1<100 &&PhCh0>8000&&PhCh0<15000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)532343
// fitgr(10000,12000)
 FCN=831.814 FROM MIGRAD    STATUS=CONVERGED      83 CALLS          84 TOTAL
                     EDM=3.86985e-08    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     7.96226e+03   1.50090e+01   1.63835e-01   1.04870e-05
   2  Mean         1.07961e+04   1.03245e+00   1.29115e-02   1.68644e-04
   3  Sigma        5.29676e+02   9.45331e-01   6.53489e-06   5.83233e-01


// cd("t-229.root")
-->      t-229.root
// T->Draw("PhCh0","PhCh0>100&&PhCh1<100 &&PhCh0>6000&&PhCh0<15000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)652168
// fitgr(9000,11000)
 FCN=543.206 FROM MIGRAD    STATUS=CONVERGED      83 CALLS          84 TOTAL
                     EDM=2.67264e-07    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     1.15955e+04   2.02270e+01   1.76726e-01  -4.68267e-05
   2  Mean         9.79696e+03   1.06289e+00   1.06308e-02  -7.82602e-05
   3  Sigma        5.70151e+02   1.02812e+00   5.35004e-06  -1.04830e+00


// cd("t-238.root")
-->      t-238.root
// T->Draw("PhCh0","PhCh0>100&&PhCh1<100 &&PhCh0>4000&&PhCh0<15000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)412254
// fitgr(7000,9000)
 FCN=261.423 FROM MIGRAD    STATUS=CONVERGED      73 CALLS          74 TOTAL
                     EDM=2.32529e-09    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   2.3 per cent
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     7.77961e+03   1.82528e+01   4.52333e-02  -3.17949e-06
   2  Mean         8.01361e+03   1.54191e+00   5.20612e-04  -4.99842e-06
   3  Sigma        6.68668e+02   1.97876e+00  -2.48704e-06   3.57533e-02


cd("t-243.root")
-->      t-243.root
// T->Draw("PhCh0","PhCh0>100&&PhCh1<100 &&PhCh0>2000&&PhCh0<15000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)534643
// fitgr(6000,8000)
 FCN=234.016 FROM MIGRAD    STATUS=CONVERGED      93 CALLS          94 TOTAL
                     EDM=2.77971e-10    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     1.08620e+04   2.27560e+01   1.25879e-01  -9.48370e-07
   2  Mean         6.82047e+03   1.75935e+00   1.16115e-02  -1.02333e-05
   3  Sigma        7.16155e+02   2.12627e+00   6.07433e-06  -9.42517e-03


cd("t-251.root")
// T->Draw("PhCh0","PhCh0>100&&PhCh1<100 &&PhCh0>100&&PhCh0<10000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)420365
// fitgr(3800,6000)
 FCN=163.949 FROM MIGRAD    STATUS=CONVERGED      76 CALLS          77 TOTAL
                     EDM=2.43584e-09    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   3.5 per cent
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     4.75921e+03   1.12366e+01   3.75082e-02  -8.83143e-06
   2  Mean         4.46160e+03   4.00379e+00   3.92331e-02   3.66551e-06
   3  Sigma        8.60451e+02   3.79164e+00  -1.68509e-05  -1.22588e-02


// cd("t-256.root")
-->      t-256.root
// T->Draw("PhCh0","PhCh0>100&&PhCh1<100 &&PhCh0>100&&PhCh0<10000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)298275
// fitgr(2000,4000)
 FCN=119.142 FROM MIGRAD    STATUS=CONVERGED     110 CALLS         111 TOTAL
                     EDM=2.38574e-08    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     3.44271e+03   9.26275e+00   3.83721e-02  -1.13697e-05
   2  Mean         2.60699e+03   7.89066e+00   2.19385e-02  -1.68281e-05
   3  Sigma        1.10264e+03   9.89662e+00   1.12298e-05   9.55836e-03


   --------------------------------------
   PhCh1
   --------------------------------------

///   // cd("t-137.root")
///   -->      t-137.root
///   // Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1_n25
///   T->Draw("PhCh1","PhCh1>100&&PhCh2<100 &&PhCh1>10000&&PhCh1<20000 &&NhitLyr0==1&&NhitLyr5==1")
///   (Long64_t)36067
///   // fitgr(15200,16000)
///    FCN=20.4711 FROM MIGRAD    STATUS=CONVERGED     104 CALLS         105 TOTAL
///                        EDM=3.44138e-09    STRATEGY= 1      ERROR MATRIX ACCURATE 
///     EXT PARAMETER                                   STEP         FIRST   
///     NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
///      1  Constant     1.45395e+02   3.30007e+00   5.47772e-03  -2.97350e-05
///      2  Mean         1.54932e+04   2.01675e+01   3.10201e-02  -3.65123e-06
///      3  Sigma        4.38387e+02   3.10315e+01   4.40762e-05  -2.21153e-03


-->      t-145.root
// T->Draw("PhCh1","PhCh1>100&&PhCh2<100 &&PhCh1>10000&&PhCh1<20000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)58535
// fitgr(14000,15000)
 FCN=167.612 FROM MIGRAD    STATUS=CONVERGED      81 CALLS          82 TOTAL
                     EDM=7.26211e-13    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   4.1 per cent
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     8.72228e+02   7.75936e+00   6.29204e-03   1.94561e-07
   2  Mean         1.43540e+04   5.57522e+00   4.04762e-03  -1.46030e-08
   3  Sigma        4.04870e+02   6.62650e+00  -5.85803e-06  -3.20611e-05

/ cd("t-150.root")
// T->Draw("PhCh1","PhCh1>100&&PhCh2<100 &&PhCh1>10000&&PhCh1<20000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)240237
// fitgr(13900,14700)
 FCN=466.85 FROM MIGRAD    STATUS=CONVERGED      84 CALLS          85 TOTAL
                     EDM=1.76874e-08    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     1.00953e+04   3.17470e+01   2.64378e-01  -1.81514e-06
   2  Mean         1.41744e+04   8.13918e-01   6.75888e-03   1.52208e-04
   3  Sigma        2.13770e+02   6.93025e-01   8.24778e-06   1.54829e-01


// cd("t-158.root")
-->      t-158.root
// T->Draw("PhCh1","PhCh1>100&&PhCh2<100 &&PhCh1>10000&&PhCh1<20000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)654665
// fitgr(12800,14000)
 FCN=225.987 FROM MIGRAD    STATUS=CONVERGED      92 CALLS          93 TOTAL
                     EDM=1.35871e-06    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     1.61200e+04   3.18582e+01   1.67215e-01   5.61603e-05
   2  Mean         1.33683e+04   9.90752e-01   7.20524e-03  -1.78014e-04
   3  Sigma        4.61716e+02   1.45451e+00   6.90641e-06   1.69139e+00


// cd("t-163.root")
// T->Draw("PhCh1","PhCh1>100&&PhCh2<100 &&PhCh1>8000&&PhCh1<16000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)684138
// fitgr(12000,13500)
 FCN=413.62 FROM MIGRAD    STATUS=CONVERGED      76 CALLS          77 TOTAL
                     EDM=2.04631e-09    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   2.2 per cent
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     1.23622e+04   2.18437e+01   2.72499e-03   2.46793e-06
   2  Mean         1.25261e+04   1.14802e+00   1.09094e-03  -5.55074e-06
   3  Sigma        4.82110e+02   1.10111e+00  -1.07079e-06   1.12624e-01

// cd("t-172.root")
-->      t-172.root
// Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1_n17
T->Draw("PhCh1","PhCh1>100&&PhCh2<100 &&PhCh1>6000&&PhCh1<16000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)692248
// fitgr(10000,12000)
 FCN=708.718 FROM MIGRAD    STATUS=CONVERGED      82 CALLS          83 TOTAL
                     EDM=6.05034e-09    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     1.33742e+04   2.30895e+01   2.23820e-01   2.27430e-06
   2  Mean         1.09511e+04   9.42768e-01   1.20761e-02   8.81345e-05
   3  Sigma        5.86303e+02   1.02854e+00   6.30187e-06   1.72117e-01

/ cd("t-177.root")
-->      t-177.root
// T->Draw("PhCh1","PhCh1>100&&PhCh2<100 &&PhCh1>6000&&PhCh1<16000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)627939
// fitgr(9000,11000)
 FCN=293.733 FROM MIGRAD    STATUS=CONVERGED      82 CALLS          83 TOTAL
                     EDM=3.59127e-08    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     1.19666e+04   2.19567e+01   1.36378e-01  -5.67443e-06
   2  Mean         9.92210e+03   1.03720e+00   8.46941e-03   1.39308e-04
   3  Sigma        5.93753e+02   1.16089e+00   4.51058e-06   2.05842e-01


// cd("t-185.root")
-->      t-185.root
// T->Draw("PhCh1","PhCh1>100&&PhCh2<100 &&PhCh1>4000&&PhCh1<14000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)698030
// fitgr(7500,9500)
 FCN=427.87 FROM MIGRAD    STATUS=CONVERGED      75 CALLS          76 TOTAL
                     EDM=8.35548e-09    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   2.3 per cent
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     1.16523e+04   2.03078e+01   1.59540e-03   5.38423e-06
   2  Mean         8.17778e+03   1.65831e+00   2.61486e-03  -1.03394e-05
   3  Sigma        6.52788e+02   1.53282e+00  -1.65870e-06   2.10465e-01

cd("t-190.root")
-->      t-190.root
// T->Draw("PhCh1","PhCh1>100&&PhCh2<100 &&PhCh1>2000&&PhCh1<12000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)662039
// fitgr(6000,8000)
 FCN=345.528 FROM MIGRAD    STATUS=CONVERGED      76 CALLS          77 TOTAL
                     EDM=6.034e-07    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   2.8 per cent
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     9.64511e+03   1.87087e+01   2.79355e-01  -4.36491e-05
   2  Mean         6.87382e+03   1.71520e+00   1.01159e-02  -1.84398e-04
   3  Sigma        7.64667e+02   2.28970e+00  -1.69709e-05   5.58231e-01
   

cd("t-198.root")
-->      t-198.root
// T->Draw("PhCh1","PhCh1>100&&PhCh2<100 &&PhCh1>100&&PhCh1<10000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)639495
// fitgr(3800,6000)
 FCN=314.977 FROM MIGRAD    STATUS=CONVERGED      93 CALLS          94 TOTAL
                     EDM=2.8434e-08    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     6.50541e+03   1.29483e+01   9.12336e-02   3.64341e-06
   2  Mean         4.38116e+03   4.15170e+00   1.80625e-02   7.78845e-06
   3  Sigma        9.06435e+02   3.70092e+00   7.25896e-06   1.55550e-01

cd("t-203.root")
-->      t-203.root
// T->Draw("PhCh1","PhCh1>100&&PhCh2<100 &&PhCh1>100&&PhCh1<10000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)498766
// fitgr(1700,4000)
 FCN=210.917 FROM MIGRAD    STATUS=CONVERGED      84 CALLS          85 TOTAL
                     EDM=2.71889e-09    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   3.7 per cent
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     4.34227e+03   1.00119e+01  -2.80509e-03   8.33713e-06
   2  Mean         2.39236e+03   6.76270e+00   9.34384e-03   1.85041e-05
   3  Sigma        1.16936e+03   7.85280e+00  -2.42460e-06   6.12986e-02

   chan1_wet_vs_ampl.dat
   ---------------------
14354     145
14174     150
13368     158
12526     163
10951     172
9922      177
8178      185
6874      190
4381      198
2392      203

   ------------>
   // g1 = new TGraph("wet_chan1.dat")
(class TGraph*)0x3b82d50
// g1->SetMarkerStyle(20)
// g1->Draw("ap")
// g1 = new TGraph("wet_chan1.dat")
(class TGraph*)0x3bdf8f0
// g1->SetMarkerStyle(20)
// g1->Draw("ap")

// fitp(0,0,4)

****************************************
Minimizer is Linear
Chi2                      =       3919.3
NDf                       =            5
p0                        = -1.33097e+06   +/-   118952      
p1                        =      31228.7   +/-   2770.43     
p2                        =     -270.451   +/-   24.0909     
p3                        =      1.03942   +/-   0.0927032   
p4                        =  -0.00150528   +/-   0.000133202 

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

ch2

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// cd("t-150.root")
// T->Draw("PhCh2","PhCh2>100&&PhCh3<100 &&PhCh2>100&&PhCh2<30000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)412161
// fitgr(1700,4500)
 FCN=135.693 FROM MIGRAD    STATUS=CONVERGED      99 CALLS         100 TOTAL
                     EDM=1.5856e-07    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     5.11632e+03   1.25881e+01   5.65469e-02  -2.71754e-05
   2  Mean         2.51113e+03   8.53638e+00   2.51085e-02  -1.08049e-04
   3  Sigma        1.39170e+03   9.44554e+00   8.62010e-06  -1.84187e-01

// cd("t-145.root")
T->Draw("PhCh2","PhCh2>100&&PhCh3<100 &&PhCh2>100&&PhCh2<30000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)682587
// fitgr(4000,7000)
 FCN=291.239 FROM MIGRAD    STATUS=CONVERGED      91 CALLS          92 TOTAL
                     EDM=3.13255e-10    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     1.04358e+04   1.93857e+01   1.31120e-01  -1.46596e-06
   2  Mean         4.86391e+03   4.08759e+00   1.92389e-02  -8.25881e-07
   3  Sigma        1.11400e+03   3.61032e+00   5.81127e-06  -1.44850e-02


// cd("t-137.root")
-->      t-137.root
// T->Draw("PhCh2","PhCh2>100&&PhCh3<100 &&PhCh2>100&&PhCh2<30000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)727610
// fitgr(7000,9000)
 FCN=106.329 FROM MIGRAD    STATUS=CONVERGED     101 CALLS         102 TOTAL
                     EDM=7.58669e-07    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     1.42737e+04   2.92662e+01   1.04563e-01  -2.76005e-05
   2  Mean         7.90752e+03   2.48702e+00   1.18966e-02  -3.41627e-04
   3  Sigma        9.33907e+02   4.22149e+00   7.15101e-06  -1.34198e-02


// cd("t-132.root")
// T->Draw("PhCh2","PhCh2>100&&PhCh3<100 &&PhCh2>2000&&PhCh2<20000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)353916
// fitgr(9000,11000)
 FCN=160.607 FROM MIGRAD    STATUS=CONVERGED      95 CALLS          96 TOTAL
                     EDM=3.85043e-08    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     4.20125e+03   1.09254e+01   5.51897e-02  -6.94172e-06
   2  Mean         9.53220e+03   4.57549e+00   1.43671e-02  -1.62220e-05
   3  Sigma        7.89666e+02   3.95622e+00   6.32527e-06   8.36440e-02


/ cd("t-124.root")
-->      t-124.root
// T->Draw("PhCh2","PhCh2>100&&PhCh3<100 &&PhCh2>6000&&PhCh2<15000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)677276
// fitgr(11000,13000)
 FCN=317.363 FROM MIGRAD    STATUS=CONVERGED      77 CALLS          78 TOTAL
                     EDM=2.98333e-08    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   3.3 per cent
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     9.35681e+03   1.71770e+01   8.38832e-02   2.32924e-05
   2  Mean         1.17117e+04   1.88329e+00   8.74986e-03   4.11930e-05
   3  Sigma        7.08096e+02   1.96677e+00  -6.68638e-06   3.45550e-01


// cd("t-119.root")
-->      t-119.root
// T->Draw("PhCh2","PhCh2>100&&PhCh3<100 &&PhCh2>8000&&PhCh2<20000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)522465
// fitgr(12500,14000)
 FCN=89.6396 FROM MIGRAD    STATUS=CONVERGED      84 CALLS          85 TOTAL
                     EDM=4.44684e-09    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   2.8 per cent
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     1.03248e+04   2.27967e+01   3.73655e-02  -4.99848e-06
   2  Mean         1.29691e+04   3.00172e+00   5.40567e-03  -5.11940e-05
   3  Sigma        6.63139e+02   3.24842e+00  -3.01090e-06  -4.87189e-02



/ cd("t-110.root")
-->      t-110.root
// T->Draw("PhCh2","PhCh2>100&&PhCh3<100 &&PhCh2>8000&&PhCh2<20000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)719169
// T->Draw("PhCh2","PhCh2>100&&PhCh3<100 &&PhCh2>12000&&PhCh2<20000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)690519
// fitgr(14500,16000)
 FCN=250.305 FROM MIGRAD    STATUS=CONVERGED      95 CALLS          96 TOTAL
                     EDM=1.47089e-09    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     1.02800e+04   1.86937e+01   1.18350e-01   2.97285e-06
   2  Mean         1.48803e+04   2.45837e+00   9.20181e-03  -1.89064e-06
   3  Sigma        5.84824e+02   2.04690e+00   5.35399e-06   1.39902e-02


// cd("t-105.root")
-->      t-105.root
// T->Draw("PhCh2","PhCh2>100&&PhCh3<100 &&PhCh2>12000&&PhCh2<20000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)692684
// ndi()
(Int_t)509
// fitgr(15500,17000)
 FCN=272.386 FROM MIGRAD    STATUS=CONVERGED      79 CALLS          80 TOTAL
                     EDM=1.63169e-10    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   2.3 per cent
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     1.07598e+04   1.93571e+01   3.47995e-02  -1.13628e-06
   2  Mean         1.59458e+04   1.93529e+00   9.44276e-03   6.29295e-06
   3  Sigma        5.67253e+02   1.79676e+00  -5.18552e-06   3.26655e-03


// cd("t-97.root")
// T->Draw("PhCh2","PhCh2>100&&PhCh3<100 &&PhCh2>15000&&PhCh2<19000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)318261
/ fitgr(16900,17600)
 FCN=273.38 FROM MIGRAD    STATUS=CONVERGED      85 CALLS          86 TOTAL
                     EDM=1.51172e-08    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     4.76298e+03   1.36719e+01   9.14441e-02   1.31165e-05
   2  Mean         1.70863e+04   1.27398e+00   8.14740e-03   1.95351e-04
   3  Sigma        2.28779e+02   9.95946e-01   7.22498e-06   1.75678e-01

   
// cd("t-92.root")
-->      t-92.root
// T->Draw("PhCh2","PhCh2>100&&PhCh3<100 &&PhCh2>15000&&PhCh2<19000 &&NhitLyr0==1&&NhitLyr5==1")
(Long64_t)52522
// fitgr(17000,17500)
 FCN=45.216 FROM MIGRAD    STATUS=CONVERGED     112 CALLS         113 TOTAL
                     EDM=3.40362e-10    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     5.40578e+02   5.55975e+00   1.26920e-02  -3.97155e-06
   2  Mean         1.72415e+04   4.03211e+00   1.33249e-02  -3.41301e-06
   3  Sigma        2.78599e+02   8.65465e+00   3.57810e-05  -2.44205e-05


chan2_ampl_vs_wet.dat
---------------------
92      17241
97      17086
105     15946
110     14880
119     12969
124     11712
132     9532
137     7908
145     4864
150     2511

// g2 = new TGraph("chan2_ampl_vs_wet.dat")
(class TGraph*)0x3e010f0
// g2->SetMarkerStyle(24)
// g2->Draw("ap")
// fitp(0,0,4)

****************************************
Minimizer is Linear
Chi2                      =        17959
NDf                       =            5
p0                        =      -335289   +/-   57917.1     
p1                        =      11974.3   +/-   1967.38     
p2                        =     -149.842   +/-   24.8295     
p3                        =     0.825229   +/-   0.138011    
p4                        =   -0.0017226   +/-   0.000285133 
// png("chan2_ampl_vs_wet.dat")

*/
