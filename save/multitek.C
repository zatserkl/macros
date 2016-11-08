// Andriy Zatserklyaniy <zatserkl@fnal.gov>

#include <TROOT.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>

#include <iostream>
#include <sstream>
#include <fstream>

using std::cout;     using std::endl;

void multitek(const char* ifname, Float_t vscale=1., Float_t tscale=1.)
{
   Int_t np = 0;

   Float_t t[10000];
   Float_t v[10000];

   std::ifstream ifile(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return;
   }

   Int_t RecordLength = 0;
   Float_t SampleInterval = 0;
   Int_t TriggerPoint = 0;
   Float_t TriggerTime = 0;
   Float_t HorizontalOffset = 0;
   Int_t FastFrameCount = 1;           //-- for multiple spectra acquisition only

   std::string s;                      // for Sample Interval, Trigger Time, and Horizontal Offset
   std::string Record, Length, Points;
   std::string Sample, Interval;
   std::string Trigger, Point, Samples;
   std::string /*Trigger,*/ Time;
   std::string quotes;
   std::string Horizontal, Offset;
   std::string FastFrame, Count, Frames;

   // NB: getline(ifile,line) after ifile >> i will read LF

   std::string line;
   std::istringstream ss;

   getline(ifile,line); ss.str(line);
   ss >> Record >> Length >> RecordLength >> Points >> t[np] >> v[np];
   t[np] *= tscale;
   v[np] *= vscale;
   ++np;

   getline(ifile,line); ss.str(line);
   ss >> Sample >> Interval >> SampleInterval >> s >> t[np] >> v[np];
   t[np] *= tscale;
   v[np] *= vscale;
   ++np;

   getline(ifile,line); ss.str(line);
   ss >> Trigger >> Point >> TriggerPoint >> Samples >> t[np] >> v[np];
   t[np] *= tscale;
   v[np] *= vscale;
   ++np;

   getline(ifile,line); ss.str(line);
   ss >> Trigger >> Time >> TriggerTime >> s >> t[np] >> v[np];
   t[np] *= tscale;
   v[np] *= vscale;
   ++np;

   getline(ifile,line); ss.str(line);
   ss >> quotes >> t[np] >> v[np];
   t[np] *= tscale;
   v[np] *= vscale;
   ++np;

   getline(ifile,line); ss.str(line);
   ss >> Horizontal >> Offset >> HorizontalOffset >> s >> t[np] >> v[np];
   t[np] *= tscale;
   v[np] *= vscale;
   ++np;

   getline(ifile,line); ss.str(line);
   if (line.find("\"FastFrame Count\"") != std::string::npos) {
      ss >> FastFrame >> Count >> FastFrameCount >> Frames >> t[np] >> v[np];
   }
   else ss >> t[np] >> v[np];
   t[np] *= tscale;
   v[np] *= vscale;
   ++np;

   TFile* ofile = new TFile(Form("%s.root",ifname), "recreate");
   TTree* tree = new TTree("tek",ifname);
   tree->SetMarkerColor(602);
   tree->Branch("t", &t, Form("t[%d]/F",RecordLength));
   tree->Branch("v", &v, Form("v[%d]/F",RecordLength));

   for (int icount=0; icount<FastFrameCount; ++icount) {
      while (np < RecordLength) {
         ifile >> t[np] >> v[np];
         t[np] *= tscale;
         v[np] *= vscale;
         ++np;
      }
      tree->Fill();
      np = 0;
   }
   np = RecordLength;      // restore value of np

   // print parameters
   cout<< Record <<" "<< Length <<" "<< RecordLength <<" "<< Points <<endl;
   cout<< Sample <<" "<< Interval <<" "<< SampleInterval <<" "<< s <<endl;
   cout<< Trigger <<" "<< Point <<" "<< TriggerPoint <<" "<< Samples <<endl;
   cout<< Trigger <<" "<< Time <<" "<< TriggerTime <<" "<< s <<endl;
   cout<< quotes <<endl;
   cout<< Horizontal <<" "<< Offset <<" "<< HorizontalOffset <<" "<< s <<endl;
   cout<< FastFrame <<" "<< Count <<" "<< FastFrameCount <<" "<< Frames <<endl;

   // // print first 10 points of the first event
   // for (int i=0; i<10; ++i) {
   //    cout<< i << "\t" << t[i] <<" "<< v[i] <<endl;
   // }

   cout<< "write " << tree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
   ofile->Write();

   // show the first event
   new TCanvas;
   cout<< "\ntek->Draw(\"v:t\", \"Entry$==0\", \"L\")" <<endl;
   tree->Draw("v:t", "Entry$==0", "L");
}

void tekPlot(Int_t event, Option_t* option="L", TTree* tree=0)
{
   if (!tree) tree = (TTree*) gDirectory->Get("tek");
   if (!tree) {
      cout<< "Could not find tree \"tek\" in the current directory" <<endl;
      return;
   }

   new TCanvas;
   tree->Draw("v:t", Form("Entry$==%d",event), option);
}

void tekPlot(const char* ifname, Int_t event)
{
   TFile* ifile = new TFile(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifile <<endl;
      return;
   }
   TTree* tree = (TTree*) ifile->Get("tek");
   if (!tree) {
      cout<< "Could not find tree \"tek\" in the file " << ifile <<endl;
      return;
   }

   tekPlot(event, "L", tree);
}

TTree* tekData(TTree* tree=0)
{
   if (!tree) tree = (TTree*) gDirectory->Get("tek");
   if (!tree) {
      cout<< "Could not find tree \"tek\"" <<endl;
      return 0;
   }

   // get the data dimension
   // Int_t RecordLength = tree->GetLeaf("t")->GetLen();
   TLeaf* tleaf = tree->GetLeaf("t");
   if (!tleaf) {
      cout<< "Could not find leaf \"t\" in the tree" <<endl;
      return 0;
   }
   Int_t RecordLength = tleaf->GetLen();
   cout<< "RecordLength = " << RecordLength <<endl;

   Float_t t[10000], v[10000];

   tree->SetBranchAddress("t", &t);
   tree->SetBranchAddress("v", &v);

   // create analysis tree
   TTree* tpulse = new TTree("tpulse", Form("analysis of %s",tree->GetTitle()));
   tpulse->SetMarkerColor(2);

   Float_t pmax;
   tpulse->Branch("pmax", &pmax, "pmax/F");

   cout<< "tree->GetEntries() = " << tree->GetEntries() <<endl;

   for (int ientry=0; ientry<tree->GetEntries(); ++ientry) {
      if (tree->LoadTree(ientry) < 0) break;
      tree->GetEntry(ientry);

      if (ientry % 100 == 0) cout<< "--> processing entry " << ientry <<endl;

      // find pulse maximum
      Float_t vmax = v[0];
      for (int i=0; i<RecordLength; ++i) {
         if (v[i] > vmax) vmax = v[i];
      }

      pmax = vmax;
      tpulse->Fill();
   }
   
   new TCanvas;
   tpulse->Draw("pmax","");

   tree->ResetBranchAddresses();
   return tpulse;
}
