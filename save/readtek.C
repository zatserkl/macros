// Andriy Zatserklyaniy <zatserkl@fnal.gov>

#include <TROOT.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TPolyMarker.h>
#include <TMarker.h>

#include <iostream>
#include <fstream>
#include <cstdio>

using std::cout;     using std::endl;

// utils.C stuff
TH1* htemp(TCanvas* can=0);
TGraph* gtemp(TCanvas* can=0);

namespace Tree {
   Double_t x[10000];
   Double_t y[10000];
   Int_t np;

   Double_t ymax;
   Double_t delta;
   Int_t nbaseline;
   Double_t baseline[10];
   Double_t time[10];
   Double_t marker_x[10];
   Double_t marker_y[10];
} // namespace Tree

void readtek(const char* ifname="data-900V/131001_124120.txt")
{
   Int_t np = 0;

   Double_t x[10000];
   Double_t y[10000];

   // std::ifstream ifile(ifname);
   // if (!ifile) {
   //    cout<< "Could not open file " << ifname <<endl;
   //    return;
   // }
   FILE* ifile;
   ifile = fopen(ifname, "r");
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return;
   }

   Int_t RecordLength;
   Double_t SampleInterval;
   Int_t TriggerPoint;
   Double_t TriggerTime;
   Double_t HorizontalOffset;

   //-- NB: \n in the format string, otherwise next stream will be screwed up
   // NB: you can add/remove spaces in the format string of fscanf even in \"Record Length\"

   int nitems = -1;
   nitems = fscanf(ifile, "\"Record Length\" %d \"Points\" %lf %lf \n", &RecordLength, &x[np],&y[np]);
   if (nitems < 0) {
      cout<< ifname << ": Problem reading Record Length and Points" <<endl;
      return;
   }
   np++;
   nitems = fscanf(ifile, "\"Sample Interval\" %lf  s  %lf %lf \n", &SampleInterval, &x[np], &y[np]);
   if (nitems < 0) {
      cout<< ifname << ": Problem reading Sample Interval" <<endl;
      return;
   }
   np++;
   nitems = fscanf(ifile, "\"Trigger Point\" %d \"Samples\" %lf %lf \n", &TriggerPoint, &x[np],&y[np]);
   if (nitems < 0) {
      cout<< ifname << ": Problem reading Trigger Point and Samples" <<endl;
      return;
   }
   np++;
   nitems = fscanf(ifile, "\"Trigger Time\" %lf s %lf %lf \n", &TriggerTime, &x[np],&y[np]);
   if (nitems < 0) {
      cout<< ifname << ": Problem reading Trigger Time" <<endl;
      return;
   }
   np++;
   nitems = fscanf(ifile, "\"\" %lf %lf \n", &x[np],&y[np]);
   if (nitems < 0) {
      cout<< ifname << ": Problem reading line" <<endl;
      return;
   }
   np++;
   nitems = fscanf(ifile, "\"Horizontal Offset\" %lf s %lf %lf \n", &HorizontalOffset, &x[np],&y[np]);
   if (nitems < 0) {
      cout<< ifname << ": Problem reading Horizontal Offset" <<endl;
      return;
   }
   np++;
   while(np < RecordLength) {
      fscanf(ifile, "%lf %lf", &x[np],&y[np]);
      np++;
   }
   fclose(ifile);

   cout<< "Record Length = " << RecordLength << " Points" <<endl;
   cout<< "Sample Interval = " << SampleInterval << " s" <<endl;
   cout<< "Trigger Point = " << TriggerPoint << " Samples" <<endl;
   cout<< "Trigger Time = " << TriggerTime << " s" <<endl;
   cout<< "Horizontal Offset = " << HorizontalOffset << " s" <<endl;

   // for (int i=0; i<np; ++i) {
   for (int i=0; i<10; ++i) {
      if (i == np) break;
      cout<< i+1 <<"\t "<< x[i] <<"\t "<< y[i] <<endl;
   }
}

void readtek(Int_t& np, Double_t* x, Double_t* y, const char* ifname)
{
   np = 0;

   // std::ifstream ifile(ifname);
   // if (!ifile) {
   //    cout<< "Could not open file " << ifname <<endl;
   //    return;
   // }
   FILE* ifile;
   ifile = fopen(ifname, "r");
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return;
   }

   Int_t RecordLength;
   Double_t SampleInterval;
   Int_t TriggerPoint;
   Double_t TriggerTime;
   Double_t HorizontalOffset;

   //-- NB: \n in the format string, otherwise next stream will be screwed up
   // NB: you can add/remove spaces in the format string of fscanf even in \"Record Length\"

   fscanf(ifile, "\"Record Length\" %d \"Points\" %lf %lf \n", &RecordLength, &x[np],&y[np]);
   np++;
   fscanf(ifile, "\"Sample Interval\" %lf  s  %lf %lf \n", &SampleInterval, &x[np], &y[np]);
   np++;
   fscanf(ifile, "\"Trigger Point\" %d \"Samples\" %lf %lf \n", &TriggerPoint, &x[np],&y[np]);
   np++;
   fscanf(ifile, "\"Trigger Time\" %lf s %lf %lf \n", &TriggerTime, &x[np],&y[np]);
   np++;
   fscanf(ifile, "\"\" %lf %lf \n", &x[np],&y[np]);
   np++;
   fscanf(ifile, "\"Horizontal Offset\" %lf s %lf %lf \n", &HorizontalOffset, &x[np],&y[np]);
   np++;
   while(np < RecordLength) {
      fscanf(ifile, "%lf %lf", &x[np],&y[np]);
      np++;
   }
   fclose(ifile);

   // cout<< "Record Length = " << RecordLength << " Points" <<endl;
   // cout<< "Sample Interval = " << SampleInterval << " s" <<endl;
   // cout<< "Trigger Point = " << TriggerPoint << " Samples" <<endl;
   // cout<< "Trigger Time = " << TriggerTime << " s" <<endl;
   // cout<< "Horizontal Offset = " << HorizontalOffset << " s" <<endl;

   // for (int i=0; i<np; ++i) {
   //    cout<< i+1 <<"\t "<< x[i] <<"\t "<< y[i] <<endl;
   // }
}

TTree* tree_list(const char* listname)
{
   // Double_t *x = new Double_t(100000);
   // Double_t *y = new Double_t(100000);
   //-- Double_t x[100000], y[100000];
   //-- Int_t np;

   ifstream listfile(listname);
   if (!listfile) {
      cout<< "Could not open list file " << listname <<endl;
      return 0;
   }

   std::string fname;
   //-- listfile >> fname;
   //-- readtek(np,x,y, fname.c_str());

   TFile* ofile = TFile::Open(Form("%s.root",listname), "recreate");
   TTree* tree = new TTree("t", listname);
   // tree->SetMarkerStyle(6);
   tree->SetMarkerColor(2);

   tree->Branch("np", &Tree::np, "np/I");
   tree->Branch("x", &Tree::x, "x[np]/D");
   tree->Branch("y", &Tree::y, "y[np]/D");
   
   //-- // fill the first entry we just read
   //-- tree->Fill();

   while (listfile >> fname) {
      cout<< fname <<endl;
      readtek(Tree::np,Tree::x,Tree::y,fname.c_str());
      for (int ipoint=0; ipoint<Tree::np; ++ipoint) {
         Tree::x[ipoint] *= 1e9;
         Tree::y[ipoint] *= -1;
      }
      tree->Fill();
   }

   cout<< "Filled " << tree->GetEntries() << " events into file " << ofile->GetName() <<endl;
   ofile->Write();

   new TCanvas;
   tree->Draw("y:x", "Entry$<10", "");

   // delete[] x;
   // delete[] y;
   return tree;
}

// void findPeaks(TH1* h)
// {
//    Int_t npeaks = h->ShowPeaks();
//    cout<< "npeaks = " << npeaks <<endl;
// }

void risetime(Int_t np, const Double_t x[], const Double_t y[])
{
   Int_t imaximum = 0;
   for (int i=1; i<np; ++i) if (y[i] > y[imaximum]) imaximum = i;
   cout<< "x[imaximum] = " << x[imaximum] << " y[imaximum] = " << y[imaximum] <<endl;

   // baseline
   Double_t w = 0;
   Double_t wy = 0;
   Double_t wyy = 0;
   for (int i=0; i<imaximum; ++i) {
      if (x[i] > -10) break;
      w += 1;
      wy += y[i];
      wyy += y[i]*y[i];
   }
   Double_t baseline = wy / w;
   Double_t sigma = TMath::Sqrt((wyy - baseline*wy)/w);
   cout<< "baseline = " << baseline << " +- " << sigma <<endl;

   // // the first peak
   // for (int i=0; i<imaximum; ++i) {
   //    if (y[i] > h->Fill(y[i]);
   // }

   TH1F htmp("htmp", "htmp", 100, 0, y[imaximum]/2);
   htmp.SetDirectory(0);
   for (int i=0; i<imaximum; ++i) if (y[i] > baseline+5*sigma) htmp.Fill(y[i]);

   Int_t ibin_max = 1;
   for (int ibin=1; ibin<=htmp.GetNbinsX(); ++ibin) {
      if (htmp.GetBinContent(ibin) > htmp.GetBinContent(ibin_max)) ibin_max = ibin;
   }

   // consider sigma as sigma of the noise, i.e. the sigma for all peaks
   Double_t xmin = htmp.GetBinCenter(ibin_max)-5*sigma;
   Double_t xmax = htmp.GetBinCenter(ibin_max)+5*sigma;
   // cout<< "Fit from " << htmp.GetBinCenter(ibin_max)-5*sigma << " to " << htmp.GetBinCenter(ibin_max)+5*sigma <<endl;
   // htmp.Fit("gaus", "R", "", htmp.GetBinCenter(ibin_max)-5*sigma, htmp.GetBinCenter(ibin_max)-5*sigma);
   cout<< "Fit from " << xmin << " to " << xmax <<endl;
   htmp.Fit("gaus", "R0", "goff", xmin, xmax);

   cout<< "htmp.GetFunction(\"gaus\")->GetParameter(\"Mean\") = " << htmp.GetFunction("gaus")->GetParameter("Mean") <<endl;

   // new TCanvas;
   // htmp.DrawCopy();  htmp.GetFunction("gaus")->Draw("same");   gPad->Modified();    gPad->Update();
}

void rise_time(Int_t evt1=0, Int_t evt2=-1, const char* ifname="data-900V.lst.root", Int_t peak_index=1)
{
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return;
   }

   TTree* tree = (TTree*) ifile->Get("t");

   //Int_t np;
   //Double_t x[100000], y[100000];

   //tree->SetBranchAddress("np", &np);
   //tree->SetBranchAddress("x", x);
   //tree->SetBranchAddress("y", y);

   TFile* ofile = TFile::Open(Form("%s-rise.root",ifname), "recreate");
   TTree* otree = new TTree("trise", "rise time tree");
   otree->SetMarkerColor(2);

   otree->Branch("x", &Tree::x, "x[10000]/D");
   otree->Branch("y", &Tree::y, "y[10000]/D");
   otree->Branch("np", &Tree::np, "np/I");
   otree->Branch("ymax", &Tree::ymax, "ymax/D");
   otree->Branch("delta", &Tree::delta, "delta/D");
   otree->Branch("nbaseline", &Tree::nbaseline, "nbaseline/I");
   otree->Branch("baseline", &Tree::baseline, "baseline[10]/D");
   otree->Branch("time", &Tree::time, "time[10]/D");
   otree->Branch("marker_x", &Tree::marker_x, "marker_x[10]/D");
   otree->Branch("marker_y", &Tree::marker_y, "marker_y[10]/D");

   //ifile->cd();

   TCanvas* can1 = new TCanvas;
   TCanvas* can2 = new TCanvas;
   // TCanvas* can1 = 0;
   // TCanvas* can2 = 0;

   TMarker* marker = new TMarker();
   marker->SetMarkerStyle(20);
   marker->SetMarkerColor(4);

   if (evt2 < evt1) evt2 = tree->GetEntries() - 1;
   for (int ientry=evt1; ientry<=evt2; ++ientry)
   {
      if (ientry < 10) cout<< "\nientry = " << ientry <<endl;

      //if (tree->LoadTree(ientry) >= 0) tree->GetEntry(ientry);
      //risetime(np, x, y);

      // can1 = new TCanvas;
      // can2 = new TCanvas;

      can1->cd();
      tree->Draw("y:x", Form("Entry$==%d", ientry));
      TGraph* g = gtemp();

      Int_t imaximum = 0;
      for (int i=0; i<g->GetN(); ++i) if (g->GetY()[i] > g->GetY()[imaximum]) imaximum = i;
      // save in the tree the complete data
      Tree::np = g->GetN();
      for (int i=0; i<g->GetN(); ++i) {
         Tree::x[i] = g->GetX()[i];
         Tree::y[i] = g->GetY()[i];
      }

      // redo the TGraph up to maximum position only
      can1->cd();
      tree->Draw("y:x", Form("Entry$==%d &&Iteration$<=%d", ientry,imaximum));
      g = gtemp();

      can2->cd();
      tree->Draw("y", Form("Entry$==%d &&Iteration$<=%d", ientry,imaximum));
      TH1* h = htemp();

      Int_t npeaks = h->ShowPeaks(1, "", 0.01);
      if (npeaks > 1) {
         TPolyMarker* pm = (TPolyMarker*) h->GetListOfFunctions()->FindObject("TPolyMarker");
         Double_t pm_x[10];
         //Double_t pm_y[10];
         Int_t index[10];
         if (pm) {
            //cout<< "pm->GetN() = " << pm->GetN() <<endl;
            for (int ipm=0; ipm<pm->GetN(); ++ipm) {
               cout<< "pm->GetX()[ipm] = " << pm->GetX()[ipm] << " pm->GetY()[ipm] = " << pm->GetY()[ipm] <<endl;
               if (ipm < 10) {
                  pm_x[ipm] = pm->GetX()[ipm];
                  //pm_y[ipm] = pm->GetY()[ipm];
               }
            }
            TMath::Sort(pm->GetN(), pm_x, index, kFALSE);

            for (int ibaseline=0; ibaseline<Tree::nbaseline && ibaseline<10; ++ibaseline) {
               Tree::baseline[ibaseline] = pm->GetX()[index[ibaseline]];
               cout<< "baseline[" << ibaseline << "] = " << Tree::baseline[ibaseline] <<endl;
            }
            Double_t* ybaseline = 0;
            if (peak_index < pm->GetN()) ybaseline = &Tree::baseline[peak_index];
            if (!ybaseline) continue;

            Double_t delta = g->GetY()[imaximum] - *ybaseline;

            Tree::ymax = g->GetY()[imaximum];
            Tree::delta = delta;

            // find the 10% point
            Int_t it10 = imaximum;
            for (int i=imaximum; i>=0; --i) {
               if (g->GetY()[i] < *ybaseline + 0.1*delta) {
                  it10 = i;
                  break;
               }
            }
            cout<< "delta = " << delta <<endl;
            
            Int_t itime[10];
            Double_t time[10];
            for (int i=0; i<10; ++i) {
               itime[i] = 0;
               time[i] = 0;
            }
            itime[1] = it10;
            cout<< "g->GetX()[itime[1]] = " << g->GetX()[itime[1]] <<endl;

            Double_t marker_x[10], marker_y[10];
            for (int i=0; i<10; ++i) marker_x[i] = marker_y[i] = 0;
            marker_x[1] = g->GetX()[itime[1]];
            marker_y[1] = g->GetY()[itime[1]];

            cout<< "ipoint " << 1 << " x = " << g->GetX()[itime[1]] << "   time[ipoint] = " << time[1] <<endl;

            for (int i=it10+1; i<g->GetN(); ++i) {
               Double_t dt = g->GetX()[i] - g->GetX()[it10];
               Double_t dy = g->GetY()[i] - *ybaseline;
               for (int ipoint=2; ipoint<10; ++ipoint) {
                  if (itime[ipoint] == 0 && dy > 0.10*ipoint*delta) {
                     itime[ipoint] = i;
                     time[ipoint] = dt;
                     cout<< "ipoint " << ipoint << " x = " << g->GetX()[itime[ipoint]] << "   time[ipoint] = " << time[ipoint] <<endl;
                     marker_x[ipoint] = g->GetX()[i];
                     marker_y[ipoint] = g->GetY()[i];
                  }
               }
            }

            for (int i=0; i<10; ++i) {
               Tree::time[i] = time[i];
               Tree::marker_x[i] = marker_x[i];
               Tree::marker_y[i] = marker_y[i];
            }

            can1->cd();
            for (int i=1; i<10; ++i) {
               marker->DrawMarker(marker_x[i], marker_y[i]);
            }
         }
      }

      otree->Fill();
   }

   cout<< "otree->GetEntries() = " << otree->GetEntries() <<endl;
   ofile->Write();

   tree->ResetBranchAddresses();
}
