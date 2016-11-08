#include "DRS4Event.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>

#include <iostream>
#include <map>

using std::cout;		using std::endl;

void removeSpikes(Int_t triggerCell, Int_t nchannels, Float_t* y[])
{
   if (nchannels == 0) return;

   Int_t cell0 = 1023 - triggerCell;

   Float_t diff_low = 5;
   Float_t diff_med = 7.5;
   Float_t diff_high = 10;

   std::map<Int_t,Int_t> spikes;

   for (int i=0; i<1024; ++i)
   {
      Int_t curr1 = i;
      Int_t prev1 = (curr1 - 1);  if (prev1 == -1) prev1 = 1023;
      Int_t next1 = (curr1 + 1) % 1024;
      Int_t nextnext1 = (curr1 + 2) % 1024;

      //
      // remove double spikes
      //

      Int_t sign = y[0][curr1] > y[0][prev1]? 1: -1;

      if (sign*(y[0][curr1] - y[0][prev1]) > diff_high && sign*(y[0][next1] - y[0][nextnext1]) > diff_high)   // diff_high to find a candidate
      {
         // look at channel symmetric wrt cell0
         Int_t dspike = cell0 - curr1;             // NB: dspike may be negative
         Int_t curr2 = (cell0 + dspike) % 1024;    if (curr2 < 0) curr2 += 1024;
         Int_t prev2 = curr2 - 1;                  if (prev2 == -1) prev2 = 1023;
         Int_t next2 = (curr2 + 1) % 1024;
         Int_t nextnext2 = (curr2 + 2) % 1024;

         if (sign*(y[0][curr2] - y[0][prev2]) > diff_low && sign*(y[0][next2] - y[0][nextnext2]) > diff_low)  // diff_low to find symmetric wrt cell0
         {
            Float_t diff = nchannels > 3? diff_med: diff_low;                                                // diff_med (if nchannels > 2) to find confirmation in the rest of the channels
            Int_t nothers = 0;
            for (int ich=1; ich<nchannels; ++ich)
            {
               if (y[ich]) {
                  if (true
                      && sign*(y[ich][curr1] - y[ich][prev1]) > diff
                      && sign*(y[ich][next1] - y[ich][nextnext1]) > diff
                      && sign*(y[ich][curr2] - y[ich][prev2]) > diff
                      && sign*(y[ich][next2] - y[ich][nextnext2]) > diff
                     )
                  {
                     nothers++;
                     break;
                  }
               }
            }

            if (nchannels < 2 || nothers > 0) {
               spikes[curr1] += sign;
               spikes[curr2] += sign;
            }
         }
      }
   }

   for (std::map<Int_t,Int_t>::const_iterator it=spikes.begin(); it!=spikes.end(); ++it) {
      Int_t curr = it->first;
      Int_t next = curr + 1;     if (next > 1023) next = 0;
      Float_t sign = it->second > 0? 1: -1;
      for (int ich=0; ich<nchannels; ++ich) {
         if (y[ich]) {
            y[ich][curr] -= sign*14.8;
            y[ich][next] -= sign*14.8;
         }
      }
   }

   //
   // remove single spike
   //
   for (int i=0; i<1024; ++i)
   {
      Int_t curr1 = i;
      Int_t prev1 = (curr1 - 1);  if (prev1 == -1) prev1 = 1023;
      Int_t next1 = (curr1 + 1) % 1024;
      // Int_t nextnext1 = (curr1 + 2) % 1024;  // not in use for single spikes
      Int_t sign = 0;

      for (int ich=0; ich<nchannels; ++ich) {
         if (y[ich]) {
            sign = y[ich][curr1] > y[ich][prev1]? 1: -1;
            if (sign*(y[ich][curr1] - y[ich][prev1]) > sign*diff_high && sign*(y[ich][curr1] - y[ich][next1]) > sign*diff_high) y[ich][curr1] = 0.5*(y[ich][prev1]+y[ich][next1]);
         }
      }
   }
}

void removeSpikes(const char* ifname)
{
   TFile* ifile = new TFile(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return;
   }
   TTree* itree = (TTree*) gDirectory->Get("pulse");
   if (!itree) {
      cout<< "Error: could not find itree \"pulse\"" <<endl;
      return;
   }

   DRS4Event* drs4Event = new DRS4Event(itree);

   TFile* ofile = TFile::Open(Form("%s.pulse.root",ifname), "recreate");

   TTree* otree = new TTree("p", "new pulse tree");
   otree->SetMarkerStyle(6);
   otree->SetMarkerColor(2);
   otree->SetLineColor(2);

   //
   // output tree uses the same buffers
   //
   drs4Event->bookSyncTree(otree);

   ///////////////////////////////////////////
   //
   //       Spike removing
   //
   ///////////////////////////////////////////

   for (Long64_t jentry=0; jentry<itree->GetEntries(); jentry++)
   {
      if (jentry % 1000 == 0) cout<< "jentry = " << jentry <<endl;
      itree->LoadTree(jentry);
      itree->GetEntry(jentry);

      for (int iboard=0; iboard<3; ++iboard)
      {
         if (drs4Event->boardT(iboard))
         {
            Float_t* ychan[4];
            Int_t nchan = 0;

            for (int ich=0; ich<4; ++ich) {
               Float_t* ycurr = drs4Event->boardV(ich,iboard);
               if (ycurr) ychan[nchan++] = ycurr;
            }

            removeSpikes(drs4Event->triggerCell(iboard), nchan, ychan);
         }
      }

      otree->Fill();
   }

   cout<< "Write " << otree->GetEntries() << " events into file " << ofile->GetName() <<endl;
   ofile->Write();

   itree->ResetBranchAddresses();
   otree->ResetBranchAddresses();
}
