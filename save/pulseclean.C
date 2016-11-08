#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>

#include <iostream>
#include <map>

using std::cout;		using std::endl;

/*
root -l tb2014_run_082.root tb2014_run_082.root.pulse.root
cd(0)
.L pulseclean.C+
pulseclean(3, 299)
pulseclean(3, 11)
pulseclean(3, 82)
pulseclean(3, 85)
*/
void pulseclean_plot(Int_t board, Int_t event, TTree* tree=0)
{
   if (!tree) tree = (TTree*) gDirectory->Get("pulse");
   if (!tree) {
      cout<< "Error: could not find tree \"pulse\"" <<endl;
      return;
   }

   Int_t tc1, tc2, tc3;
   Float_t b1_t[1024], b1_c[4096], b2_t[1024], b2_c[4096], b3_t[1024], b3_c[4096];
   //-- Int_t tc[3];
   Float_t* b_t[3];
   Float_t* b_c[3][4];
   b_t[0] = b1_t;
   b_t[1] = b2_t;
   b_t[2] = b3_t;
   b_c[0][0] = b1_c;
   b_c[0][1] = b1_c + 1024;
   b_c[0][2] = b1_c + 2048;
   b_c[0][3] = b1_c + 3072;
   b_c[1][0] = b2_c;
   b_c[1][1] = b2_c + 1024;
   b_c[1][2] = b2_c + 2048;
   b_c[1][3] = b2_c + 3072;
   b_c[2][0] = b3_c;
   b_c[2][1] = b3_c + 1024;
   b_c[2][2] = b3_c + 2048;
   b_c[2][3] = b3_c + 3072;

   Float_t* b1_c1 = b1_c;
   Float_t* b1_c2 = b1_c + 1024;
   Float_t* b1_c3 = b1_c + 2048;
   Float_t* b1_c4 = b1_c + 3072;
   Float_t* b2_c1 = b2_c;
   Float_t* b2_c2 = b2_c + 1024;
   Float_t* b2_c3 = b2_c + 2048;
   Float_t* b2_c4 = b2_c + 3072;
   Float_t* b3_c1 = b3_c;
   Float_t* b3_c2 = b3_c + 1024;
   Float_t* b3_c3 = b3_c + 2048;
   Float_t* b3_c4 = b3_c + 3072;

   tree->SetBranchAddress("event", &event);
   tree->SetBranchAddress("tc1", &tc1);
   tree->SetBranchAddress("tc2", &tc2);
   tree->SetBranchAddress("tc3", &tc3);
   tree->SetBranchAddress("b1_t", &b1_t);
   tree->SetBranchAddress("b2_t", &b2_t);
   tree->SetBranchAddress("b3_t", &b3_t);
   tree->SetBranchAddress("b1_c", &b1_c);
   tree->SetBranchAddress("b2_c", &b2_c);
   tree->SetBranchAddress("b3_c", &b3_c);

   Float_t* x = 0;
   Float_t* y1 = 0;
   Float_t* y2 = 0;
   Float_t* y3 = 0;
   Float_t* y4 = 0;

   Int_t* tc = 0;

   switch (board) {
      case 1: tc = &tc1;
              x = b1_t;
              y1 = b1_c1;
              y2 = b1_c2;
              y3 = b1_c3;
              y4 = b1_c4;
              break;
      case 2: tc = &tc2;
              x = b2_t;
              y1 = b2_c1;
              y2 = b2_c2;
              y3 = b2_c3;
              y4 = b2_c4;
              break;
      case 3: tc = &tc3;
              x = b3_t;
              y1 = b3_c1;
              y2 = b3_c2;
              y3 = b3_c3;
              y4 = b3_c4;
              break;
      default: cout<< "board " << board << " is out of range 1..3" <<endl;
               return;
   }

   if (tree->LoadTree(event) < 0) {
      cout<< "Entry was not found: " << event <<endl;
      tree->ResetBranchAddresses();
      return;
   }
   tree->GetEntry(event);

   cout<< "tc1 = " << tc1 << " tc2 = " << tc2 << " tc3 = " << tc3 <<endl;

   cout<< "*tc = " << *tc << " x[" << *tc << "] = " << x[*tc] <<endl;
   
   Int_t chan = 1;
   TGraph* g1 = new TGraph(1024,x,y1);
   g1->SetNameTitle(Form("g%d_evt%d",chan,event), Form("board %d chan %d event %d",board,chan,event));
   g1->SetMarkerStyle(6);
   g1->SetMarkerColor(2);

   new TCanvas;
   g1->Draw("ap");

   Int_t cell0 = 1023 - *tc;
   cout<< "cell0 = " << cell0 <<endl;

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
      Int_t sign = 0;

      /// //
      /// // remove single spike
      /// //
      /// sign = y1[curr1] > y1[prev1]? 1: -1;
      /// if (sign*(y1[curr1] - y1[prev1]) > sign*diff_high && sign*(y1[curr1] - y1[next1]) > sign*diff_high) y1[curr1] = 0.5*(y1[prev1]+y1[next1]);
      /// sign = y2[curr1] > y2[prev1]? 1: -1;                                      
      /// if (sign*(y2[curr1] - y2[prev1]) > sign*diff_high && sign*(y2[curr1] - y2[next1]) > sign*diff_high) y2[curr1] = 0.5*(y2[prev1]+y2[next1]);
      /// sign = y3[curr1] > y3[prev1]? 1: -1;                                      
      /// if (sign*(y3[curr1] - y3[prev1]) > sign*diff_high && sign*(y3[curr1] - y3[next1]) > sign*diff_high) y3[curr1] = 0.5*(y3[prev1]+y3[next1]);
      /// sign = y4[curr1] > y4[prev1]? 1: -1;                                      
      /// if (sign*(y4[curr1] - y4[prev1]) > sign*diff_high && sign*(y4[curr1] - y4[next1]) > sign*diff_high) y4[curr1] = 0.5*(y4[prev1]+y4[next1]);

      //
      // remove double spikes
      //
      // problem event 299 from the file tb2014_run_082.root: negative spikes survive the procedure. Solved now. 
      // also unclear with event 302

      sign = y1[curr1] > y1[prev1]? 1: -1;

      cout<< "curr1 = " << curr1 << " prev1 = " << prev1 << " next1 = " << next1 << " nextnext1 = " << nextnext1
         << "   x[curr1] = " << x[curr1]
         << " y1[curr1] = " << y1[curr1] << " y1[prev1] = " << y1[prev1]
         << " y1[curr1]-y1[prev1] = " << y1[curr1] - y1[prev1]
         << " y1[next1]-y1[nextnext1] = " << y1[next1] - y1[nextnext1]
         << " sign = " << sign <<endl;

      if (sign*(y1[curr1] - y1[prev1]) > diff_high && sign*(y1[next1] - y1[nextnext1]) > diff_high)
      {
         // cout<< "x[curr1] = " << x[curr1] << " y1[curr1] - y1[prev1] = " << y1[curr1] - y1[prev1] << " sign = " << sign <<endl;
         cout<< "x[curr1] = " << x[curr1] << " y1[curr1]-y1[prev1] = " << y1[curr1] - y1[prev1] << " sign = " << sign <<endl;
         // look at channel symmetric wrt cell0
         Int_t dspike = cell0 - curr1;
         Int_t curr2 = (cell0 + dspike) % 1024;    if (curr2 < 0) curr2 += 1024;
         Int_t prev2 = curr2 - 1;                  if (prev2 == -1) prev2 = 1023;
         Int_t next2 = (curr2 + 1) % 1024;
         Int_t nextnext2 = (curr2 + 2) % 1024;

         cout<< "sign*(y1[curr2]-y1[prev2]) = " << sign*(y1[curr2] - y1[prev2]) << " sign*(y1[next2]-y1[nextnext2]) = " << sign*(y1[next2] - y1[nextnext2]) <<endl;
         if (sign*(y1[curr2] - y1[prev2]) > diff_low && sign*(y1[next2] - y1[nextnext2]) > diff_low)
         {
            cout<< "spikes are at curr1 = " << curr1 << " curr2 = " << curr2 <<endl;
            cout<< "curr2 = " << curr2 << " prev2 = " << prev2 << " next2 = " << next2 << " nextnext2 = " << nextnext2 <<endl;

            Double_t diff = diff_med;
            //-- Double_t diff = diff_low;
            if (sign*(y2[curr1] - y2[prev1]) > diff && sign*(y2[next1] - y2[nextnext1]) > diff && sign*(y2[curr2] - y2[prev2]) > diff && sign*(y2[next2] - y2[nextnext2]) > diff) cout<< "y2" <<endl;
            cout
               << "   x[curr1] = " << x[curr1]
               << " y2[curr1] = " << y2[curr1] << " y2[prev1] = " << y2[prev1]
               << " y2[curr1]-y2[prev1] = " << y2[curr1] - y2[prev1]
               << " y2[next1]-y2[nextnext1] = " << y2[next1] - y2[nextnext1]
               << " y2[curr2]-y2[prev2] = " << y2[curr2] - y2[prev2] << " y2[next2]-y2[nextnext2] = " << y2[next2] - y2[nextnext2] <<endl;
            if (sign*(y3[curr1] - y3[prev1]) > diff && sign*(y3[next1] - y3[nextnext1]) > diff && sign*(y3[curr2] - y3[prev2]) > diff && sign*(y3[next2] - y3[nextnext2]) > diff) cout<< "y3" <<endl;
            cout
               << "   x[curr1] = " << x[curr1]
               << " y3[curr1] = " << y3[curr1] << " y3[prev1] = " << y3[prev1]
               << " y3[curr1]-y3[prev1] = " << y3[curr1] - y3[prev1]
               << " y3[next1]-y3[nextnext1] = " << y3[next1] - y3[nextnext1]
               << " y3[curr2]-y3[prev2] = " << y3[curr2] - y3[prev2] << " y3[next2]-y3[nextnext2] = " << y3[next2] - y3[nextnext2] <<endl;
            if (sign*(y4[curr1] - y4[prev1]) > diff && sign*(y4[next1] - y4[nextnext1]) > diff && sign*(y4[curr2] - y4[prev2]) > diff && sign*(y4[next2] - y4[nextnext2]) > diff) cout<< "y4" <<endl;
            cout
               << "   x[curr1] = " << x[curr1]
               << " y4[curr1] = " << y4[curr1] << " y4[prev1] = " << y4[prev1]
               << " y4[curr1]-y4[prev1] = " << y4[curr1] - y4[prev1]
               << " y4[next1]-y4[nextnext1] = " << y4[next1] - y4[nextnext1]
               << " y4[curr2]-y4[prev2] = " << y4[curr2] - y4[prev2] << " y4[next2]-y4[nextnext2] = " << y4[next2] - y4[nextnext2] <<endl;
            cout<<endl;
            Int_t nothers = 0;
            if (sign*(y2[curr1] - y2[prev1]) > diff && sign*(y2[next1] - y2[nextnext1]) > diff && sign*(y2[curr2] - y2[prev2]) > diff && sign*(y2[next2] - y2[nextnext2]) > diff) nothers++;
            if (sign*(y3[curr1] - y3[prev1]) > diff && sign*(y3[next1] - y3[nextnext1]) > diff && sign*(y3[curr2] - y3[prev2]) > diff && sign*(y3[next2] - y3[nextnext2]) > diff) nothers++;
            if (sign*(y4[curr1] - y4[prev1]) > diff && sign*(y4[next1] - y4[nextnext1]) > diff && sign*(y4[curr2] - y4[prev2]) > diff && sign*(y4[next2] - y4[nextnext2]) > diff) nothers++;

            if (nothers > 0) {
               spikes[curr1] += sign;
               spikes[curr2] += sign;
            }
         }
      }
   }

   for (std::map<Int_t,Int_t>::const_iterator it=spikes.begin(); it!=spikes.end(); ++it) {
      cout<< "it->first = " << it->first << " it->second = " << it->second <<endl;
      Int_t curr = it->first;
      Int_t next = curr + 1;     if (next > 1023) next = 0;
      Double_t sign = it->second > 0? 1: -1;
      //gc1->SetPoint(curr, gc1->GetX()[curr], gc1->GetY()[curr] - sign*14.8);
      //gc1->SetPoint(next, gc1->GetX()[next], gc1->GetY()[next] - sign*14.8);
      y1[curr] -= sign*14.8;
      y1[next] -= sign*14.8;
      y2[curr] -= sign*14.8;
      y2[next] -= sign*14.8;
      y3[curr] -= sign*14.8;
      y3[next] -= sign*14.8;
      y4[curr] -= sign*14.8;
      y4[next] -= sign*14.8;
   }

   //
   // remove single spike
   //
   for (int i=0; i<1024; ++i)
   {
      Int_t curr1 = i;
      Int_t prev1 = (curr1 - 1);  if (prev1 == -1) prev1 = 1023;
      Int_t next1 = (curr1 + 1) % 1024;
      //Int_t nextnext1 = (curr1 + 2) % 1024;
      Int_t sign = 0;

      sign = y1[curr1] > y1[prev1]? 1: -1;
      if (sign*(y1[curr1] - y1[prev1]) > sign*diff_high && sign*(y1[curr1] - y1[next1]) > sign*diff_high) y1[curr1] = 0.5*(y1[prev1]+y1[next1]);
      sign = y2[curr1] > y2[prev1]? 1: -1;                                      
      if (sign*(y2[curr1] - y2[prev1]) > sign*diff_high && sign*(y2[curr1] - y2[next1]) > sign*diff_high) y2[curr1] = 0.5*(y2[prev1]+y2[next1]);
      sign = y3[curr1] > y3[prev1]? 1: -1;                                      
      if (sign*(y3[curr1] - y3[prev1]) > sign*diff_high && sign*(y3[curr1] - y3[next1]) > sign*diff_high) y3[curr1] = 0.5*(y3[prev1]+y3[next1]);
      sign = y4[curr1] > y4[prev1]? 1: -1;                                      
      if (sign*(y4[curr1] - y4[prev1]) > sign*diff_high && sign*(y4[curr1] - y4[next1]) > sign*diff_high) y4[curr1] = 0.5*(y4[prev1]+y4[next1]);
   }

   TGraph* gc1 = new TGraph(1024, x, y1);
   gc1->SetNameTitle(Form("g%d_evt%d",chan,event), Form("Corrected board %d chan %d event %d",board,chan,event));
   gc1->SetMarkerStyle(6);
   gc1->SetMarkerColor(2);

   new TCanvas;
   gc1->Draw("ap");

   tree->ResetBranchAddresses();
}

void pulseclean(const char* ifname)
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

   Int_t event, tc1, tc2, tc3;
   Float_t b1_t[1024], b1_c[4096], b2_t[1024], b2_c[4096], b3_t[1024], b3_c[4096];
   //-- Int_t tc[3];
   Float_t* b_t[3];
   Float_t* b_c[3][4];
   b_t[0] = b1_t;
   b_t[1] = b2_t;
   b_t[2] = b3_t;
   b_c[0][0] = b1_c;
   b_c[0][1] = b1_c + 1024;
   b_c[0][2] = b1_c + 2048;
   b_c[0][3] = b1_c + 3072;
   b_c[1][0] = b2_c;
   b_c[1][1] = b2_c + 1024;
   b_c[1][2] = b2_c + 2048;
   b_c[1][3] = b2_c + 3072;
   b_c[2][0] = b3_c;
   b_c[2][1] = b3_c + 1024;
   b_c[2][2] = b3_c + 2048;
   b_c[2][3] = b3_c + 3072;

   Float_t* b1_c1 = b1_c;
   Float_t* b1_c2 = b1_c + 1024;
   Float_t* b1_c3 = b1_c + 2048;
   Float_t* b1_c4 = b1_c + 3072;
   Float_t* b2_c1 = b2_c;
   Float_t* b2_c2 = b2_c + 1024;
   Float_t* b2_c3 = b2_c + 2048;
   Float_t* b2_c4 = b2_c + 3072;
   Float_t* b3_c1 = b3_c;
   Float_t* b3_c2 = b3_c + 1024;
   Float_t* b3_c3 = b3_c + 2048;
   Float_t* b3_c4 = b3_c + 3072;

   itree->SetBranchAddress("event", &event);
   itree->SetBranchAddress("tc1", &tc1);
   itree->SetBranchAddress("tc2", &tc2);
   itree->SetBranchAddress("tc3", &tc3);
   itree->SetBranchAddress("b1_t", &b1_t);
   itree->SetBranchAddress("b2_t", &b2_t);
   itree->SetBranchAddress("b3_t", &b3_t);
   itree->SetBranchAddress("b1_c", &b1_c);
   itree->SetBranchAddress("b2_c", &b2_c);
   itree->SetBranchAddress("b3_c", &b3_c);

   //
   // output tree uses the same variables
   //

   TFile* ofile = TFile::Open(Form("%s.pulse.root",ifname), "recreate");

   TTree* otree = new TTree("p", "new pulse tree");
   otree->SetMarkerStyle(6);
   otree->SetMarkerColor(2);
   otree->SetLineColor(2);

   otree->Branch("event", &event, "event/I");
   otree->Branch("tc1", &tc1, "tc1/I");
   otree->Branch("t1",  b1_t, "t1[1024]/F");
   otree->Branch("c1", b1_c1, "c1[1024]/F");
   otree->Branch("c2", b1_c2, "c2[1024]/F");
   otree->Branch("c3", b1_c3, "c3[1024]/F");
   otree->Branch("c4", b1_c4, "c4[1024]/F");
   otree->Branch("tc2", &tc2, "tc2/I");
   otree->Branch("t2",  b2_t, "t2[1024]/F");
   otree->Branch("c5", b2_c1, "c5[1024]/F");
   otree->Branch("c6", b2_c2, "c6[1024]/F");
   otree->Branch("c7", b2_c3, "c7[1024]/F");
   otree->Branch("c8", b2_c4, "c8[1024]/F");
   otree->Branch("tc3", &tc3, "tc3/I");
   otree->Branch("t3",  b3_t, "t3[1024]/F");
   otree->Branch("c9", b3_c1, "c9[1024]/F");
   otree->Branch("c10", b3_c2, "c10[1024]/F");
   otree->Branch("c11", b3_c3, "c11[1024]/F");
   otree->Branch("c12", b3_c4, "c12[1024]/F");

   ///////////////////////////////////////////
   //
   //       Spike removing
   //
   ///////////////////////////////////////////

   Float_t* x = 0;
   Float_t* y1 = 0;
   Float_t* y2 = 0;
   Float_t* y3 = 0;
   Float_t* y4 = 0;

   Int_t* tc = 0;

   for (Long64_t jentry=0; jentry<itree->GetEntries(); jentry++)
   {
      if (jentry % 1000 == 0) cout<< "jentry = " << jentry <<endl;
      itree->LoadTree(jentry);
      itree->GetEntry(jentry);

      for (int iboard=1; iboard<=3; ++ iboard)
      {
         switch (iboard) {
            case 1: tc = &tc1;
                    x = b1_t;
                    y1 = b1_c1;
                    y2 = b1_c2;
                    y3 = b1_c3;
                    y4 = b1_c4;
                    break;
            case 2: tc = &tc2;
                    x = b2_t;
                    y1 = b2_c1;
                    y2 = b2_c2;
                    y3 = b2_c3;
                    y4 = b2_c4;
                    break;
            case 3: tc = &tc3;
                    x = b3_t;
                    y1 = b3_c1;
                    y2 = b3_c2;
                    y3 = b3_c3;
                    y4 = b3_c4;
                    break;
            default: cout<< "iboard " << iboard << " is out of range 1..3" <<endl;
                     return;
         }

         Int_t cell0 = 1023 - *tc;

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

            Int_t sign = y1[curr1] > y1[prev1]? 1: -1;

            if (sign*(y1[curr1] - y1[prev1]) > diff_high && sign*(y1[next1] - y1[nextnext1]) > diff_high)   // diff_high to find a candidate
            {
               // look at channel symmetric wrt cell0
               Int_t dspike = cell0 - curr1;             // NB: dspike may be negative
               Int_t curr2 = (cell0 + dspike) % 1024;    if (curr2 < 0) curr2 += 1024;
               Int_t prev2 = curr2 - 1;                  if (prev2 == -1) prev2 = 1023;
               Int_t next2 = (curr2 + 1) % 1024;
               Int_t nextnext2 = (curr2 + 2) % 1024;

               if (sign*(y1[curr2] - y1[prev2]) > diff_low && sign*(y1[next2] - y1[nextnext2]) > diff_low)  // diff_low to find symmetric wrt cell0
               {
                  Double_t diff = diff_med;                                                                 // diff_med to find confirmation in the rest of the channels
                  Int_t nothers = 0;
                  if (sign*(y2[curr1] - y2[prev1]) > diff && sign*(y2[next1] - y2[nextnext1]) > diff && sign*(y2[curr2] - y2[prev2]) > diff && sign*(y2[next2] - y2[nextnext2]) > diff) nothers++;
                  if (sign*(y3[curr1] - y3[prev1]) > diff && sign*(y3[next1] - y3[nextnext1]) > diff && sign*(y3[curr2] - y3[prev2]) > diff && sign*(y3[next2] - y3[nextnext2]) > diff) nothers++;
                  if (sign*(y4[curr1] - y4[prev1]) > diff && sign*(y4[next1] - y4[nextnext1]) > diff && sign*(y4[curr2] - y4[prev2]) > diff && sign*(y4[next2] - y4[nextnext2]) > diff) nothers++;

                  if (nothers > 0) {
                     spikes[curr1] += sign;
                     spikes[curr2] += sign;
                  }
               }
            }
         }

         for (std::map<Int_t,Int_t>::const_iterator it=spikes.begin(); it!=spikes.end(); ++it) {
            Int_t curr = it->first;
            Int_t next = curr + 1;     if (next > 1023) next = 0;
            Double_t sign = it->second > 0? 1: -1;
            y1[curr] -= sign*14.8;
            y1[next] -= sign*14.8;
            y2[curr] -= sign*14.8;
            y2[next] -= sign*14.8;
            y3[curr] -= sign*14.8;
            y3[next] -= sign*14.8;
            y4[curr] -= sign*14.8;
            y4[next] -= sign*14.8;
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

            sign = y1[curr1] > y1[prev1]? 1: -1;
            if (sign*(y1[curr1] - y1[prev1]) > sign*diff_high && sign*(y1[curr1] - y1[next1]) > sign*diff_high) y1[curr1] = 0.5*(y1[prev1]+y1[next1]);
            sign = y2[curr1] > y2[prev1]? 1: -1;                                      
            if (sign*(y2[curr1] - y2[prev1]) > sign*diff_high && sign*(y2[curr1] - y2[next1]) > sign*diff_high) y2[curr1] = 0.5*(y2[prev1]+y2[next1]);
            sign = y3[curr1] > y3[prev1]? 1: -1;                                      
            if (sign*(y3[curr1] - y3[prev1]) > sign*diff_high && sign*(y3[curr1] - y3[next1]) > sign*diff_high) y3[curr1] = 0.5*(y3[prev1]+y3[next1]);
            sign = y4[curr1] > y4[prev1]? 1: -1;                                      
            if (sign*(y4[curr1] - y4[prev1]) > sign*diff_high && sign*(y4[curr1] - y4[next1]) > sign*diff_high) y4[curr1] = 0.5*(y4[prev1]+y4[next1]);
         }
      }

      otree->Fill();
   }

   cout<< "Write " << otree->GetEntries() << " events into file " << ofile->GetName() <<endl;
   ofile->Write();

   itree->ResetBranchAddresses();
}
