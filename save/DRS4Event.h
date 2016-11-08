#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TList.h>

#include <iostream>

using std::cout;		using std::endl;

#ifndef DRS4Event_h
#define DRS4Event_h

class DRS4Event {
public:
   TTree* tree;
   TList otree_list;       // list of output tree(s) which use the same buffer as input tree. Must be released in destructor.

   struct BoardBranch {
      TBranch* triggerCell;
      TBranch* t;
      TBranch* c[4];
      TBranch* c4096;
   };
   struct Board {
      Int_t triggerCell;
      Float_t t[1024];
      Float_t* c[4];
      Float_t buf[4096];
      Board() {for (int ich=0; ich<4; ++ich) c[ich] = buf + ich*1024;}
   };
   // buffers for three boards
   static const Int_t nboards = 3;
   BoardBranch branch[nboards];
   Board b[nboards];
public:
   void clear() {
      for (int iboard=0; iboard<nboards; ++iboard) {
         branch[iboard].triggerCell = 0;
         branch[iboard].t = 0;
         for (int ich=0; ich<4; ++ich) branch[iboard].c[ich] = 0;

         b[iboard].triggerCell = 0;
         for (int i=0; i<1024; ++i) b[iboard].t[i] = 0;
         for (int ich=0; ich<4; ++ich) for (int i=0; i<1024; ++i) b[iboard].c[ich][i] = 0;
      }
   }
   DRS4Event(TTree*& the_tree): tree(the_tree) {
      clear();
      if (!tree) tree = (TTree*) gDirectory->Get("p");
      if (!tree) tree = (TTree*) gDirectory->Get("pulse");
      if (!tree) {
         cout<< "Could not find pulse tree in the current directory" <<endl;
         return;
      }
      the_tree = tree;
      connect();
   }

   // operations with tree
   TTree* GetTree() const {return tree;}
   // Long64_t LoadTree(Long64_t entry) {return tree->LoadTree(entry);}
   Long64_t GetEntry(Long64_t entry, Int_t getall=0) {
     if (tree->LoadTree(entry) < 0) return 0;
     return tree->GetEntry(entry, getall);
   }

   // get the data
   Int_t triggerCell(Int_t board) {
      if (branch[board].t) return b[board].triggerCell;
      else {
         cout<< "***DRS4Data::triggerCell: board not found: " << board <<endl;
         return 0;
      }
   }
   Float_t* boardT(Int_t board) {
      if (branch[board].t) return b[board].t;
      else {
         // cout<< "***DRS4Data::T: board not found: " << board <<endl;
         return 0;
      }
   }
   Float_t* chanT(Int_t chan) {                       // chan = (0,11)
      Int_t board = chan/4;
      if (branch[board].t) return b[board].t;
      else {
         // cout<< "***DRS4Data::T: board not found: " << board <<endl;
         return 0;
      }
   }
   Float_t* boardV(Int_t chan, Int_t board) {         // chan = (0,11)
      if (board < 0) board = chan/4;
      chan = chan % 4;

      if (branch[board].t == 0) {                           // check the time branch
         // cout<< "***DRS4Event::V: board not found: " << board <<endl;
         return 0;
      }
      if (branch[board].c4096);                             // no channel information
      else {
         if (branch[board].c[chan] == 0) {                  // check the channel branch
            // cout<< "DRS4Event::V: No such branch: board " << board << " channel " << chan+1 <<endl;
            return 0;
         }
      }
      return b[board].c[chan];
   }
   Float_t* chanV(Int_t chan) {                       // chan = (0,11)
      Int_t board = chan/4;
      chan = chan % 4;

      if (branch[board].t == 0) {                           // check the time branch
         // cout<< "***DRS4Event::V: board not found: " << board <<endl;
         return 0;
      }
      if (branch[board].c4096);                             // no channel information
      else {
         if (branch[board].c[chan] == 0) {                  // check the channel branch
            // cout<< "DRS4Event::V: No such branch: board " << board << " channel " << chan+1 <<endl;
            return 0;
         }
      }
      return b[board].c[chan];
   }

   void uniform(Int_t N, const Float_t* x, const Float_t* y, Double_t* xuni, Double_t* yuni)
   {
     // interpolate to the tempered (униформ) scale

     Double_t step = (x[N-1] - x[0]) / (N-1);

     xuni[0] = x[0];
     yuni[0] = y[0];

     Int_t idrs = 1;
     Int_t iuni = 1;
     while (iuni < N) {
       xuni[iuni] = xuni[iuni-1] + step;
       while (x[idrs] < xuni[iuni] && idrs < N-1) idrs++;

       Double_t dx_drs = x[idrs] - x[idrs-1];
       Double_t dy_drs = y[idrs] - y[idrs-1];
       Double_t d = xuni[iuni] - x[idrs-1];
       yuni[iuni] = y[idrs-1] + d*dy_drs/dx_drs;
       ++iuni;
     }
   }
   TH1D* unihist(Int_t N, const Float_t* x, const Float_t* y)
   {
      // interpolate to the tempered (униформ) scale

      Float_t step = (x[N-1] - x[0]) / (N-1);
      TH1D* h = new TH1D("h_unihist","unihist", N, x[0]-0.5*step, x[N-1]+0.5*step);
      h->SetDirectory(0);

      Int_t idrs = 1;
      Int_t iuni = 1;
      while (iuni <= N) {
         while (x[idrs] < h->GetBinCenter(iuni) && idrs <= N-1) idrs++;

         Double_t dx_drs = x[idrs] - x[idrs-1];
         Double_t dy_drs = y[idrs] - y[idrs-1];
         Double_t d = h->GetBinCenter(iuni) - x[idrs-1];
         h->SetBinContent(iuni, y[idrs-1] + d*dy_drs/dx_drs);
         //cout<< "iuni = " << iuni << " h->GetBinCenter(iuni) = " << h->GetBinCenter(iuni) << " h->GetBinContent(iuni) = " << h->GetBinContent(iuni) <<endl;
         ++iuni;
      }
      return h;
   }

   void release() {tree->ResetBranchAddresses();}
   ~DRS4Event() {
      // cout<< "destructor ~DRS4Event" <<endl;
      // cout<< "otree_list.GetEntries() = " << otree_list.GetEntries() <<endl;
      TIter next(&otree_list);
      TTree* otree = 0;
      while ((otree = (TTree*) next())) {
         // cout<< "otree->GetName() = " << otree->GetName() <<endl;
         otree->ResetBranchAddresses();
         // cout<< "otree->GetName() = " << otree->GetName() << " has been released" <<endl;
      }
      // cout<< "release the tree" <<endl;
      release();
      // cout<< "destructor ~DRS4Event: done" <<endl;
   }
   void connect()
   {
      // old tree "pulse"
      if (strcmp(tree->GetName(), "pulse") == 0) {
         // the branches tc1 and b1_t exist in both versions
         for (int iboard=0; iboard<nboards; ++iboard) {
            if ((branch[iboard].triggerCell = tree->GetBranch(Form("tc%d",iboard+1)))) tree->SetBranchAddress(Form("tc%d",iboard+1), &b[iboard].triggerCell);
            if ((branch[iboard].t = tree->GetBranch(Form("b%d_t",iboard+1)))) tree->SetBranchAddress(Form("b%d_t",iboard+1), &b[iboard].t);
         }

         // last version of the Heejong tree: the name b1_c is unique while the name b1_t was used in previous verison
         for (int iboard=0; iboard<nboards; ++iboard) {
            if ((branch[iboard].c4096 = tree->GetBranch(Form("b%d_c",iboard+1)))) tree->SetBranchAddress(Form("b%d_c",iboard+1), &b[iboard].buf);
         }
         // previous version of the Heejong tree
         for (int iboard=0; iboard<nboards; ++iboard) {
            for (int ich=0; ich<4; ++ich) {
               if ((branch[iboard].c[ich] = tree->GetBranch(Form("b%d_c%d",iboard+1,ich+1)))) tree->SetBranchAddress(Form("b%d_c%d",iboard+1,ich+1), b[iboard].c[ich]);
            }
         }
      }

      // new tree "p"
      if (strcmp(tree->GetName(), "p") == 0) {
         for (int iboard=0; iboard<nboards; ++iboard) {
            if ((branch[iboard].triggerCell = tree->GetBranch(Form("tc%d",iboard+1)))) tree->SetBranchAddress(Form("tc%d",iboard+1), &b[iboard].triggerCell);
            if ((branch[iboard].t = tree->GetBranch(Form("t%d",iboard+1)))) tree->SetBranchAddress(Form("t%d",iboard+1), &b[iboard].t);
         }
         for (int ich=0; ich<12; ++ich) {
            if (tree->GetBranch(Form("c%d",ich+1))) {
               int board = ich/4;
               int ich_local = ich % 4;
               //cout<< "SetBranchAddress for " << Form("c%d",ich+1) << " board = " <<  board << " ich_local = " << ich_local <<endl;
               if ((branch[board].c[ich_local] = tree->GetBranch(Form("c%d",ich+1)))) tree->SetBranchAddress(Form("c%d",ich+1), b[board].c[ich_local]);
            }
         }
      }

      // for (int iboard=0; iboard<nboards; ++iboard) {
      //    cout<< "branch[" << iboard << "].c4096 = " << branch[iboard].c4096 <<endl;
      //    for (int ich=0; ich<4; ++ich) {
      //       cout<< "branch[" << iboard << "].c[" << ich << "] = " << branch[iboard].c[ich] <<endl;
      //    }
      // }
   }
   void bookSyncTree(TTree* otree)       // cannot be const because of TTree::Branch
   {
      // book output tree with the same buffers
      // NB: the destructor must release both trees
      otree_list.Add(otree);
      for (int iboard=0; iboard<nboards; ++iboard) {
         if (branch[iboard].triggerCell) otree->Branch(Form("tc%d",iboard+1), &b[iboard].triggerCell, Form("tc%d/I",iboard+1));
         if (branch[iboard].t) otree->Branch(Form("t%d",iboard+1), &b[iboard].t, Form("t%d[1024]/F",iboard+1));
         for (int ich=0; ich<4; ++ich) {
            Int_t channel = iboard*4 + ich + 1;
            if (branch[iboard].c4096) otree->Branch(Form("c%d",channel), b[iboard].c[ich%4], Form("c%d[1024]/F",channel));   // we don't know exact channels, assume all
            else if (branch[iboard].c[ich%4]) otree->Branch(Form("c%d",channel), b[iboard].c[ich%4], Form("c%d[1024]/F",channel));
         }
      }

      // show branches
      cout<< "branches of the output tree " << tree->GetName() <<endl;
      TIter next(otree->GetListOfBranches());
      TBranch* bra;
      while ((bra = (TBranch*) next())) {
         cout<< bra->GetName() << "\t address = " << (TBranch*) bra->GetAddress() <<endl;
      }
   }
   TGraph* plot(Int_t event, Int_t chan, Int_t board=-1)
   {
      if (board < 0) board = chan/4;
      chan = chan % 4;

      if (branch[board].t == 0) {                           // check the time branch
         cout<< "DRS4Event::plot: board not found: " << board <<endl;
         return 0;
      }
      if (branch[board].c4096);                             // no channel information
      else {
         if (branch[board].c[chan] == 0) {                  // check the channel branch
            cout<< "No such branch: board " << board << " channel " << chan+1 <<endl;
            return 0;
         }
      }

      // if (branch[board].c4096 == 0 && branch[board].c[chan] == 0) {
      //    cout<< "No such branch: board " << board << " channel " << chan+1 <<endl;
      //    return 0;
      // }

      if (tree->LoadTree(event) < 0) {
         cout<< "DRS4Event::plot: event " << event << " is out of range (0, " << tree->GetEntries()-1 << ")" <<endl;
         return 0;
      }
      tree->GetEntry(event);

      TGraph* graph = new TGraph(1024, b[board].t, b[board].c[chan]);
      graph->SetNameTitle(Form("b_%d_evt_%d_ch_%d",board,event,chan+1), Form("board %d evt %d ch %d",board,event,chan+1));
      graph->SetMarkerStyle(7);
      graph->SetMarkerColor(46);
      graph->SetLineColor(46);

      new TCanvas;
      graph->Draw("apl");

      return graph;
   }
};

#endif   // DRS4Event_h
