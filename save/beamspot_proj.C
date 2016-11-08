// Andriy Zatserklyaniy, April 17, 2014

#include "Reco.h"
#include "DataFormat.h"

#include <TCanvas.h>
#include <TH2.h>
#include <TFile.h>

namespace Tree {
   Bool_t ok;
   Double_t vbeamspotIn;
   Double_t vbeamspotOut;
   Double_t vreal[4];
   Double_t vproj[4];
   Double_t tbeamspotIn;
   Double_t tbeamspotOut;
   Double_t treal[4];
   Double_t tproj[4];

   void clear() {
      ok = kFALSE;
      vbeamspotIn = 0;
      vbeamspotOut = 0;
      tbeamspotIn = 0;
      tbeamspotOut = 0;
      for (int ilayer=0; ilayer<4; ++ilayer) {
         vreal[ilayer] = 0;
         vproj[ilayer] = 0;
         treal[ilayer] = 0;
         tproj[ilayer] = 0;
      }
   }
   void book(TTree* tree) {
      tree->Branch("ok", &ok, "ok/B");
      tree->Branch("vbeamspotIn", &vbeamspotIn, "vbeamspotIn/D");
      tree->Branch("vbeamspotOut", &vbeamspotOut, "vbeamspotOut/D");
      tree->Branch("vreal", &vreal, "vreal[4]/D");
      tree->Branch("vproj", &vproj, "vproj[4]/D");
      tree->Branch("tbeamspotIn", &tbeamspotIn, "tbeamspotIn/D");
      tree->Branch("tbeamspotOut", &tbeamspotOut, "tbeamspotOut/D");
      tree->Branch("treal", &treal, "treal[4]/D");
      tree->Branch("tproj", &tproj, "tproj[4]/D");
   }
   void connect(TTree* tree)                                // need for event-by-event analysis
   {  
      // connects tree buffers with variables to use for event-by-event analysis
      tree->SetBranchAddress("ok",         &ok);
      tree->SetBranchAddress("vbeamspotIn",         &vbeamspotIn);
      tree->SetBranchAddress("vbeamspotIn",         &vbeamspotIn);
      tree->SetBranchAddress("vreal",             &vreal);
      tree->SetBranchAddress("vproj",             &vproj);
      tree->SetBranchAddress("tbeamspotIn",         &tbeamspotIn);
      tree->SetBranchAddress("tbeamspotOut",         &tbeamspotOut);
      tree->SetBranchAddress("treal",             &treal);
      tree->SetBranchAddress("tproj",             &tproj);
   }
}  // namespace Tree

void beamspot_proj(Int_t event1=0, Int_t event2=-1, TTree* tree=0)
{
   Bool_t debug = kFALSE;
   if (debug) cout<< "debug in on" <<endl;

   if (!tree) tree = (TTree*) gDirectory->Get("r");
   if (!tree) {
      cout<< "Could not find tree r" <<endl;
      return;
   }

   // use runHeader to get the run number
   RunHeader* runHeader = (RunHeader*) tree->GetUserInfo()->First();
   cout<< "runHeader->GetRun() = " << runHeader->GetRun() <<endl;
   time_t start_time = runHeader->GetTime();
   cout<< "run start time: " << std::ctime(&start_time);
   cout<< "program version is " << runHeader->GetVersion() <<endl;
   if (runHeader->GetTimeTag()) cout<< "event time tag was written out" <<endl;
   else cout<< "event time tag was not written out" <<endl;
   cout<<endl;

   const RecoEvent* recoEvent = 0;
   tree->SetBranchAddress("revent", &recoEvent);

   TFile* ofile = new TFile(Form("%s-proj.root",gDirectory->GetName()), "recreate");
   TTree* otree = new TTree("proj", "Track projections");
   otree->SetMarkerColor(602);
   otree->SetLineColor(602);

   Tree::book(otree);

   if (event2 < event1) event2 = tree->GetEntries()-1;

   for (int ientry=event1; ientry<=event2; ++ientry)
   {
      if (tree->LoadTree(ientry) < 0) {
         cout<< "Could not load event " << ientry <<endl;
         break;
      }
      tree->GetEntry(ientry);

      if (false
          || (ientry < 10)
          || (ientry < 10000 && ientry%1000 == 0)
          || (ientry%100000 == 0)
         ) cout<< "---> processing entry " << ientry <<endl;

      Tree::clear();

      if (recoEvent->track->GetLast()+1 == 1)   // require one reconstructed track
      {
         const SuperTrack* superTrack = (const SuperTrack*) recoEvent->track->At(0);

         Tree::vbeamspotIn = superTrack->Vfront(-3500.);
         Tree::vbeamspotOut = superTrack->Vrear(-3500.);
         Tree::tbeamspotIn = superTrack->Tfront(-3500.);
         Tree::tbeamspotOut = superTrack->Trear(-3500.);

         // all hits should be real
         if (true
             && superTrack->vTrack_->itrack2D_->hit1_->sensorId_ > 0
             && superTrack->vTrack_->itrack2D_->hit2_->sensorId_ > 0
             && superTrack->vTrack_->otrack2D_->hit1_->sensorId_ > 0
             && superTrack->vTrack_->otrack2D_->hit2_->sensorId_ > 0
             && superTrack->tTrack_->itrack2D_->hit1_->sensorId_ > 0
             && superTrack->tTrack_->itrack2D_->hit2_->sensorId_ > 0
             && superTrack->tTrack_->otrack2D_->hit1_->sensorId_ > 0
             && superTrack->tTrack_->otrack2D_->hit2_->sensorId_ > 0
            ) {

            Tree::ok = kTRUE;

            // u-positions of the layers
            Double_t uv[4];
            Double_t ut[4];

            uv[0] = superTrack->vTrack_->itrack2D_->hit1_->u_;
            uv[1] = superTrack->vTrack_->itrack2D_->hit2_->u_;
            uv[2] = superTrack->vTrack_->otrack2D_->hit1_->u_;
            uv[3] = superTrack->vTrack_->otrack2D_->hit2_->u_;
            ut[0] = superTrack->vTrack_->itrack2D_->hit1_->u_;
            ut[1] = superTrack->vTrack_->itrack2D_->hit2_->u_;
            ut[2] = superTrack->vTrack_->otrack2D_->hit1_->u_;
            ut[3] = superTrack->vTrack_->otrack2D_->hit2_->u_;

            Tree::vreal[0] = superTrack->vTrack_->itrack2D_->hit1_->pos_;
            Tree::vreal[1] = superTrack->vTrack_->itrack2D_->hit2_->pos_;
            Tree::vreal[2] = superTrack->vTrack_->otrack2D_->hit1_->pos_;
            Tree::vreal[3] = superTrack->vTrack_->otrack2D_->hit2_->pos_;
            Tree::treal[0] = superTrack->tTrack_->itrack2D_->hit1_->pos_;
            Tree::treal[1] = superTrack->tTrack_->itrack2D_->hit2_->pos_;
            Tree::treal[2] = superTrack->tTrack_->otrack2D_->hit1_->pos_;
            Tree::treal[3] = superTrack->tTrack_->otrack2D_->hit2_->pos_;

            // projections

            Double_t RbeamSpot = 0.5*(14.16 + 17.08);                         // Sep2014 beam test
            BeamSpot beamSpotIn(35.52,-1.07,-3500, RbeamSpot);                // Sep2014 beam test
            BeamSpot beamSpotOut(13.76,-6.96,-3500, 140.);                    // Sep2014 beam test

            // vertex for input telescope
            SensorHit vertexHitIn_t(0, 0, 0, beamSpotIn.u_, beamSpotIn.t_);          // vertex hit for recovery of the t-hits
            SensorHit vertexHitIn_v(0, 0, 0, beamSpotIn.u_, beamSpotIn.v_);          // vertex hit for recovery of the v-hits
            // vertex for output telescope
            SensorHit vertexHitOut_t(0, 0, 0, beamSpotOut.u_, beamSpotOut.t_);          // vertex hit for recovery of the t-hits
            SensorHit vertexHitOut_v(0, 0, 0, beamSpotOut.u_, beamSpotOut.v_);          // vertex hit for recovery of the v-hits

            // vertex track
            Track2D track_v0(&vertexHitIn_v, superTrack->vTrack_->itrack2D_->hit1_);  // track through v0 hit
            Tree::vproj[1] = track_v0.at(uv[1]);
            Track2D track_v1(&vertexHitIn_v, superTrack->vTrack_->itrack2D_->hit2_);  // track through v1 hit
            Tree::vproj[0] = track_v1.at(uv[0]);
            Track2D track_t0(&vertexHitIn_t, superTrack->tTrack_->itrack2D_->hit1_);  // track through t0 hit
            Tree::tproj[1] = track_t0.at(ut[1]);
            Track2D track_t1(&vertexHitIn_t, superTrack->tTrack_->itrack2D_->hit2_);  // track through t1 hit
            Tree::tproj[0] = track_t1.at(ut[0]);

            Track2D track_v2(&vertexHitOut_v, superTrack->vTrack_->otrack2D_->hit1_);  // track through v0 hit
            Tree::vproj[3] = track_v2.at(uv[3]);
            Track2D track_v3(&vertexHitOut_v, superTrack->vTrack_->otrack2D_->hit2_);  // track through v1 hit
            Tree::vproj[2] = track_v3.at(uv[2]);
            Track2D track_t2(&vertexHitOut_t, superTrack->tTrack_->otrack2D_->hit1_);  // track through t0 hit
            Tree::tproj[3] = track_t2.at(ut[3]);
            Track2D track_t3(&vertexHitOut_t, superTrack->tTrack_->otrack2D_->hit2_);  // track through t1 hit
            Tree::tproj[2] = track_t3.at(ut[2]);
         }
      }
      otree->Fill();
   }

   cout<< "Writing " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
   otree->Write();
}

void beamspot_proj(const char* ifname, Int_t event1=0, Int_t event2=-1)
{
   Bool_t debug = kFALSE;
   if (debug) cout<< "debug in on" <<endl;

   TFile* ifile = new TFile(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return;
   }

   TTree* tree = (TTree*) gDirectory->Get("r");
   if (!tree) {
      cout<< "Could not find tree r" <<endl;
      return;
   }

   beamspot_proj(event1, event2, tree);
}
