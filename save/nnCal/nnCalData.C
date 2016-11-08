// Andriy Zatserklyaniy, May 13, 2015

#include "Reco.h"
#include "DataFormat.h"

#include <TCanvas.h>
#include <TH2.h>
#include <TFile.h>
#include <TCut.h>
#include <TEventList.h>
#include <TMultiLayerPerceptron.h>
#include <TMLPAnalyzer.h>
#include <TChain.h>

#include <iomanip>

// Energy in MeV
// distance in cm

namespace NN
{
    Int_t evt;
    Double_t v;             // v and t at the steps (approx.). Important for the real data only.
    Double_t t;
    Double_t angle;         // angle between the input and output tracks: important for the real data only.
    Double_t uv;            // v-track: u-coordinate of intersection point of the input and output tracks
    Double_t xv;            // v-track: x-coordinate of intersection point of the input and output tracks
    Double_t ut;            // t-track: u-coordinate of intersection point of the input and output tracks
    Double_t xt;            // t-track: x-coordinate of intersection point of the input and output tracks
    Double_t weplReco;      // WEPL computed for the real data by Reco
    Double_t d;             // actual degrader thickness in cm
    Int_t slot;             // # slot in t from the center in the direction of the t-axis
    Int_t step;             // degrader thickness in the number of 6.35 mm steps
    Int_t brick;            // degrader thickness in the number of 6.35 mm steps
    Double_t Ed;            // energy deposited in the degrader
    Double_t Esci;          // energy deposited in scintillator
    Double_t a[5];          // raw ADC counts for the data
    Double_t E[5];          // energy deposited in each of scintillator stage (TV-corrected for the data)
    Double_t& E0 = E[0];
    Double_t& E1 = E[1];
    Double_t& E2 = E[2];
    Double_t& E3 = E[3];
    Double_t& E4 = E[4];
    Double_t Esum;

    void clear() {
        evt = -1;
        v = 0;
        t = 0;
        angle = 0;
        uv = xv = ut = xt = 0;
        weplReco = 0;
        d = 0;
        slot = -1;
        step = -1;
        brick = -1;
        Ed = 0;
        Esci = 0;
        for (int istage=0; istage<5; ++istage) {
            E[istage] = 0;
            a[istage] = 0;
        }
        Esum = 0;
    }

    void book(TTree* tree) {
        tree->Branch("evt",         &evt,       "evt/I");
        tree->Branch("v",           &v,         "v/D");
        tree->Branch("t",           &t,         "t/D");
        tree->Branch("angle",       &angle,     "angle/D");
        tree->Branch("uv",          &uv,        "uv/D");
        tree->Branch("xv",          &xv,        "xv/D");
        tree->Branch("ut",          &ut,        "ut/D");
        tree->Branch("xt",          &xt,        "xt/D");
        tree->Branch("weplReco",    &weplReco,  "weplReco/D");
        tree->Branch("d",           &d,         "d/D");
        tree->Branch("slot",        &slot,      "slot/I");
        tree->Branch("step",        &step,      "step/I");
        tree->Branch("brick",       &brick,     "brick/I");
        tree->Branch("Ed",          &Ed,        "Ed/D");
        tree->Branch("Esci",        &Esci,      "Esci/D");
        tree->Branch("a",           &a,         "a[5]/D");
        tree->Branch("E",           &E,         "E[5]/D");
        tree->Branch("E0",          &E0,        "E0/D");
        tree->Branch("E1",          &E1,        "E1/D");
        tree->Branch("E2",          &E2,        "E2/D");
        tree->Branch("E3",          &E3,        "E3/D");
        tree->Branch("E4",          &E4,        "E4/D");
        tree->Branch("Esum",        &Esum,      "Esum/D");
    }

    void connect(TTree* tree) {
        tree->SetBranchAddress("evt",        &evt);
        tree->SetBranchAddress("d",          &d);
        tree->SetBranchAddress("slot",       &slot);
        tree->SetBranchAddress("step",       &step);
        tree->SetBranchAddress("brick",      &brick);
        tree->SetBranchAddress("v",          &v);
        tree->SetBranchAddress("t",          &t);
        tree->SetBranchAddress("angle",      &angle);
        tree->SetBranchAddress("uv",         &uv);
        tree->SetBranchAddress("xv",         &xv);
        tree->SetBranchAddress("ut",         &ut);
        tree->SetBranchAddress("xt",         &xt);
        tree->SetBranchAddress("weplReco",   &weplReco);
        tree->SetBranchAddress("Ed",         &Ed);
        tree->SetBranchAddress("Esci",       &Esci);
        tree->SetBranchAddress("a",          &a);
        tree->SetBranchAddress("E",          &E);
        tree->SetBranchAddress("E0",         &E0);
        tree->SetBranchAddress("E1",         &E1);
        tree->SetBranchAddress("E2",         &E2);
        tree->SetBranchAddress("E3",         &E3);
        tree->SetBranchAddress("E4",         &E4);
        tree->SetBranchAddress("Esum",       &Esum);
    }
}  // namespace NN

TChain* chain() {
    TChain* s = new TChain("s");
    s->SetMarkerStyle(7);
    s->SetMarkerColor(602);
    s->SetLineColor(602);
    s->AddFile("Calib_0043_000.dat.root.reco.root.good.root.nn.root");          // 0 bricks
    // s->AddFile("Calib_0044_000.dat.root.reco.root.good.root.nn.root");           // 1 bricks
    // s->AddFile("Calib_0045_000.dat.root.reco.root.good.root.nn.root");           // 2 bricks
    // s->AddFile("Calib_0046_000.dat.root.reco.root.good.root.nn.root");           // 3 bricks
    // s->AddFile("Calib_0047_000.dat.root.reco.root.good.root.nn.root");           // 4 bricks

    cout<< "Connect the chain with NN buffers" <<endl;
    NN::connect(s);
    cout<< "--> The total number of entries in the chain s->GetEntries() = " << s->GetEntries() <<endl;
    return s;
}

TChain* s = chain();        //--NB: global variable

void nnCalData(TTree* tree, TTree* otree, TTree* nntree, Int_t nbricks)
{
    Bool_t debug = kFALSE;
    if (debug) cout<< "debug in on" <<endl;

    // use runHeader to get the run number
    RunHeader* runHeader = (RunHeader*) tree->GetUserInfo()->First();
    cout<< "runHeader->GetRun() = " << runHeader->GetRun() <<endl;
    time_t start_time = runHeader->GetTime();
    cout<< "run start time: " << std::ctime(&start_time);
    cout<< "program version is " << runHeader->GetVersion() <<endl;
    if (runHeader->GetTimeTag()) cout<< "event time tag was written out" <<endl;
    else cout<< "event time tag was not written out" <<endl;
    cout<<endl;

    // connect trees to make possible copying of events from tree to otree

    //-- const RecoEvent* recoEvent = 0;
    RecoEvent* recoEvent = new RecoEvent();                     // NB: create an object
    tree->SetBranchAddress("revent", &recoEvent);               // NB: address of pointer to recoEvent here
    //-- otree->Branch("revent", &recoEvent);
    otree->Branch("revent", recoEvent);                         // NB: pointer to recoEvent here

    const Int_t nsteps = 9;
    Double_t t1[nsteps];
    Double_t t2[nsteps];
    Double_t d[nsteps];
    const Double_t width = 2.;              // use +-width from the center
    const Double_t step = 6.35;
    Double_t center = 0 + step/2 + 50.8;
    Double_t depth = 0;
    for (int istep=0; istep<nsteps; ++istep) {
        t1[istep] = center - width;
        t2[istep] = center + width;
        d[istep] = nbricks*50.8 + depth;
        depth += step;                          // going up in depth
        center -= step;                         // going backward in t
    }
    // modify ends
    t2[0] += step;              // extend right boundary for the first step
    t1[nsteps-1] -= step;       // extend left boundary for the last step
    d[0] = nbricks*50.8;        // just to have exact number for the Double_t number

    for (int istep=0; istep<nsteps; ++istep) cout<< istep << "\t" << istep+nbricks*8 << "\t " << d[istep] << "\t " << t1[istep] << "\t " << t2[istep] <<endl;
    // for (int istep=0; istep<nsteps; ++istep) cout<< std::setw(10) << istep << std::setw(10) << d[istep] << std::setw(10) << t1[istep] << std::setw(10) << t2[istep] <<endl;

    for (int ientry=0; ientry<tree->GetEntries(); ++ientry) {
        if (tree->LoadTree(ientry) < 0) break;
        tree->GetEntry(ientry);

        if (false
            || (ientry < 100000 && ientry%10000 == 0)
            || ientry%100000 == 0
           )
        cout<< "processing entry " << ientry <<endl;

        if (recoEvent->track->GetLast()+1 == 0) continue;
        const SuperTrack* superTrack = (const SuperTrack*) recoEvent->track->At(0);

        NN::clear();

        NN::v = superTrack->V(-100);
        if (TMath::Abs(NN::v) > 30) continue;       // v acceptance

        NN::t = superTrack->T(-100);
        NN::Esum = 0;

        for (int istep=0; istep<nsteps; ++istep) {
            if (NN::t > t1[istep] && NN::t <= t2[istep]) {
                NN::d = d[istep];
                NN::slot = (nsteps-1) - istep;
                NN::step = nbricks*8 + istep;
                NN::brick = nbricks;
                NN::angle = superTrack->angle;
                CRay2D::intersect(superTrack->vTrack_->itrack2D_, superTrack->vTrack_->otrack2D_, NN::uv, NN::xv);
                CRay2D::intersect(superTrack->tTrack_->itrack2D_, superTrack->tTrack_->otrack2D_, NN::ut, NN::xt);
                NN::weplReco = recoEvent->wepl;
                for (int istage=0; istage<5; ++istage) NN::a[istage] = recoEvent->a[istage];    // raw ADC
                //-- for (int istage=0; istage<5; ++istage) {
                //--     NN::E[istage] = recoEvent->atv[istage];  // use TV-corrected ADC values
                //--     NN::Esum += recoEvent->atv[istage];
                //-- }
                otree->Fill();
                nntree->Fill();
                break;
            }
        }
    }

    otree->Write();
    nntree->Write();

    tree->ResetBranchAddresses();
    otree->ResetBranchAddresses();
    nntree->ResetBranchAddresses();
}

void nnCalData(const char* ifname, Int_t nbricks)
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

    //-- TFile* ofile = new TFile(Form("%s.good.root",ifname), "recreate");
    //-- TTree* otree = new TTree("r", Form("selected events from %s",ifname));
    //-- otree->SetMarkerColor(602);
    //-- otree->SetLineColor(602);

    TFile* nnfile = new TFile(Form("%s.good.root.nn.root",ifname), "recreate");
    TTree* nntree = new TTree("s", Form("selected events from %s in stage format",ifname));
    nntree->SetMarkerColor(602);
    nntree->SetLineColor(602);
    nntree->SetMarkerStyle(7);
    NN::book(nntree);

    TTree* otree = new TTree("r", Form("selected events from %s",ifname));
    otree->SetMarkerColor(602);
    otree->SetLineColor(602);

    nnCalData(tree,otree,nntree, nbricks);

    cout<< "otree->GetEntries() = " << otree->GetEntries() <<endl;
    cout<< "nntree->GetEntries() = " << nntree->GetEntries() <<endl;
    cout<< "Write simple tree in file " << nnfile->GetName() <<endl;
    //-- ofile->Write();
    nnfile->Write();
}

/// TTree* mlpCalArrayData(Int_t ntrain=100)
/// {
///     //const char* ifname="Calib_0042_000.dat.root.reco.root.good.root.nn.root";
///     const char* ifname="nn.root";                                                   //--NB: will be used to create ofname
/// 
///     // open input file with the data. It also has a type information: variable d
/// 
///     // TFile* ifile = new TFile(ifname);
///     // if (!ifile) {
///     //     cout<< "Could not open file " << ifname <<endl;
///     //     return 0;
///     // }
/// 
///     // const Double_t thres = 100.;                                             // threshold in ADC counts
///     // //-- const Double_t thres = 0;                                                   // threshold in ADC counts
///     // if (thres > 0) cout<< "--> apply threshold " << thres << " MeV" <<endl;
/// 
///     // TTree* tree = (TTree*) ifile->Get("s");
/// 
///     TChain* tree = chain();                         // get chain from the standalone function
///     NN::connect(tree);
/// 
///     const Double_t thres = 5.;                      // threshold in ~MeV
/// 
///     TMultiLayerPerceptron* mlp[5] = {0,0,0,0,0};
/// 
///     int istage = 0;
/// 
///     // TCut odd = "Entry$%2";
///     // TCut even = "(Entry$+1)%2";
///     TCut odd = "Entry$%32==0";
///     TCut even = "(Entry$+1)%32==0";
///     TCut cut = "";
///     TEventList* elist = 0;
/// 
///     TEventList* elist_train[5];
///     TEventList* elist_test[5];
/// 
///     /// for (istage=0; istage<5-1; ++istage)
///     /// {
///     ///     cut = Form("E[%d]>=%f&&E[%d]<%f",istage,thres,istage+1,thres);
///     ///     //cout<< "cut.GetTitle() = " << cut.GetTitle() <<endl;
///     ///     // train
///     ///     tree->Draw(Form(">>elist_train_%d",istage), even + cut);                  // train: even entries
///     ///     elist = (TEventList*) gDirectory->Get(Form("elist_train_%d",istage));
///     ///     elist_train[istage] = elist;
///     ///     // test
///     ///     tree->Draw(Form(">>elist_test_%d",istage), odd + cut);                   // test: odd entries
///     ///     elist = (TEventList*) gDirectory->Get(Form("elist_test_%d",istage));
///     ///     elist_test[istage] = elist;
///     /// }
///     /// istage = 4;
///     /// cut = Form("E[%d]>=%f",istage,thres);
/// 
///     // cuts
///     // for (istage=4; istage<5-1; ++istage)
///     // {
///     //     cut = "E0>3000&&E0<4500&&E4>100&&E4<5000";      // good events for the Bragg peak in the last stage
///     //     cout<< "cut.GetTitle() = " << cut.GetTitle() <<endl;
///     //     // train
///     //     tree->Draw(Form(">>elist_train_%d",istage), even + cut);                  // train: even entries
///     //     elist = (TEventList*) gDirectory->Get(Form("elist_train_%d",istage));
///     //     elist_train[istage] = elist;
///     //     // test
///     //     tree->Draw(Form(">>elist_test_%d",istage), odd + cut);                   // test: odd entries
///     //     elist = (TEventList*) gDirectory->Get(Form("elist_test_%d",istage));
///     //     elist_test[istage] = elist;
///     // }
///     // separate processing of the last stage
///     istage = 4;
///     //--ok E4-- cut = "E0>3000&&E0<4500&&E4>100&&E4<5000 &&(abs(d)<1||abs(6.35-d)<1||abs(12.7-d)<1)";      // good events for the Bragg peak in the last stage
///     std::stringstream cut_ss;
///     cut_ss << "E0>20&&E0<35&&E3>30&&E3<70&&Esum>140&&Esum<210"
///         // << " &&(abs( 0.00 - d) < 1 || abs( 6.35 - d) < 1 || abs(12.70 - d) < 1 || abs(19.05 - d) < 1)";
///         << " &&(abs( 0.00 - d) < 1 || abs( 6.35 - d) < 1 || abs(12.70 - d) < 1)";
///     cut = cut_ss.str().c_str();
///     // train
///     tree->Draw(Form(">>elist_train_%d",istage), even + cut);                      // train: even entries
///     elist = (TEventList*) gDirectory->Get(Form("elist_train_%d",istage));
///     elist_train[istage] = elist;
///     // test
///     tree->Draw(Form(">>elist_test_%d",istage), odd + cut);                       // test: odd entries
///     elist = (TEventList*) gDirectory->Get(Form("elist_test_%d",istage));
///     elist_test[istage] = elist;
/// 
///     for (istage=4; istage<5; ++istage) {        //--NB: for the last stage only
///         cout<< "elist_train[" << istage << "]->GetN() = " << elist_train[istage]->GetN() <<endl;
///         cout<< "elist_test[" << istage << "]->GetN() = " << elist_test[istage]->GetN() <<endl;
///     }
/// 
///     // mlp[0] = new TMultiLayerPerceptron("@E0:1:d", "E0", tree, elist_train[0], elist_test[0]);
///     // mlp[1] = new TMultiLayerPerceptron("@E0,@E1:2:d", "E0+E1", tree, elist_train[1], elist_test[1]);
///     // mlp[2] = new TMultiLayerPerceptron("@E0,@E1,@E2:3:2:d", tree, elist_train[2], elist_test[2]);
///     // mlp[3] = new TMultiLayerPerceptron("@E0,@E1,@E2,@E3:4:2:d", tree, elist_train[3], elist_test[3]);
///     //-- mlp[4] = new TMultiLayerPerceptron("@E0,@E1,@E2,@E3,@E4:5:3:2:d", tree, elist_train[4], elist_test[4]);
///     //--ok-- mlp[4] = new TMultiLayerPerceptron("@E0,@E1,@E2,@E3,@E4:5:2:d", tree, elist_train[4], elist_test[4]);
/// 
///     //-- mlp[4] = new TMultiLayerPerceptron("@E0,@E1,@E2,@E3,@E4,@Esum:10:5:2:d", "Esum", tree, elist_train[4], elist_test[4]);
///     // mlp[4] = new TMultiLayerPerceptron("@E0,@E1,@E2,@E3,@E4,@Esum:5:2:d", tree, elist_train[4], elist_test[4]);
///     mlp[4] = new TMultiLayerPerceptron("@E0,@E3,@E4,@Esum:5:2:d", tree, elist_train[4], elist_test[4]);
/// 
///     //-- return 0;
/// 
///     for (istage=0; istage<5; ++istage)
///     {
///         if (!mlp[istage]) continue;
/// 
///         cout<< "--> processing istage = " << istage <<endl;
/// 
///         mlp[istage]->Train(ntrain, "text,graph,update=10");
/// 
///         mlp[istage]->Export(Form("nnWeplData%d",istage),"C++");
/// 
///         /// // Use TMLPAnalyzer to see what it looks for
///         /// TCanvas* mlpa_canvas = new TCanvas(Form("mlpa_canvas_%d",istage),Form("Network analysis for stage %d",istage));
///         /// //-- mlpa_canvas->Divide(2,2);
///         /// mlpa_canvas->Divide(1,2);
///         /// TMLPAnalyzer ana(mlp[istage]);
///         /// // Initialisation
///         /// ana.GatherInformations();
///         /// // output to the console
///         /// ana.CheckNetwork();
///         /// mlpa_canvas->cd(1);
///         /// // shows how each variable influences the network
///         /// ana.DrawDInputs();
///         /// mlpa_canvas->cd(2);
///         /// // shows the network structure
///         /// mlp[istage]->Draw();
///     }
/// 
///     //
///     //  output file with wepl
///     //
/// 
///     TFile* ofile = new TFile(Form("%s.wepl.root",ifname), "recreate");
///     TTree* otree = new TTree("nn","WEPL tree");
/// 
///     cout<< "open output file " << ofile->GetName() <<endl;
/// 
///     Double_t wepl;
///     Double_t dtrue;
///     otree->Branch("wepl", &wepl, "wepl/D");
///     otree->Branch("dtrue", &dtrue, "dtrue/D");
/// 
///     cout<< "\nFill output tree. " << tree->GetEntries() << " events to process" <<endl;
/// 
///     for (int ientry=0; ientry<tree->GetEntries(); ++ientry) {
///         if (tree->LoadTree(ientry) < 0) break;
///         tree->GetEntry(ientry);
/// 
///         if (ientry % 10000 == 0) cout<< "processiong entry " << ientry <<endl;
/// 
///         dtrue = NN::d;
///         wepl = -1001.;
///         if (false
///             || TMath::Abs( 0.00 - dtrue) < 1.
///             || TMath::Abs( 6.35 - dtrue) < 1.
///             || TMath::Abs(12.70 - dtrue) < 1.
///             //-- || TMath::Abs(19.05 - dtrue) < 1.
///            )
///         {
///             TMultiLayerPerceptron* mlp_ptr = 0;
///             // for (istage=0; istage<5-1; ++istage)
///             // {
///             //     if (NN::E[istage] >= thres && NN::E[istage+1] < thres) {
///             //         mlp_ptr = mlp[istage];
///             //         break;
///             //     }
///             // }
///             istage = 4;
///             if (!mlp_ptr) if (NN::E[istage] >= thres) mlp_ptr = mlp[istage];
/// 
///             if (true
///                 && NN::E0 > 20 && NN::E0 < 35       // good for slot 8 to 0
///                 && NN::E3 > 30 && NN::E3 < 70
///                 && NN::Esum > 140 && NN::Esum < 210
///                )
///             {
///                 //wepl = mlp_ptr? mlp_ptr->Evaluate(0, NN::E): -1001.;
///                 Double_t* parameters;
///                 Int_t ipar = 0;
///                 parameters[ipar++] = NN::E0;
///                 parameters[ipar++] = NN::E3;
///                 parameters[ipar++] = NN::E4;
///                 parameters[ipar++] = NN::Esum;
///                 if (mlp_ptr) wepl = mlp_ptr->Evaluate(0, parameters);
///             }
///         }
///         otree->Fill();
///     }
/// 
///     cout<< "Write " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
///     otree->Write();
/// 
///     // for (int istage=0; istage<5; ++istage) delete mlp[istage];
/// 
///     //--NB-- otree->ResetBranchAddresses();
/// 
///     return otree;
/// }

void mlp(Int_t ntrain=100)
{
    TChain* tree = chain();
    //TFile* ifile = new TFile("Calib_0043_000.dat.root.reco.root.good.root.nn.root");
    //TTree* tree = (TTree*) ifile->Get("s");

    TCut even = "Entry$%128==0";
    TCut odd = "(Entry$+1)%128==0";

    std::stringstream cut_ss;
    cut_ss << "E0>20&&E0<35"
        << "&&E3>30&&E3<70"
        // << "&&E4>5"
        << "&&E4>10"
        << "&&Esum>140&&Esum<210"
        << " &&(abs(0.00 - d) < 1 || abs(6.35 - d) < 1 || abs(12.70 - d) < 1 || abs(19.05 - d) < 1)";
        // << " &&(abs(0.00 - d) < 1 || abs(6.35 - d) < 1 || abs(12.70 - d) < 1)"
        ;
    TCut cut = cut_ss.str().c_str();

    TCut cut_train = cut + odd;
    TCut cut_test = cut + even;

    //-- TMultiLayerPerceptron *mlp = new TMultiLayerPerceptron("@E0,@E3,@E4,@Esum:5:3:d", "Esum", tree, cut_train, cut_test);
    TMultiLayerPerceptron *mlp = new TMultiLayerPerceptron("@E0,@E1,@E2,@E3,@E4,@Esum:5:3:d", "Esum", tree, cut_train, cut_test);

    mlp->Train(ntrain, "text,graph,update=5");
    mlp->Export("nntest","C++");

    // Use TMLPAnalyzer to see what it looks for
    TCanvas* mlpa_canvas = new TCanvas("mlpa_canvas","Network analysis");
    mlpa_canvas->Divide(1,2);
    TMLPAnalyzer ana(mlp);
    // Initialisation
    ana.GatherInformations();
    // output to the console
    ana.CheckNetwork();
    mlpa_canvas->cd(1);
    // shows how each variable influences the network
    ana.DrawDInputs();
    mlpa_canvas->cd(2);
    // shows the network structure
    mlp->Draw();
    mlpa_canvas->cd();

    new TCanvas;
    // ana.DrawNetwork(0, "slot==1", "slot==2");
    ana.DrawNetwork(0, "abs(6.35-d)<1", "abs(12.7-d)<1");

    //
    //  output file with wepl
    //

    // tree = chain();         // Reset the input tree

    TFile* ofile = new TFile("nntest.root", "recreate");
    TTree* otree = new TTree("nn","WEPL tree");
    otree->SetMarkerStyle(7);

    cout<< "\nopen output file " << ofile->GetName() <<endl;

    Double_t wepl;
    Double_t dtrue;
    otree->Branch("wepl", &wepl, "wepl/D");
    otree->Branch("dtrue", &dtrue, "dtrue/D");

    cout<< "\nFill output tree. " << tree->GetEntries() << " events to process" <<endl;

    Int_t nmess = 0;
    NN::clear();
    NN::connect(tree);

    for (int ientry=0; ientry<tree->GetEntries(); ++ientry) {
        if (tree->LoadTree(ientry) < 0) break;
        tree->GetEntry(ientry);

        if (ientry % 10000 == 0) cout<< "processiong entry " << ientry <<endl;

        dtrue = NN::d;
        wepl = -1001.;

        if (false
            || TMath::Abs( 0.00 - NN::d) < 1.
            || TMath::Abs( 6.35 - NN::d) < 1.
            || TMath::Abs(12.70 - NN::d) < 1.
            || TMath::Abs(19.05 - NN::d) < 1.
           )
        {
            if (++nmess <= 10) cout<< "d = " << NN::d << " E0 = " << NN::E0 << " E3 = " << NN::E3 << " Esum = " << NN::Esum <<endl;
            if (true
                && NN::E0 > 20 && NN::E0 < 35           // good for slot 8 to 0
                && NN::E3 > 30 && NN::E3 < 70
                && NN::E4 > 5
                // && NN::E4 > 10
                && NN::Esum > 140 && NN::Esum < 210
               )
            {
                Double_t parameters[100];
                Int_t ipar = 0;
                parameters[ipar++] = NN::E0;
                parameters[ipar++] = NN::E1;
                parameters[ipar++] = NN::E2;
                parameters[ipar++] = NN::E3;
                parameters[ipar++] = NN::E4;
                parameters[ipar++] = NN::Esum;
                // cout<< "Evaluate d = " << NN::d << " for E0 = " << NN::E0 << " E3 = " << NN::E3 << " Esum = " << NN::Esum <<endl;
                wepl = mlp->Evaluate(0, parameters);
                // cout<< "after mlp->Evaluate(0, parameters) wepl = " << wepl <<endl;
            }
        }
        otree->Fill();
    }

    cout<< "Write " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
    otree->Write();
}

// 
//  content of the file NNTest.C to use the authogenerated nntest.h and nntest.cxx
//
//  NB: it has in the beginning
//
//  #include "nntest.h
//  and also
//  namespace NN and function chain()
//

// void nnTest()
// {
//     TChain* tree = chain();
// 
//     //gROOT->LoadMacro("nntest.cxx++");
// 
//     TFile* ofile = new TFile("nntest.root", "recreate");
// 
//     TTree* otree = new TTree("nn", "built with nntest.cxx");
//     otree->SetMarkerStyle(7);
// 
//     Double_t wepl;
//     Double_t dtrue;
//     otree->Branch("wepl", &wepl, "wepl/D");
//     otree->Branch("dtrue", &dtrue, "dtrue/D");
// 
//     nntest nnTest;
//     Int_t nmess = 0;
// 
//     for (int ientry=0; ientry<tree->GetEntries(); ++ientry) {
//         if (tree->LoadTree(ientry) < 0) break;
//         tree->GetEntry(ientry);
// 
//         if (ientry % 10000 == 0) cout<< "processiong entry " << ientry <<endl;
// 
//         dtrue = NN::d;
//         wepl = -1001.;
// 
//         if (false
//             || TMath::Abs( 0.00 - NN::d) < 1.
//             || TMath::Abs( 6.35 - NN::d) < 1.
//             || TMath::Abs(12.70 - NN::d) < 1.
//             || TMath::Abs(19.05 - NN::d) < 1.
//            )
//         {
//             if (++nmess <= 10) cout<< "d = " << NN::d << " E0 = " << NN::E0 << " E3 = " << NN::E3 << " Esum = " << NN::Esum <<endl;
//             if (true
//                 // && NN::E0 > 20 && NN::E0 < 35           // good for slot 8 to 0
//                 // && NN::E3 > 30 && NN::E3 < 70
//                 // // && NN::E4 > 5
//                 // && NN::E4 > 10
//                 // && NN::Esum > 140 && NN::Esum < 210
//                )
//             {
//                 // cout<< "Evaluate d = " << NN::d << " for E0 = " << NN::E0 << " E3 = " << NN::E3 << " Esum = " << NN::Esum <<endl;
//                 wepl = nnTest.Value(0, NN::E0, NN::E1, NN::E2, NN::E3, NN::E4, NN::Esum);
//                 // cout<< "after mlp->Evaluate(0, parameters) wepl = " << wepl <<endl;
//             }
//         }
//         otree->Fill();
//     }
// 
//     cout<< "Write " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
//     otree->Write();
// }
