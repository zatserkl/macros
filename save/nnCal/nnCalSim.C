#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TMultiLayerPerceptron.h>
#include <TMLPAnalyzer.h>
#include <TH1.h>
#include <TEventList.h>
#include <TCut.h>

#include <iostream>

using std::cout;     using std::endl;

// Energy in MeV
// distance (and thickness) in cm

Double_t Range(Double_t Ein, Double_t density=1.05) {
    return TMath::Power(Ein,1.7637) / (247.93*1.7637) / density;
}

Double_t Eout(Double_t Ein, Double_t x, Double_t density=1.05) {
    Double_t arg = TMath::Power(Ein,1.7637) - 247.93*1.7637*x*density;
    if (arg > 0) return TMath::Power(arg,1./1.7637);
    else return 0;
}

// Double_t EdepositedStage(Double_t Ein, Double_t thickness, deltaE)
// {
//    Double_t EoutExast = Eout(Ein,thickness);
// 
//    if (EoutExact == 0) {
//       return 0;
//    }
//    return 0;
// }

class EnergyLoss {
    //
    // Energy in MeV
    // length in cm
    //
public:
    const Double_t eps;
    const Double_t density;       // polystyrene density = 1.05 g/cm3
    const Double_t dstage;        // stage thickness
    Double_t Emin;                // minimum proton energy to leave the stage = 81.2 MeV
    TRandom3 random;
public:
    EnergyLoss(UInt_t seed=0): eps(1e-7), density(1.05), dstage(5.08)
                               , Emin(81.2)          // recalc in constructor
    {
        random.SetSeed(seed);

        // calculate the Emin and add eps to that value
        // Emin is the energy that proton should have to pass a range 5.08 cm of polysterene
        // Emin is parametrized by two constants:
        // a = 247.93
        // b = 1.7637
        Emin = TMath::Power(247.93*1.7637*5.08*density, 1/1.7637) + eps;  // = 81.2 MeV
    }
    Double_t Range(Double_t Ein) const {
        return TMath::Power(Ein,1.7637) / (247.93*1.7637) / density;
    }
    Double_t EoutX(Double_t Ein, Double_t x) const
    {
        // proton energy (MeV) after passing of x cm of water
        Double_t arg = TMath::Power(Ein,1.7637) - 247.93*1.7637*x*density;
        if (arg > 0) return TMath::Power(arg,1./1.7637);
        else return 0;                                        // got stopped in medium
    }
    Double_t EoutStage(Double_t Ein) const {return EoutX(Ein, dstage);}
    Double_t EoutStageLandau(Double_t Ein)
    {
        Double_t dE = Ein - EoutStage(Ein);       // average energy loss in the stage
        Double_t Eloss = random.Landau(dE);       // energy loss with Landau fluctuation
        if (Eloss < 0) Eloss = 0;
        if (Eloss > Ein) Eloss = Ein;             // avoid acceleration

        return Ein - Eloss;
    }
    Double_t EoutStageGaus(Double_t Ein, Double_t sigmaPercent=0.05)
    {
        Double_t dE = Ein - EoutStage(Ein);                   // average energy loss in the stage
        Double_t Eloss = random.Gaus(dE, dE*sigmaPercent);    // energy loss with Gauss fluctuation
        if (Eloss < 0) Eloss = 0;
        if (Eloss > Ein) Eloss = Ein;                         // avoid acceleration

        return Ein - Eloss;
    }
    Double_t EoutXGaus(Double_t Ein, Double_t x, Double_t sigmaPercent=0.05)
    {
        Double_t dE = Ein - EoutX(Ein,x);                  // average energy loss after thickness x
        Double_t Eloss = random.Gaus(dE, dE*sigmaPercent); // energy loss with Gauss fluctuation
        if (Eloss < 0) Eloss = 0;
        if (Eloss > Ein) Eloss = Ein;                      // avoid acceleration

        return Ein - Eloss;
    }
};

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
    Double_t E[5];          // energy deposited in each of scintillator stage
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

/*
root -l 'nnCalSim.C+(10000,0.05,25.,2)'
*/
void nnCalSim(Int_t nevents_per_step=1000       // the number event per each of degrader step
          , Double_t sigma=0.05                 // energy spread, for real data sigma ~ 0.05 cm
          , Double_t dmax=25.                   // max degrader thickness
          , Double_t dstep=1.                   // degrader step
          , Double_t E0=200.
          , const char* ofname="nnCalSim.root"
          )
{
    bool debug = false;

    const Double_t eps = 1.e-7;

    UInt_t seed = 0;
    //-- seed = 1000;         // train
    seed = 2000;         // test

    TFile* ofile = new TFile(ofname, "recreate");
    TTree* tree = new TTree("s", "stage output");
    tree->SetMarkerStyle(7);
    NN::book(tree);

    EnergyLoss energyLoss;

    //-- const Double_t thres = 1.;                                             // threshold in MeV
    const Double_t thres = 0;                                                   // threshold in MeV
    if (thres > 0) cout<< "--> apply threshold " << thres << " MeV" <<endl;

    Int_t nevents = 0;                          // the total number of events generated

    for (Double_t d=0; d<dmax+eps; d+=dstep) {
        //
        // generate nevents_per_step
        //
        for (int ievent=0; ievent<nevents_per_step; ++ievent)
        {
            NN::clear();                                          // clear buffers before new stage

            Double_t EoutD = energyLoss.EoutXGaus(E0,d,sigma);      // energy after the degrader
            //-- Double_t EoutD = energyLoss.EoutX(E0,d);           // energy after the degrader NB: no energy smearing

            Double_t Ed = E0 - EoutD;                               // energy lost in the degrader

            NN::evt = nevents++;
            NN::d = d;                    // degrader thickness
            NN::Ed = Ed;                  // energy lost in the degrader

            Double_t Ecurr = E0 - Ed;       // current proton energy
            if (Ecurr > 0)
            {
                // go through the stages
                for (int istage=0; istage<5; ++istage)
                {
                    Double_t Eout = 0;
                    //if (istage < 4) Eout = energyLoss.EoutStageLandau(Ecurr);
                    //else Eout = energyLoss.EoutStage(Ecurr);

                    //-- Eout = energyLoss.EoutStageLandau(Ecurr);

                    Eout = energyLoss.EoutStageGaus(Ecurr, sigma);      // Gauss spread
                    //-- Eout = energyLoss.EoutStage(Ecurr);            // energy after the current stage: no energy spread

                    NN::E[istage] = Ecurr - Eout;                     // energy deposited in the current stage
                    NN::Esci += NN::E[istage];                      // the total energy lost in the scintillator
                    Ecurr = Eout;                                       // update Ecurr, the current proton energy

                    if (debug) cout<< NN::E[istage] << "   ";

                    if (thres > 0) if (NN::E[istage] < thres) NN::E[istage] = 0;

                    if (Ecurr < eps) break;
                }
                if (debug) cout<<endl;
            }
            tree->Fill();
        }
    }

    //tree->Write();
    ofile->Write();
}

void mlpCalSim(const char* ifname="nnCalSim.root", Int_t ntrain=100)
{
    // open input file with the data. It also has a type information: variable d

    TFile* ifile = new TFile(ifname);
    if (!ifile) {
        cout<< "Could not open file " << ifname <<endl;
        return;
    }

    TTree* tree = (TTree*) ifile->Get("s");
    NN::connect(tree);

    // build and train NN
    //-- TMultiLayerPerceptron *mlp = new TMultiLayerPerceptron("@E0,@E1,@E2,@E3,@E4:7:3:d",tree,"Entry$%2","(Entry$+1)%2");
    //-- TMultiLayerPerceptron *mlp = new TMultiLayerPerceptron("@E0,@E1,@E2,@E3,@E4:5:d",tree,"Entry$%2","(Entry$+1)%2");
    TMultiLayerPerceptron *mlp = new TMultiLayerPerceptron("@E0,@E1,@E2,@E3,@E4:9:4:d",tree,"Entry$%2","(Entry$+1)%2");

    mlp->Train(ntrain, "text,graph,update=10");

    mlp->Export("nnWeplSim","C++");

    // Use TMLPAnalyzer to see what it looks for
    TCanvas* mlpa_canvas = new TCanvas("mlpa_canvas","Network analysis");
    //-- mlpa_canvas->Divide(2,2);
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

    // mlpa_canvas->cd(3);
    //
    //      AZ: skip this part for now
    //
    //-- // draws the resulting network
    //-- ana.DrawNetwork(0,"type==1","type==0");

    // mlpa_canvas->cd(4);
    //
    //      AZ: skip this part for now
    //

    // test
    // Double_t E[] = {26.4947, 30.6918, 37.9094, 62.6504, 23.1756};       // entry 49000, d = 4 cm
    // cout<< "\nmlp->Evaluate(0,E) = " << mlp->Evaluate(0,E) <<endl;

    //
    //  output file with wepl
    //

    TFile* ofile = new TFile(Form("%s.wepl.root",ifname), "recreate");
    TTree* otree = new TTree("nn","WEPL tree");

    Double_t wepl;
    Double_t dtrue;
    otree->Branch("wepl", &wepl, "wepl/D");
    otree->Branch("dtrue", &dtrue, "dtrue/D");

    cout<< "\nCreate output tree. " << tree->GetEntries() << " events to process" <<endl;

    for (int ientry=0; ientry<tree->GetEntries(); ++ientry) {
        if (tree->LoadTree(ientry) < 0) break;
        tree->GetEntry(ientry);

        if (ientry % 10000 == 0) cout<< "processiong entry " << ientry <<endl;

        dtrue = NN::d;
        wepl = mlp->Evaluate(0, NN::E);

        otree->Fill();
    }

    cout<< "Write " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
    otree->Write();
}

void mlpCalArraySim(const char* ifname="nnCalSim.root", Int_t ntrain=100)
{
    // open input file with the data. It also has a type information: variable d

    TFile* ifile = new TFile(ifname);
    if (!ifile) {
        cout<< "Could not open file " << ifname <<endl;
        return;
    }

    const Double_t thres = 1.;                                             // threshold in MeV
    //-- const Double_t thres = 0;                                                   // threshold in MeV
    if (thres > 0) cout<< "--> apply threshold " << thres << " MeV" <<endl;

    TTree* tree = (TTree*) ifile->Get("s");
    NN::connect(tree);

    TMultiLayerPerceptron* mlp[5] = {0,0,0,0,0};

    int istage = 0;

    TCut even = "Entry$%2";
    TCut odd = "(Entry$+1)%2";
    TCut cut = "";
    TEventList* elist = 0;

    TEventList* elist_train[5];
    TEventList* elist_test[5];
    for (istage=0; istage<5-1; ++istage)
    {
        cut = Form("E[%d]>=%f&&E[%d]<%f",istage,thres,istage+1,thres);
        //cout<< "cut.GetTitle() = " << cut.GetTitle() <<endl;
        // train
        tree->Draw(Form(">>elist_train_%d",istage), even + cut);                  // train: even entries
        elist = (TEventList*) gDirectory->Get(Form("elist_train_%d",istage));
        elist_train[istage] = elist;
        // test
        tree->Draw(Form(">>elist_test_%d",istage), odd + cut);                   // test: odd entries
        elist = (TEventList*) gDirectory->Get(Form("elist_test_%d",istage));
        elist_test[istage] = elist;
    }
    istage = 4;
    cut = Form("E[%d]>=%f",istage,thres);
    // train
    tree->Draw(Form(">>elist_train_%d",istage), even + cut);                      // train: even entries
    elist = (TEventList*) gDirectory->Get(Form("elist_train_%d",istage));
    elist_train[istage] = elist;
    // test
    tree->Draw(Form(">>elist_test_%d",istage), odd + cut);                       // test: odd entries
    elist = (TEventList*) gDirectory->Get(Form("elist_test_%d",istage));
    elist_test[istage] = elist;

    for (istage=0; istage<5; ++istage) {
        cout<< "elist_train[" << istage << "]->GetN() = " << elist_train[istage]->GetN() <<endl;
        cout<< "elist_test[" << istage << "]->GetN() = " << elist_test[istage]->GetN() <<endl;
    }

    mlp[0] = new TMultiLayerPerceptron("@E0:1:d", "E0", tree, elist_train[0], elist_test[0]);
    mlp[1] = new TMultiLayerPerceptron("@E0,@E1:2:d", "E0+E1", tree, elist_train[1], elist_test[1]);
    mlp[2] = new TMultiLayerPerceptron("@E0,@E1,@E2:3:2:d", tree, elist_train[2], elist_test[2]);
    mlp[3] = new TMultiLayerPerceptron("@E0,@E1,@E2,@E3:4:2:d", tree, elist_train[3], elist_test[3]);
    mlp[4] = new TMultiLayerPerceptron("@E0,@E1,@E2,@E3,@E4:5:3:2:d", tree, elist_train[4], elist_test[4]);

    //-- return;

    for (istage=0; istage<5; ++istage)
    {
        cout<< "--> processing istage = " << istage <<endl;

        if (!mlp[istage]) continue;

        mlp[istage]->Train(ntrain, "text,graph,update=10");

        mlp[istage]->Export(Form("nnWeplSim%d",istage),"C++");

        // Use TMLPAnalyzer to see what it looks for
        TCanvas* mlpa_canvas = new TCanvas(Form("mlpa_canvas_%d",istage),Form("Network analysis for stage %d",istage));
        //-- mlpa_canvas->Divide(2,2);
        mlpa_canvas->Divide(1,2);
        TMLPAnalyzer ana(mlp[istage]);
        // Initialisation
        ana.GatherInformations();
        // output to the console
        ana.CheckNetwork();
        mlpa_canvas->cd(1);
        // shows how each variable influences the network
        ana.DrawDInputs();
        mlpa_canvas->cd(2);
        // shows the network structure
        mlp[istage]->Draw();
    }

    //
    //  output file with wepl
    //

    TFile* ofile = new TFile(Form("%s.wepl.root",ifname), "recreate");
    TTree* otree = new TTree("nn","WEPL tree");

    Double_t wepl;
    Double_t dtrue;
    otree->Branch("wepl", &wepl, "wepl/D");
    otree->Branch("dtrue", &dtrue, "dtrue/D");

    cout<< "\nCreate output tree. " << tree->GetEntries() << " events to process" <<endl;

    for (int ientry=0; ientry<tree->GetEntries(); ++ientry) {
        if (tree->LoadTree(ientry) < 0) break;
        tree->GetEntry(ientry);

        if (ientry % 10000 == 0) cout<< "processiong entry " << ientry <<endl;

        dtrue = NN::d;
        TMultiLayerPerceptron* mlp_ptr = 0;
        for (istage=0; istage<5-1; ++istage)
        {
            if (NN::E[istage] >= thres && NN::E[istage+1] < thres) {
                mlp_ptr = mlp[istage];
                break;
            }
        }
        istage = 4;
        if (!mlp_ptr) if (NN::E[istage] >= thres) mlp_ptr = mlp[istage];

        wepl = mlp_ptr? mlp_ptr->Evaluate(0, NN::E): -1001.;
        otree->Fill();
    }

    cout<< "Write " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
    otree->Write();
}

/// void mlpD(const char* ifname="esim.root")
/// {
///     Int_t ntrain = 100;                    // the number of epoch
/// 
///     TFile* ifile = new TFile(ifname);
///     if (!ifile) {
///         cout<< "Could not open file " << ifname <<endl;
///         return;
///     }
/// 
///     TTree* tree = (TTree*) ifile->Get("s");
///     NN::connect(tree);
/// 
///     //-- Double_t dmax = 18.;    // cm
/// 
///     // TMultiLayerPerceptron *mlp = new TMultiLayerPerceptron("@E0,@E1,@E2,@E3,@E4:9:3:d",tree,"Entry$%2","(Entry$+1)%2");
///     //-- TMultiLayerPerceptron *mlp = new TMultiLayerPerceptron("@Esci,@E0,@E1,@E2,@E3,@E4:8:d","Esci",tree,"Entry$%2","(Entry$+1)%2");
///     //-- TMultiLayerPerceptron *mlp = new TMultiLayerPerceptron("@E0,@E1,@E2,@E3,@E4:8:d",tree,"Entry$%2","(Entry$+1)%2");
///     //----- TMultiLayerPerceptron *mlp = new TMultiLayerPerceptron("@E0,@E1,@E2,@E3,@E4:8:3:d",tree,"Entry$%2","(Entry$+1)%2");
///     //--not-bad-- TMultiLayerPerceptron *mlp = new TMultiLayerPerceptron("@E0,@E1,@E2,@E3,@E4:7:d",tree,"Entry$%2","(Entry$+1)%2");
///     //--25-and-32-not-resolved-- TMultiLayerPerceptron *mlp = new TMultiLayerPerceptron("@E0,@E1,@E2,@E3,@E4:6:d",tree,"Entry$%2","(Entry$+1)%2");
///     TMultiLayerPerceptron *mlp = new TMultiLayerPerceptron("@E0,@E1,@E2,@E3,@E4:7:d",tree,"Entry$%2","(Entry$+1)%2");
///     mlp->Train(ntrain, "text,graph,update=10");
///     mlp->Export("nnD","C++");
/// }
