#include <TTree.h>

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
