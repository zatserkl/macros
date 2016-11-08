#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>

#include <vector>
#include <string>
#include <iostream>

using std::cout;    using std::endl;

//
//  Converts tree filled by std::vector<double> to plane array of double
//
//  See example https://root.cern.ch/root/html/tutorials/math/mathcoreVectorCollection.C.html
//

void convertVector(const char* ifname, const char* ofname="")
{
    TFile* ifile = new TFile(ifname);
    if (!ifile) {
        cout<< "Could not find file " << ifname <<endl;
        return;
    }
    TTree* itree = (TTree*) ifile->Get("Data");

    std::vector<double>* vt1 = new std::vector<double>(20000);      // time, ns
    std::vector<double>* vt2 = new std::vector<double>(20000);
    std::vector<double>* vt3 = new std::vector<double>(20000);
    std::vector<double>* vy1 = new std::vector<double>(20000);      // amplitude, V
    std::vector<double>* vy2 = new std::vector<double>(20000);
    std::vector<double>* vy3 = new std::vector<double>(20000);

    itree->SetBranchAddress("vt1",  &vt1);
    itree->SetBranchAddress("vt2",  &vt2);
    itree->SetBranchAddress("vt3",  &vt3);
    itree->SetBranchAddress("vy1",  &vy1);
    itree->SetBranchAddress("vy2",  &vy2);
    itree->SetBranchAddress("vy3",  &vy3);

    // load some event to get the number of points

    Int_t event = 4;
    Int_t chan = 1;

    if (itree->LoadTree(event) < 0) return;
    itree->GetEntry(event);

    Int_t npoints = vt1->size();
    cout<< "npoints = " << npoints <<endl;

    TGraph* g = new TGraph(102, vt1->data(), vy2->data());
    g->SetNameTitle(Form("evt_%d_ch_%d",event,chan), Form("evt_%d_ch_%d",event,chan));
    g->SetMarkerStyle(6);
    new TCanvas;
    g->Draw("ap");

    std::string ofnameStr;
    if (ofname && *ofname) ofnameStr = ofname;
    else ofnameStr = Form("%s.wfm.root",ifname);

    TFile* ofile = new TFile(ofnameStr.c_str(), "recreate");
    TTree* otree = new TTree("wfm", Form("wfm from %s",ifname));

    otree->SetMarkerStyle(6);
    otree->SetMarkerColor(602);
    otree->SetLineColor(602);

    otree->Branch("t1", vt1->data(), Form("t1[%d]/D",npoints));
    otree->Branch("t2", vt2->data(), Form("t2[%d]/D",npoints));
    otree->Branch("t3", vt3->data(), Form("t3[%d]/D",npoints));
    otree->Branch("w1", vy1->data(), Form("w1[%d]/D",npoints));
    otree->Branch("w2", vy2->data(), Form("w2[%d]/D",npoints));
    otree->Branch("w3", vy3->data(), Form("w3[%d]/D",npoints));

    for (int ientry=0; ientry<itree->GetEntries(); ++ientry) {
        if (itree->LoadTree(ientry) < 0) break;
        itree->GetEntry(ientry);
        if (ientry % 10000 == 0) cout<< "processing entry " << ientry <<endl;
        otree->Fill();
    }

    cout<< "Write " << otree->GetEntries() << " into output wfm file " << ofile->GetName() <<endl;
    ofile->Write();

    itree->ResetBranchAddresses();
}
