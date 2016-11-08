#include "NN.h"

#include <TROOT.h>
#include <TFile.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TCut.h>

#include <iostream>

using std::cout;    using std::endl;

/*
Read about the covariance, e.g. at
http://mathworld.wolfram.com/Covariance.html
*/

void covariance(const char* ifname="Calib_0043_000.dat.root.reco.root.good.root.nn.root")
{
    TFile* ifile = new TFile(ifname);
    if (!ifile) {
        cout<< "Could not open file " << ifname <<endl;
        return;
    }
    TTree* tree = (TTree*) ifile->Get("s");
    cout<< "tree->GetEntries() = " << tree->GetEntries() <<endl;

    NN::connect(tree);

    tree->SetMarkerStyle(6); 

    Double_t d = 0;         // mm
    // Double_t d = 6.35;      // mm
    // Double_t d = 12.7;      // mm

    cout<< "d = " << d << " mm" <<endl;

    const Int_t NCAL = 5;                   // 5 stages of the energy detector
    Double_t sum_x[NCAL];
    Double_t sum_xx[NCAL][NCAL];
    //
    // sum_xx 10 elements in use (r/c shows the row and column numbers):
    //
    //  r/c 0   1   2   3   4
    //  0   00  01  02  03  04
    //  1   .   11  12  13  14
    //  2   .   .   22  23  24
    //  3   .   .   .   33  34
    //  4   .   .   .   .   44
    //
    for (int i=0; i<NCAL; ++i) for (int j=0; j<NCAL; ++j) sum_xx[i][j] = 0;

    // cuts that work for 0 <= d <= 12.70 mm
    Double_t amin[5] = {3000,   3000,   2000,   4500,   1000};
    Double_t amax[5] = {4500,   4000,   3500,   7000,   5000};

    Int_t np = 0;
    for (int ientry=0; ientry<tree->GetEntries(); ++ientry) {
        if (tree->LoadTree(ientry) < 0) break;
        tree->GetEntry(ientry);

        if (ientry % 100000 == 0) cout<< "processing entry " << ientry <<endl;

        bool selection = false                  // selection should be the same as in cuts used in Draw below
            || TMath::Abs(d - NN::d) < 1.
            && NN::a[0] > amin[0] && NN::a[0] < amax[0]
            && NN::a[1] > amin[1] && NN::a[1] < amax[1]
            && NN::a[2] > amin[2] && NN::a[2] < amax[2]
            && NN::a[3] > amin[3] && NN::a[3] < amax[3]
            && NN::a[4] > amin[4] && NN::a[4] < amax[4]
            ;

        if (selection) {
            ++np;
            for (int ical=0; ical<NCAL; ++ical) {
                sum_x[ical] += NN::a[ical];
                for (int next=ical; next<NCAL; ++next) sum_xx[ical][next] += NN::a[ical] * NN::a[next];
            }
        }
    }

    for (int ical=0; ical<NCAL; ++ical) {
        sum_x[ical] /= np;
        for (int next=ical; next<NCAL; ++next) sum_xx[ical][next] /= np;
    }

    Double_t cov[NCAL][NCAL];               // cov is covariance
    Double_t cor[NCAL][NCAL];               // cor is correlation coefficient
    for (int i=0; i<NCAL; ++i) for (int j=0; j<NCAL; ++j) {
        cov[i][j] = 0;
        cor[i][j] = 0;
    }
    cout<< "\nCovariance coefficients: will be used to get the sigmas of sum of the stages" <<endl<<endl;
    for (int ical=0; ical<NCAL; ++ical) {
        for (int next=ical; next<NCAL; ++next) {
            cov[ical][next] = sum_xx[ical][next] - sum_x[ical]*sum_x[next];
            cout<< "cov[" << ical << "][" << next << "] = " << cov[ical][next] << " ";
        }
        cout<<endl;
    }

    cout<< "\nCorrelation coefficients: just to show correlation of the stages" <<endl<<endl;
    for (int ical=0; ical<NCAL; ++ical) {
        for (int next=ical; next<NCAL; ++next) {
            cor[ical][next] = cov[ical][next] / TMath::Sqrt(cov[ical][ical] * cov[next][next]); // divide by sigmas
            cout<< "cor[" << ical << "][" << next << "] = " << cor[ical][next] << " ";
        }
        cout<<endl;
    }

    cout<< "\nSigmas of the channels: sqrt of diagonal covariance coefficients: sqrt(cov[i][i]):" <<endl;
    for (int ical=0; ical<NCAL; ++ical) {
        cout<< "sigma[" << ical << "] = " << TMath::Sqrt(cov[ical][ical]) << " ";
    }
    cout<<endl<<endl;

    Double_t VarienceOfSum[NCAL];               // varience of sum of stages
    Double_t sum_var = 0;                       // running sum of channel variances: diagonal elements of cov
    Double_t sum_cov = 0;                       // running sum of channel covariances
    for (int ical=0; ical<NCAL; ++ical) {
        sum_var += cov[ical][ical];
        for (int prev=0; prev<ical; ++prev) {
            sum_cov += cov[prev][ical];
        }
        VarienceOfSum[ical] = sum_var + 2.*sum_cov;
    }

    cout<< "Sigma of the channels sums: sqrt of varience of sum of the channels:" <<endl;
    for (int ical=0; ical<NCAL; ++ical) {
        cout<< "sigma(sum 0.." << ical << ") = " << TMath::Sqrt(VarienceOfSum[ical]) << "   ";
    }
    cout<<endl;

    //
    //  Draw
    //

    // s->Draw("a[4]","abs(0-d)<1 &&a[0]>3000&&a[0]<4500&&a[1]>3000&&a[1]<4000&&a[2]>2000&&a[2]<3500&&a[3]>4500&&a[3]<7000&&a[4]>1000&&a[4]<5000")

    TCut cutd = Form("abs(%0.2f-d)<1",d);

    TCut cut[5];
    TCut cutall = cutd;
    for (int ical=0; ical<5; ++ical) {
        TCut cuta = Form("a[%d]>%0.0f&&a[%d]<%0.0f",ical,amin[ical],ical,amax[ical]);
        cut[ical] = cutd + cuta;
        cutall += cuta;
    }

    cout<< "\ncutall.GetTitle() = " << cutall.GetTitle() <<endl;

    TCanvas* can;

    int ncanvas = 1;
    while (gROOT->GetListOfCanvases()->FindObject(Form("covariance_%d",ncanvas))) ncanvas++;

    can = new TCanvas(Form("covariance_%d",ncanvas),Form("covariance_%d",ncanvas), 1.10*700, 0.55*500);
    ++ncanvas;
    can->Divide(2,1);
    can->cd(1); tree->Draw("a[0]", cutd + cut[0]);
    can->cd(2); tree->Draw("a[0]", cutd + cutall);
    can = new TCanvas(Form("covariance_%d",ncanvas),Form("covariance_%d",ncanvas), 1.10*700, 0.55*500);
    ++ncanvas;
    can->Divide(2,1);
    can->cd(1); tree->Draw("a[1]", cutd + cut[1]);
    can->cd(2); tree->Draw("a[1]", cutd + cutall);
    can = new TCanvas(Form("covariance_%d",ncanvas),Form("covariance_%d",ncanvas), 1.10*700, 0.55*500);
    ++ncanvas;
    can->Divide(2,1);
    can->cd(1); tree->Draw("a[2]", cutd + cut[2]);
    can->cd(2); tree->Draw("a[2]", cutd + cutall);
    can = new TCanvas(Form("covariance_%d",ncanvas),Form("covariance_%d",ncanvas), 1.10*700, 0.55*500);
    ++ncanvas;
    can->Divide(2,1);
    can->cd(1); tree->Draw("a[3]", cutd + cut[3]);
    can->cd(2); tree->Draw("a[3]", cutd + cutall);
    can = new TCanvas(Form("covariance_%d",ncanvas),Form("covariance_%d",ncanvas), 1.10*700, 0.55*500);
    ++ncanvas;
    can->Divide(2,1);
    can->cd(1); tree->Draw("a[4]", cutd + cut[4]);
    can->cd(2); tree->Draw("a[4]", cutd + cutall);

    new TCanvas;
    tree->Draw("a[1]:a[0]",                 cutd + cutall);

    new TCanvas;
    tree->Draw("a[0]+a[1]",                 cutd + cutall);

    new TCanvas;
    tree->Draw("a[0]+a[1]+a[2]",            cutd + cutall);

    new TCanvas;
    tree->Draw("a[0]+a[1]+a[2]+a[3]",       cutd + cutall);

    new TCanvas;
    tree->Draw("a[0]+a[1]+a[2]+a[3]+a[4]",  cutd + cutall);
}
