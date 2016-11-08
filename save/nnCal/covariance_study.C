#include "NN.h"

#include <TROOT.h>
#include <TFile.h>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>

using std::cout;    using std::endl;

void covariance_study(Double_t d=0, const char* ifname="Calib_0043_000.dat.root.reco.root.good.root.nn.root")
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
    new TCanvas;
    tree->Draw("a[0]","abs(0-d)<1&&a[0]>3000&&a[0]<4500&&a[1]>3000&&a[1]<4500");
    new TCanvas;
    tree->Draw("a[1]","abs(0-d)<1&&a[0]>3000&&a[0]<4500&&a[1]>3000&&a[1]<4500");
    new TCanvas;
    tree->Draw("a[1]:a[0]","abs(0-d)<1 &&a[0]>3000&&a[0]<4500&&a[1]>3000&&a[1]<4500");
    new TCanvas;
    tree->Draw("a[0]+a[1]","abs(0-d)<1&&a[0]>3000&&a[0]<4500&&a[1]>3000&&a[1]<4500");

    Double_t amean0 = 0;
    Double_t amean1 = 0;
    Int_t np = 0;

    cout<< "get the mean value" <<endl;

    for (int ientry=0; ientry<tree->GetEntries(); ++ientry) {
        if (tree->LoadTree(ientry) < 0) break;
        tree->GetEntry(ientry);

        if (ientry % 100000 == 0) cout<< "processing entry " << ientry <<endl;

        bool selection = false
            || TMath::Abs(d - NN::d) < 1.
            && NN::a[0]>3000 && NN::a[0]<4500
            && NN::a[1]>3000 && NN::a[1]<4500
            ;

        if (selection) {
            ++np;
            amean0 += NN::a[0];
            amean1 += NN::a[1];
        }
    }
    cout<< "np = " << np <<endl;
    if (np > 0) {
        amean0 /= np;
        amean1 /= np;
    }

    cout<< "calculate covariance" <<endl;

    Double_t var0 = 0;
    Double_t var1 = 0;
    Double_t cov01 = 0;
    np = 0;

    for (int ientry=0; ientry<tree->GetEntries(); ++ientry) {
        if (tree->LoadTree(ientry) < 0) break;
        tree->GetEntry(ientry);

        if (ientry % 100000 == 0) cout<< "processing entry " << ientry <<endl;

        bool selection = false
            || TMath::Abs(d - NN::d) < 1.
            && NN::a[0]>3000 && NN::a[0]<4500
            && NN::a[1]>3000 && NN::a[1]<4500
            ;

        if (selection) {
            ++np;
            var0 += (NN::a[0] - amean0) * (NN::a[0] - amean0);
            var1 += (NN::a[1] - amean1) * (NN::a[1] - amean1);
            cov01 += (NN::a[0] - amean0) * (NN::a[1] - amean1);
        }
    }
    var0 /= np;
    var1 /= np;
    cov01 /= np;

    Double_t sigma0 = TMath::Sqrt(var0);
    Double_t sigma1 = TMath::Sqrt(var1);
    Double_t variance = var0 + var1 + 2.*cov01;
    Double_t sigma01 = TMath::Sqrt(variance);
    Double_t correlation = cov01/TMath::Sqrt(var0)/TMath::Sqrt(var1);

    cout<< "amean0 = " << amean0 << " amean1 = " << amean1 << " var0 = " << var0 << " var1 = " << var1 << " cov01 = " << cov01 << " variance = " << variance << " sigma0 = " << sigma0 << " sigma1 = " << sigma1 << " sigma01 = " << sigma01 << " correlation coefficient = " << correlation <<endl;

    cout<< "calculate in the fast way" <<endl;

    Double_t x0 = 0;
    Double_t x1 = 0;
    Double_t sqrx0 = 0;
    Double_t sqrx1 = 0;
    Double_t x0x1 = 0;
    np = 0;

    for (int ientry=0; ientry<tree->GetEntries(); ++ientry) {
        if (tree->LoadTree(ientry) < 0) break;
        tree->GetEntry(ientry);

        if (ientry % 100000 == 0) cout<< "processing entry " << ientry <<endl;

        bool selection = false
            || TMath::Abs(d - NN::d) < 1.
            && NN::a[0]>3000 && NN::a[0]<4500
            && NN::a[1]>3000 && NN::a[1]<4500
            ;

        if (selection) {
            ++np;
            x0 += NN::a[0];
            x1 += NN::a[1];
            sqrx0 += NN::a[0]*NN::a[0];
            sqrx1 += NN::a[1]*NN::a[1];
            x0x1 += NN::a[0]*NN::a[1];
        }
    }
    x0 /= np;
    x1 /= np;
    sqrx0 /= np;
    sqrx1 /= np;
    x0x1 /= np;

    Double_t fast_var0 = sqrx0 - x0*x0;
    Double_t fast_var1 = sqrx1 - x1*x1;
    Double_t fast_cov01 = x0x1 - x0*x1;
    Double_t fast_variance = fast_var0 + fast_var1 + 2.*fast_cov01;
    Double_t fast_sigma0 = TMath::Sqrt(fast_var0);
    Double_t fast_sigma1 = TMath::Sqrt(fast_var1);
    Double_t fast_sigma01 = TMath::Sqrt(fast_variance);
    Double_t fast_correlation = fast_cov01/TMath::Sqrt(fast_var0)/TMath::Sqrt(fast_var1);

    cout<< "x0 = " << x0 << " x1 = " << x1 << " fast_var0 = " << fast_var0 << " fast_var1 = " << fast_var1 << " fast_cov01 = " << fast_cov01 << " fast_variance = " << fast_variance << " fast_sigma0 = " << fast_sigma0 << " fast_sigma1 = " << fast_sigma1 << " fast_sigma01 = " << fast_sigma01 << " fast correlation coefficient = " << fast_correlation <<endl;
}
