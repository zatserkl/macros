#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TMath.h>

#include <iostream>
#include <cassert>

using std::cout;     using std::endl;

void fft_inverse(const Int_t N, Double_t reX[], Double_t imX[], Double_t x[])
{
    assert(N%2 == 0);

    //-- Double_t x[N];                   // amplitudes in the time domain == time samples

    Double_t reXX[N/2+1];            // scaled Re amplitudes in the frequency domain = frequency samples
    Double_t imXX[N/2+1];            // scaled Im amplitudes in the frequency domain = frequency samples
    for (int k=0; k<N/2+1; ++k) {
        reXX[k] = reX[k]/(N/2);
        imXX[k] = -1.*imX[k]/(N/2);
    }
    reXX[0]    /= 2;
    reXX[N/2]  /= 2;

    for (int isample=0; isample<N; ++isample) x[isample] = 0;

    for (int k=0; k<=N/2; ++k) {
        for (int isample=0; isample<N; ++isample) {
            x[isample] += reXX[k]*TMath::Cos(2.*TMath::Pi()*k*isample/N);
            x[isample] += imXX[k]*TMath::Sin(2.*TMath::Pi()*k*isample/N);
        }
    }
}

void fft_correlation(const Int_t N, Double_t* x, Double_t* reX, Double_t* imX)
{
    //
    //  Correlation is basically a projection, like projection of a vector to basis. 
    //  In statistics it's the same: each variable play role of a basis vector.
    //

    assert(N%2 == 0);

    for (int k=0; k<=N/2; ++k) {
        reX[k] = 0;
        imX[k] = 0;
    }

    for (int k=0; k<=N/2; ++k) {
        for (int isample=0; isample<N; ++isample) {
            reX[k] += x[isample]*TMath::Cos(2.*TMath::Pi()*k*isample/N);
            imX[k] -= x[isample]*TMath::Sin(2.*TMath::Pi()*k*isample/N);
        }
    }
}

void Polar(const Int_t N, Double_t reX[], Double_t imX[], Double_t magX[], Double_t phaX[])
{
    // NB: N is the number of samples in the time domain
    // The number of samples in the frequecy domain is N/2+1
    // Consider the N as even number

    const Double_t eps = 1e-7;    // calculate phase only for points with mag > eps

    for (int k=0; k<=N/2; ++k) {
        magX[k] = TMath::Sqrt(reX[k]*reX[k]+imX[k]*imX[k]);
        phaX[k] = magX[k] > eps? TMath::ATan2(imX[k],reX[k]): 0;
    }
}

void fft_correlation()
{
    const char* ifname = "mppc-72.0V-dark.dat.root";

    Int_t event = 0;

    TFile* ifile = new TFile(ifname);
    if (!ifile) {
        cout<< "Could not open file " << ifname <<endl;
        return;
    }

    TTree* tree = (TTree*) ifile->Get("p");
    if (!tree) {
        cout<< "Could not find tree \"p\" in " << ifname <<endl;
        return;
    }

    tree->Draw("c1:t1",Form("Entry$==%d",event),"goff");

    Double_t t[1024], v[1024];
    Int_t np = tree->GetSelectedRows(); 
    for (int i=0; i<np; ++i) {
        t[i] = tree->GetV2()[i];
        v[i] = -1.*tree->GetV1()[i];
    }

    //--64-- const Int_t N = 64;           // the number of samples
    const Int_t N = 512;           // the number of samples

    Int_t tmin = 0;
    Int_t i1 = 0;
    while (i1 < np && t[i1] < tmin) ++i1;

    cout<< "i1 = " << i1 << "\t t[i1] = " << t[i1] <<endl;

    Double_t baseline = 0;
    for (int isample=0; isample<N; ++isample) baseline += v[isample];
    baseline /= N;
    for (int i=0; i<1024; ++i) v[i] -= baseline;

    Double_t* x = &v[i1];         // pointer to the first point

    TGraph* g = new TGraph(np, t, v);
    g->SetNameTitle(Form("g%d",event), Form("Event %d", event));
    g->SetMarkerStyle(7);
    g->SetMarkerColor(2);
    g->SetLineColor(2);

    new TCanvas;
    g->Draw("apl");

    Double_t reX[N/2+1];
    Double_t imX[N/2+1];
    for (int k=0; k<=N/2; ++k) {
        reX[k] = 0;
        imX[k] = 0;
    }

    // get the reX and imX

    fft_correlation(N, x, reX, imX);

    Double_t sample[1024];
    for (int isample=0; isample<1024; ++isample) sample[isample] = isample;

    TGraph* gsample = new TGraph(N, sample, x);
    gsample->SetNameTitle("gsample", "Time Domain;sample number;V, mV");
    gsample->SetMarkerStyle(7);
    gsample->SetMarkerColor(1);
    gsample->SetLineColor(1);

    new TCanvas;
    gsample->Draw("apl");

    TGraph* gRe = new TGraph(N/2+1, sample, reX);
    gRe->SetNameTitle("gRe", "ReX");
    gRe->SetMarkerStyle(2);
    gRe->SetMarkerColor(2);

    new TCanvas;
    gRe->Draw("ap");

    TGraph* gIm = new TGraph(N/2+1, sample, imX);
    gIm->SetNameTitle("gIm", "ImX");
    gIm->SetMarkerStyle(2);
    gIm->SetMarkerColor(4);

    new TCanvas;
    gIm->Draw("ap");

    // polar coordinates

    Double_t magX[N];
    Double_t phaX[N];
    Polar(N, reX, imX, magX, phaX);

    // clean up the frequency domain

    for (int k=0; k<=N/2; ++k) {
        // if (k == 8 || k == 9 || k == 10 || k == 11) continue;
        if (k > 100) {
            reX[k] = 0;
            imX[k] = 0;
        }
        // if (magX[k] < 150) reX[k] = imX[k] = 0;
        if (magX[k] < 40) reX[k] = imX[k] = 0;
    }

    Double_t xinv[N];
    fft_inverse(N, reX, imX, xinv);

    TGraph* ginv = new TGraph(N, sample, xinv);
    ginv->SetNameTitle("ginv", "Time Domain from inverse;sample number;amplitude");
    ginv->SetMarkerStyle(7);
    ginv->SetMarkerColor(46);
    ginv->SetLineColor(46);

    new TCanvas;
    ginv->Draw("apl");

    TGraph* gmag = new TGraph(N/2+1, sample, magX);
    gmag->SetNameTitle("gmag", "Magnitude;frequency sample number;magnitude");
    gmag->SetMarkerStyle(2);
    gmag->SetMarkerColor(2);
    gmag->SetLineColor(2);

    new TCanvas;
    gmag->Draw("ap");

    TGraph* gpha = new TGraph(N/2+1, sample, phaX);
    gpha->SetNameTitle("gpha", "Phase;frequency sample number;phase");
    gpha->SetMarkerStyle(2);
    gpha->SetMarkerColor(4);
    gpha->SetLineColor(4);

    new TCanvas;
    gpha->Draw("ap");
}
