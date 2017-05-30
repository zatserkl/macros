// Andriy Zatserklyaniy <zatserkl@gmail.com> May 30, 2017

#include <TROOT.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TMath.h>
#include <TVirtualFFT.h>
#include <TGraph.h>

#include <iostream>
#include <sstream>

using std::cout;     using std::endl;

#ifndef macro_utils_C
const char* nextname(const char* base) {
    // if name 'base' already exists increaments it like h_1, h_2, etc.
    // returns pointer to circular buffer owned by ROOT function Form
    bool found = gDirectory->Get(base);
    if (!found) return base;
    Int_t icycle = 0;
    while (found) {
        std::stringstream ss;
        ss.str("");
        ss << base << "_" << ++icycle;
        //cout<< "         current hname: " << ss.str() <<endl;
        found = gDirectory->Get(ss.str().c_str());
    }
    //cout<< "new hname: " << Form("%s_%d",base,icycle) <<endl;
    return Form("%s_%d",base,icycle);
}
#endif

void fftturbine(const char *fname="T0005CH1.dat", Int_t Nset=0)
{
    TGraph *g = new TGraph(fname);

    if (!g && g->GetN() == 0) {
        cout<< "File not found: " << fname <<endl;
        return;
    }
    if (g->GetN() % 2 == 0) g->Set(g->GetN()-1);    // if the number points in the graph is odd set it to even
    g->SetName("g");
    g->SetMarkerStyle(6);

    new TCanvas;
    g->Draw("AP");

    Double_t deltat = (g->GetX()[g->GetN()-1] - g->GetX()[0]) / g->GetN();  // sampling step in s

    //
    // the number of samples
    //

    Int_t Ndata = g->GetN();        // the number of points in the real data
    Int_t N = Ndata;                // the number of samples
    if (Nset > 0) N = Nset;         // if Nset > N, the rest will be padded by zeros

    //
    // We will not change the deltat, but we may change the number of samples.
    // Two cases:
    // if the required number of samples Nset < Ndata, we will use Nset samples for transform
    // if the required number of samples Nset > Ndata, we will pad the rest with zeroes
    //

    cout<< "The number of samples N = " << N << " sampling step = " << deltat << " s" <<endl;

    Int_t NDIM = 2*(N/2 + 1);   // more than double of the maximum number of samples: a room for positive and negative frequencies

    // Make our own TVirtualFFT object (using option "K") for forward and inverse transforms.
    // They can be reused for the same size and type of transform.
    //
    // Third parameter (option) consists of 3 parts:
    // - transform type:
    //     real input/complex output in our case
    // - transform flag: 
    //     the amount of time spent in planning
    //     the transform (see TVirtualFFT class description)
    // - to create a new TVirtualFFT object (option "K") or use the global (default)

    TVirtualFFT *fft = TVirtualFFT::FFT(1, &N, "R2C ES K");     // forward transform
    TVirtualFFT* fft_inv = TVirtualFFT::FFT(1, &N, "C2R M K");  // inverse trasform

    // prepeare the data for the transform

    cout<< "The sampling frequency is 1/deltat = " << 1/deltat << " Hz. The Nyquist frequency is " << 0.5/deltat << " Hz" <<endl;

    TH1D* hInput = new TH1D(nextname("hinput"),Form("Input signal for %s;time, s",fname), N, 0, N*deltat);    // histo to visualize the input data

    Double_t* input = new Double_t[NDIM];
    for (int i=0; i<NDIM; ++i) input[i] = 0;            // fill the whole input array with zeros

    // fill the real data
    for (Int_t isample=0; isample<Ndata; isample++)
    {
        if (isample == N) break;                        // break if N < Ndata

        input[isample] = g->GetY()[isample];

        Int_t bin = isample+1;
        hInput->SetBinContent(bin,input[isample]);
    }

    new TCanvas;
    hInput->Draw();

    // set the input data and do the trasform

    fft->SetPoints(input);
    fft->Transform();

    // Allocate an array big enough to hold the transform output
    // Transform output in 1d contains, for a transform of size N, 
    // N/2+1 complex numbers, i.e. 2*(N/2+1) real numbers
    Double_t* re = new Double_t[NDIM];
    Double_t* im = new Double_t[NDIM];
    fft->GetPointsComplex(re, im);

    // for (int i=0; i<N/2+1; ++i) {cout<< i << "\tre = " << re[i] << "   " << im[i] <<endl;}   // array size is N/2 + 1

    Double_t* mag = new Double_t[NDIM];
    Double_t* phase = new Double_t[NDIM];
    for (int isample=0; isample<N/2+1; ++isample) {
        mag[isample] = TMath::Sqrt(re[isample]*re[isample] + im[isample]*im[isample]);
        phase[isample] = TMath::ATan2(im[isample],re[isample]);
        // cout<< "mag[" << isample << "] = " << mag[isample] <<endl;
    }

    // plot magnitude as samples in the frequency domain

    TH1D* hMagSamples = new TH1D(nextname("hMagSamples"),Form("Magnitude in samples for %s;samples",fname), N/2+1, 0, N/2+1);
    for (int isample=0; isample<N/2+1; ++isample) {
        hMagSamples->SetBinContent(isample+1, mag[isample]);
    }
    new TCanvas;
    hMagSamples->Draw();

    // plot magnitude as fraction of the sampling frequency

    // NB comment from the FFT.C
    // NOTE: for "real" frequencies you have to divide the x-axes range with the range of your function 
    // y-axes has to be rescaled by a factor of 1/SQRT(N) to be right: this is not done automatically!

    TH1D* hMagFraction = new TH1D(nextname("hMagFraction"),
                                  Form("Magnitude in the fraction of the sampling frequency for %s;fraction of the sampling frequency",fname),
                                  N/2+1, 0, Double_t(N/2+1)/N);
    for (int isample=0; isample<N/2+1; ++isample) {
        hMagFraction->SetBinContent(isample+1, mag[isample]/TMath::Sqrt(N));    // rescale by sqrt(N)
    }
    new TCanvas;
    hMagFraction->Draw();

    // plot magnitude as real world frequency

    TH1D* hMagFrequency = new TH1D(nextname("hMagFrequency"),
                                   Form("Magnitude in the real world frequency for %s;real world frequency, Hz",fname),
                                   N/2+1, 0, Double_t(N/2+1)/N/deltat);
    for (int isample=0; isample<N/2+1; ++isample) {
        hMagFrequency->SetBinContent(isample+1, mag[isample]/TMath::Sqrt(N));   // same values as in nMagFraction, just different x-axis scale
    }
    new TCanvas;
    hMagFrequency->Draw();

    // fill up the rest of the Re and Im array using Hermitian symmetry
    for (int isample=N/2+1; isample<N; ++isample) {
        re[isample] = re[N-isample];
        im[isample] = -im[N-isample];
    }

    //
    // inverse transform
    //

    fft_inv->SetPointsComplex(re, im);
    fft_inv->Transform();

    // get output of the inverse tranform and scale it to 1/N

    Double_t* output = fft_inv->GetPointsReal();
    for (int isample=0; isample<N; ++isample) {
        output[isample] /= Double_t(N);
        // cout<< "output[" << isample << "] = " << output[isample]/Double_t(N) <<endl;
    }

    TH1D* hOutput = new TH1D(nextname("hOutput"), Form("Output of the inverse transform for %s;time, s",fname), N, 0, Ndata*deltat);
    for (int i=0; i<N; ++i) hOutput->SetBinContent(i+1, output[i]);
    new TCanvas;
    hOutput->Draw();

    delete[] input;
    delete[] re;
    delete[] im;
    delete[] mag;
    delete[] phase;
}
