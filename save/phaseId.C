#include <TROOT.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMath.h>
#include <TVirtualFFT.h>
#include <TRandom3.h>

#include <vector>
#include <iostream>

using std::cout;    using std::endl;

// utils.C stuff
void zoom(Axis_t xmin=0, Axis_t xmax=0, Axis_t ymin=0, Axis_t ymax=0);
void resize(Int_t width=0, Int_t height=0, Bool_t internal=kFALSE);  // resizes current canvas

void phaseId_view()
{
    Double_t t[20000];
    Double_t w[20000];
    //-- Int_t N = 8002;
    Int_t N = 256;

    Double_t tmin = -3.2;            // ns
    Double_t tmax = +3.2;            // ns
    Double_t dt = (tmax - tmin)/(N-1);
    cout<< "tmin = " << tmin << " tmax = " << tmax << " dt = " << dt <<endl;

    //-- Double_t baseline = 0.020;
    Double_t baseline = 0.;

    for (int i=0; i<N; ++i) {
        t[i] = tmin + i*dt;     // set time axis
        w[i] = baseline;
    }

    cout<< "N = " << N << " t[N-1] - t[0] = " << t[N-1] << " - " << t[0] << " = " << t[N-1] - t[0] << " dt = " << dt <<endl;

    /// // signal: Gaussian
    /// Double_t signal_pmax = 0.200;           // V
    /// Double_t signal_sigma = 0.250;          // width at the base in ns
    /// Double_t signal_pos = 0;                // position of the signal (signal maximum)
    /// for (int i=0; i<N; ++i) w[i] += signal_pmax*TMath::Gaus(t[i], signal_pos, signal_sigma);

    /// // low freq signal
    /// Double_t freq = 0.100;                  // GHz for time in ns
    /// Double_t phase0 = 0;
    /// Double_t ampl = 0.020;
    /// for (int i=0; i<N; ++i) w[i] += ampl*TMath::Sin(2.*TMath::Pi()*freq*t[i] + phase0);

    Double_t freq = 0.0;        // GHz for time in ns
    Double_t phase0 = 0;
    Double_t ampl = 0.0;

    freq = 1.0;        // GHz for time in ns
    phase0 = 0;
    ampl = 0.010;
    for (int i=0; i<N; ++i) w[i] += ampl*TMath::Sin(2.*TMath::Pi()*freq*t[i] + phase0);

    freq = 1.5;        // GHz for time in ns
    phase0 = 0;
    ampl = 0.020;
    for (int i=0; i<N; ++i) w[i] += ampl*TMath::Sin(2.*TMath::Pi()*freq*t[i] + phase0);

    freq = 2.0;                // GHz for time in ns
    phase0 = 0;
    ampl = 0.020;
    for (int i=0; i<N; ++i) w[i] += ampl*TMath::Cos(2.*TMath::Pi()*freq*t[i] + phase0);

    freq = 2.5;                // GHz for time in ns
    phase0 = 0;
    ampl = 0.030;
    //-- for (int i=0; i<N; ++i) w[i] += ampl*TMath::Sin(2.*TMath::Pi()*freq*t[i] + phase0);
    for (int i=0; i<N; ++i) w[i] += ampl*TMath::Cos(2.*TMath::Pi()*freq*t[i] + phase0);          // Cos instead of Sin, phase0 = 0

    /// // Gaussian spread 1 mV
  
    /// TRandom3 rand;           // default parameter is UInt_t seed = 4357
    /// Double_t sig_mean = 0;          // mV
    /// Double_t sig_sigma = 0.001;     // mV
    /// for (int i=0; i<N; ++i) w[i] += rand.Gaus(sig_mean,sig_sigma);

    TGraph* graw = new TGraph(N, t, w);
    graw->SetNameTitle("graw", "Raw Data");
    graw->SetMarkerStyle(6);
    graw->SetMarkerColor(9);
    graw->SetLineColor(9);

    new TCanvas;
    graw->Draw("apl");

    /// // set minimum and maximum for the graph because of visual aliasing
    /// // graw->SetMinimum(-0.10);
    /// // graw->SetMaximum(0.20);
    /// // resize and zoom the region of interest
    /// resize(1200,300); zoom(-20,5);

    //
    // analyse with FFT
    //

    //
    // get Re and Im parts with the standard routine and syntese the original graph using sin and cos
    //

    TVirtualFFT* fft = TVirtualFFT::FFT(1, &N, "R2C K");         // Option K: keep own TVirtualFFT object. Need to be deleted at the end.

    fft->SetPoints(w);
    fft->Transform();

    Double_t re[20000];
    Double_t im[20000];

    const Double_t pattern = 123456789.;
    for (int i=0; i<20000; ++i) re[i] = im[i] = pattern;

    fft->GetPointsComplex(re, im);

    // find the number of non-zero terms
    Int_t nre = sizeof(re)/sizeof(Double_t) - 1;
    while (re[nre] == pattern && nre >= 0) --nre;
    cout<< "nre = " << nre <<endl;
    Int_t nim = sizeof(im)/sizeof(Double_t) - 1;
    while (im[nim] == pattern && nim >= 0) --nim;
    cout<< "nim = " << nim <<endl;

    cout<< "re[0] = " << re[0] << " re[N/2] = " << re[N/2] <<endl;
    cout<< "im[0] = " << im[0] << " im[N/2] = " << im[N/2] <<endl;
    cout<< "re[0]/N = " << re[0]/N << " re[N/2]/N = " << re[N/2]/N <<endl;
    cout<< "im[0]/N = " << im[0]/N << " im[N/2]/N = " << im[N/2]/N <<endl;

    //
    // plot coeffs
    //

    Double_t index[20000];
    for (int i=0; i<20000; ++i) index[i] = i;

    TGraph* gRe = new TGraph(N/2+1, index, re);
    gRe->SetNameTitle("gRe","Real part");
    gRe->SetMarkerStyle(6);
    gRe->SetMarkerColor(8);
    gRe->SetLineColor(8);
    new TCanvas;
    gRe->Draw("apl");
    TGraph* gIm = new TGraph(N/2+1, index, im);
    gIm->SetNameTitle("gIm","Imaginary part");
    gIm->SetMarkerStyle(6);
    gIm->SetMarkerColor(6);
    gIm->SetLineColor(6);
    new TCanvas;
    gIm->Draw("apl");

    Double_t out[20000];
    for (int i=0; i<N; ++i) {
        out[i] = (re[0] + re[N/2]) / N;
        for (int k=1; k<N/2; ++k) {
            out[i] += (re[k]*TMath::Cos(TMath::TwoPi()*k*i/N) - im[k]*TMath::Sin(TMath::TwoPi()*k*i/N)) / (N/2.);
        }
    }
    TGraph* gOut = new TGraph(N, index, out);
    //gOut->SetNameTitle("gOut","#sum_{k=1}^{#frac{N}{2}-1}[ReX[k]cos(2#pik#frac{n}{N}) + ImX[k]sin(2#pik#frac{n}{N})]");
    gOut->SetNameTitle("gOut","Reconstructed from re and im parts (using cos and sin)");
    gOut->SetMarkerStyle(6);
    gOut->SetMarkerColor(46);
    gOut->SetLineColor(46);
    new TCanvas;
    gOut->Draw("apl");

    // // I don't understand that: I thought it should be zero.
    // Double_t outIm[20000];
    // for (int i=0; i<N; ++i) {
    //     outIm[i] = (re[0] + re[N/2]) / N;
    //     for (int k=1; k<N/2; ++k) {
    //         outIm[i] += (re[k]*TMath::Sin(TMath::TwoPi()*k*i/N) + im[k]*TMath::Cos(TMath::TwoPi()*k*i/N)) / (N/2.);
    //     }
    // }
    // TGraph* gOutIm = new TGraph(N, index, outIm);
    // gOutIm->SetNameTitle("gOutIm","Imaginary part of reconstruction");
    // gOutIm->SetMarkerStyle(6);
    // gOutIm->SetMarkerColor(30);
    // gOutIm->SetLineColor(30);
    // new TCanvas;
    // gOutIm->Draw("apl");

    //
    // phase and amplitude
    //

    Double_t phase[20000];
    Double_t mag[20000];
    phase[0] = 0;
    phase[N/2] = 0;
    mag[0] = re[0]/N;
    mag[N/2] = re[N/2]/N;
    for (int k=1; k<N/2; ++k) {
        phase[k] = TMath::ATan2(im[k],re[k]);
        mag[k] = TMath::Sqrt((re[k]*re[k]+im[k]*im[k])) / (N/2);
    }

    Double_t outPhase[20000];
    for (int i=0; i<N; ++i) {
        outPhase[i] = 0;
        for (int k=0; k<N/2+1; ++k) {
            outPhase[i] += mag[k]*TMath::Cos(TMath::TwoPi()*k*i/N + phase[k]);
        }
    }
    TGraph* gOutPhase = new TGraph(N, index, outPhase);
    gOutPhase->SetNameTitle("gOutPhase","Reconstructed from mag and phase");
    gOutPhase->SetMarkerStyle(6);
    gOutPhase->SetMarkerColor(49);
    gOutPhase->SetLineColor(49);
    new TCanvas;
    gOutPhase->Draw("apl");

    TGraph* gmag = new TGraph(N/2+1, index, mag);
    gmag->SetNameTitle("gmag","cosine's mag");
    gmag->SetMarkerStyle(6);
    gmag->SetMarkerColor(38);
    gmag->SetLineColor(38);
    new TCanvas;
    gmag->Draw("apl");

    TGraph* gphase = new TGraph(N/2+1, index, phase);
    gphase->SetNameTitle("gphase","cosine's phase");
    gphase->SetMarkerStyle(6);
    gphase->SetMarkerColor(28);
    gphase->SetLineColor(28);
    new TCanvas;
    gphase->Draw("apl");

    //
    // inverse transform
    //

    TVirtualFFT* fft_inv = TVirtualFFT::FFT(1, &N, "C2R K");  // inverse trasform. Option K: keep own TVirtualFFT object
    fft_inv->SetPointsComplex(re, im);
    fft_inv->Transform();

    // get output of the inverse tranform and scale it to 1/N

    Double_t y[20000];                                          // reconstructed back

    Double_t* output_inv = fft_inv->GetPointsReal();
    for (int isample=0; isample<N; ++isample) {
        y[isample] = output_inv[isample]/Double_t(N);           // NB: divide to N
    }

    TGraph* gback = new TGraph(N, t, y);
    gback->SetNameTitle("gback", "Raw Data reconstructed back");
    gback->SetMarkerStyle(6);
    gback->SetMarkerColor(2);
    gback->SetLineColor(2);

    new TCanvas;
    gback->Draw("apl");

    delete fft;         // delete owned (option 'K') TVirtualFFT object
    delete fft_inv;     // delete owned (option 'K') TVirtualFFT object
}

void MagPhase(Int_t N, Double_t re[], Double_t im[], Double_t mag[], Double_t phase[])
{
    phase[0] = 0;
    phase[N/2] = 0;
    mag[0] = re[0]/N;
    mag[N/2] = re[N/2]/N;
    for (int k=1; k<N/2; ++k) {
        phase[k] = TMath::ATan2(im[k],re[k]);
        mag[k] = TMath::Sqrt((re[k]*re[k]+im[k]*im[k])) / (N/2);
    }
}

void AddFrequency(Int_t k, Int_t N, Double_t mag[], Double_t phase[], Double_t ampl[])
{
    for (int i=0; i<N; ++i) {
        ampl[i] += mag[k]*TMath::Cos(TMath::TwoPi()*k*i/N + phase[k]);
    }
}

void RecoMagPhase(Int_t N, Double_t mag[], Double_t phase[], Double_t ampl[])
{
    for (int i=0; i<N; ++i) {
        for (int k=0; k<N/2+1; ++k) {
            ampl[i] += mag[k]*TMath::Cos(TMath::TwoPi()*k*i/N + phase[k]);
        }
    }
}

Double_t DerivPhase(Int_t N, Double_t phaseSet[], Int_t k, Int_t np, Double_t deriv[])
{
    // N is the number of samples in time domain. The number of points in freq domain is N/2+1.
    // k is No. of sample in freq domain
    // np is the number of offsets with N/2+1 elements each
    // dim: phaseSet[N/2+1], deriv[np]

    Double_t diffAverage = 0;
    for (int ioffset=0; ioffset<np-1; ++ioffset) {
        Double_t diff = phaseSet[(ioffset+1)*(N/2+1) + k] - phaseSet[ioffset*(N/2+1) + k];
        if (diff < -5) diff += TMath::TwoPi();
        if (diff >  5) diff -= TMath::TwoPi();
        deriv[ioffset] = diff;
        diffAverage += diff;
    }
    diffAverage /= np-1;
    return diffAverage;
}

void phaseId()
{
    Double_t t[20000];
    Double_t w[20000];
    //-- Int_t Ndata = 8002;
    //-- Int_t Ndata = 256;
    Int_t Ndata = 512;

    // Double_t tmin = -3.2;            // ns
    // Double_t tmax = +3.2;            // ns
    Double_t tmin = -6.4;            // ns
    Double_t tmax = +6.4;            // ns
    Double_t dt = (tmax - tmin)/(Ndata-1);
    cout<< "tmin = " << tmin << " tmax = " << tmax << " dt = " << dt <<endl;

    //-- Double_t baseline = 0.020;
    Double_t baseline = 0.;

    for (int i=0; i<Ndata; ++i) {
        t[i] = tmin + i*dt;     // set time axis
        w[i] = baseline;
    }

    cout<< "Ndata = " << Ndata << " t[Ndata-1] - t[0] = " << t[Ndata-1] << " - " << t[0] << " = " << t[Ndata-1] - t[0] << " dt = " << dt <<endl;

    // signal: Gaussian
    Double_t signal_pmax = 0.200;           // V
    Double_t signal_sigma = 0.250;          // width at the base in ns
    //-- Double_t signal_pos = 0;                // position of the signal (signal maximum)
    Double_t signal_pos = -2.0;                // position of the signal (signal maximum)
    for (int i=0; i<Ndata; ++i) w[i] += signal_pmax*TMath::Gaus(t[i], signal_pos, signal_sigma);

    Double_t freq = 0.0;        // GHz for time in ns
    Double_t phase0 = 0;
    Double_t ampl = 0.0;

    /// // low freq signal
    /// freq = 0.100;                  // GHz for time in ns
    /// phase0 = 0;
    /// ampl = 0.020;
    /// for (int i=0; i<Ndata; ++i) w[i] += ampl*TMath::Sin(2.*TMath::Pi()*freq*t[i] + phase0);

    freq = 1.0;        // GHz for time in ns
    phase0 = 0;
    ampl = 0.010;
    for (int i=0; i<Ndata; ++i) w[i] += ampl*TMath::Sin(2.*TMath::Pi()*freq*t[i] + phase0);

    freq = 1.5;        // GHz for time in ns
    phase0 = 0;
    ampl = 0.020;
    for (int i=0; i<Ndata; ++i) w[i] += ampl*TMath::Sin(2.*TMath::Pi()*freq*t[i] + phase0);

    freq = 2.0;                // GHz for time in ns
    phase0 = 0;
    ampl = 0.020;
    for (int i=0; i<Ndata; ++i) w[i] += ampl*TMath::Cos(2.*TMath::Pi()*freq*t[i] + phase0);

    freq = 2.5;                // GHz for time in ns
    phase0 = 0;
    ampl = 0.030;
    //-- for (int i=0; i<Ndata; ++i) w[i] += ampl*TMath::Sin(2.*TMath::Pi()*freq*t[i] + phase0);
    for (int i=0; i<Ndata; ++i) w[i] += ampl*TMath::Cos(2.*TMath::Pi()*freq*t[i] + phase0);          // Cos instead of Sin, phase0 = 0

    /// // Gaussian spread 1 mV
  
    /// TRandom3 rand;           // default parameter is UInt_t seed = 4357
    /// Double_t sig_mean = 0;          // mV
    /// Double_t sig_sigma = 0.001;     // mV
    /// for (int i=0; i<Ndata; ++i) w[i] += rand.Gaus(sig_mean,sig_sigma);

    TGraph* graw = new TGraph(Ndata, t, w);
    //-- TGraph* graw = new TGraph(256, t, w);
    graw->SetNameTitle("graw", "Raw Data");
    graw->SetMarkerStyle(6);
    graw->SetMarkerColor(9);
    graw->SetLineColor(9);

    new TCanvas;
    graw->Draw("apl");

    /// // set minimum and maximum for the graph because of visual aliasing
    /// // graw->SetMinimum(-0.10);
    /// // graw->SetMaximum(0.20);
    /// // resize and zoom the region of interest
    /// resize(1200,300); zoom(-20,5);

    //
    // analyse with FFT
    //

    //
    // get Re and Im parts with the standard routine and syntese the original graph using sin and cos
    //

    //
    // Set N to 256
    //
    Int_t N = 256;

    Int_t np1ns = 40;           // 40 points correspond to 1 ns

    TVirtualFFT* fft = TVirtualFFT::FFT(1, &N, "R2C K");         // Option K: keep own TVirtualFFT object. Need to be deleted at the end.

    Double_t* phase1ns = new Double_t[(N/2+1)*np1ns];
    Double_t* mag1ns = new Double_t[(N/2+1)*np1ns];

    Double_t re[20000];
    Double_t im[20000];

    for (int ioffset=0; ioffset<np1ns; ++ioffset)
    {
        fft->SetPoints(&w[ioffset]);
        fft->Transform();
        fft->GetPointsComplex(re, im);

        MagPhase(N, re, im, &mag1ns[ioffset*(N/2+1)], &phase1ns[ioffset*(N/2+1)]);
    }

    Double_t index[20000];
    for (int i=0; i<20000; ++i) index[i] = i;

    // plot the initial point

    Double_t* magPtr = &mag1ns[(N/2+1)*0];          // magnitude for the no offset
    Double_t* phasePtr = &phase1ns[(N/2+1)*0];      // phase for the no offset;

    TGraph* gmag = new TGraph(N/2+1, index, magPtr);
    gmag->SetNameTitle("gmag","cosine's mag");
    gmag->SetMarkerStyle(6);
    gmag->SetMarkerColor(38);
    gmag->SetLineColor(38);
    new TCanvas;
    gmag->Draw("apl");

    TGraph* gphase = new TGraph(N/2+1, index, phasePtr);
    gphase->SetNameTitle("gphase","cosine's phase");
    gphase->SetMarkerStyle(6);
    gphase->SetMarkerColor(28);
    gphase->SetLineColor(28);
    new TCanvas;
    gphase->Draw("apl");

    // reconstruct amplitude from the mag and phase

    Double_t outPhase[20000];
    for (unsigned i=0; i<sizeof(outPhase)/sizeof(Double_t); ++i) outPhase[i] = 0;
    RecoMagPhase(N, magPtr, phasePtr, outPhase);
    for (int i=0; i<10; ++i) {
        cout<< i << "\t" << outPhase[i] <<endl;
    }

    TGraph* gOutPhase = new TGraph(N, index, outPhase);
    gOutPhase->SetNameTitle("gOutPhase","Reconstructed from mag and phase");
    gOutPhase->SetMarkerStyle(6);
    gOutPhase->SetMarkerColor(49);
    gOutPhase->SetLineColor(49);
    new TCanvas;
    gOutPhase->Draw("apl");

    cout<< "gOutPhase->GetN() = " << gOutPhase->GetN() <<endl;

    // phase vs offset for phase point #16

    Double_t mag[20000];
    Double_t phase[20000];

    Int_t kFreq = 0;
    Double_t diffAverage = 0;
    Double_t deriv[20000];
    
    kFreq = 16;

    diffAverage = DerivPhase(N, phase1ns, kFreq, np1ns, deriv);
    cout<< "diffAverage for kFreq = " << kFreq << " is " << diffAverage <<endl;

    for (int i=0; i<np1ns; ++i) {
        mag[i] = mag1ns[i*(N/2+1) + kFreq];
        phase[i] = phase1ns[i*(N/2+1) + kFreq];
    }

    TGraph* gMag16 = new TGraph(np1ns, index, mag);
    gMag16->SetNameTitle("gMag16","Magnitude vs offset for kFreq = 16");
    gMag16->SetMarkerStyle(6);
    gMag16->SetMarkerColor(9);
    gMag16->SetLineColor(9);
    new TCanvas;
    gMag16->Draw("apl");

    TGraph* gPhase16 = new TGraph(np1ns, index, phase);
    gPhase16->SetNameTitle("gPhase16","Phase vs offset for kFreq = 16");
    gPhase16->SetMarkerStyle(6);
    gPhase16->SetMarkerColor(2);
    gPhase16->SetLineColor(2);
    new TCanvas;
    gPhase16->Draw("apl");
    
    kFreq = 20;

    diffAverage = DerivPhase(N, phase1ns, kFreq, np1ns, deriv);
    cout<< "diffAverage for kFreq = " << kFreq << " is " << diffAverage <<endl;

    for (int i=0; i<np1ns; ++i) {
        mag[i] = mag1ns[i*(N/2+1) + kFreq];
        phase[i] = phase1ns[i*(N/2+1) + kFreq];
    }

    TGraph* gMag20 = new TGraph(np1ns, index, mag);
    gMag20->SetNameTitle("gMag20","Magnitude vs offset for kFreq = 20");
    gMag20->SetMarkerStyle(6);
    gMag20->SetMarkerColor(9);
    gMag20->SetLineColor(9);
    new TCanvas;
    gMag20->Draw("apl");

    TGraph* gPhase20 = new TGraph(np1ns, index, phase);
    gPhase20->SetNameTitle("gPhase20","Phase vs offset for kFreq = 20");
    gPhase20->SetMarkerStyle(6);
    gPhase20->SetMarkerColor(2);
    gPhase20->SetLineColor(2);
    new TCanvas;
    gPhase20->Draw("apl");

    Double_t diff[20000];
    for (int k=0; k<N/2+1; ++k) {
        diff[k] = DerivPhase(N, phase1ns, k, np1ns, deriv);
    }

    TGraph* gDiff = new TGraph(N/2+1, index, diff);
    gDiff->SetNameTitle("gDiff","Diff vs freq sample;k, frequency sample");
    gDiff->SetMarkerStyle(7);
    gDiff->SetMarkerColor(6);
    gDiff->SetLineColor(6);
    new TCanvas;
    gDiff->Draw("apl");

    // phase space (phase derivative)x(mag)

    Double_t mag_diff[20000];
    for (int k=0; k<N/2+1; ++k) {
        mag_diff[k] = magPtr[k]*diff[k];
    }

    TGraph* gMagDiff = new TGraph(N/2+1, index, mag_diff);
    gMagDiff->SetNameTitle("gMagDiff","phase space (phase derivative)x(mag) vs freq sample;k, frequency sample");
    gMagDiff->SetMarkerStyle(7);
    gMagDiff->SetMarkerColor(36);
    gMagDiff->SetLineColor(36);
    new TCanvas;
    gMagDiff->Draw("apl");

    //
    // threshold value for average derivative of the phase
    //
    Double_t thres_deriv = 0.2;

    std::vector<Int_t> vfreq;
    for (int k=0; k<np1ns; ++k) if (diff[k] > thres_deriv) vfreq.push_back(k); 

    // output array for reconstructed amplitude
    Double_t out[20000];
    for (unsigned i=0; i<sizeof(out)/sizeof(Double_t); ++i) out[i] = 0;     // NB: all arrays in ROOT are global

    for (unsigned k=0; k<vfreq.size(); ++k) {
        AddFrequency(vfreq[k], N, mag1ns, phase1ns, out);
    }

    TGraph* gReco = new TGraph(N, index, out);
    gReco->SetNameTitle("gReco",Form("Reco from highest %lu frequencies",vfreq.size()));
    gReco->SetMarkerStyle(7);
    gReco->SetMarkerColor(6);
    gReco->SetLineColor(2);
    new TCanvas;
    gReco->Draw("apl");

    // subrtract the noise from the signal

    TGraph* gsub = (TGraph*) gOutPhase->Clone("gsub");
    gsub->SetTitle("Signal with noise subtracted");
    for (int i=0; i<gsub->GetN(); ++i) {
        gsub->SetPoint(i, gsub->GetX()[i], gsub->GetY()[i]-gReco->GetY()[i]);
    }
    new TCanvas;
    gsub->Draw("apl");

    delete fft;         // delete owned (option 'K') TVirtualFFT object
}
