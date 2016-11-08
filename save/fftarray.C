#include "TROOT.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "TVirtualFFT.h"

#include <iostream>

using std::cout;     using std::endl;

void fftarray()
{
   // The number of samples
   Int_t N = 32;
   // Int_t N = 64;
   // Int_t N = 256;
   // Int_t N = 1024;
   // Int_t N = 8192;
   
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

   // create data

   // L is arbitrary length of the signal, seconds
   // Given the N and L the sampling frequency is N/L
   // The Nyquist frequency is 0.5*N/L

   Double_t L = 1.;   // length of the time domain signal in seconds

   cout<< "The sampling frequency is N/L = " << N << "/" << L << " = " << N/L << " Hz. The Nyquist frequency is " << 0.5*N/L << " Hz" <<endl;

   //
   // Use signal in form sin(2*pi*f*x) = sin(2*pi*x/T) where f is frequency in Hz (if x in seconds)
   //
   TF1 *fsin = new TF1("fsin", "sin(TMath::TwoPi()*1*x)+sin(TMath::TwoPi()*3*x)+sin(TMath::TwoPi()*5*x)+1", 0, L);
   fsin->SetNpx(20000);

   TH1D* hInput = new TH1D("hinput","Input signal;time, s", N, 0, L);   // histo to visualize the input data

   Double_t input[20000];
   for (Int_t isample=0; isample<N; isample++)
   {
      Double_t x = L*isample/N;
      Double_t val = fsin->Eval(x);
      input[isample] = val;

      Int_t bin = isample+1;
      hInput->SetBinContent(bin,val);
   }

   new TCanvas;
   hInput->Draw();

   // set the input data and do the trasform

   fft->SetPoints(input);
   fft->Transform();

   // Allocate an array big enough to hold the transform output
   // Transform output in 1d contains, for a transform of size N, 
   // N/2+1 complex numbers, i.e. 2*(N/2+1) real numbers
	Double_t re[20000];
	Double_t im[20000];
	fft->GetPointsComplex(re, im);

	// for (int i=0; i<N/2+1; ++i) {cout<< i << "\tre = " << re[i] << "   " << im[i] <<endl;}   // array size is N/2 + 1

	Double_t mag[20000];
	Double_t phase[20000];
	for (int isample=0; isample<N/2+1; ++isample) {
      mag[isample] = TMath::Sqrt(re[isample]*re[isample] + im[isample]*im[isample]);
      phase[isample] = TMath::ATan2(im[isample],re[isample]);
      // cout<< "mag[" << isample << "] = " << mag[isample] <<endl;
   }

   // plot magnitude as samples in the frequency domain

   TH1D* hMagSamples = new TH1D("hMagSamples","Magnitude in samples;samples", N/2+1, 0, N/2+1);
   for (int isample=0; isample<N/2+1; ++isample) hMagSamples->SetBinContent(isample+1, mag[isample]);
   new TCanvas;
   hMagSamples->Draw();

   // plot magnitude as fraction of the sampling frequency

   // NB comment from the FFT.C
   // NOTE: for "real" frequencies you have to divide the x-axes range with the range of your function 
   // (in this case L); y-axes has to be rescaled by a factor of 1/SQRT(N) to be right: this is not done automatically!

   TH1D* hMagFraction = new TH1D("hMagFraction","Magnitude in the fraction of the sampling frequency;fraction of the sampling frequency", N/2+1, 0, Double_t(N/2+1)/N);
   for (int isample=0; isample<N/2+1; ++isample) hMagFraction->SetBinContent(isample+1, mag[isample]/TMath::Sqrt(N));
   new TCanvas;
   hMagFraction->Draw();

   // plot magnitude as real world frequency

   TH1D* hMagFrequency = new TH1D("hMagFrequency","Magnitude in the real world frequency;real world frequency, Hz", N/2+1, 0, (N/2+1)/L);
   for (int isample=0; isample<N/2+1; ++isample) hMagFrequency->SetBinContent(isample+1, mag[isample]/TMath::Sqrt(N));
   new TCanvas;
   hMagFrequency->Draw();

   // cut the frequency spectrum above 4 Hz

   Double_t real_world_freq_cut = 4.;           // Hz
   Int_t sample_cut = real_world_freq_cut*L;    // sample #
   cout<< "real_world_freq_cut = " << real_world_freq_cut << " Hz, sample_cut = real_world_freq_cut*L = " << real_world_freq_cut << "*" << L << " = " << sample_cut << " samples" <<endl;

   for (int isample=0; isample<N/2+1; ++isample) {
      re[isample] = mag[isample]*TMath::Cos(phase[isample]);
      im[isample] = mag[isample]*TMath::Sin(phase[isample]);
      // cout<< "re[" << isample << "] = " << re[isample] << "\t im[" << isample << "] = " << im[isample] <<endl;

      if (isample > sample_cut) {
         re[isample] = 0;
         im[isample] = 0;
      }
   }
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

	TH1D* hOutput = new TH1D("hOutput", "Output of the inverse transform;time, s", N, 0, L);
	for (int i=0; i<N; ++i) hOutput->SetBinContent(i+1, output[i]);
	new TCanvas;
	hOutput->Draw();
}
