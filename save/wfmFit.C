// Andriy Zatserklyaniy <zatserkl@fnal.gov> Feb 27, 2016

#ifndef macro_wfmFit_C
#define macro_wfmFit_C

#include "FunRC.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TVirtualFFT.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TObject.h>

#include <TEnv.h>
#include <TList.h>
#include <TNamed.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

using std::cout;     using std::endl;

// utils.C stuff
TGraph* gtemp(TCanvas* can=0, Bool_t updateNameTitle=kTRUE);
TGraph* gtemp(Int_t index, TCanvas* can, Bool_t updateNameTitle=kTRUE);
TH1* htemp(TCanvas* can=0);
TGraph* gtemp(TCanvas* can, Bool_t updateNameTitle) {return gtemp(-1, can, updateNameTitle);}
TGraph* gtemp(Int_t index, TCanvas* can, Bool_t updateNameTitle)
{
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
        cout<< "No canvas found" <<endl;
        return 0;
    }
    if (can == 0) can = (TCanvas*) gPad;

    TObjArray graphs;
    TIter next(can->GetListOfPrimitives());
    TObject* obj = 0;
    while ((obj = next())) {
        if (obj->IsA()->InheritsFrom(TGraph::Class())) graphs.Add(obj);
    }
    // if (graphs.GetEntries() == 0) {
    //    cout<< "Could not find TGraph object. ";
    //    lstemp(can);
    //    return 0;
    // }

    TGraph* gr = 0;

    if (index < 0) gr = (TGraph*) graphs.Last();
    else if (index < graphs.GetEntries()) gr = (TGraph*) graphs[index];

    if (gr && updateNameTitle) {
        //
        // Try to correct the graph's name and title in case the graph was created by TTree::Draw
        //
        // if the gtemp was created by TTree::Draw
        // 1) it has name and title "Graph"
        // 2) there is an empty histogram with hame htemp and meaninful title
        // If so, assign the htemp title to the graph
        if (strcmp(gr->GetName(), "Graph")==0 && strcmp(gr->GetTitle(), "Graph")==0) {
            TH1* h = htemp(0);
            if (h && strcmp(h->GetName(),"htemp")==0 && h->GetEntries() == 0) {
                gr->SetName("gtemp");
                gr->SetTitle(h->GetTitle());
            }
        }
    }

    return gr;
}
TGraph* gtemp(Int_t index, TCanvas* can, const char* name, const char* title)
{
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
        cout<< "No canvas found" <<endl;
        return 0;
    }
    if (can == 0) can = (TCanvas*) gPad;

    TObjArray graphs;
    TIter next(can->GetListOfPrimitives());
    TObject* obj = 0;
    while ((obj = next())) {
        if (obj->IsA()->InheritsFrom(TGraph::Class())) graphs.Add(obj);
    }

    TGraph* gr = 0;

    if (index < 0) gr = (TGraph*) graphs.Last();
    else if (index < graphs.GetEntries()) gr = (TGraph*) graphs[index];

    if (gr) {
        //
        // Try to correct the graph's name and title in case the graph was created by TTree::Draw
        //
        // if the gtemp was created by TTree::Draw
        // 1) it has name and title "Graph"
        // 2) there is an empty histogram with hame htemp and meaninful title
        // If so, assign the htemp title to the graph
        if (strcmp(gr->GetName(), "Graph")==0 && strcmp(gr->GetTitle(), "Graph")==0) {
            TH1* h = htemp(0);
            if (h && strcmp(h->GetName(),"htemp")==0 && h->GetEntries() == 0) {
                gr->SetName("gtemp");
                gr->SetTitle(h->GetTitle());
            }
        }
        if (name && *name) gr->SetName(name);
        if (title && *title) gr->SetTitle(title);
    }

    return gr;
}
TH1* htemp(Int_t index, TCanvas* can)
{
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
        cout<< "No canvas found" <<endl;
        return 0;
    }
    if (can == 0) can = (TCanvas*) gPad;

    TObjArray histos;
    TIter next(can->GetListOfPrimitives());
    TObject* obj = 0;
    while ((obj = next())) {
        if (obj->IsA()->InheritsFrom(TH1::Class())) histos.Add(obj);
    }
    // if (histos.GetEntries() == 0) {
    //    cout<< "Could not find TH1 object. ";
    //    lstemp(can);
    //    return 0;
    // }

    if (index < 0) return (TH1*) histos.Last();

    if (index < histos.GetEntries()) return (TH1*) histos[index];
    else return 0;
}
TH1* htemp(TCanvas* can) {return htemp(-1,can);}

class TextObject {
    //
    // To write source code into TTree UserInfo Object.
    // Also use static method TextObject::print(list) to print the content of the list.
    //
public:
    TList* list;
    TextObject(TList* list_0): list(list_0) {
        // write command line
        const char* fname = gEnv->GetValue("Rint.History", "");
        // cout<< "Rint.History file is " << fname <<endl;
        std::ifstream file(fname);
        std::string line = "N/A";
        if (file) {
            std::string line_read;
            while (std::getline(file, line_read)) line = line_read;
        }
        else cout<< "\n***Warning TextObject::TextObject: could not find history file " << fname <<endl<<endl;
        // cout<< "Command line is: " << line <<endl;

        // create an object and write it into the tree
        TNamed* textObject = new TNamed("Command line", line.c_str());
        list->AddLast(textObject);
        file.close();
    }
    void writeFile(const char* fname) {
        std::ifstream file(fname);
        if (file) {
            std::stringstream ss;
            std::string line;
            while (std::getline(file, line)) ss << line << endl;

            // create an object and write it into the tree
            TNamed* textObject = new TNamed(fname, ss.str());
            list->AddLast(textObject);
            file.close();
        }
        else cout<< "\n***Warning TextObject::TextObject: could not find file " << fname <<endl<<endl;
    }
    static void print(TList* the_list)
    {
        cout<< "the number of entries in the list is " << the_list->GetEntries() <<endl;

        TIter next(the_list);
        TNamed* object;
        while ((object = (TNamed*) next())) {
            cout<< "\nobject->GetName() = " << object->GetName() <<endl;
            cout<< object->GetTitle() <<endl;
        }
    }
};


namespace TreeRes {
    Int_t npar[4];
    Double_t par0[4];
    Double_t par1[4];
    Double_t par2[4];
    Double_t par3[4];
    Double_t par4[4];
    Double_t par5[4];
    Double_t ts[4];
    Double_t pmax[4];
    Double_t pmaxx[4];
    Double_t pmax_raw[4];
    Double_t bkg_raw[4];
    Double_t ebkg_raw[4];
    Double_t bkg[4];
    Double_t ebkg[4];
    Double_t fwhm[4];
    Double_t fwqm[4];
    Double_t q[4];
    Double_t xm[4];
    Double_t ym[4];
    Double_t x5mV[4];
    Double_t x10mV[4];
    Double_t tx[4];
    Double_t cf10[4];
    Double_t cf15[4];
    Double_t cf20[4];
    Double_t cf25[4];
    Double_t cf30[4];
    Double_t cf35[4];
    Double_t cf40[4];
    Double_t cf45[4];
    Double_t cf50[4];
    Double_t cf55[4];
    Double_t cf60[4];
    Double_t dwdt10[4];
    Double_t dwdt15[4];
    Double_t dwdt20[4];
    Double_t dwdt25[4];
    Double_t dwdt30[4];
    Double_t dwdt35[4];
    Double_t dwdt40[4];
    Double_t dwdt45[4];
    Double_t dwdt50[4];
    Double_t dwdt55[4];
    Double_t dwdt60[4];
    Double_t dwdtfrac[4];

    void clear() {
        for (int ich=0; ich<4; ++ich) {
            npar[ich] = 0;
            ts[ich] = 0;
            pmax[ich] = 0;
            pmaxx[ich] = 0;
            pmax_raw[ich] = 0;
            bkg_raw[ich] = 0;
            ebkg_raw[ich] = 0;
            bkg[ich] = 0;
            ebkg[ich] = 0;
            fwhm[ich] = 0;
            fwqm[ich] = 0;
            q[ich] = 0;
            xm[ich] = 0;
            ym[ich] = 0;
            x5mV[ich] = 0;
            x10mV[ich] = 0;
            tx[ich] = 0;
            par0[ich] = 0;
            par1[ich] = 0;
            par2[ich] = 0;
            par3[ich] = 0;
            par4[ich] = 0;
            par5[ich] = 0;
            cf10[ich] = 0;
            cf15[ich] = 0;
            cf20[ich] = 0;
            cf25[ich] = 0;
            cf30[ich] = 0;
            cf35[ich] = 0;
            cf40[ich] = 0;
            cf45[ich] = 0;
            cf50[ich] = 0;
            cf55[ich] = 0;
            cf60[ich] = 0;
            dwdt10[ich] = 0;
            dwdt15[ich] = 0;
            dwdt20[ich] = 0;
            dwdt25[ich] = 0;
            dwdt30[ich] = 0;
            dwdt35[ich] = 0;
            dwdt40[ich] = 0;
            dwdt45[ich] = 0;
            dwdt50[ich] = 0;
            dwdt55[ich] = 0;
            dwdt60[ich] = 0;
            dwdtfrac[ich] = 0;
        }
    }

    void book(TTree* tree) {
        tree->Branch("npar",          &npar,         "npar[4]/I");
        tree->Branch("par0",          &par0,         "par0[4]/D");
        tree->Branch("par1",          &par1,         "par1[4]/D");
        tree->Branch("par2",          &par2,         "par2[4]/D");
        tree->Branch("par3",          &par3,         "par3[4]/D");
        tree->Branch("par4",          &par4,         "par4[4]/D");
        tree->Branch("par5",          &par5,         "par5[4]/D");
        tree->Branch("ts",          &ts,         "ts[4]/D");
        tree->Branch("pmax",          &pmax,         "pmax[4]/D");
        tree->Branch("pmaxx",          &pmaxx,         "pmaxx[4]/D");
        tree->Branch("pmax_raw",          &pmax_raw,         "pmax_raw[4]/D");
        tree->Branch("bkg_raw",          &bkg_raw,         "bkg_raw[4]/D");
        tree->Branch("ebkg_raw",          &ebkg_raw,         "ebkg_raw[4]/D");
        tree->Branch("bkg",          &bkg,         "bkg[4]/D");
        tree->Branch("ebkg",          &ebkg,         "ebkg[4]/D");
        tree->Branch("fwhm",          &fwhm,         "fwhm[4]/D");
        tree->Branch("fwqm",          &fwqm,         "fwqm[4]/D");
        tree->Branch("q",          &q,         "q[4]/D");
        tree->Branch("xm",          &xm,         "xm[4]/D");
        tree->Branch("ym",          &ym,         "ym[4]/D");
        tree->Branch("x5mV",          &x5mV,         "x5mV[4]/D");
        tree->Branch("x10mV",          &x10mV,         "x10mV[4]/D");
        tree->Branch("tx",          &tx,         "tx[4]/D");
        tree->Branch("cf10",          &cf10,         "cf10[4]/D");
        tree->Branch("cf15",          &cf15,         "cf15[4]/D");
        tree->Branch("cf20",          &cf20,         "cf20[4]/D");
        tree->Branch("cf25",          &cf25,         "cf25[4]/D");
        tree->Branch("cf30",          &cf30,         "cf30[4]/D");
        tree->Branch("cf35",          &cf35,         "cf35[4]/D");
        tree->Branch("cf40",          &cf40,         "cf40[4]/D");
        tree->Branch("cf45",          &cf45,         "cf45[4]/D");
        tree->Branch("cf50",          &cf50,         "cf50[4]/D");
        tree->Branch("cf55",          &cf55,         "cf55[4]/D");
        tree->Branch("cf60",          &cf60,         "cf60[4]/D");
        tree->Branch("dwdt10",          &dwdt10,         "dwdt10[4]/D");
        tree->Branch("dwdt15",          &dwdt15,         "dwdt15[4]/D");
        tree->Branch("dwdt20",          &dwdt20,         "dwdt20[4]/D");
        tree->Branch("dwdt25",          &dwdt25,         "dwdt25[4]/D");
        tree->Branch("dwdt30",          &dwdt30,         "dwdt30[4]/D");
        tree->Branch("dwdt35",          &dwdt35,         "dwdt35[4]/D");
        tree->Branch("dwdt40",          &dwdt40,         "dwdt40[4]/D");
        tree->Branch("dwdt45",          &dwdt45,         "dwdt45[4]/D");
        tree->Branch("dwdt50",          &dwdt50,         "dwdt50[4]/D");
        tree->Branch("dwdt55",          &dwdt55,         "dwdt55[4]/D");
        tree->Branch("dwdt60",          &dwdt60,         "dwdt60[4]/D");
        tree->Branch("dwdtfrac",          &dwdtfrac,         "dwdtfrac[4]/D");
    }

    void connect(TTree* tree) {
        tree->SetBranchAddress("npar",          &npar);
        tree->SetBranchAddress("par0",          &par0);
        tree->SetBranchAddress("par1",          &par1);
        tree->SetBranchAddress("par2",          &par2);
        tree->SetBranchAddress("par3",          &par3);
        tree->SetBranchAddress("par4",          &par4);
        tree->SetBranchAddress("par5",          &par5);
        tree->SetBranchAddress("ts",          &ts);
        tree->SetBranchAddress("pmax",          &pmax);
        tree->SetBranchAddress("pmaxx",          &pmaxx);
        tree->SetBranchAddress("pmax_raw",          &pmax_raw);
        tree->SetBranchAddress("bkg_raw",          &bkg_raw);
        tree->SetBranchAddress("ebkg_raw",          &ebkg_raw);
        tree->SetBranchAddress("bkg",          &bkg);
        tree->SetBranchAddress("ebkg",          &ebkg);
        tree->SetBranchAddress("fwhm",          &fwhm);
        tree->SetBranchAddress("fwqm",          &fwqm);
        tree->SetBranchAddress("q",          &q);
        tree->SetBranchAddress("xm",          &xm);
        tree->SetBranchAddress("ym",          &ym);
        tree->SetBranchAddress("x5mV",          &x5mV);
        tree->SetBranchAddress("x10mV",          &x10mV);
        tree->SetBranchAddress("tx",          &tx);
        tree->SetBranchAddress("cf10",          &cf10);
        tree->SetBranchAddress("cf15",          &cf15);
        tree->SetBranchAddress("cf20",          &cf20);
        tree->SetBranchAddress("cf25",          &cf25);
        tree->SetBranchAddress("cf30",          &cf30);
        tree->SetBranchAddress("cf35",          &cf35);
        tree->SetBranchAddress("cf40",          &cf40);
        tree->SetBranchAddress("cf45",          &cf45);
        tree->SetBranchAddress("cf50",          &cf50);
        tree->SetBranchAddress("cf55",          &cf55);
        tree->SetBranchAddress("cf60",          &cf60);
        tree->SetBranchAddress("dwdt10",          &dwdt10);
        tree->SetBranchAddress("dwdt15",          &dwdt15);
        tree->SetBranchAddress("dwdt20",          &dwdt20);
        tree->SetBranchAddress("dwdt25",          &dwdt25);
        tree->SetBranchAddress("dwdt30",          &dwdt30);
        tree->SetBranchAddress("dwdt35",          &dwdt35);
        tree->SetBranchAddress("dwdt40",          &dwdt40);
        tree->SetBranchAddress("dwdt45",          &dwdt45);
        tree->SetBranchAddress("dwdt50",          &dwdt50);
        tree->SetBranchAddress("dwdt55",          &dwdt55);
        tree->SetBranchAddress("dwdt60",          &dwdt60);
        tree->SetBranchAddress("dwdtfrac",          &dwdtfrac);
    }
}  // namespace TreeRes

//--------- Fast Fourier Transform -----------

Double_t fft_filter(Int_t N, Double_t* x, Double_t* y, Double_t fcut=0.2)  // fcut in GHz
{
    // NB: array y is in use for both input and output

    TVirtualFFT* fft = TVirtualFFT::FFT(1, &N, "R2C ES K");

    fft->SetPoints(y);
    fft->Transform();

    Double_t re[20000];
    Double_t im[20000];
    fft->GetPointsComplex(re, im);

    Double_t baseline = re[0]/N;
    // cout<< "re[0]/N = " << re[0]/N <<endl;

    // cut frequency band between the real world frequencies fcut1 and fcut2
    Double_t L = x[N-1] - x[0];
    Int_t scut = fcut*L;             // convert from real world frequency to frequency sample number

    for (int isample=0; isample<N/2+1; ++isample) {
        if (isample > scut) {
            re[isample] = 0;
            im[isample] = 0;
        }
    }

    // inverse transform
    TVirtualFFT* fft_inv = TVirtualFFT::FFT(1, &N, "C2R M K");  // inverse trasform
    fft_inv->SetPointsComplex(re, im);
    fft_inv->Transform();

    // get output of the inverse tranform and scale it to 1/N

    Double_t* output_inv = fft_inv->GetPointsReal();
    for (int isample=0; isample<N; ++isample) {
        y[isample] = output_inv[isample]/Double_t(N);
    }

    return baseline;     // re[0]/N
}

//--------------- moving average -------------------

void moving_average(Int_t naver, Int_t np, const Double_t* y, Double_t* yaver)   // Double_t version
{
    // assume naver odd
    if (naver % 2 == 0) {
        naver += 1;
        cout<< "\n--> set naver to odd number " << naver <<endl;
    }

    Double_t acc = 0;
    for (int i=0; i<naver; ++i) acc += y[i];

    yaver[0] = acc/naver;
    Int_t icurr = naver/2;
    for (int i=0; i<=icurr; ++i) yaver[i] = y[0];

    Int_t ifirst = 0;
    Int_t inext = naver;

    while (inext<np) {
        acc += y[inext++] - y[ifirst++];
        yaver[++icurr] = acc/naver;
    }

    for (int i=icurr+1; i<np; ++i) yaver[i] = yaver[icurr];
}

TGraph* filter_graph(Double_t fcut=0.1, TGraph* g=0)  // fcut in GHz
{
    if (!g) g = gtemp();

    Int_t N = g->GetN();

    TVirtualFFT* fft = TVirtualFFT::FFT(1, &N, "R2C ES K");

    fft->SetPoints(g->GetY());
    fft->Transform();

    Double_t re[20000];
    Double_t im[20000];
    fft->GetPointsComplex(re, im);

    Double_t baseline = re[0]/N;
    //-- cout<< "re[0]/N = " << re[0]/N <<endl;

    Double_t mag[20000];
    // Double_t phase[20000];

    // cut frequency band between the real world frequencies fcut1 and fcut2
    Double_t L = g->GetX()[g->GetN()-1] - g->GetX()[0];
    Int_t scut = fcut*L;             // convert from real world frequency to frequency sample number

    for (int isample=0; isample<N/2+1; ++isample) {
        // cut the frequency band
        if (isample > scut) {
            re[isample] = 0;
            im[isample] = 0;
        }
        mag[isample] = TMath::Sqrt(re[isample]*re[isample] + im[isample]*im[isample]);
        // phase[isample] = TMath::ATan2(im[isample],re[isample]);
    }

    TObject* hobj = gROOT->FindObjectAny("hmag");
    if (hobj && hobj->IsA()->InheritsFrom(TH1::Class())) delete hobj;
    TH1D* hmag = new TH1D("hmag", Form("Magnitude for %s",g->GetName()), N/2+1, 0, (N/2+1)/L);
    for (int i=0; i<N/2+1; ++i) hmag->SetBinContent(i+1, mag[i]/TMath::Sqrt(N));
    // new TCanvas;
    // hmag->Draw();

    // inverse transform
    TVirtualFFT* fft_inv = TVirtualFFT::FFT(1, &N, "C2R M K");  // inverse trasform
    fft_inv->SetPointsComplex(re, im);
    fft_inv->Transform();

    // get output of the inverse tranform and scale it to 1/N

    Double_t output[20000];

    Double_t* output_inv = fft_inv->GetPointsReal();
    for (int isample=0; isample<N; ++isample) {
        output[isample] = output_inv[isample]/Double_t(N);
        // cout<< "output[" << isample << "] = " << output[isample]/Double_t(N) <<endl;
    }

    TGraph* gout = new TGraph(g->GetN(), g->GetX(), output);
    gout->SetNameTitle(Form("%s_f1_%0.2f",g->GetName(),fcut), Form("%s f1=%0.2f",g->GetTitle(),fcut));
    gout->SetMarkerColor(46);
    gout->SetLineColor(46);

    // new TCanvas;
    // gout->Draw("apl");
    gout->Draw("pl");

    Double_t xbkg_max = -10.;
    Int_t nbkg = 0;
    Double_t baseline_diff = 0;
    for (int isample=0; isample<gout->GetN(); ++isample) {
        if (gout->GetX()[isample] > xbkg_max) break;
        Double_t diff = gout->GetY()[isample] - baseline;
        baseline_diff += diff*diff;
        ++nbkg;
    }
    baseline_diff = TMath::Sqrt(baseline_diff/nbkg);
    //-- cout<< "average diff with baseline = " << baseline_diff <<endl;

    return gout;
}

Double_t analyseFFT(Double_t xmin, Double_t xmax, TGraph* g=0, TH1D* hmag=0, TH1D* hphase=0)
{
    if (!g) g = gtemp();

    Int_t xmin_i = -1;
    Int_t xmax_i = -1;
    Int_t N = 0;

    if (!g) g = gtemp();
    if (!g) {
        cout<< "Could not find initial graph" <<endl;
        return 0;
    }

    for (int i=0; i<g->GetN(); ++i) {
        if (g->GetX()[i] < xmin) continue;
        if (xmin_i < 0) xmin_i = i;
        if (g->GetX()[i] > xmax) break;
        xmax_i = i;
        ++N;
    }

    TVirtualFFT* fft = TVirtualFFT::FFT(1, &N, "R2C ES K");

    fft->SetPoints(&g->GetY()[xmin_i]);
    fft->Transform();

    Double_t re[20000];
    Double_t im[20000];
    fft->GetPointsComplex(re, im);

    Double_t baseline = re[0]/N;
    cout<< "re[0]/N = " << re[0]/N <<endl;

    Double_t L = g->GetX()[xmax_i] - g->GetX()[xmin_i];

    // delete histogram hmag if it exists
    TObject* hobj = gROOT->FindObjectAny("hmag");
    if (hobj && hobj->IsA()->InheritsFrom(TH1::Class())) delete hobj;
    // create a new histogram hmag
    hmag = new TH1D("hmag", Form("Magnitude for %s;frequency, GHz",g->GetName()), N/2+1, 0, (N/2+1)/L);

    // delete histogram hphase if it exists
    hobj = gROOT->FindObjectAny("hphase");
    if (hobj && hobj->IsA()->InheritsFrom(TH1::Class())) delete hobj;
    // create a new histogram hphase
    hphase = new TH1D("hphase", Form("Phase for %s;frequency, GHz",g->GetName()), N/2+1, 0, (N/2+1)/L);

    for (int isample=0; isample<N/2+1; ++isample) {
        Double_t mag = TMath::Sqrt(re[isample]*re[isample] + im[isample]*im[isample]);
        Double_t phase = TMath::ATan2(im[isample],re[isample]);
        hmag->SetBinContent(isample+1, mag/TMath::Sqrt(N));
        hphase->SetBinContent(isample+1, phase);
    }

    return baseline;
}

//-------------------------- end of tools --------------------------

/// Double_t peakFinder(Int_t N, Double_t *t, Double_t *w, Double_t thres, Double_t trig, Double_t gate, Double_t width, bool debug)    // returns end point
/// {
///     // baseline = 0 for now
/// 
///     Int_t thres_i = 0;
///     while (w[thres_i] < thres && thres_i < N) ++thres_i;
///     if (debug) cout<< "thres_i = " << thres_i << " t[thres_i] = " << t[thres_i] <<endl;
///     if (thres_i == N) return t[thres_i];            // no signal
/// 
///     Int_t maximum_i = thres_i;
///     for (int i=thres_i+1; i<N; ++i) {
///         if (w[i] > w[maximum_i]) maximum_i = i;
///         if (w[i] < 0.9*w[maximum_i]) break;              // stop at the first maximum
///     }
///     if (debug) cout<< "maximum_i = " << maximum_i << " t[maximum_i] = " << t[maximum_i] <<endl;
/// 
///     Int_t halfmax_i1 = maximum_i - 1;
///     if (halfmax_i1 < 0) halfmax_i1 = 0;
///     while (halfmax_i1 > 0 && w[halfmax_i1] > 0.5*w[maximum_i]) --halfmax_i1;
/// 
///     Int_t halfmax_i2 = maximum_i + 1;
///     if (halfmax_i2 >= N-1) halfmax_i2 = N-1;
///     while (halfmax_i2 < N-2 && w[halfmax_i2] > 0.5*w[maximum_i]) ++halfmax_i2;
/// 
///     Double_t fwhm = t[halfmax_i2] - t[halfmax_i1];
/// 
///     if (debug) cout<< "t[halfmax_i1] = " << t[halfmax_i1] << " t[halfmax_i2] = " << t[halfmax_i2] << " fwhm = " << fwhm <<endl;
/// 
///     // find the signal region
/// 
///     Int_t signal_i1 = maximum_i - 3*(maximum_i - halfmax_i1);
///     if (signal_i1 < 0) signal_i1 = 0;
///     Int_t signal_i2 = maximum_i + 3*(halfmax_i2 - maximum_i);
///     if (signal_i2 > N-1) signal_i1 = N-1;
///     if (debug) cout<< "t[signal_i1] = " << t[signal_i1] << " t[signal_i2] = " << t[signal_i2] <<endl;
/// 
///     if (fwhm < width) {
///         // perform a recursive call
///     }
/// }

void channel_fit(Int_t N, Double_t* t, Double_t* w
                 , Double_t thres
                 , Double_t trig
                 , Double_t cutoff_GHz
                 , Double_t yfit1=0.30
                 , Double_t yfit2=0.70
                 , Int_t chan=0
                 , Int_t event=-1
                 , bool debug=false)
{
    if (debug) cout<< "channel_fit event = " << event << " channel = " << chan <<endl;

    if (chan < 1) {
        cout<< "invalid channel number" <<endl;
        return;
    }

    if (event < 0) {
        cout<< "invalid event number" <<endl;
        return;
    }

    TGraph* g = new TGraph(N,t,w);
    g->SetNameTitle(Form("g_evt_%d_ch_%d",event,chan), Form("g_evt_%d_ch_%d",event,chan));
    g->SetMarkerStyle(7);
    g->SetMarkerColor(602);
    g->SetLineColor(602);

    if (debug) {
        new TCanvas;
        g->Draw("ap");
    }

    // pmax in the raw data
    Int_t maximum_raw_i = 0;
    for (int i=0; i<N; ++i) {
        if (w[i] > w[maximum_raw_i]) maximum_raw_i = i;
    }
    TreeRes::pmax_raw[chan-1] = w[maximum_raw_i];
    if (debug) cout<< "maximum_raw_i = " << maximum_raw_i << " t[maximum_raw_i] = " << t[maximum_raw_i] << " w[maximum_raw_i] = " << w[maximum_raw_i] <<endl;

    // apply low pass filter
    Double_t re0 = 0;
    if (cutoff_GHz > 0) re0 = fft_filter(N, t, w, cutoff_GHz);
    if (cutoff_GHz > 0 && debug) cout<< "re0 = " << re0 <<endl;

    Double_t ey[20000];
    for (int i=0; i<N; ++i) ey[i] = 0.001;

    //TGraph* gout = new TGraph(N,t,w);
    TGraphErrors* gout = new TGraphErrors(N,t,w,0,ey);
    gout->SetNameTitle(Form("gout_evt_%d_ch_%d",event,chan), Form("gout_evt_%d_ch_%d",event,chan));
    gout->SetMarkerStyle(7);
    gout->SetMarkerColor(46);
    gout->SetLineColor(46);

    if (debug) {
        gout->Draw("p");
    }

    // find the trigger time point: the point to start

    Int_t trig_i;
    for (trig_i=0; trig_i<N; ++trig_i) if (t[trig_i] > trig) break;

    Int_t thres_i;
    for (thres_i=trig_i; thres_i<N; ++thres_i) if (w[thres_i] > thres) break;
    if (thres_i == N) thres_i = N-1;
    if (debug) cout<< "thres_i = " << thres_i << " t[thres_i] = " << t[thres_i] <<endl;
    if (thres_i == N-1) return;                          // no signal above the threshold

    Int_t maximum_i = thres_i;
    for (int i=thres_i+1; i<N; ++i) {
        if (w[i] > w[maximum_i]) maximum_i = i;
        //-- if (w[i] < 0.90*w[maximum_i]) break;              // stop at the first maximum
        if (w[i] < 0.75*w[maximum_i]) break;              // stop at the first maximum
    }
    if (debug) cout<< "maximum_i = " << maximum_i << " t[maximum_i] = " << t[maximum_i] << " w[maximum_i] = " << w[maximum_i] <<endl;

    Int_t halfmax_i1;
    for (halfmax_i1=maximum_i; halfmax_i1>=trig_i; --halfmax_i1) if (w[halfmax_i1] < 0.5*w[maximum_i]) break;

    Int_t halfmax_i2;
    for (halfmax_i2=maximum_i; halfmax_i2<N; ++halfmax_i2) if (w[halfmax_i2] < 0.5*w[maximum_i]) break;
    if (halfmax_i2 == N) halfmax_i2 = N-1;

    Double_t fwhm = t[halfmax_i2] - t[halfmax_i1];

    Int_t quartermax_i1;
    for (quartermax_i1=maximum_i; quartermax_i1>=trig_i; --quartermax_i1) if (w[quartermax_i1] < 0.25*w[maximum_i]) break;

    Int_t quartermax_i2;
    for (quartermax_i2=maximum_i; quartermax_i2<N; ++quartermax_i2) if (w[quartermax_i2] < 0.25*w[maximum_i]) break;
    if (quartermax_i2 == N) quartermax_i2 = N-1;

    Double_t fwqm = t[quartermax_i2] - t[quartermax_i1];

    if (debug) cout<< "t[halfmax_i1] = " << t[halfmax_i1] << " t[halfmax_i2] = " << t[halfmax_i2] << " fwhm = " << fwhm <<endl;

    // signal region

    Int_t signal_i1 = maximum_i - 3*(maximum_i - halfmax_i1);
    if (signal_i1 < 0) signal_i1 = 0;
    Int_t signal_i2 = maximum_i + 3*(halfmax_i2 - maximum_i);
    if (signal_i2 > N-1) signal_i2 = N-1;
    Double_t charge = 0;
    for (int i=signal_i1; i<=signal_i2; ++i) charge += w[i];
    if (debug) cout<< "t[signal_i1] = " << t[signal_i1] << " t[signal_i2] = " << t[signal_i2] << " charge = " << charge <<endl;

    // fill constant fraction thresholds
    Double_t pmax = w[maximum_i];
    Int_t cf_i = maximum_i;
    Double_t frac;
    Double_t wfrac;
    Double_t dwdt;
    Double_t dwdtmax = 0;
    Double_t dwdtfrac = 0;
    Double_t cf;
    // go from 60% to 5% of the maximum value
    frac = 0.60;
    wfrac = frac*pmax;
    while (cf_i >= 0 && w[cf_i] > wfrac) --cf_i;
    dwdt = (w[cf_i+1] - w[cf_i]) / (t[cf_i+1] - t[cf_i]);
    cf = TMath::Abs(dwdt) > 1e-7? t[cf_i] + (wfrac - w[cf_i])/dwdt: t[cf_i];
    TreeRes::cf60[chan-1] = cf;
    TreeRes::dwdt60[chan-1] = dwdt;
    if (dwdt > dwdtmax) {dwdtmax = dwdt; dwdtfrac = frac;}

    frac = 0.55;
    wfrac = frac*pmax;
    while (cf_i >= 0 && w[cf_i] > wfrac) --cf_i;
    dwdt = (w[cf_i+1] - w[cf_i]) / (t[cf_i+1] - t[cf_i]);
    cf = TMath::Abs(dwdt) > 1e-7? t[cf_i] + (wfrac - w[cf_i])/dwdt: t[cf_i];
    TreeRes::cf55[chan-1] = cf;
    TreeRes::dwdt55[chan-1] = dwdt;
    if (dwdt > dwdtmax) {dwdtmax = dwdt; dwdtfrac = frac;}

    frac = 0.50;
    wfrac = frac*pmax;
    while (cf_i >= 0 && w[cf_i] > wfrac) --cf_i;
    dwdt = (w[cf_i+1] - w[cf_i]) / (t[cf_i+1] - t[cf_i]);
    cf = TMath::Abs(dwdt) > 1e-7? t[cf_i] + (wfrac - w[cf_i])/dwdt: t[cf_i];
    TreeRes::cf50[chan-1] = cf;
    TreeRes::dwdt50[chan-1] = dwdt;
    if (dwdt > dwdtmax) {dwdtmax = dwdt; dwdtfrac = frac;}

    frac = 0.45;
    wfrac = frac*pmax;
    while (cf_i >= 0 && w[cf_i] > wfrac) --cf_i;
    dwdt = (w[cf_i+1] - w[cf_i]) / (t[cf_i+1] - t[cf_i]);
    cf = TMath::Abs(dwdt) > 1e-7? t[cf_i] + (wfrac - w[cf_i])/dwdt: t[cf_i];
    TreeRes::cf45[chan-1] = cf;
    TreeRes::dwdt45[chan-1] = dwdt;
    if (dwdt > dwdtmax) {dwdtmax = dwdt; dwdtfrac = frac;}

    frac = 0.40;
    wfrac = frac*pmax;
    while (cf_i >= 0 && w[cf_i] > wfrac) --cf_i;
    dwdt = (w[cf_i+1] - w[cf_i]) / (t[cf_i+1] - t[cf_i]);
    cf = TMath::Abs(dwdt) > 1e-7? t[cf_i] + (wfrac - w[cf_i])/dwdt: t[cf_i];
    TreeRes::cf40[chan-1] = cf;
    TreeRes::dwdt40[chan-1] = dwdt;
    if (dwdt > dwdtmax) {dwdtmax = dwdt; dwdtfrac = frac;}

    frac = 0.35;
    wfrac = frac*pmax;
    while (cf_i >= 0 && w[cf_i] > wfrac) --cf_i;
    dwdt = (w[cf_i+1] - w[cf_i]) / (t[cf_i+1] - t[cf_i]);
    cf = TMath::Abs(dwdt) > 1e-7? t[cf_i] + (wfrac - w[cf_i])/dwdt: t[cf_i];
    TreeRes::cf35[chan-1] = cf;
    TreeRes::dwdt35[chan-1] = dwdt;
    if (dwdt > dwdtmax) {dwdtmax = dwdt; dwdtfrac = frac;}

    frac = 0.30;
    wfrac = frac*pmax;
    while (cf_i >= 0 && w[cf_i] > wfrac) --cf_i;
    dwdt = (w[cf_i+1] - w[cf_i]) / (t[cf_i+1] - t[cf_i]);
    cf = TMath::Abs(dwdt) > 1e-7? t[cf_i] + (wfrac - w[cf_i])/dwdt: t[cf_i];
    TreeRes::cf30[chan-1] = cf;
    TreeRes::dwdt30[chan-1] = dwdt;
    if (dwdt > dwdtmax) {dwdtmax = dwdt; dwdtfrac = frac;}

    frac = 0.25;
    wfrac = frac*pmax;
    while (cf_i >= 0 && w[cf_i] > wfrac) --cf_i;
    dwdt = (w[cf_i+1] - w[cf_i]) / (t[cf_i+1] - t[cf_i]);
    cf = TMath::Abs(dwdt) > 1e-7? t[cf_i] + (wfrac - w[cf_i])/dwdt: t[cf_i];
    TreeRes::cf25[chan-1] = cf;
    TreeRes::dwdt25[chan-1] = dwdt;
    if (dwdt > dwdtmax) {dwdtmax = dwdt; dwdtfrac = frac;}

    frac = 0.20;
    wfrac = frac*pmax;
    while (cf_i >= 0 && w[cf_i] > wfrac) --cf_i;
    dwdt = (w[cf_i+1] - w[cf_i]) / (t[cf_i+1] - t[cf_i]);
    cf = TMath::Abs(dwdt) > 1e-7? t[cf_i] + (wfrac - w[cf_i])/dwdt: t[cf_i];
    TreeRes::cf20[chan-1] = cf;
    TreeRes::dwdt20[chan-1] = dwdt;
    if (dwdt > dwdtmax) {dwdtmax = dwdt; dwdtfrac = frac;}

    frac = 0.15;
    wfrac = frac*pmax;
    while (cf_i >= 0 && w[cf_i] > wfrac) --cf_i;
    dwdt = (w[cf_i+1] - w[cf_i]) / (t[cf_i+1] - t[cf_i]);
    cf = TMath::Abs(dwdt) > 1e-7? t[cf_i] + (wfrac - w[cf_i])/dwdt: t[cf_i];
    TreeRes::cf15[chan-1] = cf;
    TreeRes::dwdt15[chan-1] = dwdt;
    if (dwdt > dwdtmax) {dwdtmax = dwdt; dwdtfrac = frac;}

    frac = 0.10;
    wfrac = frac*pmax;
    while (cf_i >= 0 && w[cf_i] > wfrac) --cf_i;
    dwdt = (w[cf_i+1] - w[cf_i]) / (t[cf_i+1] - t[cf_i]);
    cf = TMath::Abs(dwdt) > 1e-7? t[cf_i] + (wfrac - w[cf_i])/dwdt: t[cf_i];
    TreeRes::cf10[chan-1] = cf;
    TreeRes::dwdt10[chan-1] = dwdt;
    if (dwdt > dwdtmax) {dwdtmax = dwdt; dwdtfrac = frac;}

    TreeRes::dwdtfrac[chan-1] = dwdtfrac;

    if (debug) cout<< "dwdtfrac = " << dwdtfrac << " dwdt: 0.10: " << TreeRes::dwdt10[chan-1] << " 0.15: " << TreeRes::dwdt15[chan-1] << " 0.20: " << TreeRes::dwdt20[chan-1]  << " 0.25: " << TreeRes::dwdt25[chan-1] << " 0.30: " << TreeRes::dwdt30[chan-1] << " 0.35: " << TreeRes::dwdt35[chan-1] << " 0.40: " << TreeRes::dwdt40[chan-1] << " 0.45: " << TreeRes::dwdt45[chan-1] << " 0.50: " << TreeRes::dwdt50[chan-1] << " 0.55: " << TreeRes::dwdt55[chan-1] << " 0.60: " << TreeRes::dwdt60[chan-1] <<endl;

    Double_t baseline = 0;
    if (signal_i1 > 7) {
        // Double_t W = 0;      // weight sum
        // Double_t mean = 0;
        // for (int i=0; i<signal_i; ++i) {
        //    baseline += w[i];
        //    Double_t w = 1.;                    // equal weights
        //    mean = mean + w*(w[i] - mean)/(W+w);
        //    W = W + w;
        // }
        // baseline /= signal_i;
        // TreeRes::bkg[chan-1] = mean;
        // TreeRes::ebkg[chan-1] = mean/TMath::Sqrt(W);
        TreeRes::bkg_raw[chan-1] = TMath::Mean(signal_i1, w);
        TreeRes::ebkg_raw[chan-1] = TMath::RMS(signal_i1, w);
        TreeRes::bkg[chan-1] = TMath::Mean(signal_i1, w);
        TreeRes::ebkg[chan-1] = TMath::RMS(signal_i1, w);
        baseline = TreeRes::bkg[chan-1];
    }
    if (debug) cout<< "baseline = " << baseline <<endl;
    baseline = 0;

    for (int i=0; i<N; ++i) {
        w[i] -= baseline;
        g->SetPoint(i, g->GetX()[i], g->GetY()[i]-baseline);
        gout->SetPoint(i, gout->GetX()[i], gout->GetY()[i]-baseline);
    }

    /// Double_t ey[20000];
    /// //-----------------------------------------------------for (int i=0; i<N; ++i) ey[i] = TreeRes::ebkg[chan-1];
    /// for (int i=0; i<N; ++i) ey[i] = 0.001;

    /// //TGraph* gout = new TGraph(N,t,w);
    /// TGraphErrors* gout = new TGraphErrors(N,t,w,0,ey);
    /// gout->SetNameTitle(Form("gout_evt_%d_ch_%d",event,chan), Form("gout_evt_%d_ch_%d",event,chan));
    /// gout->SetMarkerStyle(7);
    /// gout->SetMarkerColor(46);
    /// gout->SetLineColor(46);

    /// if (debug) {
    ///     gout->Draw("p");
    /// }

    // Int_t i88 = maximum_i;
    // while (w[i88] > 0.88*w[maximum_i] && i88 > 0) --i88;

    // Int_t i12 = i88;
    // while (w[i12] > 0.12*w[maximum_i] && i12 > 0) --i12;

    // if (debug) cout<< "i88 = " << i88 << " t[i88] = " << t[i88] << " i12 = " << i12 << " t[i12] = " << t[i12] <<endl;

    Int_t i80 = maximum_i;
    while (w[i80] > yfit2*w[maximum_i] && i80 > 0) --i80;

    Int_t i20 = i80;
    while (w[i20] > yfit1*w[maximum_i] && i20 > 0) --i20;

    if (debug) cout<< "i80 = " << i80 << " t[i80] = " << t[i80] << " i20 = " << i20 << " t[i20] = " << t[i20] <<endl;

    Double_t fit_x1 = 0;
    Double_t fit_x2 = 0;
    Double_t eps = 1e-7;
    // fit_x1 = t[i12] - eps;
    // fit_x2 = t[i88] + eps;
    fit_x1 = t[i20] - eps;
    fit_x2 = t[i80] + eps;
    //-------------------------------------------------------------------------fit_x1 -= 6.;  // ns
    if (debug) cout<< "fit_x1 = " << fit_x1 << " fit_x2 = " << fit_x2 <<endl;

    /// // Double_t tau = 0.5*(t[maximum_i] - t[i12]);
    /// // Double_t x0 = t[i12];
    /// // Double_t A = w[maximum_i]*(tau/2.)*TMath::Exp(2);
    /// Double_t tau = 0.5*(t[maximum_i] - t[i20]);
    /// Double_t x0 = t[i20];
    /// Double_t A = w[maximum_i]*(tau/2.)*TMath::Exp(2);

    /// //-- TF1* f = FunRC::RC2(fit_x1,fit_x2, A,x0,tau);
    /// Double_t sigma = 0.5;
    /// // TF1* f = FunRC::RC1gb(fit_x1,fit_x2, A,x0,tau,sigma);        // NB: RC1gb
    /// TF1* f = FunRC::RC2gb(fit_x1,fit_x2, A,x0,tau,sigma);        // NB: RC2gb

    //
    // NB: ROOT plots a result of the first fit.
    //

    //
    //  Fit gaus fit first (will be shown on the plot) then straight line
    //
    // Double_t dmax = 0.300;  // ps
    // if (debug) gout->Fit("gaus","","",t[maximum_i]-dmax,t[maximum_i]+dmax);
    // else gout->Fit("gaus","q","goff",t[maximum_i]-dmax,t[maximum_i]+dmax);
    // TreeRes::pmaxx[chan-1] = gout->GetFunction("gaus")->GetParameter(1);
    //
    // if (debug) gout->Fit("pol1","+","",fit_x1,fit_x2);
    // else gout->Fit("pol1","+q","goff",fit_x1,fit_x2);

    //
    //  Fit straight line fit first (will be shown on the plot) then gaus
    //
    if (debug) gout->Fit("pol1","","",fit_x1,fit_x2);
    else gout->Fit("pol1","q","goff",fit_x1,fit_x2);

    Double_t dmax = 0.300;  // ps
    if (debug) gout->Fit("gaus","+","",t[maximum_i]-dmax,t[maximum_i]+dmax);
    else gout->Fit("gaus","q+","goff",t[maximum_i]-dmax,t[maximum_i]+dmax);
    gout->GetFunction("gaus")->SetLineColor(1);
    gout->GetFunction("gaus")->SetLineWidth(1);
    TreeRes::pmaxx[chan-1] = gout->GetFunction("gaus")->GetParameter(1);

    Double_t delta = 0.1;
    TF1* fun = gout->GetFunction("pol1");
    if (!fun) {
        cout<< "--> No fit function found!" <<endl;
        //-- return;
    }

    if (fun) {
        //-- fun->SetRange(fit_x1-delta, fit_x2+delta);      // extend the fit function range
        TreeRes::npar[chan-1] = fun->GetNpar();
        if (fun->GetNpar() > 0) TreeRes::par0[chan-1] = fun->GetParameter(0);
        if (fun->GetNpar() > 1) TreeRes::par1[chan-1] = fun->GetParameter(1);
        if (fun->GetNpar() > 2) TreeRes::par2[chan-1] = fun->GetParameter(2);
        if (fun->GetNpar() > 3) TreeRes::par3[chan-1] = fun->GetParameter(3);
        if (fun->GetNpar() > 4) TreeRes::par4[chan-1] = fun->GetParameter(4);
        if (fun->GetNpar() > 5) TreeRes::par5[chan-1] = fun->GetParameter(5);
    }

    // get baseline from the fit
    //-- Double_t bkg = fun->GetParameter("bkg");
    //-- Double_t bkg = fun->GetParameter(4);
    Double_t bkg = 0;

    TreeRes::fwhm[chan-1] = fwhm;
    TreeRes::fwqm[chan-1] = fwqm;
    TreeRes::q[chan-1] = charge;
    TreeRes::pmax[chan-1] = w[maximum_i];
    //TreeRes::pmaxx[chan-1] = t[maximum_i];
    if (fun) {
        TreeRes::ym[chan-1] = bkg + 0.5*(w[maximum_i] - bkg);
        TreeRes::xm[chan-1] = fun->GetX(TreeRes::ym[chan-1], t[i20], t[i80]);
        //TreeRes::ts = fun->GetX(0.5*w[maximum_i], t[i20], t[i80]);
        //-- Double_t derivative = (2./(TreeRes::xm[chan-1]-TreeRes::par[chan-1][1]) - 1/TreeRes::par[chan-1][2]) * TreeRes::ym[chan-1];
        Double_t derivative = fun->Derivative(TreeRes::xm[chan-1]);
        TreeRes::ts[chan-1] = TreeRes::xm[chan-1] - (TreeRes::ym[chan-1] - bkg)/derivative;
        if (debug) cout<< "TreeRes::ts[chan-1] = " << TreeRes::ts[chan-1] << " derivative = " << derivative <<endl;

        TreeRes::tx[chan-1] = fit_x1;
        for (int i=i20; i<i80; ++i) {
            derivative = fun->Derivative(t[i]);
            Double_t tx = t[i] - (fun->Eval(t[i]) - bkg)/derivative;
            if (tx > TreeRes::tx[chan-1]) TreeRes::tx[chan-1] = tx;
        }

        TreeRes::x5mV[chan-1] = fun->GetX(5+bkg, fit_x1, t[i80]);
        TreeRes::x10mV[chan-1] = fun->GetX(10+bkg, fit_x1, t[i80]);
    }
}

/// void channel_fit(Int_t N, Double_t* t, Double_t* w
///                  , Double_t thres
///                  , Double_t trig
///                  , Double_t cutoff_GHz
///                  , Double_t yfit1=0.30
///                  , Double_t yfit2=0.70
///                  , Int_t chan
///                  , Int_t event
///                  , bool debug=false)

void channel1_fit(Int_t N, Double_t* t, Double_t* w, Double_t thres, Double_t trig, Double_t cutoff_GHz, Double_t yfit1, Double_t yfit2, Int_t chan=1, Int_t event=-1, bool debug=false)
{
    if (debug) cout<< "channel1_fit event = " << event <<endl;
    if (debug) cout<< "thres = " << thres << " trig = " << trig << " cutoff_GHz = " << cutoff_GHz << " yfit1 = " << yfit1 << " yfit2 = " << yfit2 << " chan = " << chan << " event = " << event << " debug = " << debug <<endl;
    channel_fit(N,t,w,thres,trig,cutoff_GHz,yfit1,yfit2,chan,event,debug);
}

void channel2_fit(Int_t N, Double_t* t, Double_t* w, Double_t thres, Double_t trig, Double_t cutoff_GHz, Double_t yfit1, Double_t yfit2, Int_t chan=2, Int_t event=-1, bool debug=false)
{
    if (debug) cout<< "channel2_fit event = " << event <<endl;
    if (debug) cout<< "thres = " << thres << " trig = " << trig << " cutoff_GHz = " << cutoff_GHz << " yfit1 = " << yfit1 << " yfit2 = " << yfit2 << " chan = " << chan << " event = " << event << " debug = " << debug <<endl;
    channel_fit(N,t,w,thres,trig,cutoff_GHz,yfit1,yfit2,chan,event,debug);
}

void channel3_fit(Int_t N, Double_t* t, Double_t* w, Double_t thres, Double_t trig, Double_t cutoff_GHz, Double_t yfit1, Double_t yfit2, Int_t chan=3, Int_t event=-1, bool debug=false)
{
    if (debug) cout<< "channel3_fit event = " << event <<endl;
    if (debug) cout<< "thres = " << thres << " trig = " << trig << " cutoff_GHz = " << cutoff_GHz << " yfit1 = " << yfit1 << " yfit2 = " << yfit2 << " chan = " << chan << " event = " << event << " debug = " << debug <<endl;
    channel_fit(N,t,w,thres,trig,cutoff_GHz,yfit1,yfit2,chan,event,debug);
}

void channel4_fit(Int_t N, Double_t* t, Double_t* w, Double_t thres, Double_t trig, Double_t cutoff_GHz, Double_t yfit1, Double_t yfit2, Int_t chan=4, Int_t event=-1, bool debug=false)
{
    if (debug) cout<< "channel4_fit event = " << event <<endl;
    if (debug) cout<< "thres = " << thres << " trig = " << trig << " cutoff_GHz = " << cutoff_GHz << " yfit1 = " << yfit1 << " yfit2 = " << yfit2 << " chan = " << chan << " event = " << event << " debug = " << debug <<endl;
    channel_fit(N,t,w,thres,trig,cutoff_GHz,yfit1,yfit2,chan,event,debug);
}

/// void channel2_fit(Int_t N, Double_t* t, Double_t* w, Double_t thres, Double_t trig, Double_t cutoff_GHz, Double_t yfit1, Double_t yfit2, Int_t chan=2, Int_t event=-1, bool debug=false)
/// {
///     if (debug) cout<< "channel2_fitRC event = " << event <<endl;
///     if (debug) cout<< "thres = " << thres << " trig = " << trig << " cutoff_GHz = " << cutoff_GHz << " yfit1 = " << yfit1 << " yfit2 = " << yfit2 << " chan = " << chan << " event = " << event << " debug = " << debug <<endl;
/// 
///     //--     if (debug) cout<< "channel_fit event = " << event << " channel = " << chan <<endl;
/// 
///     //-- channel_fit(N,t,w,thres,trig,cutoff_GHz,yfit1,yfit2,chan,event,debug);
/// 
///     if (chan < 1) {
///         cout<< "invalid channel number" <<endl;
///         return;
///     }
/// 
///     if (event < 0) {
///         cout<< "invalid event number" <<endl;
///         return;
///     }
/// 
///     TGraph* g = new TGraph(N,t,w);
///     g->SetNameTitle(Form("g_evt_%d_ch_%d",event,chan), Form("g_evt_%d_ch_%d",event,chan));
///     g->SetMarkerStyle(7);
///     g->SetMarkerColor(602);
///     g->SetLineColor(602);
/// 
///     if (debug) {
///         new TCanvas;
///         g->Draw("ap");
///     }
/// 
///     // pmax in the raw data
///     Int_t maximum_raw_i = 0;
///     for (int i=0; i<N; ++i) {
///         if (w[i] > w[maximum_raw_i]) maximum_raw_i = i;
///     }
///     TreeRes::pmax_raw[chan-1] = w[maximum_raw_i];
///     if (debug) cout<< "maximum_raw_i = " << maximum_raw_i << " t[maximum_raw_i] = " << t[maximum_raw_i] << " w[maximum_raw_i] = " << w[maximum_raw_i] <<endl;
/// 
///     // apply low pass filter
///     Double_t re0 = 0;
///     if (cutoff_GHz > 0) re0 = fft_filter(N, t, w, cutoff_GHz);
///     if (cutoff_GHz > 0 && debug) cout<< "re0 = " << re0 <<endl;
/// 
///     // find the trigger time point: the point to start
/// 
///     Int_t trig_i;
///     for (trig_i=0; trig_i<N; ++trig_i) if (t[trig_i] > trig) break;
/// 
///     Int_t thres_i;
///     for (thres_i=trig_i; thres_i<N; ++thres_i) if (w[thres_i] > thres) break;
///     if (thres_i == N) thres_i = N-1;
///     if (debug) cout<< "thres_i = " << thres_i << " t[thres_i] = " << t[thres_i] <<endl;
///     if (thres_i == N-1) return;                          // no signal above the threshold
/// 
///     Int_t maximum_i = thres_i;
///     for (int i=thres_i+1; i<N; ++i) {
///         if (w[i] > w[maximum_i]) maximum_i = i;
///         if (w[i] < 0.9*w[maximum_i]) break;              // stop at the first maximum
///     }
///     if (debug) cout<< "maximum_i = " << maximum_i << " t[maximum_i] = " << t[maximum_i] <<endl;
/// 
///     Int_t halfmax_i1;
///     for (halfmax_i1=maximum_i; halfmax_i1>=trig_i; --halfmax_i1) if (w[halfmax_i1] < 0.5*w[maximum_i]) break;
/// 
///     Int_t halfmax_i2;
///     for (halfmax_i2=maximum_i; halfmax_i2<N; ++halfmax_i2) if (w[halfmax_i2] < 0.5*w[maximum_i]) break;
///     if (halfmax_i2 == N) halfmax_i2 = N-1;
/// 
///     Double_t fwhm = t[halfmax_i2] - t[halfmax_i1];
/// 
///     if (debug) cout<< "t[halfmax_i1] = " << t[halfmax_i1] << " t[halfmax_i2] = " << t[halfmax_i2] << " fwhm = " << fwhm <<endl;
/// 
///     // signal region
/// 
///     Int_t signal_i1 = maximum_i - 3*(maximum_i - halfmax_i1);
///     if (signal_i1 < 0) signal_i1 = 0;
///     Int_t signal_i2 = maximum_i + 3*(halfmax_i2 - maximum_i);
///     if (signal_i2 > N-1) signal_i2 = N-1;
///     Double_t charge = 0;
///     for (int i=signal_i1; i<=signal_i2; ++i) charge += w[i];
///     if (debug) cout<< "t[signal_i1] = " << t[signal_i1] << " t[signal_i2] = " << t[signal_i2] << " charge = " << charge <<endl;
/// 
///     Double_t baseline = 0;
///     if (signal_i1 > 7) {
///         // Double_t W = 0;      // weight sum
///         // Double_t mean = 0;
///         // for (int i=0; i<signal_i; ++i) {
///         //    baseline += w[i];
///         //    Double_t w = 1.;                    // equal weights
///         //    mean = mean + w*(w[i] - mean)/(W+w);
///         //    W = W + w;
///         // }
///         // baseline /= signal_i;
///         // TreeRes::bkg[chan-1] = mean;
///         // TreeRes::ebkg[chan-1] = mean/TMath::Sqrt(W);
///         TreeRes::bkg_raw[chan-1] = TMath::Mean(signal_i1, w);
///         TreeRes::ebkg_raw[chan-1] = TMath::RMS(signal_i1, w);
///         TreeRes::bkg[chan-1] = TMath::Mean(signal_i1, w);
///         TreeRes::ebkg[chan-1] = TMath::RMS(signal_i1, w);
///         baseline = TreeRes::bkg[chan-1];
///     }
/// 
///     for (int i=0; i<N; ++i) {
///         w[i] -= baseline;
///         g->SetPoint(i, g->GetX()[i], g->GetY()[i]-baseline);
///     }
/// 
///     Double_t ey[20000];
///     //-----------------------------------------------------for (int i=0; i<N; ++i) ey[i] = TreeRes::ebkg[chan-1];
///     for (int i=0; i<N; ++i) ey[i] = 0.001;
/// 
///     //TGraph* gout = new TGraph(N,t,w);
///     TGraphErrors* gout = new TGraphErrors(N,t,w,0,ey);
///     gout->SetNameTitle(Form("gout_evt_%d_ch_%d",event,chan), Form("gout_evt_%d_ch_%d",event,chan));
///     gout->SetMarkerStyle(7);
///     gout->SetMarkerColor(46);
///     gout->SetLineColor(46);
/// 
///     if (debug) {
///         gout->Draw("p");
///     }
/// 
///     // Int_t i88 = maximum_i;
///     // while (w[i88] > 0.88*w[maximum_i] && i88 > 0) --i88;
/// 
///     // Int_t i12 = i88;
///     // while (w[i12] > 0.12*w[maximum_i] && i12 > 0) --i12;
/// 
///     // if (debug) cout<< "i88 = " << i88 << " t[i88] = " << t[i88] << " i12 = " << i12 << " t[i12] = " << t[i12] <<endl;
/// 
///     Int_t i80 = maximum_i;
///     while (w[i80] > yfit2*w[maximum_i] && i80 > 0) --i80;
/// 
///     Int_t i20 = i80;
///     while (w[i20] > yfit1*w[maximum_i] && i20 > 0) --i20;
/// 
///     if (debug) cout<< "i80 = " << i80 << " t[i80] = " << t[i80] << " i20 = " << i20 << " t[i20] = " << t[i20] <<endl;
/// 
///     Double_t fit_x1 = 0;
///     Double_t fit_x2 = 0;
///     Double_t eps = 1e-7;
///     // fit_x1 = t[i12] - eps;
///     // fit_x2 = t[i88] + eps;
/// 
///     //-- fit_x1 = t[i20] - eps;
///     //-- fit_x2 = t[i80] + eps;
/// 
///     fit_x1 = t[signal_i1] - eps;
///     fit_x2 = t[halfmax_i2] + eps;
///     //-------------------------------------------------------------------------fit_x1 -= 6.;  // ns
///     if (debug) cout<< "fit_x1 = " << fit_x1 << " fit_x2 = " << fit_x2 <<endl;
/// 
///     // Double_t tau = 0.5*(t[maximum_i] - t[i12]);
///     // Double_t x0 = t[i12];
///     // Double_t A = w[maximum_i]*(tau/2.)*TMath::Exp(2);
///     Double_t tau = 0.5*(t[maximum_i] - t[i20]);
///     Double_t x0 = t[i20];
///     Double_t A = w[maximum_i]*(tau/2.)*TMath::Exp(2);
/// 
///     tau = 0.2538;
///     Double_t sigma = 0.1256;
/// 
///     Double_t AvsQ_par0 = 0.00651518;
///     Double_t AvsQ_par1 = 0.028516;
///     A = AvsQ_par0 + AvsQ_par1*charge;
/// 
///     //-- TF1* f = FunRC::RC2(fit_x1,fit_x2, A,x0,tau);
///     //Double_t sigma = 0.5;
///     TF1* f = FunRC::RC1g(fit_x1,fit_x2, A,x0,tau,sigma);        // NB: RC1gb
///     // TF1* f = FunRC::RC2gb(fit_x1,fit_x2, A,x0,tau,sigma);        // NB: RC2gb
///     // TF1* f = FunRC::RC2g(fit_x1,fit_x2, A,x0,tau,sigma);        // NB: RC2g
/// 
///     // TF1* f = new TF1("pol1", "pol1", fit_x1, fit_x2);          // straight line
/// 
///     // f->FixParameter(0, A);
///     // f->FixParameter(2, tau);
///     // f->FixParameter(3, sigma);
/// 
///     if (debug) gout->Fit(f,"","",fit_x1,fit_x2);
///     else gout->Fit(f,"q","goff",fit_x1,fit_x2);
/// 
///     Double_t delta = 0.1;
///     //--line-- TF1* fun = gout->GetFunction("pol1");
///     TF1* fun = gout->GetFunction(f->GetName());
///     if (!fun) {
///         cout<< "--> No fit function found!" <<endl;
///         //-- return;
///     }
/// 
///     if (fun) {
///         //-- fun->SetRange(fit_x1-delta, fit_x2+delta);      // extend the fit function range
///         TreeRes::npar[chan-1] = fun->GetNpar();
///         if (fun->GetNpar() > 0) TreeRes::par0[chan-1] = fun->GetParameter(0);
///         if (fun->GetNpar() > 1) TreeRes::par1[chan-1] = fun->GetParameter(1);
///         if (fun->GetNpar() > 2) TreeRes::par2[chan-1] = fun->GetParameter(2);
///         if (fun->GetNpar() > 3) TreeRes::par3[chan-1] = fun->GetParameter(3);
///         if (fun->GetNpar() > 4) TreeRes::par4[chan-1] = fun->GetParameter(4);
///         if (fun->GetNpar() > 5) TreeRes::par5[chan-1] = fun->GetParameter(5);
///     }
/// 
///     // get baseline from the fit
///     //-- Double_t bkg = fun->GetParameter("bkg");
///     //-- Double_t bkg = fun->GetParameter(4);
///     Double_t bkg = 0;
/// 
///     TreeRes::fwhm[chan-1] = fwhm;
///     TreeRes::q[chan-1] = charge;
///     TreeRes::pmax[chan-1] = w[maximum_i];
///     if (fun) {
///         TreeRes::ym[chan-1] = bkg + 0.5*(w[maximum_i] - bkg);
///         TreeRes::xm[chan-1] = fun->GetX(TreeRes::ym[chan-1], t[i20], t[i80]);
///         //TreeRes::ts = fun->GetX(0.5*w[maximum_i], t[i20], t[i80]);
///         //-- Double_t derivative = (2./(TreeRes::xm[chan-1]-TreeRes::par[chan-1][1]) - 1/TreeRes::par[chan-1][2]) * TreeRes::ym[chan-1];
///         Double_t derivative = fun->Derivative(TreeRes::xm[chan-1]);
///         TreeRes::ts[chan-1] = TreeRes::xm[chan-1] - (TreeRes::ym[chan-1] - bkg)/derivative;
///         if (debug) cout<< "TreeRes::ts[chan-1] = " << TreeRes::ts[chan-1] << " derivative = " << derivative <<endl;
/// 
///         TreeRes::tx[chan-1] = fit_x1;
///         for (int i=i20; i<i80; ++i) {
///             derivative = fun->Derivative(t[i]);
///             Double_t tx = t[i] - (fun->Eval(t[i]) - bkg)/derivative;
///             if (tx > TreeRes::tx[chan-1]) TreeRes::tx[chan-1] = tx;
///         }
/// 
///         TreeRes::x5mV[chan-1] = fun->GetX(5+bkg, fit_x1, t[i80]);
///         TreeRes::x10mV[chan-1] = fun->GetX(10+bkg, fit_x1, t[i80]);
///     }
/// }

void channel2_run(const char* ifname
                  , Double_t thres=0.020
                  , Double_t trig=-150.
                  , Double_t cutoff_GHz=0.350
                  , Double_t yfit1=0.30
                  , Double_t yfit2=0.70
                  , Int_t event1=0, Int_t event2=-1, bool debug=false)
{
#if defined(__CINT__) && !defined(__MAKECINT__)
    cout<< "\n***Warning: This script needs to be compiled. Load it with the plus sign after the name, e.g.\n" <<endl;
    cout<< ".L " << __FILE__ << "+" <<endl;
    return;
#endif

    TFile* ifile = new TFile(ifname);
    if (!ifile) {
        cout<< "Could not find file " << ifname <<endl;
        return;
    }
    TTree* wfm = (TTree*) ifile->Get("wfm");
    if (!wfm) {
        cout<< "Could not find tree 'wfm'" <<endl;
        return;
    }

    Int_t chan = 2;

    TFile* ofile = new TFile(Form("%s-chan%d.time.root",ifname,chan), "recreate");
    TTree* otree = new TTree("res", Form("Resulting tree for %s",ifname));
    otree->SetMarkerColor(46);
    otree->SetLineColor(46);

    TreeRes::book(otree);

    Double_t t[20000];
    Double_t w[20000];

    wfm->SetBranchAddress(Form("t%d",chan), &t);
    wfm->SetBranchAddress(Form("w%d",chan), &w);
    Int_t N = wfm->GetLeaf(Form("t%d",chan))->GetLen();

    if (event2 < event1) event2 = wfm->GetEntries() - 1;

    for (int ientry=event1; ientry<=event2; ++ientry)
    {
        if (wfm->LoadTree(ientry) < 0) break;
        wfm->GetEntry(ientry);

        if (false
            || (ientry-event1) < 10
            || (ientry-event1)%1000 == 0
           ) cout<< "--- processing entry " << ientry <<endl;
        if (ientry-event1+1 >= 10) debug = false;

        TreeRes::clear();
        channel2_fit(N, t, w, thres, trig, cutoff_GHz, yfit1, yfit2, chan, ientry, debug);
        otree->Fill();
    }

    cout<< "Writing " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
    otree->Write();

    wfm->ResetBranchAddresses();
}

void channel3_run(const char* ifname
                  , Double_t thres=0.030
                  , Double_t trig=-10.
                  , Double_t cutoff_GHz=0.500
                  , Double_t yfit1=0.30
                  , Double_t yfit2=0.70
                  , Int_t event1=0, Int_t event2=-1, bool debug=false)
{
#if defined(__CINT__) && !defined(__MAKECINT__)
    cout<< "\n***Warning: This script needs to be compiled. Load it with the plus sign after the name, e.g.\n" <<endl;
    cout<< ".L " << __FILE__ << "+" <<endl;
    return;
#endif

    TFile* ifile = new TFile(ifname);
    if (!ifile) {
        cout<< "Could not find file " << ifname <<endl;
        return;
    }
    TTree* wfm = (TTree*) ifile->Get("wfm");
    if (!wfm) {
        cout<< "Could not find tree 'wfm'" <<endl;
        return;
    }

    Int_t chan = 3;

    TFile* ofile = new TFile(Form("%s-chan%d.time.root",ifname,chan), "recreate");
    TTree* otree = new TTree("res", Form("Resulting tree for %s",ifname));
    otree->SetMarkerColor(46);
    otree->SetLineColor(46);

    TreeRes::book(otree);

    Double_t t[20000];
    Double_t w[20000];

    wfm->SetBranchAddress(Form("t%d",chan), &t);
    wfm->SetBranchAddress(Form("w%d",chan), &w);
    Int_t N = wfm->GetLeaf(Form("t%d",chan))->GetLen();

    if (event2 < event1) event2 = wfm->GetEntries() - 1;

    for (int ientry=event1; ientry<=event2; ++ientry)
    {
        if (wfm->LoadTree(ientry) < 0) break;
        wfm->GetEntry(ientry);

        if (false
            || (ientry-event1) < 10
            || (ientry-event1)%1000 == 0
           ) cout<< "--- processing entry " << ientry <<endl;
        if (ientry-event1+1 >= 10) debug = false;

        TreeRes::clear();
        channel3_fit(N, t, w, thres, trig, cutoff_GHz, yfit1, yfit2, chan, ientry, debug);
        otree->Fill();
    }

    cout<< "Writing " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
    otree->Write();

    wfm->ResetBranchAddresses();
}

void channel4_run(const char* ifname
                  , Double_t thres=0.030
                  , Double_t trig=-10.
                  , Double_t cutoff_GHz=0.500
                  , Double_t yfit1=0.30
                  , Double_t yfit2=0.70
                  , Int_t event1=0, Int_t event2=-1, bool debug=false)
{
#if defined(__CINT__) && !defined(__MAKECINT__)
    cout<< "\n***Warning: This script needs to be compiled. Load it with the plus sign after the name, e.g.\n" <<endl;
    cout<< ".L " << __FILE__ << "+" <<endl;
    return;
#endif

    TFile* ifile = new TFile(ifname);
    if (!ifile) {
        cout<< "Could not find file " << ifname <<endl;
        return;
    }
    TTree* wfm = (TTree*) ifile->Get("wfm");
    if (!wfm) {
        cout<< "Could not find tree 'wfm'" <<endl;
        return;
    }

    Int_t chan = 4;

    TFile* ofile = new TFile(Form("%s-chan%d.time.root",ifname,chan), "recreate");
    TTree* otree = new TTree("res", Form("Resulting tree for %s",ifname));
    otree->SetMarkerColor(46);
    otree->SetLineColor(46);

    TreeRes::book(otree);

    Double_t t[20000];
    Double_t w[20000];

    wfm->SetBranchAddress(Form("t%d",chan), &t);
    wfm->SetBranchAddress(Form("w%d",chan), &w);
    Int_t N = wfm->GetLeaf(Form("t%d",chan))->GetLen();

    if (event2 < event1) event2 = wfm->GetEntries() - 1;

    for (int ientry=event1; ientry<=event2; ++ientry)
    {
        if (wfm->LoadTree(ientry) < 0) break;
        wfm->GetEntry(ientry);

        if (false
            || (ientry-event1) < 10
            || (ientry-event1)%1000 == 0
           ) cout<< "--- processing entry " << ientry <<endl;
        if (ientry-event1+1 >= 10) debug = false;

        TreeRes::clear();
        channel3_fit(N, t, w, thres, trig, cutoff_GHz, yfit1, yfit2, chan, ientry, debug);
        otree->Fill();
    }

    cout<< "Writing " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
    otree->Write();

    wfm->ResetBranchAddresses();
}

void channel23_run(const char* ifname, Int_t event1=0, Int_t event2=-1, bool debug=false)
{
#if defined(__CINT__) && !defined(__MAKECINT__)
    cout<< "\n***Warning: This script needs to be compiled. Load it with the plus sign after the name, e.g.\n" <<endl;
    cout<< ".L " << __FILE__ << "+" <<endl;
    return;
#endif

    TFile* ifile = new TFile(ifname);
    if (!ifile) {
        cout<< "Could not find file " << ifname <<endl;
        return;
    }
    TTree* wfm = (TTree*) ifile->Get("wfm");
    if (!wfm) {
        cout<< "Could not find tree 'wfm'" <<endl;
        return;
    }

    // channels to run

    Int_t chan[4];
    Int_t nchan = 0;
    chan[nchan++] = 2;
    chan[nchan++] = 3;

    // thresholds
    Double_t thres[4];
    //---------------------------------------thres[0] = 0.060;    // mV
    thres[0] = 0.020;    // mV
    //-- thres[0] = 0.003;    // mV          // <-- for single pe peak
    //-- thres[0] = 0.010;    // mV
    //-- thres[0] = 0.100;    // mV
    thres[1] = 0.003;    // mV          // <-- for single pe peak
    //-- thres[1] = 0.020;    // mV

    // trigger time point == starting point
    Double_t trig[4];
    trig[0] = -5.;    // ns
    trig[1] = -5.;     // ns
    //-- trig[0] = -100.;    // ns
    //-- trig[1] = -100.;     // ns

    // cutoff frequency
    Double_t cutoff_GHz[4];
    //-- cutoff_GHz[0] = 1.500;
    //-- cutoff_GHz[0] = 1.750;

    //-- cutoff_GHz[0] = 1.500;
    //-- cutoff_GHz[1] = 0.750;

    //-- cutoff_GHz[0] = 1.500;
    //-- cutoff_GHz[1] = 0.750;

    //-- cutoff_GHz[1] = 0.350;

    //-- cutoff_GHz[0] = 1.0;
    //-- cutoff_GHz[1] = 1.0;
    cutoff_GHz[0] = 10.;
    cutoff_GHz[1] = 10.;
    // cutoff_GHz[0] = 0.750;
    // cutoff_GHz[1] = 0.750;
    //-- cutoff_GHz[0] = 0.500;
    //-- cutoff_GHz[1] = 0.500;

    // fit range
    Double_t yfit1[4];
    Double_t yfit2[4];
    yfit1[0] = 0.3;
    yfit2[0] = 0.7;
    // yfit1[0] = 0.2;
    // yfit2[0] = 0.6;
    yfit1[1] = 0.3;
    yfit2[1] = 0.7;

    void (*channel_fit[4])(Int_t N, Double_t* t, Double_t* w, Double_t thres, Double_t trig, Double_t cutoff_GHz, Double_t yfit1, Double_t yfit2, Int_t chan, Int_t event, bool debug);
    // void (*channel_fit[4])(Int_t N, Double_t* t, Double_t* w, Double_t thres, Double_t trig, Double_t cutoff_GHz, Double_t, Double_t, Int_t chan, Int_t event, bool debug);
    channel_fit[0] = channel2_fit;
    channel_fit[1] = channel3_fit;

    if (event2 < event1) event2 = wfm->GetEntries() - 1;

    //TFile* ofile = new TFile(Form("%s.time.root",ifname), "recreate");
    TFile* ofile;
    if (event1 != 0 || event2 != wfm->GetEntries()-1) ofile = new TFile(Form("%s.time-%d-%d.root",ifname,event1,event2), "recreate");
    else ofile = new TFile(Form("%s.time.root",ifname), "recreate");
    TTree* otree = new TTree("res", Form("Resulting tree for %s",ifname));
    otree->SetMarkerColor(46);
    otree->SetLineColor(46);

    TreeRes::book(otree);

    cout<< "Add to the res->GetUserInfo() sources FunRC.h and " << __FILE__ <<endl;
    TextObject textObject(otree->GetUserInfo());
    textObject.writeFile("FunRC.h");
    textObject.writeFile(__FILE__);

    Double_t t[4][20000];
    Double_t w[4][20000];

    for (int ich=0; ich<nchan; ++ich) {
        wfm->SetBranchAddress(Form("t%d",chan[ich]), t[ich]);
        wfm->SetBranchAddress(Form("w%d",chan[ich]), w[ich]);
    }
    Int_t N = wfm->GetLeaf(Form("t%d",chan[0]))->GetLen();

    for (int ientry=event1; ientry<=event2; ++ientry)
    {
        if (wfm->LoadTree(ientry) < 0) break;
        wfm->GetEntry(ientry);

        if (false
            || (ientry-event1) < 10
            || (ientry-event1)%1000 == 0
           ) cout<< "--- processing entry " << ientry <<endl;
        if (ientry-event1+1 > 30) debug = false;
        //cout<< "--> entry = " << ientry <<endl;

        // if (debug) for (int i=0; i<10; ++i) {
        //     cout<< i << "\t";
        //     for (int ich=0; ich<nchan; ++ich) cout<< " t[" << ich << "][" << i << "] = " << t[ich][i] << " w[" << ich << "][" << i << "] = " << w[ich][i];
        //     cout<<endl;
        // }

        TreeRes::clear();

        for (int ich=0; ich<nchan; ++ich)
        {
            // if (debug && ich == 1) continue;        // to look at LGAD only
            // if (debug && ich == 2) continue;        // to look at Cherenkov only

            if (ich == 0) for (int i=0; i<N; ++i) w[ich][i] = -w[ich][i];   // invert SiPM data
            if (ich == 1) for (int i=0; i<N; ++i) w[ich][i] = -w[ich][i];   // invert SiPM data
            channel_fit[ich](N, t[ich], w[ich], thres[ich], trig[ich], cutoff_GHz[ich], yfit1[ich],yfit2[ich], chan[ich], ientry, debug);
        }

        otree->Fill();
    }

    cout<< "Writing " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
    otree->Write();

    wfm->ResetBranchAddresses();
}

void channel234_run(const char* ifname, Int_t event1=0, Int_t event2=-1, bool debug=false)
{
#if defined(__CINT__) && !defined(__MAKECINT__)
    cout<< "\n***Warning: This script needs to be compiled. Load it with the plus sign after the name, e.g.\n" <<endl;
    cout<< ".L " << __FILE__ << "+" <<endl;
    return;
#endif

    TFile* ifile = new TFile(ifname);
    if (!ifile) {
        cout<< "Could not find file " << ifname <<endl;
        return;
    }
    TTree* wfm = (TTree*) ifile->Get("wfm");
    if (!wfm) {
        cout<< "Could not find tree 'wfm'" <<endl;
        return;
    }

    // channels to run

    Int_t chan[4];
    Int_t nchan = 0;
    chan[nchan++] = 2;
    chan[nchan++] = 3;
    chan[nchan++] = 4;

    // thresholds
    Double_t thres[4];
    //---------------------------------------thres[0] = 0.060;    // mV
    thres[0] = 0.010;    // mV
    //-- thres[0] = 0.010;    // mV
    //-- thres[0] = 0.100;    // mV
    thres[1] = 0.010;    // mV
    thres[2] = 0.010;    // mV

    // trigger time point == starting point
    Double_t trig[4];
    trig[0] = 23.;    // ns
    trig[1] = 23.;     // ns
    trig[2] = 23.;     // ns

    // cutoff frequency
    Double_t cutoff_GHz[4];
    //-- cutoff_GHz[0] = 1.500;
    //-- cutoff_GHz[0] = 1.750;

    //-- cutoff_GHz[0] = 1.500;
    //-- cutoff_GHz[1] = 0.750;

    //-- cutoff_GHz[0] = 1.500;
    //-- cutoff_GHz[1] = 0.750;

    cutoff_GHz[0] = 10.0;
    cutoff_GHz[1] = 10.0;
    cutoff_GHz[2] = 10.0;

    // fit range
    Double_t yfit1[4];
    Double_t yfit2[4];
    yfit1[0] = 0.3;
    yfit2[0] = 0.7;
    // yfit1[0] = 0.2;
    // yfit2[0] = 0.6;
    yfit1[1] = 0.3;
    yfit2[1] = 0.7;
    yfit1[2] = 0.3;
    yfit2[2] = 0.7;

    void (*channel_fit[4])(Int_t N, Double_t* t, Double_t* w, Double_t thres, Double_t trig, Double_t cutoff_GHz, Double_t yfit1, Double_t yfit2, Int_t chan, Int_t event, bool debug);
    // void (*channel_fit[4])(Int_t N, Double_t* t, Double_t* w, Double_t thres, Double_t trig, Double_t cutoff_GHz, Double_t, Double_t, Int_t chan, Int_t event, bool debug);
    channel_fit[0] = channel2_fit;
    channel_fit[1] = channel3_fit;
    channel_fit[2] = channel4_fit;

    if (event2 < event1) event2 = wfm->GetEntries() - 1;

    //TFile* ofile = new TFile(Form("%s.time.root",ifname), "recreate");
    TFile* ofile;
    if (event1 != 0 || event2 != wfm->GetEntries()-1) ofile = new TFile(Form("%s.time-%d-%d.root",ifname,event1,event2), "recreate");
    else ofile = new TFile(Form("%s.time.root",ifname), "recreate");
    TTree* otree = new TTree("res", Form("Resulting tree for %s",ifname));
    otree->SetMarkerColor(46);
    otree->SetLineColor(46);

    TreeRes::book(otree);

    Double_t t[4][20000];
    Double_t w[4][20000];

    for (int ich=0; ich<nchan; ++ich) {
        wfm->SetBranchAddress(Form("t%d",chan[ich]), t[ich]);
        wfm->SetBranchAddress(Form("w%d",chan[ich]), w[ich]);
    }
    Int_t N = wfm->GetLeaf(Form("t%d",chan[0]))->GetLen();

    for (int ientry=event1; ientry<=event2; ++ientry)
    {
        if (wfm->LoadTree(ientry) < 0) break;
        wfm->GetEntry(ientry);

        if (false
            || (ientry-event1) < 10
            || (ientry-event1)%1000 == 0
           ) cout<< "--- processing entry " << ientry <<endl;
        if (ientry-event1+1 > 30) debug = false;
        //cout<< "--> entry = " << ientry <<endl;

        // if (debug) for (int i=0; i<10; ++i) {
        //     cout<< i << "\t";
        //     for (int ich=0; ich<nchan; ++ich) cout<< " t[" << ich << "][" << i << "] = " << t[ich][i] << " w[" << ich << "][" << i << "] = " << w[ich][i];
        //     cout<<endl;
        // }

        TreeRes::clear();

        for (int ich=0; ich<nchan; ++ich)
        {
            // if (debug && ich == 1) continue;        // to look at LGAD only
            // if (debug && ich == 2) continue;        // to look at LGAD only
            // if (debug && ich == 3) continue;        // to look at Cherenkov only

            if (ich == 0) for (int i=0; i<N; ++i) w[ich][i] = -w[ich][i];   // invert LGAD data
            if (ich == 1) for (int i=0; i<N; ++i) w[ich][i] = -w[ich][i];   // invert LGAD data
            if (ich == 2) for (int i=0; i<N; ++i) w[ich][i] = -w[ich][i];   // invert LGAD data
            channel_fit[ich](N, t[ich], w[ich], thres[ich], trig[ich], cutoff_GHz[ich], yfit1[ich],yfit2[ich], chan[ich], ientry, debug);
        }

        otree->Fill();
    }

    cout<< "Writing " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
    otree->Write();

    wfm->ResetBranchAddresses();
}

void channel123_run(const char* ifname, Int_t event1=0, Int_t event2=-1, bool debug=false)
{
#if defined(__CINT__) && !defined(__MAKECINT__)
    cout<< "\n***Warning: This script needs to be compiled. Load it with the plus sign after the name, e.g.\n" <<endl;
    cout<< ".L " << __FILE__ << "+" <<endl;
    return;
#endif

    TFile* ifile = new TFile(ifname);
    if (!ifile) {
        cout<< "Could not find file " << ifname <<endl;
        return;
    }
    TTree* wfm = (TTree*) ifile->Get("wfm");
    if (!wfm) {
        cout<< "Could not find tree 'wfm'" <<endl;
        return;
    }

    // channels to run

    Int_t chan[4];
    Int_t nchan = 0;
    chan[nchan++] = 1;
    chan[nchan++] = 2;
    chan[nchan++] = 3;

    // thresholds
    Double_t thres[4];
    //---------------------------------------thres[0] = 0.060;    // mV
    thres[0] = 0.100;    // mV
    //-- thres[0] = 0.010;    // mV
    //-- thres[0] = 0.100;    // mV
    thres[1] = 0.010;    // mV
    //-- thres[2] = 0.040;    // mV
    thres[2] = 0.010;    // mV

    // trigger time point == starting point
    Double_t trig[4];
    trig[0] = 20.;  // ns
    trig[1] = 20.;  // ns
    trig[2] = 20.;  // ns

    // cutoff frequency
    Double_t cutoff_GHz[4];
    //-- cutoff_GHz[0] = 1.500;
    //-- cutoff_GHz[0] = 1.750;

    //-- cutoff_GHz[0] = 1.500;
    //-- cutoff_GHz[1] = 0.750;

    //-- cutoff_GHz[0] = 1.500;
    //-- cutoff_GHz[1] = 0.750;

    // std
    // cutoff_GHz[0] = 2.000;
    // cutoff_GHz[1] = 0.750;
    // cutoff_GHz[2] = 2.000;

    // increase
    // cutoff_GHz[0] = 2.4;
    // cutoff_GHz[1] = 2.4;
    // cutoff_GHz[2] = 2.4;

    // eliminate
    cutoff_GHz[0] = 10.0;
    cutoff_GHz[1] = 10.0;
    cutoff_GHz[2] = 10.0;

    // fit range
    Double_t yfit1[4];
    Double_t yfit2[4];
    yfit1[0] = 0.3;
    yfit2[0] = 0.7;
    yfit1[1] = 0.3;
    yfit2[1] = 0.7;
    yfit1[2] = 0.3;
    yfit2[2] = 0.7;

    void (*channel_fit[4])(Int_t N, Double_t* t, Double_t* w, Double_t thres, Double_t trig, Double_t cutoff_GHz, Double_t yfit1, Double_t yfit2, Int_t chan, Int_t event, bool debug);
    // void (*channel_fit[4])(Int_t N, Double_t* t, Double_t* w, Double_t thres, Double_t trig, Double_t cutoff_GHz, Double_t, Double_t, Int_t chan, Int_t event, bool debug);
    channel_fit[0] = channel1_fit;
    channel_fit[1] = channel2_fit;
    channel_fit[2] = channel3_fit;

    if (event2 < event1) event2 = wfm->GetEntries() - 1;

    //TFile* ofile = new TFile(Form("%s.time.root",ifname), "recreate");
    TFile* ofile;
    if (event1 != 0 || event2 != wfm->GetEntries()-1) ofile = new TFile(Form("%s.time-%d-%d.root",ifname,event1,event2), "recreate");
    else ofile = new TFile(Form("%s.time.root",ifname), "recreate");
    TTree* otree = new TTree("res", Form("Resulting tree for %s",ifname));
    otree->SetMarkerColor(46);
    otree->SetLineColor(46);

    TreeRes::book(otree);

    Double_t t[4][20000];
    Double_t w[4][20000];

    for (int ich=0; ich<nchan; ++ich) {
        wfm->SetBranchAddress(Form("t%d",chan[ich]), t[ich]);
        wfm->SetBranchAddress(Form("w%d",chan[ich]), w[ich]);
    }
    Int_t N = wfm->GetLeaf(Form("t%d",chan[0]))->GetLen();

    for (int ientry=event1; ientry<=event2; ++ientry)
    {
        if (wfm->LoadTree(ientry) < 0) break;
        wfm->GetEntry(ientry);

        if (false
            || (ientry-event1) < 10
            || (ientry-event1)%1000 == 0
           ) cout<< "--- processing entry " << ientry <<endl;
        if (ientry-event1+1 > 30) debug = false;
        //cout<< "--> entry = " << ientry <<endl;

        // if (debug) for (int i=0; i<10; ++i) {
        //     cout<< i << "\t";
        //     for (int ich=0; ich<nchan; ++ich) cout<< " t[" << ich << "][" << i << "] = " << t[ich][i] << " w[" << ich << "][" << i << "] = " << w[ich][i];
        //     cout<<endl;
        // }

        TreeRes::clear();

        for (int ich=0; ich<nchan; ++ich)
        {
            // if (debug && ich == 1) continue;        // to look at LGAD only
            // if (debug && ich == 2) continue;        // to look at LGAD only
            // if (debug && ich == 3) continue;        // to look at Cherenkov only

            if (ich == 0) for (int i=0; i<N; ++i) w[ich][i] = -w[ich][i];   // invert SiPM data
            if (ich == 1) for (int i=0; i<N; ++i) w[ich][i] = -w[ich][i];   // invert LGAD data
            if (ich == 2) for (int i=0; i<N; ++i) w[ich][i] = -w[ich][i];   // invert LGAD data
            channel_fit[ich](N, t[ich], w[ich], thres[ich], trig[ich], cutoff_GHz[ich], yfit1[ich],yfit2[ich], chan[ich], ientry, debug);
        }

        otree->Fill();
    }

    cout<< "Writing " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
    otree->Write();

    wfm->ResetBranchAddresses();
}

void channel14_run(const char* ifname, Int_t event1=0, Int_t event2=-1, bool debug=false)
{
#if defined(__CINT__) && !defined(__MAKECINT__)
    cout<< "\n***Warning: This script needs to be compiled. Load it with the plus sign after the name, e.g.\n" <<endl;
    cout<< ".L " << __FILE__ << "+" <<endl;
    return;
#endif

    TFile* ifile = new TFile(ifname);
    if (!ifile) {
        cout<< "Could not find file " << ifname <<endl;
        return;
    }
    TTree* wfm = (TTree*) ifile->Get("wfm");
    if (!wfm) {
        cout<< "Could not find tree 'wfm'" <<endl;
        return;
    }

    // channels to run

    Int_t chan[4];
    Int_t nchan = 0;
    chan[nchan++] = 1;
    chan[nchan++] = 4;

    // thresholds
    Double_t thres[4];
    //---------------------------------------thres[0] = 0.060;    // mV
    thres[0] = 0.040;    // mV
    thres[1] = 0.005;    // mV

    // trigger time point == starting point
    Double_t trig[4];
    trig[0] = 20.;  // ns
    trig[1] = 20.;  // ns

    // cutoff frequency
    Double_t cutoff_GHz[4];
    // eliminate frequency cutoff: set it to 10 GHz
    cutoff_GHz[0] = 10.0;
    cutoff_GHz[1] = 10.0;

    // fit range
    Double_t yfit1[4];
    Double_t yfit2[4];
    yfit1[0] = 0.3;
    yfit2[0] = 0.7;
    yfit1[1] = 0.3;
    yfit2[1] = 0.7;

    void (*channel_fit[4])(Int_t N, Double_t* t, Double_t* w, Double_t thres, Double_t trig, Double_t cutoff_GHz, Double_t yfit1, Double_t yfit2, Int_t chan, Int_t event, bool debug);
    // void (*channel_fit[4])(Int_t N, Double_t* t, Double_t* w, Double_t thres, Double_t trig, Double_t cutoff_GHz, Double_t, Double_t, Int_t chan, Int_t event, bool debug);
    channel_fit[0] = channel1_fit;
    channel_fit[1] = channel4_fit;

    if (event2 < event1) event2 = wfm->GetEntries() - 1;

    //TFile* ofile = new TFile(Form("%s.time.root",ifname), "recreate");
    TFile* ofile;
    if (event1 != 0 || event2 != wfm->GetEntries()-1) ofile = new TFile(Form("%s.time-%d-%d.root",ifname,event1,event2), "recreate");
    else ofile = new TFile(Form("%s.time.root",ifname), "recreate");
    TTree* otree = new TTree("res", Form("Resulting tree for %s",ifname));
    otree->SetMarkerColor(46);
    otree->SetLineColor(46);

    TreeRes::book(otree);

    Double_t t[4][20000];
    Double_t w[4][20000];

    for (int ich=0; ich<nchan; ++ich) {
        wfm->SetBranchAddress(Form("t%d",chan[ich]), t[ich]);
        wfm->SetBranchAddress(Form("w%d",chan[ich]), w[ich]);
    }
    Int_t N = wfm->GetLeaf(Form("t%d",chan[0]))->GetLen();

    for (int ientry=event1; ientry<=event2; ++ientry)
    {
        if (wfm->LoadTree(ientry) < 0) break;
        wfm->GetEntry(ientry);

        if (false
            || (ientry-event1) < 10
            || (ientry-event1)%1000 == 0
           ) cout<< "--- processing entry " << ientry <<endl;
        if (ientry-event1+1 > 30) debug = false;
        //cout<< "--> entry = " << ientry <<endl;

        // if (debug) for (int i=0; i<10; ++i) {
        //     cout<< i << "\t";
        //     for (int ich=0; ich<nchan; ++ich) cout<< " t[" << ich << "][" << i << "] = " << t[ich][i] << " w[" << ich << "][" << i << "] = " << w[ich][i];
        //     cout<<endl;
        // }

        TreeRes::clear();

        for (int ich=0; ich<nchan; ++ich)
        {
            // if (debug && ich == 1) continue;        // to look at LGAD only
            // if (debug && ich == 2) continue;        // to look at LGAD only
            // if (debug && ich == 3) continue;        // to look at Cherenkov only

            if (ich == 0) for (int i=0; i<N; ++i) w[ich][i] = -w[ich][i];   // invert SiPM data
            if (ich == 1) for (int i=0; i<N; ++i) w[ich][i] = -w[ich][i];   // invert LGAD data
            channel_fit[ich](N, t[ich], w[ich], thres[ich], trig[ich], cutoff_GHz[ich], yfit1[ich],yfit2[ich], chan[ich], ientry, debug);
        }

        otree->Fill();
    }

    cout<< "Writing " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
    otree->Write();

    wfm->ResetBranchAddresses();
}

void channel1234_run(const char* ifname, Int_t event1=0, Int_t event2=-1, bool debug=false)
{
#if defined(__CINT__) && !defined(__MAKECINT__)
    cout<< "\n***Warning: This script needs to be compiled. Load it with the plus sign after the name, e.g.\n" <<endl;
    cout<< ".L " << __FILE__ << "+" <<endl;
    return;
#endif

    TFile* ifile = new TFile(ifname);
    if (!ifile) {
        cout<< "Could not find file " << ifname <<endl;
        return;
    }
    TTree* wfm = (TTree*) ifile->Get("wfm");
    if (!wfm) {
        cout<< "Could not find tree 'wfm'" <<endl;
        return;
    }

    // channels to run

    Int_t chan[4];
    Int_t nchan = 0;
    chan[nchan++] = 1;
    chan[nchan++] = 2;
    chan[nchan++] = 3;
    chan[nchan++] = 4;

    // thresholds
    Double_t thres[4];
    //---------------------------------------thres[0] = 0.060;    // mV
    thres[0] = 0.100;    // mV
    //-- thres[0] = 0.010;    // mV
    //-- thres[0] = 0.100;    // mV
    thres[1] = 0.100;    // mV
    thres[2] = 0.100;    // mV
    thres[3] = 0.040;    // mV

    // trigger time point == starting point
    Double_t trig[4];
    trig[0] = -3.;    // ns
    trig[1] = -3.;     // ns
    trig[2] = -3.;     // ns
    trig[3] = -3.;     // ns

    // cutoff frequency
    Double_t cutoff_GHz[4];
    //-- cutoff_GHz[0] = 1.500;
    //-- cutoff_GHz[0] = 1.750;

    //-- cutoff_GHz[0] = 1.500;
    //-- cutoff_GHz[1] = 0.750;

    //-- cutoff_GHz[0] = 1.500;
    //-- cutoff_GHz[1] = 0.750;

    // std
    // cutoff_GHz[0] = 2.000;
    // cutoff_GHz[1] = 0.750;
    // cutoff_GHz[2] = 2.000;

    // increase
    // cutoff_GHz[0] = 2.4;
    // cutoff_GHz[1] = 2.4;
    // cutoff_GHz[2] = 2.4;

    // cutoff_GHz[0] = 1.0;
    // cutoff_GHz[1] = 1.0;
    // cutoff_GHz[2] = 1.0;
    // cutoff_GHz[3] = 1.0;

    // eliminate
    cutoff_GHz[0] = 10.0;
    cutoff_GHz[1] = 10.0;
    cutoff_GHz[2] = 10.0;
    cutoff_GHz[3] = 10.0;

    // fit range
    Double_t yfit1[4];
    Double_t yfit2[4];
    yfit1[0] = 0.3;
    yfit2[0] = 0.7;
    yfit1[1] = 0.3;
    yfit2[1] = 0.7;
    yfit1[2] = 0.30;
    yfit2[2] = 0.70;
    yfit1[3] = 0.30;
    yfit2[3] = 0.70;
    // yfit1[0] = 0.1;
    // yfit2[0] = 0.5;
    // yfit1[1] = 0.1;
    // yfit2[1] = 0.5;
    // yfit1[2] = 0.10;
    // yfit2[2] = 0.50;
    // yfit1[3] = 0.10;
    // yfit2[3] = 0.50;

    void (*channel_fit[4])(Int_t N, Double_t* t, Double_t* w, Double_t thres, Double_t trig, Double_t cutoff_GHz, Double_t yfit1, Double_t yfit2, Int_t chan, Int_t event, bool debug);
    // void (*channel_fit[4])(Int_t N, Double_t* t, Double_t* w, Double_t thres, Double_t trig, Double_t cutoff_GHz, Double_t, Double_t, Int_t chan, Int_t event, bool debug);
    channel_fit[0] = channel1_fit;
    channel_fit[1] = channel2_fit;
    channel_fit[2] = channel3_fit;
    channel_fit[3] = channel4_fit;

    if (event2 < event1) event2 = wfm->GetEntries() - 1;

    //TFile* ofile = new TFile(Form("%s.time.root",ifname), "recreate");
    TFile* ofile;
    if (event1 != 0 || event2 != wfm->GetEntries()-1) ofile = new TFile(Form("%s.time-%d-%d.root",ifname,event1,event2), "recreate");
    else ofile = new TFile(Form("%s.time.root",ifname), "recreate");
    TTree* otree = new TTree("res", Form("Resulting tree for %s",ifname));
    otree->SetMarkerColor(46);
    otree->SetLineColor(46);

    TreeRes::book(otree);

    Double_t t[4][20000];
    Double_t w[4][20000];

    for (int ich=0; ich<nchan; ++ich) {
        wfm->SetBranchAddress(Form("t%d",chan[ich]), t[ich]);
        wfm->SetBranchAddress(Form("w%d",chan[ich]), w[ich]);
    }
    Int_t N = wfm->GetLeaf(Form("t%d",chan[0]))->GetLen();

    for (int ientry=event1; ientry<=event2; ++ientry)
    {
        if (wfm->LoadTree(ientry) < 0) break;
        wfm->GetEntry(ientry);

        if (false
            || (ientry-event1) < 10
            || (ientry-event1)%1000 == 0
           ) cout<< "--- processing entry " << ientry <<endl;
        if (ientry-event1+1 > 30) debug = false;
        //cout<< "--> entry = " << ientry <<endl;

        // if (debug) for (int i=0; i<10; ++i) {
        //     cout<< i << "\t";
        //     for (int ich=0; ich<nchan; ++ich) cout<< " t[" << ich << "][" << i << "] = " << t[ich][i] << " w[" << ich << "][" << i << "] = " << w[ich][i];
        //     cout<<endl;
        // }

        TreeRes::clear();

        for (int ich=0; ich<nchan; ++ich)
        {
            // if (debug && ich == 1) continue;        // to look at LGAD only
            // if (debug && ich == 2) continue;        // to look at LGAD only
            // if (debug && ich == 3) continue;        // to look at Cherenkov only

            if (ich == 3) for (int i=0; i<N; ++i) w[ich][i] = -w[ich][i];   // invert SiPM data
            // if (ich == 3) for (int i=0; i<N; ++i) w[ich][i] = -w[ich][i];   // invert Cardinelli data
            channel_fit[ich](N, t[ich], w[ich], thres[ich], trig[ich], cutoff_GHz[ich], yfit1[ich],yfit2[ich], chan[ich], ientry, debug);
        }

        otree->Fill();
    }

    cout<< "Writing " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
    otree->Write();

    wfm->ResetBranchAddresses();
}

#endif  // macro_wfmFit_C
