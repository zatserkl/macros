// Andriy Zatserklyaniy <zatserkl@fnal.gov> Feb 27, 2016

#ifndef macro_wfmLoop_C
#define macro_wfmLoop_C

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TTimer.h>
#include <TVirtualFFT.h>
#include <TMath.h>

#include <iostream>
#include <sstream>
#include <string>
#include <cstdio>
#include <vector>

using std::cout;     using std::endl;

// utils.C stuff
TGraph* gtemp(TCanvas* can=0, Bool_t updateNameTitle=kTRUE);
Int_t countpads(TVirtualPad *pad);

#ifndef macro_utils_C

// functions gtemp and countpads (with dependencies) from utils.C

Int_t countpads(TVirtualPad *pad)
{  // Rene Brun, http://root.cern.ch/root/roottalk/roottalk02/0654.html
    // count the number of pads in pad
    if (!pad) return 0;
    Int_t npads = 0;
    TObject *obj;
    TIter next(pad->GetListOfPrimitives());
    while ((obj = next())) {
        if (obj->InheritsFrom(TVirtualPad::Class())) npads++;
    }
    if (npads == 0) npads = 1;                                   //--AZ Feb 27, 2016
    return npads;
}

TH1* htemp(TCanvas* can=0);

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

TH1* htemp(Int_t index) {return htemp(index,0);}

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

TGraph* gtemp(TCanvas* can, Bool_t updateNameTitle) {return gtemp(-1, can, updateNameTitle);}

#endif  // macro_utils_C

namespace namespace_wfmLoop
{
    //-- Double_t cutoff_GHz[8] = {0.750, 0.500, 0};
    Double_t cutoff_GHz[8] = {0.750, 0.350, 0};
}   // namespace namespace_wfmLoop

TTree* getTree()
{
    TTree* tree = (TTree*) gDirectory->Get("wfm");
    if (!tree) {
        // the wfm tree may be a friend of the res tree
        TTree* res = (TTree*) gDirectory->Get("res");
        if (res) {
            tree = res->GetFriend("wfm");
        }
    }
    return tree;
}

//--------- Fast Fourier Transform -----------

Double_t fft_filter_wfmLoop(Int_t N, Double_t* x, Double_t* y, Double_t fcut=0.2);  // fcut in GHz
void moving_average_wfmLoop(Int_t naver, Int_t np, const Double_t* y, Double_t* yaver);   // Double_t version
TGraph* filter_graph_wfmLoop(Double_t fcut=0.1, TGraph* g=0);  // fcut in GHz

void compare(Int_t evt
             , Int_t chan1=1
             , Int_t chan2=0
             , Int_t chan3=0
             , Int_t chan4=0
             , Int_t chan5=0
             , Int_t chan6=0
             , Int_t chan7=0
             , Int_t chan8=0
             //, Double_t cutoff_GHz=0
             , bool cutoff_GHz=false
             , TTree* tree=0
             , const char* wname="compare"
            )
{
    if (!tree) tree = getTree();
    if (!tree) {
        cout<< "Could not find wfm tree" <<endl;
        return;
    }

    Int_t channels[8];
    Double_t sign[8];
    Int_t nchannels = 0;

    if (chan1 != 0) {channels[nchannels] = TMath::Abs(chan1); sign[nchannels] = chan1 > 0? 1.: -1; ++nchannels;}
    if (chan2 != 0) {channels[nchannels] = TMath::Abs(chan2); sign[nchannels] = chan2 > 0? 1.: -1; ++nchannels;}
    if (chan3 != 0) {channels[nchannels] = TMath::Abs(chan3); sign[nchannels] = chan3 > 0? 1.: -1; ++nchannels;}
    if (chan4 != 0) {channels[nchannels] = TMath::Abs(chan4); sign[nchannels] = chan4 > 0? 1.: -1; ++nchannels;}
    if (chan5 != 0) {channels[nchannels] = TMath::Abs(chan5); sign[nchannels] = chan5 > 0? 1.: -1; ++nchannels;}
    if (chan6 != 0) {channels[nchannels] = TMath::Abs(chan6); sign[nchannels] = chan6 > 0? 1.: -1; ++nchannels;}
    if (chan7 != 0) {channels[nchannels] = TMath::Abs(chan7); sign[nchannels] = chan7 > 0? 1.: -1; ++nchannels;}
    if (chan8 != 0) {channels[nchannels] = TMath::Abs(chan8); sign[nchannels] = chan8 > 0? 1.: -1; ++nchannels;}

    for (int ich=0; ich<nchannels; ++ich) if (!tree->GetBranch(Form("w%d",channels[ich]))) {
        cout<< "Could not find branch " << Form("w%d",channels[ich]) <<endl;
        return;
    }

    //cout<< "nchannels = " << nchannels <<endl;
    //for (int ich=0; ich<nchannels; ++ich) cout<< "channels[" << ich << "] = " << channels[ich] <<endl;

    Int_t ndivy = 2;                          // 2 divisions along the y-axis
    Int_t ndivx = (nchannels + 1) / ndivy;
    if (nchannels == 1) {
        ndivx = 1;
        ndivy = 1;
    }
    TCanvas* can = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(wname);
    if (can) {
        // cout<< "countpads(can) = " << countpads(can) << " ndivx = " << ndivx << " ndivy = " << ndivy <<endl;
        if (countpads(can) == ndivx*ndivy) gPad = can;
        else {
            delete can;
            can = 0;
        }
    }
    if (!can) {
        // cout<< "create new canvas" <<endl;
        can = new TCanvas(wname,wname, 0,0, ndivx*gStyle->GetCanvasDefW(), ndivy*gStyle->GetCanvasDefH());
        if (nchannels > 1) can->Divide(ndivx,ndivy);
    }

    for (int ix=0; ix<ndivx; ++ix) {
        for (int iy=0; iy<ndivy; ++iy)
        {
            can->cd(ix + ndivx*iy + 1);

            Int_t ich = ix*ndivy + iy;
            if (ich < nchannels) {
                // Int_t board = (channels[ich]-1)/4 + 1;

                if (sign[ich] > 0) tree->Draw(Form("w%d:t%d",channels[ich],channels[ich]), Form("Entry$==%d",evt), "lp");    //--- Draw command
                else tree->Draw(Form("-w%d:t%d",channels[ich],channels[ich]), Form("Entry$==%d",evt), "lp");    //--- Draw command

                gtemp()->SetTitle(tree->GetHistogram()->GetTitle());
                gtemp()->SetName(Form("evt_%d_chan_%d",evt,channels[ich]));
                // cout<< "ix = " << ix << " iy = " << iy << " ich = " << ich << " old channel = " << (channels[ich]-1) % 4 + 1 << " board = " << board <<endl;
                //-- if (cutoff_GHz > 0) filter_graph_wfmLoop(cutoff_GHz);
                if (cutoff_GHz) filter_graph_wfmLoop(namespace_wfmLoop::cutoff_GHz[ich]);
                //-- if (cutoff_GHz) {
                //--     cout<< "namespace_wfmLoop::cutoff_GHz[" << ich << "] = " << namespace_wfmLoop::cutoff_GHz[ich] <<endl;
                //--     filter_graph_wfmLoop(namespace_wfmLoop::cutoff_GHz[ich]);
                //-- }
            }
        }
    }
    can->cd();
    can->Update();
}

void wfmLoop(Int_t evtNo=0
             , Int_t chan1=1
             , Int_t chan2=0
             , Int_t chan3=0
             , Int_t chan4=0
             , Int_t chan5=0
             , Int_t chan6=0
             , Int_t chan7=0
             , Int_t chan8=0
             //, Double_t cutoff_GHz=0
             , bool cutoff_GHz=false
             , TTree* tree=0
             , const char* wname="compare"
            )
{
#if defined(__CINT__) && !defined(__MAKECINT__)
    cout<< "\n***Warning: This script needs to be compiled. Load it with the plus sign after the name, e.g.\n" <<endl;
    cout<< ".L " << __FILE__ << "+" <<endl;
    return;
#endif

    TTimer timer("gSystem->ProcessEvents();",50,kFALSE);  //-- process mouse events every 50 ms

    if (cutoff_GHz) for (int ich=0; ich<8; ++ich) if (namespace_wfmLoop::cutoff_GHz[ich] > 0) cout<< "namespace_wfmLoop::cutoff_GHz[" << ich << "] = " << namespace_wfmLoop::cutoff_GHz[ich] <<endl;

    std::vector<std::string> commands;

    bool plot = true;

    while (true)
    {
        //--wrong-- timer.Start();    //-- start processing of mouse events. NB: Wrong place for timer.Start before the compare!!!!!

        std::string line;

        if (plot) compare(evtNo,chan1,chan2,chan3,chan4,chan5,chan6,chan7,chan8,cutoff_GHz,tree,wname);

        plot = true;

        cout<< "<CR>: " << evtNo+1 << " -, Clone, Save, Quit, command: ";
        timer.Start();                                                    // start processing of mouse events
        std::getline(cin, line);                                          //    enclose the interactive part (getline) in timer.Start() -- timer.Stop()
        timer.Stop();                                                     // disable processing of mouse events

        // Possible inputs
        //
        // 0) <CR>: line.size() == 0:
        //    show next event
        // 1) number:
        //    interpret as event number: show this event
        // 2) not a number:
        //    can be one character or more than one character
        //       2.1) one character, line.size() == 1:
        //            interpret as menu command: C or S or Q
        //       2.2) more than one character, line.size() > 1:
        //            interpret as a ROOT command, try to execute

        // 0) <CR>
        if (line.size() == 0)
        {
            ++evtNo;
            continue;
        }

        // 1) number
        std::stringstream ss(line);
        Int_t number = -1;
        if (ss >> number && ss.eof()) {
            evtNo = number;
            continue;
        }

        // 2) not a number
        if (line.size() == 1)
        {
            // 2.1) input was just one character: interpret as a menu command
            switch (toupper(line[0])) {
                case '-':
                    if (evtNo > 0) evtNo -= 1;          // previous event
                    break;
                case 'C':
                    if (gROOT->GetListOfCanvases()->FindObject(wname)) {
                        ((TCanvas*) gROOT->GetListOfCanvases()->FindObject(wname))->cd();    // clone the whole canvas
                        gPad->DrawClone();
                    }
                    else cout<< "Could not find TCanvas " << wname <<endl;
                    plot = false;
                    break;
                case 'S':
                    if (getTree()) gPad->SaveAs(Form("%s_evt_%d.png", getTree()->GetCurrentFile()->GetName(), evtNo));
                    else gPad->SaveAs(Form("evt_%d.png", evtNo));
                    plot = false;
                    break;
                case 'Q':
                    return;
                default: cout<< "No such key" <<endl;
            }
            continue;
        }
        else {
            // 2.2) input was more than one character: interpret as a ROOT command
            if (line == std::string(".q")) {
                cout<< "To terminate the ROOT session exit the macro first" <<endl;
                break;
            }
            else {
                if (unsigned(line[0]) == 27 and unsigned(line[1]) == 91) {
                    // input was some arrow: ESC + '[' + letter A or B or C or D
                    for (unsigned i=0; i<commands.size(); ++i) cout<< commands[i] <<endl;
                    continue;
                }

                commands.push_back(line);
                gROOT->ProcessLine(line.c_str());
                plot = false;
                continue;
            }
        }
    }
}

//
//  FFT tools
//

Double_t fft_filter_wfmLoop(Int_t N, Double_t* x, Double_t* y, Double_t fcut)  // fcut in GHz
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

void moving_average_wfmLoop(Int_t naver, Int_t np, const Double_t* y, Double_t* yaver)   // Double_t version
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

TGraph* filter_graph_wfmLoop(Double_t fcut, TGraph* g)  // fcut in GHz
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
        mag[isample] = TMath::Sqrt(re[isample]*re[isample] + im[isample]*im[isample]);
        // phase[isample] = TMath::ATan2(im[isample],re[isample]);
        // cut the frequency band
        if (isample > scut) {
            re[isample] = 0;
            im[isample] = 0;
        }
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

#endif  // macro_wfmLoop_C
