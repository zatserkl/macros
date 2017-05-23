// Andriy Zatserklyaniy <zatserkl@fnal.gov>

#ifndef macro_utils_C
#define macro_utils_C

#include <TROOT.h>
#include <TRint.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TEnv.h>
#include <TMath.h>
#include <TF1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TEventList.h>
#include <TBrowser.h>
#include <TPaveStats.h>
#include <TCanvas.h>
#include <TList.h>
#include <TColor.h>
#include <TTimer.h>
#include <TLatex.h>
#include <TKey.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <cstdio>

using std::cout;                using std::endl;

/* 
   to use utils.C function in the root macro declare them like
   TH1* htemp(TCanvas* can=0);
*/ 

/*
class HTree: public TH1 {
private:
   Int_t Fill(Double_t, Double_t) {return 0;}
   Int_t Fill(const char*, Double_t) {return 0;}
public:
   TTree* tree;
   Float_t x;
public:
   HTree(const char* name="ht", const char* title="ht", const char* varname="x"): TH1() {
      if (name && *name) SetName(name);
      if (title && *title) SetTitle(title);
      tree = new TTree(Form("%s_tree",GetName()),Form("%s tree",GetName()));
      tree->Branch(varname, &x, Form("%s/F",GetName()));
      tree->SetDirectory(0);
      tree->SetEstimate(-1);           // use all entries to estimate the range

      gDirectory->Append(this);
   }
   ~HTree() {
      tree->SetDirectory(0);
      delete tree;
   }
   void SetBranchName(const char* varname) {
      TBranch* branch = (TBranch*) tree->GetListOfBranches()->First();
      branch->SetName(varname);
   }
   Int_t Fill(Double_t value) {
      x = value;
      tree->Fill();
      return tree->GetEntries();
   }
   virtual void Draw(Option_t* option="") {
      TBranch* branch = (TBranch*) tree->GetListOfBranches()->First();
      const char* bname = branch->GetName();
      tree->Draw(bname, option);
      if (*GetTitle()) tree->GetHistogram()->SetTitle(GetTitle());
      if (*GetXaxis()->GetTitle()) tree->GetHistogram()->GetXaxis()->SetTitle(GetXaxis()->GetTitle());
      if (*GetYaxis()->GetTitle()) tree->GetHistogram()->GetYaxis()->SetTitle(GetYaxis()->GetTitle());
   }
   void SetLineColor(Color_t color) {tree->SetLineColor(color);}
   void SetFillColor(Color_t color) {tree->SetFillColor(color);}
   void SetFillStyle(Style_t style) {tree->SetFillStyle(style);}
   ClassDef(HTree, 1);
};
ClassImp(HTree);
*/

// HTree* hTree(const char* name="ht", const char* title="ht", const char* varname="x");
// HTree* hTree(const char* name, const char* title, const char* varname)
// {
//    return new HTree(name, title, varname);
// }

/// class TextObject {
///     //
///     // To write source code into TTree UserInfo Object.
///     // Also use static method TextObject::print(list) to print the content of the list.
///     //
/// public:
///     TList* list;
///     TextObject(TList* list_0): list(list_0) {
///         // write command line
///         const char* fname = gEnv->GetValue("Rint.History", "");
///         // cout<< "Rint.History file is " << fname <<endl;
///         std::ifstream file(fname);
///         std::string line = "N/A";
///         if (file) {
///             std::string line_read;
///             while (std::getline(file, line_read)) line = line_read;
///         }
///         else cout<< "\n***Warning TextObject::TextObject: could not find history file " << fname <<endl<<endl;
///         // cout<< "Command line is: " << line <<endl;
/// 
///         // create an object and write it into the tree
///         TNamed* textObject = new TNamed("Command line", line.c_str());
///         list->AddLast(textObject);
///         file.close();
///     }
///     void writeFile(const char* fname) {
///         std::ifstream file(fname);
///         if (file) {
///             std::stringstream ss;
///             std::string line;
///             while (std::getline(file, line)) ss << line << endl;
/// 
///             // create an object and write it into the tree
///             TNamed* textObject = new TNamed(fname, ss.str());
///             list->AddLast(textObject);
///             file.close();
///         }
///         else cout<< "\n***Warning TextObject::TextObject: could not find file " << fname <<endl<<endl;
///     }
///     static void print(TList* the_list)
///     {
///         cout<< "the number of entries in the list is " << the_list->GetEntries() <<endl;
/// 
///         TIter next(the_list);
///         TNamed* object;
///         while ((object = (TNamed*) next())) {
///             cout<< "\nobject->GetName() = " << object->GetName() <<endl;
///             cout<< object->GetTitle() <<endl;
///         }
///     }
/// };

TGraph* differentiate(const TGraph* gr);
TGraph* differentiate(const TGraph* gr)
{
    TGraph* gr_diff = (TGraph*) gr->Clone(Form("diff_%s",gr->GetName()));
    gr_diff->SetTitle(Form("diff %s",gr->GetName()));
    for (int i=0; i<gr->GetN()-1; ++i) {
        Double_t dx = i > 0? gr->GetX()[i] - gr->GetX()[i-1]: gr->GetX()[i+1] - gr->GetX()[i];
        Double_t dy = i > 0? gr->GetY()[i] - gr->GetY()[i-1]: 0;
        const Double_t eps = 1e-7;
        Double_t dydx = TMath::Abs(dx) > eps? dy/dx: dy/eps;
        gr_diff->SetPoint(i, gr->GetX()[i], dydx);
    }
    gr_diff->Set(gr->GetN()-1);
    return gr_diff;
}

// 
// modified functor from /srv/opt/root/root/tutorials/math/exampleFunctor.C
// Comments from there:
//
// in order to work with interpreter the function object must be created and lived all time for all time 
// of the TF1. In compiled mode, the function object can be passed by value (reccomended) and there 
// is also no need to specify the type of the function class. Example is as follow: 
// TF1 * f2 = new TF1("f2",MyDerivFunc(f1), xmin, xmax,0); // only for C++ compiled mode

/*
DerivativeFunctor* dfunctor = new DerivativeFunctor(fun);
TF1* dfun = new TF1("dfun", dfunctor, fun->GetXmin(), fun->GetXmax(), 0, "DerivativeFunctor");
 */
class DerivativeFunctor { 
    TF1* tf1_; 
public:
    DerivativeFunctor(TF1* tf1): tf1_(tf1) {}
    double operator() (double* x, double*) const { 
        return tf1_->Derivative(*x);
    }
};
/*
IntegralFunctor* ifunctor = new IntegralFunctor(fun);
TF1* ifun = new TF1("ifun", ifunctor, fun->GetXmin(), fun->GetXmax(), 0, "IntegralFunctor");
 */
class IntegralFunctor { 
    TF1* tf1_; 
public:
    IntegralFunctor(TF1* tf1): tf1_(tf1) {}
    double operator() (double* x, double*) const { 
        double a = tf1_->GetXmin();
        return tf1_->Integral(a, *x);
    }
};
TF1* derfun(TF1* fun);
TF1* derfun(TF1* fun) {
    DerivativeFunctor* dfunctor = new DerivativeFunctor(fun);
    TF1* dfun = new TF1(Form("d%s",fun->GetName()), dfunctor, fun->GetXmin(), fun->GetXmax(), 0, "DerivativeFunctor");
    return dfun;
}
TF1* intfun(TF1* fun);
TF1* intfun(TF1* fun) {
    IntegralFunctor* ifunctor = new IntegralFunctor(fun);
    TF1* ifun = new TF1(Form("i%s",fun->GetName()), ifunctor, fun->GetXmin(), fun->GetXmax(), 0, "IntegralFunctor");
    return ifun;
}

std::ostream& endn(std::ostream& os);
std::ostream& endn(std::ostream& os)
{ 
    // example of usage
    //
    // bool debug = true;
    // std::ostream& (*end)(std::ostream&);      // end is a pointer to function
    // if (debug) end = std::endl;
    // else end = endn;

    os << "\n";
    return os;
}

std::ostream& exitl(std::ostream&);
std::ostream& exitl(std::ostream&) {exit(0);}
// std::ostream& exitl(std::ostream& os) {
//    os << endl;
//    exit(0);
// }

// void GetCanvasDefWH(Int_t& w, Int_t& h);
// void GetCanvasDefW();
// void GetCanvasDefH();

void time();
void notime();
bool fexists(const char* ifname);
const char* substitute(const std::string& str0, const std::string what, const std::string by);
int nlen(const std::string str);                      // string length without leading and trailing spaces
void tcol(Int_t color=602, TTree* tree=0);            // set tree color
void gcol(Int_t color=602, TGraph* gr=0);             // set graph color
void hcol(Int_t color=602, TH1* h=0, Int_t fill=-1);  // set histogram color
void showbra(TTree* tree, const char* pattern="");    // show tree branches with pattern
void bra(TTree* tree=0, const char* pattern="");
void history(Int_t nlines=10);                        // show history
void history(const char* pattern, Int_t nlines=10);
void hi(Int_t nlines=10);                             // abbriviation of the history
void hi(const char* pattern, Int_t nlines=10);
void hiDraw(Int_t nlines=10);                         // show history of the Draw commands
void hiDraw(const char* pattern, Int_t nlines=10);
void ref();                                           // refresh current canvas
void log();                                           // toggle log/lin
void logy();
void logz();
void logx();
void gri();                                           // set grid
void nogri();
Int_t nbi(Int_t nx=0, Int_t ny=0);                    // gEnv->SetValue("Hist.Binning.1D.x", n)
void prompt(const char* prompt="// ");
void rex();                                           // redraw axis
Int_t axisNdivisions(Int_t Ndivisions=0);
Int_t ndi(Int_t npad=0);                              // change x-axis grid numbers
void clone();                                         // clone current canvas
void statover();                                      // include underflows and overflows in calculation of statistics
void nostatover();
void zoomx(Axis_t xmin=0, Axis_t xmax=0);             // zoom x-axis region
void zoomy(Axis_t ymin=0, Axis_t ymax=0);
void zoom(Axis_t xmin=0, Axis_t xmax=0, Axis_t ymin=0, Axis_t ymax=0);
void zoomall(Axis_t xmin=0, Axis_t xmax=0, Axis_t ymin=0, Axis_t ymax=0);
void unzoom();                                        // same as unzoomx()
void unzoomx();
void unzoomy();
void unzoomxy();
void minmax(Axis_t ymin, Axis_t ymax);                // set minimum and maximum for y-axis
TCanvas* fonts();                                     // show all fonts
void colors();                                        // show all colors
Int_t Color(Int_t number);                            // select nice consecutive colors
void fwlite();                                        // CMS specific
void tb(const char* fname=0);                         // loads old (better than new) TBrowser
Int_t canvas_max_number(const char* pattern="c1_n");
const char* nextname(const char* base="h");
const char* nextcan(const char* base="can");
const char* defcan();
// TCanvas* newcan(Int_t width=525, Int_t height=375, const char* name="", const char* title="");
// TCanvas* newcan(Int_t width=522, Int_t height=354, const char* name="", const char* title="");
TCanvas* newcan(Int_t width=0, Int_t height=0, const char* name="", const char* title="");
TCanvas* can(Int_t width=0, Int_t height=0, const char* name="", const char* title="");
TString NextName(TString base="h");

Double_t qsum(Double_t e1, Double_t e2);              // sqrt(e1*e1 + e2*e2)
void sqroot(double a, double b, double c, double& x1, double& x2, bool verb=false); // ax^2+bx+c = 0
void sqroot(double a, double b, double c);

TPaveStats* movestat(Double_t x1NDC=0, Double_t y1NDC=0, Double_t xlen=0, Double_t ylen=0);
TPaveStats* moveleft(Double_t xleft=0.13, TCanvas* can=0);
TPaveStats* moveright(Double_t xright=0.98, TCanvas* can=0);
void leftgaus();                                         // move stat box to the left (optimized for gaus fit)
void rightgaus();
void left();                                             // short for leftgaus
void right();
TPaveText* titsize(Double_t x1NDC, Double_t x2NDC=0, Double_t y1NDC=0, Double_t y2NDC=0);
//TPaveText* titmax(Double_t x1NDC=0, Double_t x2NDC=1., Double_t y1NDC=0.920, Double_t y2NDC=0.999);
//TPaveText* titmax(Double_t x1NDC=0, Double_t x2NDC=1., Double_t y1NDC=0.901, Double_t y2NDC=0.999);
//-- TPaveText* titmax(Double_t x1NDC=0, Double_t x2NDC=1., Double_t y1NDC=0.91, Double_t y2NDC=1.);
TPaveText* titmax(Double_t x1NDC=0, Double_t x2NDC=1., Double_t y1NDC=0.92, Double_t y2NDC=1.);
TPaveText* titlemax(Double_t x1NDC=0, Double_t x2NDC=1., Double_t y1NDC=0.92, Double_t y2NDC=1.);
void cd(const char* filename="");                        // provides default parameter
void cd(Int_t fnumber);                                  // NB: no default parameter
void resize(Int_t width=0, Int_t height=0, Bool_t internal=kFALSE);  // resizes current canvas
void square(Int_t size=0, Bool_t internal=kFALSE);

// subrtact function from the graph or from the histogram
TGraphErrors* gsubf(const TGraph* gr=0, const TF1* f=0, const char* name="", const char* title="");   // subtracts function from graph
TH1* hsubf(const TH1* h=0, const TF1* f=0, const char* name="", const char* title="");

void gsubg(TGraph* g, const TGraph* gsub);

Int_t countpads(TVirtualPad *pad);
TCanvas* gpad(TCanvas* can=0, bool verb=true);                                         // returns pointer to the current canvas
void lslist(TList* list);
TF1* lsfun(const char* fname, Int_t nfun=-1);
TF1* lsfun(Int_t nfun=-1, bool verb=true);
void lstemp(TCanvas* can=0);                             // shows content of the current canvas
TGraph* gtemp(Int_t index, TCanvas* can, Bool_t updateNameTitle=kTRUE);
TGraph* gtemp(Int_t index, TCanvas* can, const char* name, const char* title="");   // change name, title
TGraph* gtemp(TCanvas* can=0, Bool_t updateNameTitle=kTRUE);
TGraph* gtemp(Int_t index, Bool_t updateNameTitle=kTRUE);
TGraph* gtempClone(Int_t index, const char* name, const char* title, TCanvas* can);
TGraph* gtempClone(const char* name="", const char* title="", TCanvas* can=0);
TGraph* gtempClone(Int_t index, const char* name="", const char* title="");
//-- TGraphErrors* getemp(const char* name=0, const char* title=0);
//-- TGraphErrors* getempClone(const char* name=0, const char* title=0);
TH1* htemp(Int_t index, TCanvas* can);
TH1* htemp(TCanvas* can=0);
TH1* htemp(Int_t index);
/// TH1* htempClone(const char* name, const char* title="", TCanvas* can=0);
/// TH1* htempClone(const char* name, const char* title, Int_t index);
/// TH1* htempClone(const char* name, const char* title, TCanvas* can, Int_t index);
TH1* htempClone(const char* name="", const char* title="", TCanvas* can=0);
TH1* htempClone(Int_t index, const char* name, const char* title);
TH1* htempClone(Int_t index, const char* name, const char* title, TCanvas* can);
TH2* htemp2(Int_t index, TCanvas* can);
TH2* htemp2(TCanvas* can=0);
TH2* htemp2(Int_t index);
TH2* htemp2Clone(const char* name="", const char* title="", TCanvas* can=0);
TH2* htemp2Clone(Int_t index, const char* name, const char* title);
TH2* htemp2Clone(Int_t index, const char* name, const char* title, TCanvas* can);
TH2* htemp2Create(const char* name="", const char* title="");
Double_t wtemp();                                        // returns width of the bin #1 of htemp()
TF1* ftemp(TCanvas* can, Int_t index);
TF1* ftemp(TCanvas* can=0);
TF1* ftemp(Int_t index);
TMultiGraph* mgtemp(const char* name=0, const char* title=0);
TPaveStats* stemp();                                            // return pointer to the stat box
TPaveStats* stemp(Int_t index);
void statcol();                         // color stabox into (last) hist or graph line color
void statcol(Int_t index);              // color stabox into (index) hist or graph line color

void drawpan(TH1* h=0);
void fitpan(TH1* h=0);
void InitPredefinedFunctions();
void fitgaus(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="");
void fitg(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="");
void fitgl(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="");
void fitgr(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="");
void fitexpo(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="");
void fite(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="");
void fitel(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="");
void fiter(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="");
void fitpol(Double_t xmin=0, Double_t xmax=0, Int_t power=1, const char* opt="", const char* gopt="");
void fitp(Double_t xmin=0, Double_t xmax=0, Int_t power=1, const char* opt="", const char* gopt="");
Double_t gaus_ampl();
Double_t gaus_mean();
Double_t gaus_sigma();
Double_t gaus_eampl();
Double_t gaus_emean();
Double_t gaus_esigma();
Double_t gaus_area();                                    // divide to wtemp() to get Nevents in histogram
Double_t ngaus(TH1* h=0);                                // the number of events in gaus fit to histogram
void pargaus(Double_t& a, Double_t& mean, Double_t& sigma);
void npe(TH1* h=0, Double_t* p_Npe=0, Double_t* p_eNpe=0);
void hnpe(TH1* h=0, Double_t* p_Npe=0, Double_t* p_eNpe=0);
void pop(TPad* pad);
void pic(TString pathname="", TString suffix="");
void eps(TString pathname="", TString suffix="");
void png(TString pathname="", TString suffix="", bool check_fexists=true);
void fpng(TString pathname="", TString suffix="", bool check_fexists=true);
void pngall();
Double_t intgaus(TH1* h=0);
void print_htemp(TH1* h=0);
void hprint(TH1* h=0);
void tax(const char* format="");
void nul(Axis_t min=0);
void ct();
void ctadd();
void ctlist();
const char* addtit(const char* title, TCanvas* can=0);
TChain* chlist(const char* tree_name, const char* filelist);
// Stat_t hint(char* hname, Double_t x1=0, Double_t x2=0);
// Stat_t hint(TH1* h, Double_t x1=0, Double_t x2=0);
Stat_t hint(Double_t x1=0, Double_t x2=0, TH1* h=0);
void h2subtract(Double_t pedestal, TH2* h2=0);           // subtracts pedestal from TH2 cells. Set negative values to 0
void cloop(const char* command);

// void GetCanvasDefWH(Int_t& w, Int_t& h)
// {
//    // w = 700;    h = 500;    // default
//    // w = 526;    h = 382;    // 3/4 of internal area of the default: 526x382 --> 522x354 as newcan in utils.C
//    w = 333;    h = 266;    // 2/3 of internal area of the default 700x500
// }
// 
// Int_t GetCanvasDefW()
// {
//    Int_t w, h;
//    GetCanvasDefWH(w,h);
//    return w;
// }
// Int_t GetCanvasDefH()
// {
//    Int_t w, h;
//    GetCanvasDefWH(w,h);
//    return w;
// }

void tcol(Int_t color, TTree* tree)
{
    if (color == 0) color = 602;
    if (!tree) {
        TObjArray trees;
        TIter next(gDirectory->GetList());
        for (TObject* obj=0; (obj = next());) {
            if (obj->IsA()->InheritsFrom(TTree::Class())) trees.Add(obj);
        }
        if (trees.GetEntries() == 0) {
            TIter nextkey(gDirectory->GetListOfKeys());
            for (TObject* obj=0; (obj = nextkey());) {
                TKey* key = (TKey*) obj;
                if (strcmp(key->GetClassName(),"TTree")==0) {
                    // add only one instance (highest cycle)
                    TObject* key_obj = gDirectory->Get(key->GetName());
                    if (!trees.FindObject(key->GetName())) trees.Add(key_obj);
                }
            }
            if (trees.GetEntries() == 0) {
                cout<< "No TTree objects found" <<endl;
                return;
            }
        }
        if (trees.GetEntries() == 1) tree = (TTree*) trees.At(0);
        else {
            for (int i=0; i<trees.GetEntries(); ++i) {
                cout<< i << "   " << trees.At(i)->GetName() << "   " << trees.At(i)->GetTitle() <<endl;
            }
            cout<< "Enter tree number (<CR>=Quit): ";

            std::string input;
            std::getline(std::cin,input);

            if (input.size() == 0) return;   //-- user hit <CR>

            std::stringstream ss;
            ss << input;

            int number;
            ss >> number;
            if (ss.eof())
            {
                //-- user entered a number

                if (number >= 0 && number < trees.GetEntries()) {
                    tree = (TTree*) trees.At(number);
                }
                else {
                    cout<< "Tree number is out of range" <<endl;
                    return;
                }
            }
            else {
                cout<< "Input error" <<endl;
                return;
            }
        }
    }
    tree->SetMarkerColor(color);
    tree->SetLineColor(color);
}

void gcol(Int_t color, TGraph* gr) {
    if (!gr) gr = gtemp();
    if (!gr) return;
    gr->SetMarkerColor(color);
    gr->SetLineColor(color);
}

void hcol(Int_t color, TH1* h, Int_t fill) {
    if (!h) h = htemp();
    if (!h) return;
    if (fill >= 0) h->SetFillStyle(fill);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->SetFillColor(color);
}

void time() {gROOT->Time();}
void notime() {gROOT->Time(0);}

bool fexists(const char* ifname) {
    bool fexists = false;
    std::ifstream ifile(ifname, std::ifstream::in|std::ifstream::binary);
    if (ifile.is_open()) {
        ifile.close();
        fexists = true;
    }
    return fexists;
}

const char* substitute(const std::string& str0, const std::string what, const std::string by)
{
    if (what == by) return str0.c_str();
    std::string str = str0;
    for (std::string::size_type loc = str.find(what,0); loc != std::string::npos; ) {
        str.replace(loc, what.length(), by);
        loc = str.find(what,0);
    }
    //cout<< "str0 = " << str0 << " str = " << str <<endl;
    return Form("%s",str.c_str());       // put in Form circular buffer
}
int nlen(const std::string str)
{       // string length without leading spaces
    unsigned i = 0;
    for (; i<str.length(); ++i) {
        if (isspace(str[i])) continue;
        else break;
    }
    return strlen(str.c_str()) - i;
}

void showbra(TTree* tree, const char* pattern)
{
    // prints tree's list of branches
    TIter next(tree->GetListOfBranches());
    TBranch* b;
    while ((b = (TBranch*) next())) {
        if (pattern && !strstr(b->GetName(), pattern)) continue;
        cout<< b->GetName() << "   address = " << (TBranch*) b->GetAddress() <<endl;
    }
}

// prints tree's list of branches
void bra(TTree* tree, const char* pattern)
{
    if (!tree) {
        cout<< "Usage: bra(TTree* tree, const char* pattern=\"\")" <<endl;
        return;
    }

    TIter next(tree->GetListOfBranches());
    TBranch* b;
    cout<< "----- TTree " << tree->GetName() << " -----" <<endl;
    while ((b = (TBranch*) next())) {
        if (pattern && !strstr(b->GetName(), pattern)) continue;
        cout<< b->GetName() <<endl;
    }
}

void history(Int_t nlines)
{
    gSystem->Exec(Form("tail -%d root_hist",nlines));
}
void history(const char* pattern, Int_t nlines)
{
    if (nlines > 0) gSystem->Exec(Form("fgrep %s root_hist | tail -%d",pattern,nlines));
    else gSystem->Exec(Form("fgrep %s root_hist",pattern));
}
void hi(Int_t nlines) {history(nlines);}
void hi(const char* pattern, Int_t nlines) {history(pattern,nlines);}

void hiDraw(Int_t nlines)
{
    gSystem->Exec(Form("fgrep Draw root_hist | tail -%d",nlines));
}
void hiDraw(const char* pattern, Int_t nlines)
{
    if (nlines > 0) gSystem->Exec(Form("fgrep Draw root_hist | fgrep %s | tail -%d",pattern,nlines));
    else gSystem->Exec(Form("fgrep %s root_hist",pattern));
}

void ref() {
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) return;
    gPad->Modified();
    gPad->Update();
}

void log() {
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) return;
    gPad->SetLogy(!gPad->GetLogy());
    gPad->Update();
}

void logy() {log();}

void logz() {
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) return;
    gPad->SetLogz(!gPad->GetLogz());
    gPad->Update();
}

void logx() {
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) return;
    gPad->SetLogx(!gPad->GetLogx());
    gPad->Update();
}

void statover() {
    //-- used by lq
    // underflows and overflows are used by the Fill functions in the computation of statistics (mean value, RMS). 
    // By default, underflows or overflows are not used.
    TH1::StatOverflows(kTRUE);
}
void nostatover() {
    //-- default in ROOT
    // underflows and overflows are used by the Fill functions in the computation of statistics (mean value, RMS). 
    // By default, underflows or overflows are not used.
    TH1::StatOverflows(kFALSE);
}

void gri() {
    gStyle->SetPadGridY(1);
    gStyle->SetPadGridX(1);
}
void nogri() {
    gStyle->SetPadGridY(0);
    gStyle->SetPadGridX(0);
}

Int_t nbi(Int_t nx, Int_t ny)
{
    // There are two different sets of bin numbers: for 1D and for 2D historgrams
    if (nx == 0 && ny == 0) cout<< "1D histo: nx = " << gEnv->GetValue("Hist.Binning.1D.x", 0) << "\t 2D histo: nx = " << gEnv->GetValue("Hist.Binning.2D.x", 0) << " ny = " << gEnv->GetValue("Hist.Binning.2D.y", 0) <<endl;
    if (ny > 0) {
        // set binning for 2D histo only
        gEnv->SetValue("Hist.Binning.2D.x", nx);
        gEnv->SetValue("Hist.Binning.2D.y", ny);
        return gEnv->GetValue("Hist.Binning.2D.x", 0);
    }
    if (nx > 0) gEnv->SetValue("Hist.Binning.1D.x", nx);
    return gEnv->GetValue("Hist.Binning.1D.x", 0);
}
void prompt(const char* prompt) {((TRint*) gROOT->GetApplication())->SetPrompt(prompt);}
void rex() {gPad->RedrawAxis("g");}

Int_t axisNdivisions(Int_t Ndivisions) {
    TH1* h = htemp();
    if (h) {
        // if (Ndivisions == 0) cout<< "h->GetName() = " << h->GetName() << " Ndivisions = " << h->GetXaxis()->GetNdivisions() <<endl;
        if (Ndivisions == 0) Ndivisions = h->GetXaxis()->GetNdivisions();
        else h->GetXaxis()->SetNdivisions(Ndivisions);
        return Ndivisions;
    }
    TGraph* gr = gtemp();
    if (gr) {
        // if (Ndivisions == 0) cout<< "gr->GetName() = " << h->GetName() << " Ndivisions = " << gr->GetXaxis()->GetNdivisions() <<endl;
        if (Ndivisions == 0) Ndivisions = gr->GetXaxis()->GetNdivisions();
        else gr->GetXaxis()->SetNdivisions(Ndivisions);
        return Ndivisions;
    }
    TF1* f = ftemp();
    if (f) {
        // if (Ndivisions == 0) cout<< "f->GetName() = " << f->GetName() << " Ndivisions = " << f->GetXaxis()->GetNdivisions() <<endl;
        if (Ndivisions == 0) Ndivisions = f->GetXaxis()->GetNdivisions();
        else f->GetXaxis()->SetNdivisions(Ndivisions);
        return Ndivisions;
    }
    return 0;
}
Int_t ndi(Int_t npad) {
    TPad* current = (TPad*) gPad;
    TPad* subpad = (TPad*) gPad->GetMother()->GetPad(npad);
    Int_t Ndivisions = 0;
    if (subpad) {
        gPad = subpad;
        Ndivisions = axisNdivisions(0);
        if (Ndivisions == 0) {
            cout<< "subpad number is out of range" <<endl;
            return Ndivisions;
        }
        axisNdivisions(--Ndivisions);
        gPad = current;
    }
    else cout<< "subpad number is out of range" <<endl;
    return Ndivisions;
}

void clone() {gPad->DrawClone();}                  // clone current canvas

void zoomx(Axis_t xmin, Axis_t xmax) {
    TAxis *xaxis = 0;
    if (htemp()) xaxis = htemp()->GetXaxis();
    else if (gtemp()) xaxis = gtemp()->GetXaxis();
    else if (ftemp()) xaxis = ftemp()->GetXaxis();
    if (xaxis) {
        if (xmax == 0 && xmax < xmin) xmax = xaxis->GetXmax();
        xaxis->SetRangeUser(xmin,xmax);
        gPad->Modified();
        gPad->Update();
    }
}

void zoomy(Axis_t ymin, Axis_t ymax) {
    TAxis *yaxis = 0;
    if (htemp()) yaxis = htemp()->GetYaxis();
    else if (gtemp()) yaxis = gtemp()->GetYaxis();
    else if (ftemp()) yaxis = ftemp()->GetYaxis();
    if (yaxis) {
        if (ymax == 0 && ymax < ymin) ymax = yaxis->GetXmax();
        yaxis->SetRangeUser(ymin,ymax);
        gPad->Modified();
        gPad->Update();
    }
}

void zoom(Axis_t xmin, Axis_t xmax, Axis_t ymin, Axis_t ymax) {
    TAxis *xaxis = 0;
    TAxis *yaxis = 0;
    if (htemp()) {
        xaxis = htemp()->GetXaxis();
        yaxis = htemp()->GetYaxis();
    }
    else if (gtemp()) {
        xaxis = gtemp()->GetXaxis();
        yaxis = gtemp()->GetYaxis();
    }
    else if (ftemp()) {
        xaxis = ftemp()->GetXaxis();
        yaxis = ftemp()->GetYaxis();
    }
    if (xmax > xmin && xaxis) {
        if (xmax == 0 && xmax < xmin) xmax = xaxis->GetXmax();
        xaxis->SetRangeUser(xmin,xmax);
        gPad->Modified();
        gPad->Update();
    }
    if (ymax > ymin && yaxis) {
        if (ymax == 0 && ymax < ymin) ymax = yaxis->GetXmax();
        yaxis->SetRangeUser(ymin,ymax);
        gPad->Modified();
        gPad->Update();
    }
}

void zoomall(Axis_t xmin, Axis_t xmax, Axis_t ymin, Axis_t ymax)
{
    TPad* pad = gpad(0, false);
    Int_t npads = countpads(gPad);
    for (int ipad=0; ipad<npads; ++ipad) {
	pad->cd(ipad+1);
	zoom(xmin, xmax, ymin, ymax);
    }
    pad->cd();
}

void unzoom() {                                 // same as unzoomx()
    TAxis *xaxis = 0;
    if (htemp()) xaxis = htemp()->GetXaxis();
    else if (gtemp()) xaxis = gtemp()->GetXaxis();
    else if (ftemp()) xaxis = ftemp()->GetXaxis();
    if (xaxis) {
        xaxis->UnZoom();
        gPad->Modified();
        gPad->Update();
    }
}
void unzoomx() {
    TAxis *xaxis = 0;
    if (htemp()) xaxis = htemp()->GetXaxis();
    else if (gtemp()) xaxis = gtemp()->GetXaxis();
    else if (ftemp()) xaxis = ftemp()->GetXaxis();
    if (xaxis) {
        xaxis->UnZoom();
        gPad->Modified();
        gPad->Update();
    }
}
void unzoomy() {
    TAxis *yaxis = 0;
    if (htemp()) yaxis = htemp()->GetYaxis();
    else if (gtemp()) yaxis = gtemp()->GetYaxis();
    else if (ftemp()) yaxis = ftemp()->GetYaxis();
    if (yaxis) {
        yaxis->UnZoom();
        gPad->Modified();
        gPad->Update();
    }
}
void unzoomxy() {
    TAxis *xaxis = 0;
    TAxis *yaxis = 0;
    if (htemp()) {
        xaxis = htemp()->GetXaxis();
        yaxis = htemp()->GetYaxis();
    }
    else if (gtemp()) {
        xaxis = gtemp()->GetXaxis();
        yaxis = gtemp()->GetYaxis();
    }
    else if (ftemp()) {
        xaxis = ftemp()->GetXaxis();
        yaxis = ftemp()->GetYaxis();
    }
    if (xaxis && yaxis) {
        xaxis->UnZoom();
        yaxis->UnZoom();
        gPad->Modified();
        gPad->Update();
    }
}

void minmax(Axis_t ymin, Axis_t ymax) {
    if (htemp()) {
        htemp()->SetMinimum(ymin);
        htemp()->SetMaximum(ymax);
    }
    else if (gtemp()) {
        gtemp()->SetMinimum(ymin);
        gtemp()->SetMaximum(ymax);
    }
    else if (ftemp()) {
        ftemp()->SetMinimum(ymin);
        ftemp()->SetMaximum(ymax);
    }
    gPad->Modified();
    gPad->Update();
}

void fontsdrawtext(double x, double y, int f, const char *s)
{
    TLatex* t = new TLatex(x,y,Form("#font[41]{%d :} %s",f,s));
    t->SetTextFont(f);
    t->SetTextAlign(12);
    t->SetTextSize(0.048);
    t->Draw();
}
TCanvas* fonts()
{
    TCanvas* Tf = new TCanvas("Tf", "Tf",0,0,800,400);
    Tf->Range(0,0,1,1);
    Tf->SetBorderSize(2);
    Tf->SetFrameFillColor(0);

    double y = 0.95;
    for (int f = 12; f<=152; f+=10) {
        if (f!=142) fontsdrawtext(0.02,y, f,"ABCDEFGHIJKLMNOPQRSTUVWXYZ abcdefghijklmnopqrstuvwxyz 0123456789 @#$");
        else fontsdrawtext(0.02,y, f,"ABCD efgh 01234 @#$");
        y -= 0.065;
    }
    return Tf;
}

void colors()
{
    TCanvas *colors = new TCanvas("colors","Fill Area colors",0,0,500,200);
    colors->DrawColorTable();
}

Int_t Color(Int_t number)
{
    if (number == 0) return 0;
    if (number == 1) return 602;
    if (number == 2) return 2;
    if (number == 3) return 30;
    if (number == 4) return 4;
    if (number == 5) return 41;
    if (number == 6) return 6;
    if (number == 7) return 38;
    if (number == 8) return 8;
    if (number == 9) return 9;
    if (number == 10) return 29;
    if (number == 11) return 12;
    if (number == 12) return 46;
    if (number == 13) return 3;
    if (number == 14) return 39;
    if (number == 15) return 5;
    if (number == 16) return 48;
    if (number == 17) return 36;
    if (number == 18) return 29;
    if (number == 19) return 33;
    if (number == 20) return 1;
    if (number == 21) return 37;
    if (number == 22) return 28;
    if (number == 23) return 32;
    if (number == 24) return 40;
    if (number == 25) return 44;
    if (number == 26) return 47;
    if (number == 27) return 7;
    if (number == 28) return 34;
    if (number == 29) return 49;
    if (number == 30) return 11;
    return number;
}

void fwlite() {
    gSystem->Load("libFWCoreFWLite.so");
    //-- AutoLibraryLoader::enable();
}

void tb(const char* fname) {
    // loads old (better than new) TBrowser
    if (fname) TFile::Open(fname);
    new TBrowser("TBrowser", "TBrowser", 600, 0, 670, 500);
}

Int_t canvas_max_number(const char* pattern)
{
    TSeqCollection* list = gROOT->GetListOfCanvases();

    Int_t number_max = 0;
    for (int i=0; i<list->GetEntries(); ++i)
    {
        Int_t number = 0;
        std::string format_string = std::string(pattern) + std::string("%d");
        std::sscanf(list->At(i)->GetName(), format_string.c_str(), &number);
        if (number > number_max) number_max = number;
    }
    return number_max;
}

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

const char* nextcan(const char* base) {
    const char* cname = base;
    // if canvas 'cname' already exists increaments it like can_1, can_2, etc.
    // returns pointer to circular buffer owned by ROOT function Form
    int icycle = 0;
    while (gROOT->GetListOfCanvases()->FindObject(cname)) {
        cname = Form("%s_%d", base, ++icycle);
    }
    //cout<< "cname = " << cname <<endl;
    return cname;
}

const char* defcan() {
    TList *lc = (TList*)gROOT->GetListOfCanvases();
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) return Form("c1");

    const char *defcanvas = gROOT->GetDefCanvasName();
    char *cdef = 0;
    if (lc->FindObject(defcanvas)) {
        // taken from my correction to TCanvas::TCanvas(Bool_t build)
        //-- original: cdef = StrDup(Form("%s_n%d",defcanvas,lc->GetSize()+1));
        //-- Andriy Zatserklyaniy, Jan 26, 2012
        Int_t n = lc->GetSize()+1;
        while (lc->FindObject(Form("%s_n%d",defcanvas,n))) n++;
        cdef = Form("%s_n%d",defcanvas,n);
    }
    return cdef;
}

TCanvas* newcan(Int_t width, Int_t height, const char* name, const char* title)
{
    TString name_str = name;
    // if (!name || strlen(name)==0) name_str = Form("c_%dx%d",width,height);
    //-- if (!name || !*name) name_str = Form("c_1");
    if (!name || !*name) name_str = Form("c");
    name_str = nextcan(name_str.Data());
    TString title_str = title;
    if (!title || !*title) title_str = name_str;
    if (width == 0) width = gStyle->GetCanvasDefW();
    if (height == 0) height = gStyle->GetCanvasDefH();
    TCanvas * can = new TCanvas(name_str.Data(), title_str.Data(), width, height);
    //-- can->SetWindowSize(width + (width - can->GetWw()), height + (height - can->GetWh()));
    return can;
}

TCanvas* can(Int_t width, Int_t height, const char* name, const char* title)
{
    return newcan(width,height,name,title);
}

TString NextName(TString base)
{
    // if name 'base' already exists increaments it like h_1, h_2, etc.
    // returns TString which holds buffer with name
    TString out = base;
    Int_t i = 0;
    bool found = gDirectory->Get(out.Data());
    while (found) {
        std::stringstream ss;
        ss.str("");
        ss << base << "_" << ++i;
        //cout<< "         current hname: " << out <<endl;
        out = ss.str().c_str();
        TObject* obj = gDirectory->Get(out.Data());
        //cout<< "obj = " << obj <<endl;
        //found = gDirectory->Get(out.Data());
        found = (obj != 0);
    }
    //cout<< "new hname: " << out <<endl;
    return out;
}

Double_t qsum(Double_t e1, Double_t e2) {return TMath::Sqrt(e1*e1 + e2*e2);}

void sqroot(double a, double b, double c, double& x1, double& x2, bool verb)
{
    const double eps = 1e-7;
    x1 = x2 = 0;

    double D = b*b - 4.*a*c;
    if (D < -eps) {
        if (verb) cout<< "negative Discriminant = " << D <<endl;
        return;
    }
    else {
        if (D < 0)     // but D > -eps
        {
            x1 = x2 = -b / (2.*a);
            if (verb) cout<< "x1 = x2 = " << x1 << "   D = " << D <<endl;
            return;
        }

        double bsign = b >= 0? 1.: -1;
        x1 = (-b - bsign*TMath::Sqrt(D)) / (2.*a);
        x2 = c / (a*x1);                                      // use Vieta's theorem
        if (verb) cout<< "x1 = " << x1 << " x2 = " << x2 << "   D = " << D <<endl;
        return;
    }
}
void sqroot(double a, double b, double c)
{
    double x1 = 0, x2 = 0;
    sqroot(a,b,c,x1,x2,true);
}

//
// movestat with leftgaus and rightgaus
//
TPaveStats* movestat(Double_t x1NDC, Double_t y1NDC, Double_t xlen, Double_t ylen)
{
    //--AZ obsolete: consider moveleft and moveright

    // Good settings to move box at the left upper corner is movestat(0.15,0.5, 0.3,0.4)
    // i.e. y1NDC + ylen = 0.9

    if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
        cout<< "No canvas found" <<endl;
        return 0;
    }
    TPaveStats* stat = 0;

    gPad->Modified();    // make sure that the stat box appeared
    gPad->Update();

    if (htemp()) {
        stat = (TPaveStats*) htemp()->GetListOfFunctions()->FindObject("stats");
    }
    else if (gtemp()) {
        stat = (TPaveStats*) gtemp()->GetListOfFunctions()->FindObject("stats");
    }

    if (!stat) {
        cout<< "movestat: No stat box found. Try \ngPad->Modified();\ngPad->Update();" <<endl;
        return 0;
    }

    if (x1NDC == 0 && y1NDC == 0 && xlen == 0 && ylen == 0) {
        cout<< "Current: stat->GetX1NDC() = " << stat->GetX1NDC() << " stat->GetY1NDC() = " << stat->GetY1NDC() << " xlen " << stat->GetX2NDC() - stat->GetX1NDC() << " ylen " << stat->GetY2NDC() - stat->GetY1NDC() <<endl;
        return stat;
    }

    if (x1NDC == 0) x1NDC = stat->GetX1NDC();
    if (y1NDC == 0) y1NDC = stat->GetY1NDC();
    if (xlen == 0) xlen = stat->GetX2NDC() - stat->GetX1NDC();
    if (ylen == 0) ylen = stat->GetY2NDC() - stat->GetY1NDC();

    stat->SetX1NDC(x1NDC);
    stat->SetX2NDC(x1NDC + xlen);
    stat->SetY1NDC(y1NDC);
    stat->SetY2NDC(y1NDC + ylen);

    gPad->Modified();
    gPad->GetMother()->Update();

    return stat;
}

TPaveStats* moveleft(Double_t xleft, TCanvas* can)
{
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) return 0;
    if (!can) can = (TCanvas*) gPad;

    gPad->Modified();    // make sure that the stat box appeared
    gPad->Update();

    TPaveStats* ps = (TPaveStats*) can->GetPrimitive("stats");
    if (!ps) {
        cout<< "Could not find stat box in " << can->GetName() <<endl;
        return 0;
    }

    Double_t xlen = ps->GetX2NDC() - ps->GetX1NDC();
    ps->SetX1NDC(xleft);
    ps->SetX2NDC(ps->GetX1NDC() + xlen);

    can->Modified();
    can->Update();
    return ps;
}

TPaveStats* moveright(Double_t xright, TCanvas* can)
{
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) return 0;
    if (!can) can = (TCanvas*) gPad;

    gPad->Modified();    // make sure that the stat box appeared
    gPad->Update();

    TPaveStats* ps = (TPaveStats*) can->GetPrimitive("stats");
    if (!ps) {
        cout<< "Could not find stat box in " << can->GetName() <<endl;
        return 0;
    }

    Double_t xlen = ps->GetX2NDC() - ps->GetX1NDC();
    ps->SetX2NDC(xright);
    ps->SetX1NDC(ps->GetX2NDC() - xlen);

    can->Modified();
    can->Update();
    return ps;
}

// void leftgaus()  { movestat(0.15,0.5, 0.35,0.4); }
// void rightgaus() { movestat(0.64,0.5, 0.35,0.4); }
// void left()      { movestat(0.15,0.5, 0.35,0.4); }
// void right()     { movestat(0.64,0.6, 0.35,0.4); }
void leftgaus()  { moveleft(0.13); }
void rightgaus() { moveright(0.98); }
void left()      { moveleft(0.13); }
void right()     { moveright(0.98); }

TPaveText* titsize(Double_t x1NDC, Double_t x2NDC, Double_t y1NDC, Double_t y2NDC)
{
    // default values seem to be
    // GetX1NDC = 0.15, GetX2NDC = 0.85, GetY1NDC = 0.934, GetY2NDC = 0.995

    // to create a stat box (if it does not exists)
    gPad->Modified();
    gPad->Update();

    if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
        cout<< "No canvas found" <<endl;
        return 0;
    }

    TPaveText *tit = (TPaveText*)gPad->GetPrimitive("title");

    tit->SetX1NDC(x1NDC);
    if (x2NDC > 0) tit->SetX2NDC(x2NDC);
    if (y1NDC > 0) tit->SetY1NDC(y1NDC);
    if (y2NDC > 0) tit->SetY2NDC(y2NDC);

    gPad->Modified();
    gPad->Update();

    return tit;
}

TPaveText* titmax(Double_t x1NDC, Double_t x2NDC, Double_t y1NDC, Double_t y2NDC) {
    return titsize(x1NDC,x2NDC, y1NDC,y2NDC);
}

TPaveText* titlemax(Double_t x1NDC, Double_t x2NDC, Double_t y1NDC, Double_t y2NDC) {
    return titsize(x1NDC,x2NDC, y1NDC,y2NDC);
}

void cd(const char* filename)
{
    if (gROOT->GetListOfFiles()->GetEntries() == 0 && (!filename || *filename==0)) {
        cout<< "No file open" <<endl;
        return;
    }

    TFile* file = (TFile*) gDirectory;

    if (filename && *filename)
    {
        file = (TFile*) gROOT->FindObject(filename);
        if (file) {
            file->cd();
            cout<< "--> " << "     " << gDirectory->GetName() <<endl;
        }
        else {
            if (fexists(filename)) {
                cout<< "Open file " << filename <<endl;
                file = TFile::Open(filename);
                if (!file) {
                    cout<< "Could not open file " << filename <<endl;
                    cout<< "Current dir is still " << gDirectory->GetName() <<endl;
                }
            }
            else {
                cout<< "No such ROOT file: " << filename <<endl;
                cout<< "Current dir is still " << gDirectory->GetName() <<endl;
            }
            //-- cout<< "--> " << "     " << gDirectory->GetName() <<endl;
            // print names of all opened files
            if (gROOT->GetListOfFiles()->GetEntries() == 0) return;
            TIter next(gROOT->GetListOfFiles());
            int i = 0;
            for (TObject* obj=0; (obj = next());) {
                //-- if (strcmp(obj->GetName(), gDirectory->GetName())==0) cout<< "--> " << std::setw(2) << i << " " << obj->GetName() <<endl;
                if ((TFile*) obj == file) cout<< "--> " << std::setw(2) << i << " " << obj->GetName() <<endl;
                else cout<< std::setw(6) << i << " " << obj->GetName() <<endl;
                ++i;
            }
        }
        return;
    }

    // print names of all opened files
    if (gROOT->GetListOfFiles()->GetEntries() == 0) return;
    TIter next(gROOT->GetListOfFiles());
    int i = 0;
    for (TObject* obj=0; (obj = next());) {
        //-- if (strcmp(obj->GetName(), gDirectory->GetName())==0) cout<< "--> " << std::setw(2) << i << " " << obj->GetName() <<endl;
        if ((TFile*) obj == file) cout<< "--> " << std::setw(2) << i << " " << obj->GetName() <<endl;
        else cout<< std::setw(6) << i << " " << obj->GetName() <<endl;
        ++i;
    }

    cout<< "Enter file/dir number (<CR>=Quit): ";
    std::string input;
    std::getline(std::cin,input);

    if (input.size() == 0) return;   //-- user hit <CR>

    std::stringstream ss;
    ss << input;

    int number;
    ss >> number;
    if (ss.eof())
    {
        //-- user entered a number

        if (number >= 0 && number < gROOT->GetListOfFiles()->GetEntries()) {
            file = (TFile*) gROOT->GetListOfFiles()->At(number);
            file->cd();
            cout<< "--> " << setw(3) << " " << gDirectory->GetName() <<endl;
        }
        else {
            cout<< "Out of range. Current dir is still " << gDirectory->GetName() <<endl;
        }
        return;     //-- User input was processes. Return.
    }

    //-- user entered a string

    // restore the stream
    ss.clear();                // clear EOF
    ss.seekg(0);               // back to the beginning of the stream

    std::string str;
    std::getline(ss, str);     // capable to process filenames with spaces, not the case with to ss >> str;

    // strip optional quotes around the file name: if the filename is just number, user can enter it with quotes to distinguish from integer
    if (str[0] == '"' && str[str.size()-1] == '"') {
        str.erase(0, 1);
        str.erase(str.size()-1);
    }

    file = (TFile*) gROOT->FindObject(str.c_str());
    if (file) {
        file->cd();
        cout<< "--> " << "     " << gDirectory->GetName() <<endl;
    }
    else {
        if (fexists(str.c_str())) {
            cout<< "Open file " << str.c_str() <<endl;
            file = TFile::Open(str.c_str());
            if (!file) {
                cout<< "Could not open file " << str.c_str() <<endl;
                cout<< "Current dir is still " << gDirectory->GetName() <<endl;
            }
        }
        else {
            cout<< "No such ROOT file: " << str.c_str() <<endl;
            cout<< "Current dir is still " << gDirectory->GetName() <<endl;
        }
        cout<< "--> " << "     " << gDirectory->GetName() <<endl;
    }
    return;
}

void cd(Int_t fnumber)
{
    if (gROOT->GetListOfFiles()->GetEntries() == 0) {
        cout<< "No file open" <<endl;
        return;
    }

    TFile* file = (TFile*) gDirectory;

    // print names of all opened files
    if (gROOT->GetListOfFiles()->GetEntries() == 0) return;
    TIter next(gROOT->GetListOfFiles());
    int i = 0;
    for (TObject* obj=0; (obj = next());) {
        //-- if (strcmp(obj->GetName(), gDirectory->GetName())==0) cout<< "--> " << std::setw(2) << i << " " << obj->GetName() <<endl;
        if ((TFile*) obj == file) cout<< "--> " << std::setw(2) << i << " " << obj->GetName() <<endl;
        else cout<< std::setw(6) << i << " " << obj->GetName() <<endl;
        ++i;
    }

    if (fnumber >= 0 && fnumber < gROOT->GetListOfFiles()->GetEntries()) {
        file = (TFile*) gROOT->GetListOfFiles()->At(fnumber);
        file->cd();
        cout<< "--> " << setw(3) << " " << gDirectory->GetName() <<endl;
    }
    else {
        cout<< "Out of range. Current dir is still " << gDirectory->GetName() <<endl;
    }
}

void resize(Int_t width, Int_t height, Bool_t internal)
{
    if (width == 0) width = gStyle->GetCanvasDefW();
    if (height == 0) height = gStyle->GetCanvasDefH();
    // if (gROOT->GetListOfCanvases()->GetEntries() > 0) static_cast<TCanvas*>(gROOT->GetListOfCanvases()->FindObject(gPad))->SetWindowSize(width,height);
    if (gROOT->GetListOfCanvases()->GetEntries() > 0) {
        TCanvas* can = static_cast<TCanvas*>(gROOT->GetListOfCanvases()->FindObject(gPad));
        if (internal) can->SetWindowSize(width + (width - can->GetWw()), height + (height - can->GetWh()));
        else can->SetWindowSize(width,height);
    }
}

void square(Int_t size, Bool_t internal)
{
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
        cout<< "No canvas found" <<endl;
        return;
    }
    if (size == 0) size = gPad->GetWw();
    resize(size, size, internal);
}

//------------- gtemp(), htemp() and ftemp() functions -------------------

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

TCanvas* gpad(TCanvas* can, bool verb)
{
    // returns pointer to the current canvas
    // can be used for e.g. gpad()->SetWindowSize(300, 300);      // NB: does not work: gPad->SetWindowSize(300, 300);
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) return 0;

    if (can && gROOT->GetListOfCanvases()->FindObject(can)) {
        gPad = can;
        if (verb) cout<< "gPad->GetName() = " << gPad->GetName() << "   gPad->GetTitle() = " << gPad->GetTitle() <<endl;
        return (TCanvas*) gPad;
    }

    if (verb) cout<< "gPad->GetName() = " << gPad->GetName() << "   gPad->GetTitle() = " << gPad->GetTitle() <<endl;
    return (TCanvas*) gPad;
}

void lslist(TList* list)
{
    for (int i=0; i<list->GetEntries(); ++i) {
        TObject* obj = list->At(i);
        cout<< i << "\t" << obj->GetName() << "\t " << obj->GetTitle() <<endl;
    }
}

TF1* lsfun(const char* fname, Int_t nfun)
{
    TList* list = 0;

    if (htemp() && htemp()->GetEntries() > 0) list = htemp()->GetListOfFunctions();
    // if (list) cout<< "htemp()->GetTitle(): list->GetEntries() = " << list->GetEntries() <<endl;
    if (!list || (list && list->GetEntries() == 0)) {
        if (gtemp()) {
            list = gtemp()->GetListOfFunctions();
            if (!list || (list && list->GetEntries() == 0)) return 0;
        }
    }

    if (!list) return 0;

    TList* flist = new TList;
    for (int i=0; i<list->GetEntries(); ++i) {
        TF1* f = (TF1*) list->At(i);
        if (f->IsA()->InheritsFrom(TF1::Class())) {
            if (fname && *fname) {
                if (strcmp(fname, f->GetName())==0) flist->Add(f);
            }
            else flist->Add(f);
        }
    }
    if (flist->GetEntries() == 0) {
        delete flist;
        return 0;
    }

    if (nfun < 0) nfun = flist->GetEntries() - 1;

    TF1* f = 0;
    if (nfun < flist->GetEntries()) f = (TF1*) flist->At(nfun);
    delete flist;
    return f;
}

TF1* lsfun(Int_t nfun, bool verb)
{
    TList* list = 0;

    if (htemp()) list = htemp()->GetListOfFunctions();
    if ((!list || list->GetEntries()==0) && gtemp()) {list = gtemp()->GetListOfFunctions();}
    if (!list) return 0;

    TList* flist = new TList;
    for (int i=0; i<list->GetEntries(); ++i) {
        TF1* f = (TF1*) list->At(i);
        if (f->IsA()->InheritsFrom(TF1::Class())) {
            if (verb) {
                cout<< flist->GetEntries() << " " << f << "   " << f->GetName() << " " << f->GetTitle();
                for (int ipar=0; ipar<f->GetNpar(); ++ipar) {
                    cout<< "   " << f->GetParName(ipar) << " = " << f->GetParameter(ipar);
                }
                cout<<endl;
            }
            flist->Add(f);
        }
    }
    if (flist->GetEntries() == 0) {
        delete flist;
        return 0;
    }

    if (nfun < 0) nfun = flist->GetEntries() - 1;

    TF1* f = 0;
    if (nfun < flist->GetEntries()) f = (TF1*) flist->At(nfun);
    delete flist;
    return f;
}

void lstemp(TCanvas* can)           // shows content of the current canvas
{
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
        cout<< "No canvas found" <<endl;
        return;
    }
    if (can == 0) can = (TCanvas*) gPad;

    TObjArray pads;
    TObjArray histos;
    TObjArray histos2d;
    TObjArray graphs;
    TObjArray functions;

    TIter next(can->GetListOfPrimitives());
    TObject* obj = 0;
    while ((obj = next())) {
        if (obj->IsA()->InheritsFrom(TPad::Class())) {
            pads.Add(obj);
            continue;
        }
        if (obj->IsA()->InheritsFrom(TH2::Class())) {
            histos2d.Add(obj);
            continue;
        }
        if (obj->IsA()->InheritsFrom(TH1::Class())) {
            histos.Add(obj);
            continue;
        }
        if (obj->IsA()->InheritsFrom(TGraph::Class())) {
            graphs.Add(obj);
            continue;
        }
        if (obj->IsA()->InheritsFrom(TF1::Class())) {
            functions.Add(obj);
            continue;
        }
    }

    cout<< "Contents of the canvas \"" << can->GetName() << "\" \"" << can->GetTitle() << "\"" <<endl;
    for (int i=0; i<pads.GetEntries(); ++i) cout<< "pad\t " << i << " " << pads[i]->GetName() << " " << pads[i]->GetTitle() <<endl;
    for (int i=0; i<histos.GetEntries(); ++i) cout<< "histo\t " << i << " " << histos[i]->GetName() << " " << histos[i]->GetTitle() <<endl;
    for (int i=0; i<histos2d.GetEntries(); ++i) cout<< "histo2d\t " << i << " " << histos2d[i]->GetName() << " " << histos2d[i]->GetTitle() <<endl;
    for (int i=0; i<graphs.GetEntries(); ++i) cout<< "graph\t " << i << " " << graphs[i]->GetName() << " " << graphs[i]->GetTitle() <<endl;
    for (int i=0; i<functions.GetEntries(); ++i) cout<< "function\t " << i << " " << functions[i]->GetName() << " " << functions[i]->GetTitle() <<endl;
}

// --------------------------------------------------------------------
//    gtemp. NB: Declaration of the gtemp should precede htemp2Clone
// --------------------------------------------------------------------

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

TGraph* gtemp(TCanvas* can, Bool_t updateNameTitle) {return gtemp(-1, can, updateNameTitle);}

TGraph* gtemp(Int_t index, Bool_t updateNameTitle) {return gtemp(index, 0, updateNameTitle);}

TGraph* gtempClone(Int_t index, const char* name, const char* title, TCanvas* can)
{
    // TGraph* gr = gtemp(index, can, kFALSE);   // not practical
    TGraph* gr = gtemp(index, can, kTRUE);
    if (gr) {
        gr = (TGraph*) gr->Clone();
        if (name && *name) gr->SetName(name);
        if (title && *title) gr->SetTitle(title);
    }
    return gr;
}

TGraph* gtempClone(const char* name, const char* title, TCanvas* can) {
    return gtempClone(-1, name, title, can);
}

TGraph* gtempClone(Int_t index, const char* name, const char* title) {
    return gtempClone(index, name, title, 0);
}

// -----------------------------------------------------------------------
//    htemp. NB: htemp2Clone does not depends on gtemp in this new version
// -----------------------------------------------------------------------

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

TH1* htempClone(const char* name, const char* title, TCanvas* can) {
    TH1* h = htemp(can);
    if (!h) return 0;

    h = (TH1*) h->Clone();
    if (h && name && *name) h->SetName(name);
    if (h && title && *title) h->SetTitle(title);
    return h;
}

TH1* htempClone(Int_t index, const char* name, const char* title) {
    TH1* h = htemp(index);
    if (!h) return 0;

    h = (TH1*) h->Clone();
    if (h && name && *name) h->SetName(name);
    if (h && title && *title) h->SetTitle(title);
    return h;
}

TH1* htempClone(Int_t index, const char* name, const char* title, TCanvas* can) {
    TH1* h = htemp(index,can);
    if (!h) return 0;

    h = (TH1*) h->Clone();
    if (h && name && *name) h->SetName(name);
    if (h && title && *title) h->SetTitle(title);
    return h;
}

TH2* htemp2(Int_t index, TCanvas* can)
{
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
        cout<< "No canvas found" <<endl;
        return 0;
    }
    if (can == 0) can = (TCanvas*) gPad;

    TObjArray histos2d;
    TIter next(can->GetListOfPrimitives());
    TObject* obj = 0;
    while ((obj = next())) {
        if (obj->IsA()->InheritsFrom(TH2::Class())) histos2d.Add(obj);
    }
    if (histos2d.GetEntries() == 0) return 0;

    if (index < 0) return (TH2*) histos2d.Last();

    if (index < histos2d.GetEntries()) return (TH2*) histos2d[index];
    else return 0;
}

TH2* htemp2(TCanvas* can) {return htemp2(-1, can);}

TH2* htemp2(Int_t index) {return htemp2(index, 0);}

TH2* htemp2Clone(const char* name, const char* title, TCanvas* can) {
    TH2* h2 = htemp2(can);
    if (!h2) return 0;

    h2 = (TH2*) h2->Clone();
    if (h2 && name && *name) h2->SetName(name);
    if (h2 && title && *title) h2->SetTitle(title);
    return h2;
}

TH2* htemp2Clone(Int_t index, const char* name, const char* title) {
    TH2* h2 = htemp2(index);
    if (!h2) return 0;

    h2 = (TH2*) h2->Clone();
    if (h2 && name && *name) h2->SetName(name);
    if (h2 && title && *title) h2->SetTitle(title);
    return h2;
}

TH2* htemp2Clone(Int_t index, const char* name, const char* title, TCanvas* can) {
    TH2* h2 = htemp2(index,can);
    if (!h2) return 0;

    h2 = (TH2*) h2->Clone();
    if (h2 && name && *name) h2->SetName(name);
    if (h2 && title && *title) h2->SetTitle(title);
    return h2;
}

TH2* htemp2Create(const char* name, const char* title)
{
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
        cout<< "No canvas found" <<endl;
        return 0;
    }
    // there are 2 possible pad which appears as 2-dim histogram:
    // 1) in case of TH2::Draw the pad contains 2-dim histogram indeed
    // 2) in case of TTree::Draw the pad contains empty 2-dim histogram and TGraph
    TH2* h2 = 0;
    TGraph* gr = 0;
    //cout<< "gPad->GetName() = " << gPad->GetName() << endl;
    TIter next(gPad->GetListOfPrimitives());
    TObject* obj;
    while ((obj = next())) {
        //cout<< "obj->GetName() = " << obj->GetName() <<endl;
        if (obj->IsA()->InheritsFrom(TH2::Class()) && !h2) {
            h2 = (TH2*) obj->Clone();      // clone the histo
        }
        if (obj->IsA()->InheritsFrom(TGraph::Class()) && !gr) {
            gr = (TGraph*) obj;
        }
    }
    if (!h2) return 0;
    if (h2 && !gr) {
        // looks like pure 2-dim histogram
        return h2;
    }
    if (h2 && gr) {
        // case of TTree::Draw
        // constract "normal" histogram filling up the h2
        // Note that we do not have information about the number of entries!
        for (int i=0; i<gr->GetN(); ++i) {
            h2->Fill(*(gr->GetX()+i), *(gr->GetY()+i));
        }
    }
    if (h2 && name && *name) h2->SetName(name);
    if (h2 && title && *title) h2->SetTitle(title);
    return h2;
}

Double_t wtemp() {
    TH1* h = htemp();
    if (h) return h->GetBinWidth(1);
    else return 0;
}

// --------------------------------------------
//    ftemp
// --------------------------------------------

TF1* ftemp(TCanvas* can, Int_t index)
{
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
        cout<< "No canvas found" <<endl;
        return 0;
    }
    if (can == 0) can = (TCanvas*) gPad;

    TObjArray functions;
    TIter next(can->GetListOfPrimitives());
    TObject* obj = 0;
    while ((obj = next())) {
        if (obj->IsA()->InheritsFrom(TF1::Class())) functions.Add(obj);
    }

    if (functions.GetEntries() > 0) {
        if (index < 0) return (TF1*) functions.Last();

        if (index < functions.GetEntries()) return (TF1*) functions[index];
        else return 0;
    }

    // no standalone functions in the canvas. Look at htemp() and gtemp() (in this order)
    return lsfun("", index);
}

TF1* ftemp(TCanvas* can) {return ftemp(can, -1);}

TF1* ftemp(Int_t index) {return ftemp(0, index);}

TMultiGraph* mgtemp(const char* name, const char* title)
{
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
        cout<< "No canvas found" <<endl;
        return 0;
    }
    TIter next(gPad->GetListOfPrimitives());
    TObject* obj = 0;
    while ((obj = next())) if (obj->IsA()->InheritsFrom(TMultiGraph::Class())) break; else obj=0;
    if (!obj) return 0;
    TMultiGraph* mg = (TMultiGraph*) obj;
    if (mg && name && *name) mg->SetName(name);
    if (mg && title && *title) mg->SetTitle(title);
    return mg;
}

TPaveStats* stemp()     // return pointer to the stat box
{
    TPaveStats* stats = 0;

    if (htemp()) stats = (TPaveStats*) htemp()->GetFunction("stats");
    if (!stats && gtemp()) stats = (TPaveStats*) gtemp()->GetFunction("stats");
    if (!stats) {
        // try to access the whole canvas
        if (gROOT->GetListOfCanvases()->GetEntries()) {
            stats = (TPaveStats*) gPad->GetPrimitive("stats");
        }
    }

    return stats;
}

TPaveStats* stemp(Int_t index)
{
    TPaveStats* stats = 0;

    if (htemp(index)) stats = (TPaveStats*) htemp(index)->GetFunction("stats");
    if (!stats && gtemp(index)) stats = (TPaveStats*) gtemp(index)->GetFunction("stats");
    if (!stats) {
        // try to access the whole canvas
        if (gROOT->GetListOfCanvases()->GetEntries()) {
            stats = (TPaveStats*) gPad->GetPrimitive("stats");
        }
    }

    return stats;
}

void statcol()
{
    TPaveStats* stats = 0;
    Int_t color = -1;

    if (htemp()) {
        stats = (TPaveStats*) htemp()->GetFunction("stats");
        if (stats) color = htemp()->GetLineColor();
    }
    if (!stats && gtemp()) {
        stats = (TPaveStats*) gtemp()->GetFunction("stats");
        if (stats) color = gtemp()->GetLineColor();
    }

    if (stats) {
        stats->SetLineColor(color);
        stats->SetTextColor(color);
        gPad->Modified();
        gPad->Update();
    }
}

void statcol(Int_t index)       // index of historgram or graph
{
    TPaveStats* stats = 0;
    Int_t color = -1;

    if (htemp(index)) {
        stats = (TPaveStats*) htemp(index)->GetFunction("stats");
        if (stats) color = htemp(index)->GetLineColor();
    }
    if (!stats && gtemp(index)) {
        stats = (TPaveStats*) gtemp(index)->GetFunction("stats");
        if (stats) color = gtemp(index)->GetLineColor();
    }

    if (stats) {
        stats->SetLineColor(color);
        stats->SetTextColor(color);
        gPad->Modified();
        gPad->Update();
    }
}

//--------------- end of gtemp, htemp, ftemp ------------------------

TGraphErrors* gsubf(const TGraph* gr, const TF1* f, const char* name, const char* title)
{
    if (!gr) gr = gtemp();
    if (!gr) {
        cout<< "No TGraph found" <<endl;
        return 0;
    }
    if (!f) f = ftemp();
    if (!f) {
        cout<< "No TF1 found" <<endl;
        return 0;
    }

    Double_t* x = gr->GetX();
    Double_t* y = gr->GetY();
    Double_t* ex = 0;
    Double_t* ey = 0;
    if (gr->InheritsFrom(TGraphErrors::Class())) {
        ex = gr->GetEX();
        ey = gr->GetEY();
    }
    TGraphErrors* grsub = new TGraphErrors(gr->GetN(),x,y,ex,ey);
    if (name && *name) grsub->SetName(name);
    else grsub->SetName(Form("%s-%s",gr->GetName(),f->GetName()));
    if (title && *title) grsub->SetTitle(title);
    else grsub->SetTitle(Form("%s - %s",gr->GetTitle(),f->GetName()));
    grsub->SetMarkerStyle(gr->GetMarkerStyle());
    grsub->SetMarkerColor(gr->GetMarkerColor());
    grsub->SetLineStyle(gr->GetLineStyle());
    grsub->SetLineColor(gr->GetLineColor());

    for (int i=0; i<grsub->GetN(); ++i) {
        Double_t val_old = grsub->GetY()[i];
        Double_t val_fun = f->Eval(grsub->GetX()[i]);
        Double_t val_new = val_old - val_fun;
        grsub->SetPoint(i, grsub->GetX()[i], val_new);
    }
    return grsub;
}

TH1* hsubf(const TH1* h, const TF1* f, const char* name, const char* title)
{
    if (!h) h = htemp();
    if (!h) {
        cout<< "No histogram found" <<endl;
        return 0;
    }
    if (!f) f = ftemp();
    if (!f) {
        cout<< "No TF1 found" <<endl;
        return 0;
    }

    TH1* hsub = (TH1*) h->Clone();
    if (name && *name) hsub->SetName(name);
    else hsub->SetName(Form("%s-%s",h->GetName(),f->GetName()));
    if (title && *title) hsub->SetTitle(title);
    else hsub->SetTitle(Form("%s - %s",h->GetTitle(),f->GetName()));
    hsub->SetMarkerStyle(h->GetMarkerStyle());
    hsub->SetMarkerColor(h->GetMarkerColor());
    hsub->SetLineStyle(h->GetLineStyle());
    hsub->SetLineColor(h->GetLineColor());

    for (int i=1; i<hsub->GetNbinsX(); ++i) {
        Double_t val_old = hsub->GetBinContent(i);
        Double_t val_fun = f->Eval(hsub->GetBinCenter(i));
        Double_t val_new = val_old - val_fun;
        hsub->SetBinContent(i,val_new);
    }
    return hsub;
}

void gsubg(TGraph* g, const TGraph* gsub)
{
    if (g->GetN() != gsub->GetN()) {
        cout<< "gsubg: different dimensions: g->GetN() = " << g->GetN() << " gsub->GetN() = " << gsub->GetN() <<endl;
        return;
    }

    for (int i=0; i<g->GetN(); ++i)
    {
        g->SetPoint(i, g->GetX()[i], g->GetY()[i] - gsub->GetY()[i]);
    }
}

void drawpan(TH1* h) {
    if (h == 0) h = htemp();
    if (h) {
        h->DrawPanel();
        return;
    }
    TGraph* gr = gtemp();
    if (gr) {
        gr->DrawPanel();
        return;
    }
}

void fitpan(TH1* h) {
    if (h == 0) h = htemp();
    if (h) {
        h->FitPanel();
        return;
    }
    TGraph* gr = gtemp();
    if (gr) {
        gr->FitPanel();
        return;
    }
}

void InitPredefinedFunctions()
{
    TF1 *f1;
    if (!gROOT->GetListOfFunctions()->FindObject("gaus"))    {f1 = new TF1("gaus","gaus",-1,1);       f1->SetParameters(1,0,1);}
    if (!gROOT->GetListOfFunctions()->FindObject("gausn"))   {f1 = new TF1("gausn","gausn",-1,1);     f1->SetParameters(1,0,1);}
    if (!gROOT->GetListOfFunctions()->FindObject("landau"))  {f1 = new TF1("landau","landau",-1,1);   f1->SetParameters(1,0,1);}
    if (!gROOT->GetListOfFunctions()->FindObject("landaun")) {f1 = new TF1("landaun","landaun",-1,1); f1->SetParameters(1,0,1);}
    if (!gROOT->GetListOfFunctions()->FindObject("expo"))    {f1 = new TF1("expo","expo",-1,1);       f1->SetParameters(1,1);}
    if (!gROOT->GetListOfFunctions()->FindObject("pol0")) {
        for (Int_t i=0;i<10;i++) {
            f1 = new TF1(Form("pol%d",i),Form("pol%d",i),-1,1);
            f1->SetParameters(1,1,1,1,1,1,1,1,1,1);
        }
    }
}

void fitgaus(Double_t xmin, Double_t xmax, const char* opt, const char* gopt)
{
    // histogram part
    TH1* h = htemp();
    if (h && h->GetEntries() == 0) h = 0;        // seems, this is just 2-dim frame for the plot
    if (h)
    {
        if (xmax == 0 && xmax < xmin) xmax = h->GetXaxis()->GetXmax();
        if (xmin == 0 && xmax == 0) h->Fit("gaus", opt, gopt);
        else {
            TF1* gaus = new TF1("gaus", "gaus", xmin, xmax);
            h->Fit(gaus, Form("%sR",opt), gopt);
        }
        h->SetStats(0); // the following line is needed to avoid that the automatic redrawing of stats ($ROOTSYS/tutorial/statsEditing.C)
        h->SetStats(1); //--AZ: trick to recalculate text size in the stat box
        gPad->Modified();
        gPad->Update();
        return;
    }

    // graph part
    TGraph* gr = gtemp();
    if (gr)
    {
        if (xmax == 0 && xmax < xmin) xmax = gr->GetXaxis()->GetXmax();
        if (xmin == 0 && xmax == 0) gr->Fit("gaus", opt, gopt);
        else {
            TF1* gaus = new TF1("gaus", "gaus", xmin, xmax);
            gr->Fit(gaus, Form("%sR",opt), gopt);
            //--no TGraph::SetStats: gr->SetStats(0); // the following line is needed to avoid that the automatic redrawing of stats
        }
        gPad->Modified();
        gPad->Update();
        return;
    }
} 
void fitg(Double_t xmin, Double_t xmax, const char* opt, const char* gopt) {
    fitgaus(xmin,xmax,opt,gopt);
}
void fitgl(Double_t xmin, Double_t xmax, const char* opt, const char* gopt) {
    fitgaus(xmin,xmax,opt,gopt);
    leftgaus();
}
void fitgr(Double_t xmin, Double_t xmax, const char* opt, const char* gopt) {
    fitgaus(xmin,xmax,opt,gopt);
    rightgaus();
}

void fitexpo(Double_t xmin, Double_t xmax, const char* opt, const char* gopt)
{
    // histogram part
    TH1* h = htemp();
    if (h && h->GetEntries() == 0) h = 0;        // seems, this is just 2-dim frame for the plot
    if (h)
    {
        if (xmax == 0 && xmax < xmin) xmax = h->GetXaxis()->GetXmax();
        if (xmin == 0 && xmax == 0) h->Fit("expo", opt, gopt);
        else {
            TF1* expo = new TF1("expo", "expo", xmin, xmax);
            h->Fit(expo, Form("%sR",opt), gopt);
        }
        gPad->Modified();
        gPad->Update();
        return;
    }

    // graph part
    TGraph* gr = gtemp();
    if (gr)
    {
        if (xmax == 0 && xmax < xmin) xmax = gr->GetXaxis()->GetXmax();
        if (xmin == 0 && xmax == 0) gr->Fit("expo", opt, gopt);
        else {
            TF1* expo = new TF1("expo", "expo", xmin, xmax);
            gr->Fit(expo, Form("%sR",opt), gopt);
        }
        gPad->Modified();
        gPad->Update();
        return;
    }
} 
void fite(Double_t xmin, Double_t xmax, const char* opt, const char* gopt) {
    fitexpo(xmin,xmax,opt,gopt);
}
void fitel(Double_t xmin, Double_t xmax, const char* opt, const char* gopt) {
    fitexpo(xmin,xmax,opt,gopt);
    leftgaus();
}
void fiter(Double_t xmin, Double_t xmax, const char* opt, const char* gopt) {
    fitexpo(xmin,xmax,opt,gopt);
    rightgaus();
}

void fitpol(Double_t xmin, Double_t xmax, Int_t power, const char* opt, const char* gopt)
{
    if (!gROOT->GetListOfFunctions()->FindObject("pol0")) InitPredefinedFunctions();

    TH1* h = htemp();
    if (h && h->GetEntries() == 0) h = 0;                // case of Tree::Draw for 2-dim histo
    if (h) {
        if (xmax == 0 && xmax < xmin) xmax = h->GetXaxis()->GetXmax();
        if (xmin == 0 && xmax == 0) h->Fit(Form("pol%d",power), opt, gopt);
        else {
            // TF1* pol = new TF1(Form("pol%d",power),Form("pol%d",power), xmin, xmax);
            // h->Fit(pol, Form("%sR",opt), gopt);
            h->Fit(Form("pol%d",power), opt, gopt, xmin, xmax);
        }
    }
    else {
        TGraph* gr = gtemp();
        if (gr) {
            if (xmax == 0 && xmax < xmin) xmax = gr->GetXaxis()->GetXmax();
            if (xmin == 0 && xmax == 0) gr->Fit(Form("pol%d",power), opt, gopt);
            else {
                // TF1* pol = new TF1(Form("pol%d",power),Form("pol%d",power), xmin, xmax);
                // gr->Fit(pol, Form("%sR",opt), gopt);
                gr->Fit(Form("pol%d",power), opt, gopt, xmin, xmax);
            }
        }
    }
    gPad->Modified();
    gPad->Update();
}
void fitp(Double_t xmin, Double_t xmax, Int_t power, const char* opt, const char* gopt) {
    fitpol(xmin,xmax,power,opt,gopt);
}

Double_t fit_parameter(const char* fun_name, const char* par_name)
{
    TF1* f = lsfun(fun_name);
    if (!f) return 0;

    if (par_name && *par_name) return f->GetParameter(par_name);
    else return 0;
}
Double_t fit_parameter(const char* fun_name, Int_t number)
{
    TF1* f = lsfun(fun_name);
    if (!f) return 0;

    if (number >= 0 && number < f->GetNpar()) return f->GetParameter(number);
    else return 0;
}
Double_t fit_parameter_error(const char* fun_name, Int_t number)
{
    TF1* f = lsfun(fun_name);
    if (!f) return 0;

    if (number >= 0 && number < f->GetNpar()) return f->GetParError(number);
    else return 0;
}
// gaus parameters
Double_t gaus_ampl() {return fit_parameter("gaus", "Constant");}
Double_t gaus_mean() {return fit_parameter("gaus", "Mean");}
Double_t gaus_sigma() {return fit_parameter("gaus", "Sigma");}
// gaus parameter errors
Double_t gaus_eampl() {return fit_parameter_error("gaus", 0);}
Double_t gaus_emean() {return fit_parameter_error("gaus", 1);}
Double_t gaus_esigma() {return fit_parameter_error("gaus", 2);}

Double_t gaus_area()       // divide to wtemp() to get Nevents in histogram
{
    Double_t Nevents = 0;

    TF1* fgaus = ftemp();
    if (!fgaus) return 0;
    if (strcmp(fgaus->GetName(),"gaus")==0) {
        Double_t constant = fgaus->GetParameter("Constant");
        Double_t sigma = fgaus->GetParameter("Sigma");
        Nevents = TMath::Sqrt(2.*TMath::Pi()) * constant * sigma;
        printf("Nevents = %f\n", Nevents);
    }
    else cout<< "Function ftemp() does not have name \"gaus\"" <<endl;
    return Nevents;
}

Double_t ngaus(TH1* h)
{
    // the number of events in gaus fit to histogram
    if (!h) h = htemp();
    if (!h) return 0;

    TObjArray functions;
    TF1* f = 0;
    TObject* obj = 0;
    TIter next(h->GetListOfFunctions());
    while ((obj = next())) {
        if (obj->IsA()->InheritsFrom(TF1::Class())) {
            f = dynamic_cast<TF1*>(obj);
            if (f && strcmp(f->GetName(),"gaus")==0) functions.Add(f);
        }
    }
    if (functions.GetEntries() > 0) f = (TF1*) functions.Last();
    if (!f) return 0;

    // access histogram to obtain a bin width (assumed const)
    Double_t w = h->GetBinWidth(1);

    Double_t A = f->GetParameter(0);
    //Double_t mean = f->GetParameter(1);
    Double_t sigma = f->GetParameter(2);

    Double_t sum = TMath::Sqrt(2.*TMath::Pi()) * A * sigma / w;

    cout<< "The number of events in the gaussian fit is " << sum <<endl;

    return sum;
}
// Double_t ngaus(TH1* h)
// {
//    if (!h) h = htemp();
//    if (!h) return 0;
//    TF1* gaus = h->GetFunction("gaus");
//    Double_t w = h->GetBinWidth(1);
//    Double_t nevents = 0;
//    if (gaus) {
//       Double_t anorm = gaus->GetParameter("Constant");
//       Double_t sigma = gaus->GetParameter("Sigma");
//       nevents = TMath::Sqrt(2.*TMath::Pi()) * anorm * sigma / w;
//    }
//    return nevents;
// }

void pargaus(Double_t& a, Double_t& mean, Double_t& sigma)
{
    a = mean = sigma = 0;

    TH1* h = htemp();
    if (h) {
        TObjArray functions;
        TF1* gaus = 0;
        TObject* obj = 0;
        TIter next(h->GetListOfFunctions());
        while ((obj = next())) {
            if (obj->IsA()->InheritsFrom(TF1::Class())) {
                gaus = dynamic_cast<TF1*>(obj);
                if (gaus && strcmp(gaus->GetName(),"gaus")==0) functions.Add(gaus);
            }
        }
        if (functions.GetEntries() > 0) gaus = (TF1*) functions.Last();
        if (gaus) {
            a = gaus->GetParameter(1);
            mean = gaus->GetParameter(1);
            sigma = gaus->GetParameter(2);
        }
    }
}

void npe(TH1* h, Double_t* p_Npe, Double_t* p_eNpe)
{
    if (p_Npe) *p_Npe = 0;
    if (p_eNpe) *p_eNpe = 0;

    if (!h) h = htemp();
    if (h) {
        TObjArray functions;
        TF1* gaus = 0;
        TObject* obj = 0;
        TIter next(h->GetListOfFunctions());
        while ((obj = next())) {
            if (obj->IsA()->InheritsFrom(TF1::Class())) {
                gaus = dynamic_cast<TF1*>(obj);
                if (gaus && strcmp(gaus->GetName(),"gaus")==0) functions.Add(gaus);
            }
        }
        if (functions.GetEntries() > 0) gaus = (TF1*) functions.Last();
        if (gaus) {
            Double_t mean = gaus->GetParameter(1);
            Double_t sigma = gaus->GetParameter(2);
            Double_t emean = gaus->GetParError(1);
            Double_t esigma = gaus->GetParError(2);
            Double_t ratio = mean/sigma;
            Double_t Npe = ratio*ratio;
            Double_t ra = emean/mean;
            Double_t rb = esigma/sigma;
            Double_t eNpe = 2.*Npe*TMath::Sqrt(ra*ra + rb*rb);
            // cout<< "Npe = " << Npe << " +- " << eNpe << "\t mean = " << mean << " sigma = " << sigma <<endl;
            // cout<< setprecision(3) << "Npe = " << Npe << " +- " << eNpe << "\t mean = " << mean << " sigma = " << sigma <<endl;
            printf("Npe = %0.3f +- %0.3f mean = %0.3f sigma = %0.3f\n", Npe,eNpe, mean,sigma);
            if (p_Npe) *p_Npe = Npe;
            if (p_eNpe) *p_eNpe = eNpe;
        }
    }
}

void hnpe(TH1* h, Double_t* p_Npe, Double_t* p_eNpe)
{
    if (!h) h = htemp();
    if (!h) {
        cout<< "hnpe: h = 0" <<endl;
        return;
    }

    if (p_Npe) *p_Npe = 0;
    if (p_eNpe) *p_eNpe = 0;

    const Double_t range_est = 3.;  // sigma
    const Double_t range = 2.;      // sigma

    Int_t axisBinFirst = h->GetXaxis()->GetFirst();
    Int_t axisBinLast = h->GetXaxis()->GetLast();

    Double_t x1 = h->GetMean() - range_est*h->GetRMS();
    Double_t x2 = h->GetMean() + range_est*h->GetRMS();

    h->SetAxisRange(x1, x2);

    Double_t mean = h->GetMean();
    Double_t rms = h->GetRMS();

    // subrange +- 1.5 sigma

    x1 = mean - range*rms;
    x2 = h->GetMean() + range*h->GetRMS();

    cout<< "subrange x1,x2 = " << x1 << ", " << x2 <<endl;

    h->SetAxisRange(x1, x2);

    mean = h->GetMean();
    rms = h->GetRMS();

    // restore axis range
    h->GetXaxis()->SetRange(axisBinFirst, axisBinLast);

    Double_t Npe = TMath::Power(mean/rms, 2);

    Double_t emean = h->GetMeanError();
    Double_t erms = h->GetRMSError();

    Double_t ra = emean/mean;
    Double_t rb = erms/rms;
    Double_t eNpe = 2.*Npe*TMath::Sqrt(ra*ra + rb*rb);

    printf("Histogram Npe = %0.3f +- %0.3f mean = %0.3f RMS = %0.3f\n", Npe,eNpe, mean,rms);
    if (p_Npe) *p_Npe = Npe;
    if (p_eNpe) *p_eNpe = eNpe;

    h->Fit("gaus", "", "", x1,x2);
    npe();
}

void pop(TPad* pad)
{
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
        cout<< "No canvas found" <<endl;
        return;
    }

    if (pad) gPad = pad;
    else cout<< "Could not find TPad " << pad <<endl;
}

void pic(TString pathname, TString suffix)
{
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
        cout<< "No canvas found" <<endl;
        return;
    }
    if (pathname == "") {
        if (htemp()) pathname = htemp()->GetName();
        else if (gtemp()) pathname = gtemp()->GetName();
        else if (ftemp()) pathname = ftemp()->GetName();
        else pathname = gPad->GetName();
    }
    gPad->SaveAs(pathname + suffix + ".eps");
    gPad->SaveAs(pathname + suffix + ".png");
}
// just eps
void eps(TString pathname, TString suffix)
{
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
        cout<< "No canvas found" <<endl;
        return;
    }
    if (pathname == "") {
        if (htemp()) pathname = htemp()->GetName();
        else if (gtemp()) pathname = gtemp()->GetName();
        else if (ftemp()) pathname = ftemp()->GetName();
        else pathname = gPad->GetName();
    }
    gPad->SaveAs(pathname + suffix + ".eps");
}
// just png
void png(TString pathname, TString suffix, bool check_fexists)
{
    if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
        cout<< "No canvas found" <<endl;
        return;
    }
    if (pathname == "") {
        if (htemp()) pathname = htemp()->GetName();
        else if (gtemp()) pathname = gtemp()->GetName();
        else if (ftemp()) pathname = ftemp()->GetName();
        else pathname = gPad->GetName();
    }
    TString ofname = pathname + suffix + ".png";
    if (check_fexists && fexists(ofname.Data())) {
        cout<< gPad->GetName() << " --> " << "File exists: \"" << ofname << "\"" <<endl;
        cout<< "Enter <CR> to overwrite, any other key to quit: ";
        std::string str;
        std::getline(std::cin, str);    // read complete line to avoid problems with the next call
        if (str.size() > 0) {
            cout<< "*** quit without saving" <<endl;
            return;
        }
    }
    cout<< gPad->GetName() << " --> " <<std::flush;
    gPad->SaveAs(pathname + suffix + ".png");
}
void fpng(TString pathname, TString suffix, bool check_fexists)
{
    TString fname = gDirectory->GetName();
    if (fname == "Rint") fname = "";
    pathname = fname + " " + pathname;
    png(pathname, suffix, check_fexists);
}
void pngall() {
    for (int i=0; i<gROOT->GetListOfCanvases()->GetEntries(); i++) {
        TCanvas* can=(TCanvas*)gROOT->GetListOfCanvases()->At(i);
        can->SaveAs(Form("%s.png",can->GetName()));
    }
}

Double_t intgaus(TH1* h)
{
    if (h == 0) h = htemp();

    Double_t gampl = h->GetFunction("gaus")->GetParameter(0);
    // Double_t gmean = h->GetFunction("gaus")->GetParameter(1);
    Double_t gsigma = h->GetFunction("gaus")->GetParameter(2);
    Double_t nevt_w = TMath::Sqrt(2.*TMath::Pi()) * gampl*gsigma / h->GetBinWidth(1);
    cout<< "TMath::Sqrt(2.*TMath::Pi()) * gampl*gsigma / h->GetBinWidth(1) = " << nevt_w <<endl;

    return nevt_w;
}

void print_htemp(TH1* h)
{
    if (!h) h = htemp();
    if (!h) return;

    cout<< "h->GetEntries() = " << h->GetEntries() <<endl;
    cout<< "sum of bin contents: h->Integral() = " << h->Integral() <<endl;

    Double_t integral = 0;
    Double_t integral_count = 0;
    for (int i=1; i<=h->GetNbinsX(); ++i) {
        Double_t con = h->GetBinContent(i);
        integral += con*h->GetBinCenter(i);
        integral_count += con*(i-1);
        if (con) {
            cout<< i << "\t " << h->GetBinCenter(i) << "\t " << con <<endl;
        }
    }
    cout<< "integer count integral con*(i-1) = " << integral_count <<endl;
    cout<< "integral con*h->GetBinCenter(i)  = " << integral <<endl;
}
// void hpri(TH1* h=0)
void hprint(TH1* h)
{
    if (!h) h = htemp();
    if (!h) return;

    cout<< "h->GetEntries() = " << h->GetEntries() <<endl;
    //cout<< "sum of bin contents: h->Integral() = " << h->Integral() <<endl;

    Double_t integral = 0;
    //Double_t integral_count = 0;
    for (int i=1; i<=h->GetNbinsX(); ++i) {
        Double_t con = h->GetBinContent(i);
        integral += con*h->GetBinLowEdge(i);
        //integral_count += con*(i-1);
        if (con) {
            // cout<< i << "\t " << h->GetBinCenter(i) << "\t " << con <<endl;
            cout<< i << "\t " << h->GetBinLowEdge(i) << "\t " << con <<endl;
        }
    }
    //cout<< "integer count integral con*(i-1) = " << integral_count <<endl;
    cout<< "integral con*h->GetBinLowEdge(i)     = " << integral <<endl;
}

//void tax(const char* format="%d-%H")
void tax(const char* format)
{
    // special option format = 0 (NULL) to reset(remove) time axis
    /*
    Defines the format of the labels along the time axis. 
    It can be changed using the TAxis method SetTimeFormat. 
    The time format is the one used by the C function strftime(). 
    It's a string containing the following formatting characters :
     */
    if (strcmp(format, "help")==0)
    {
        cout<< "   * for date :                                                                                   " <<endl;
        cout<< "                    o %a abbreviated weekday name         " <<endl;
        cout<< "                    o %b abbreviated month name           " <<endl;
        cout<< "                    o %d day of the month (01-31)         " <<endl;
        cout<< "                    o %m month (01-12)                    " <<endl;
        cout<< "                    o %y year without century             " <<endl;
        cout<< "                    o %Y year with century                " <<endl;
        cout<< "   * for time :                                           " <<endl;
        cout<< "                    o %H hour (24-hour clock)             " <<endl;
        cout<< "                    o %I hour (12-hour clock)             " <<endl;
        cout<< "                    o %p local equivalent of AM or PM     " <<endl;
        cout<< "                    o %M minute (00-59)                   " <<endl;
        cout<< "                    o %S seconds (00-61)                  " <<endl;
        cout<< "                    o %% %                                " <<endl;
        return;
    }
    /*
    The other characters are output as is. For example to have a format like dd/mm/yyyy one should do:
    h->GetXaxis()->SetTimeFormat("%d\/%m\/%Y");
    If the time format is not defined, a default one will be computed automatically.
     */
    TAxis* axis = 0;
    if (htemp()) {
        axis = htemp()->GetXaxis();
    }
    if (!axis) {
        if (gtemp()) {
            axis = gtemp()->GetXaxis();
        }
    }
    if (axis) {
        if (format) {
            axis->SetTimeDisplay(1);
            axis->SetTimeFormat(format);
        }
        else {
            // special option format = 0 (NULL) to reset time axis
            axis->SetTimeDisplay(0);
        }
        gPad->Modified();
        gPad->Update();
    }
    else cout<< "***error: no hist or graph found" <<endl;
}

void nul(Axis_t min)
{
    TIter next(gPad->GetListOfPrimitives());
    TObject* obj = 0;
    while ((obj = next())) if (obj->IsA()->InheritsFrom(TH1::Class())) break; else obj=0;
    if (obj) {
        TH1* h = (TH1*) obj;
        h->SetMinimum(min);
        h->Draw();
    }
}

void ct()
{
    // substitute canvas window title by object title
    TIter next(gROOT->GetListOfCanvases());
    TCanvas* c;
    while ((c = (TCanvas*) next())) {
        //cout<< c->GetName() <<endl;
        TIter next_obj(c->GetListOfPrimitives());
        TObject* obj = 0;
        while ((obj = next_obj())) {
            if ( false
                 || obj->IsA()->InheritsFrom(TH1::Class())
                 || obj->IsA()->InheritsFrom(TGraph::Class())
                 || obj->IsA()->InheritsFrom(TMultiGraph::Class())
                 || obj->IsA()->InheritsFrom(TF1::Class())
               ) break;
            else obj=0;
        }
        if (obj) {
            const char* tit = obj->GetTitle();
            //cout<< "      " << tit <<endl;
            c->SetTitle(tit);
        }
    }
}
void ctadd()
{
    // add to canvas window title an object title
    TIter next(gROOT->GetListOfCanvases());
    TCanvas* c;
    while ((c = (TCanvas*) next())) {
        //cout<< c->GetName() <<endl;
        TIter next_obj(c->GetListOfPrimitives());
        TObject* obj = 0;
        while ((obj = next_obj())) {
            if ( false
                 || obj->IsA()->InheritsFrom(TH1::Class())
                 || obj->IsA()->InheritsFrom(TGraph::Class())
                 || obj->IsA()->InheritsFrom(TMultiGraph::Class())
                 || obj->IsA()->InheritsFrom(TF1::Class())
               ) break;
            else obj=0;
        }
        if (obj) {
            const char* tit = obj->GetTitle();
            //cout<< "      " << tit <<endl;
            c->SetTitle(Form("%s %s",c->GetName(),tit));
        }
    }
}
void ctlist()
{
    // list canvas window title and object title
    TIter next(gROOT->GetListOfCanvases());
    TCanvas* c;
    while ((c = (TCanvas*) next())) {
        //cout<< c->GetName() <<endl;
        TIter next_obj(c->GetListOfPrimitives());
        TObject* obj = 0;
        while ((obj = next_obj())) {
            if ( false
                 || obj->IsA()->InheritsFrom(TH1::Class())
                 || obj->IsA()->InheritsFrom(TGraph::Class())
                 || obj->IsA()->InheritsFrom(TMultiGraph::Class())
                 || obj->IsA()->InheritsFrom(TF1::Class())
               ) break;
            else obj=0;
        }
        if (obj) {
            //const char* tit = obj->GetTitle();
            //cout<< "   " << tit <<endl;
            //c->SetTitle(Form("%s %s",c->GetName(),tit));
            cout<< c->GetTitle() << "\t " << obj->GetTitle() <<endl;
        }
    }
}
const char* addtit(const char* title, TCanvas* can)
{
    // adds text to existing canvas title to have e.g. "c1_n3 before the fit"
    if (can == 0) can = (TCanvas*) gPad; // conversion from ‘TVirtualPad*’ to ‘TCanvas*
    if (can == 0) return 0;

    std::stringstream ss;
    ss.str("");
    ss << can->GetTitle() << " " << title;
    can->SetTitle(ss.str().c_str());
    return can->GetTitle();
}

TChain* chlist(const char* tree_name, const char* filelist) {
    TChain* ch = new TChain(tree_name);
    ifstream f(filelist);
    if (!f.is_open()) {
        cout<< "Filelist not found: " << filelist <<endl;
        return 0;
    }
    std::string str;
    while (f >> str) {
        // // check if file exist (TChain::Add do not check it)
        // ifstream ftmp(str.c_str());
        // if (!ftmp.is_open()) {
        //    cout<< "File not found: " << str <<endl;
        //    return 0;
        // }
        // else ftmp.close();
        //  add to chain
        ch->Add(str.c_str(),-1);
    }
    cout<< "Total events: " << ch->GetEntries() <<endl;
    return ch;
}

// Stat_t hint(char* hname, Double_t x1, Double_t x2)
// {
//    TObject* obj = (TObject*) gDirectory->Get(hname);
//    if (obj && strcmp(obj->IsA()->GetName(),TH1F::Class()->GetName())==0)
//    {
//       TH1* h = (TH1*) obj;
//       Int_t bin1 = h->GetXaxis()->FindBin(x1);
//       if (x2 == 0) x2 = h->GetXaxis()->GetXmax();
//       Int_t bin2 = h->GetXaxis()->FindBin(x2);
//       Stat_t sum = h->Integral(bin1,bin2);
//       Stat_t tot = h->Integral(h->GetXaxis()->GetFirst(),h->GetXaxis()->GetLast());
//       //cout<< "Intergal is " << sum << " ==> " << 100.*Double_t(sum)/Double_t(tot) << " %" <<endl;
//       cout<< "Intergal from " << x1 << " to " << x2 << " is " << sum << " ==> " << 100.*Double_t(sum)/Double_t(tot) << " %" <<endl;
//       return sum;
//    }
//    else {
//       cout<< "Could not found histogram " << hname <<endl;
//       return 0;
//    }
// }
// 
// Stat_t hint(TH1* h, Double_t x1, Double_t x2)
// {
//    Int_t bin1 = h->GetXaxis()->FindBin(x1);
//    if (x2 == 0) x2 = h->GetXaxis()->GetXmax();
//    Int_t bin2 = h->GetXaxis()->FindBin(x2);
//    Stat_t sum = h->Integral(bin1,bin2);
//    Stat_t tot = h->Integral(h->GetXaxis()->GetFirst(),h->GetXaxis()->GetLast());
//    cout<< "Intergal is " << sum << " ==> " << 100.*Double_t(sum)/Double_t(tot) << " %" <<endl;
//    return sum;
// }

Stat_t hint(Double_t x1, Double_t x2, TH1* h)
{
    if (!h) h = htemp();
    if (!h) return 0;

    Int_t bin1 = 1;
    Int_t bin2 = h->GetNbinsX();
    if (x1 == 0 && x2 == 0) {
        bin1 = 1;
        bin2 = h->GetNbinsX();
    }
    else {
        bin1 = h->GetXaxis()->FindBin(x1);
        if (x2 <= x1) x2 = h->GetXaxis()->GetXmax();
        bin2 = h->GetXaxis()->FindBin(x2);
    }
    Stat_t sum = h->Integral(bin1,bin2);
    Stat_t tot = h->Integral(h->GetXaxis()->GetFirst(),h->GetXaxis()->GetLast());
    cout<< "Intergal is " << sum << " ==> " << 100.*Double_t(sum)/Double_t(tot) << " %" <<endl;
    return sum;
}

void h2subtract(Double_t pedestal, TH2* h2)
{
    // subtracts pedestal from each cell. Set content to 0 if it's negative.
    bool current_canvas = false;
    if (!h2) {
        h2 = htemp2();
        current_canvas = true;
    }
    if (!h2) return;

    for (int ix=1; ix<=h2->GetNbinsX(); ++ix) {
        for (int iy=1; iy<=h2->GetNbinsY(); ++iy) {
            Int_t bin = h2->GetBin(ix,iy);
            Double_t val = h2->GetBinContent(bin);
            if (val < pedestal) val = 0;
            else val -= pedestal;
            h2->SetBinContent(bin, val);
        }
    }

    if (current_canvas) {
        gPad->Modified();
        gPad->Update();
    }
}

/*
t->Draw(">>elist","np==1&&Nt(0.1)>0&&Nt(0.1)<=3"); TEventList* elist=(TEventList*) gDirectory->Get("elist"); for(int i=0; i<10; ++i) {cout<< elist->GetEntry(i) <<endl;}
Int_t n = -1                           // to have n to be equal to the current event number
cloop("FillDrawCells(elist->GetEntry(++n))")

----> to see prev event issue command n-=2
 */
void cloop(const char* command)
{
    TTimer timer("gSystem->ProcessEvents();",50,kFALSE);  // process mouse events every 50 ms

    while (true)
    {
        gROOT->ProcessLine(command);

        timer.Start();    // start processing of mouse events on TCanvas
        cout<< "<CR>=next, Clone, Quit, command ";
        std::string line;
        std::getline(std::cin, line);
        timer.Stop();     // disable timer after <CR>

        if (line.size() == 0) continue;           // <CR>

        if (line == "Q" || line == "q") break;    // Q or q

        if (line == "C" || line == "c") {         // C or c
            gPad->DrawClone();
            continue;
        }

        if (line == ".q") {
            cout<< "--> Confirm command to quit ROOT: ";
            std::getline(std::cin, line);
        }
        gROOT->ProcessLine(line.c_str());         // interpret as a command
    }
}

#endif  // define macro_utils_C
