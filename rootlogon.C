#include <TROOT.h>
#include <TRint.h>
#include <TSystem.h>
#include <TStyle.h>

#include <iostream>

//--#define utils_cxx

using std::cout;     using std::endl;

void rootlogon()
{
   //-- cout<< "*-- Default rootlogon" <<endl;

   gSystem->Exec("echo //-- `date` >> root_hist");


   // gSystem->SetIncludePath("-I$HOME/macros -I.");
   // gROOT->LoadMacro("utils.C+");
   gROOT->LoadMacro("utils.C");     // do not compile in ROOT 6.13/02
   
   // gStyle->SetCanvasPreferGL(kTRUE);   // anti-aliasing with OpenGL -- does not work now?!

   /* Add to include path directory with utils.C to use it for ACLiC like
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <utils.C>
#endif
   */
   // gSystem->AddIncludePath(gSystem->ExpandPathName("-I$(HOME)/macros"));
   //cout<< "gSystem->GetIncludePath() = " << gSystem->GetIncludePath() <<endl;

   ((TRint*)gROOT->GetApplication())->SetPrompt("// ");

   gStyle->SetOptFit(1);
      
   gStyle->SetFillStyle(3001);

   //-- my long used sizes, 80% of default size
   //gStyle->SetCanvasDefW(540); gStyle->SetCanvasDefH(400);
   //-- 3/4 of internal area of the default: 526x382 --> 522x354 as newcan in utils.C
   // gStyle->SetCanvasDefW(526); gStyle->SetCanvasDefH(382);

   // gStyle->SetCanvasDefW(700*2/3); gStyle->SetCanvasDefH(500*2/3);   // 2/3 of default
   // gStyle->SetCanvasDefW(700/2); gStyle->SetCanvasDefH(500/2);   // 1/2 of default
   // gStyle->SetCanvasDefW(420); gStyle->SetCanvasDefH(300);   // 2/3 - 10% of default
   //-- gStyle->SetCanvasDefW(0.55*700); gStyle->SetCanvasDefH(0.55*500);    // 385x275
   //-- gStyle->SetCanvasDefW(0.54*700); gStyle->SetCanvasDefH(0.54*500);    // 378x270
   gStyle->SetCanvasDefW(0.53*700); gStyle->SetCanvasDefH(0.53*500);    // 371x265

   gStyle->SetCanvasDefX(0);
   gStyle->SetCanvasDefY(0);
   
   gStyle->SetPadGridX(kTRUE);
   gStyle->SetPadGridY(kTRUE);

   // no box around title
   // gStyle->SetTitleX(0);
   // gStyle->SetTitleW(0);
   // gStyle->SetTitleBorderSize(0);

   // To make title box transparent https://docs.google.com/Doc?docid=0AaltnJcgAgafZGNoOThnaDZfMzEzZjlxOTZ4Zmg&hl=en
   // gStyle->SetStatStyle(0);
   // gStyle->SetTitleStyle(0);
   // gROOT->ForceStyle();

   // gStyle->SetStatStyle(3001);

   ///-- // text sizes
   ///-- /// gStyle->SetTitleFontSize(.07);
   ///-- gStyle->SetLabelSize(.05,"xyz");
   ///-- gStyle->SetTitleSize(.05,"xyz");
   ///-- // and fonts
   ///-- // changed 42 --> 43 for retina
   ///-- gStyle->SetTitleFont(42,"xyz");  // if "xyz" set all 3 axes, any other value of axis will set the pad title font
   ///-- gStyle->SetTitleFont(42,"a");
   ///-- gStyle->SetLabelFont(42,"xyz");
   ///-- gStyle->SetLabelFont(42,"a");
   ///-- gStyle->SetStatFont(42);
   ///-- gStyle->SetTextFont(42);

   // // only font/size settings for retina and 2/3 - 10% canvas
   // gStyle->SetLabelSize(.048,"xyz");
   // gStyle->SetTitleSize(.048,"xyz");
   gStyle->SetLabelSize(.050,"xyz");
   gStyle->SetTitleSize(.050,"xyz");

   gStyle->SetTitleOffset(0.90);    // to fit into canvas subscript in capital, t_{UFSD}-t_{SiPM}
}
