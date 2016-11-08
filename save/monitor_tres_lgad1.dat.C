{
TGraph* g = new TGraph(__FILE__);
TH1D* h_mon_lgad1 = new TH1D("h_mon_lgad1","Time resolution of SiPM-LGAD1 @150V;time resolution, ps;runs", 45,40,55);  
for (int i=0; i<g->GetN(); ++i) h_mon_lgad1->Fill(g->GetY()[i]);
new TCanvas;
h_mon_lgad1->Draw();
fitgl(0,0,"L");
titmax();
}

run     ps
31      48.7
32      47.5
33      50.8
34      50.2
35      49.8
36      49.1
37      48.3
38      48.8
39      50.2
40      51.3
41      49.3
42      49.5
43      50.5
44      50.0
