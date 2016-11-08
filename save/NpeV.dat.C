{
TGraph* g = new TGraph(__FILE__);
TH1D* h32V = new TH1D("h32V","Npe for bias voltage 32 V;Npe;runs", 20,40,100);
h32V->SetLineColor(4);
TH1D* h28V = new TH1D("h28V","Npe for bias voltage 28 V;Npe;runs", 20,40,100);
h28V->SetLineColor(2);
for (int i=0; i<g->GetN(); ++i) {
    if (g->GetX()[i] == 32) h32V->Fill(g->GetY()[i]);
    if (g->GetX()[i] == 28) h28V->Fill(g->GetY()[i]);
}

new TCanvas;

TH1* h = h28V;
h->SetTitle("Npe for bias voltages 28V and 32V");

h->Draw();
fitgr(0,0,"L");
ftemp()->SetLineColor(h->GetLineColor());

h32V->Draw("sames");
fitgl(0,0,"L");
ftemp()->SetLineColor(h32V->GetLineColor());

titmax();
}

V       Npe
-----------
32      57.9
32      57.6
32      60.0
32      65.1
32      70.7
32      63.9
32      65.3
32      60.2

28      49.9
28      53.4
28      53.7
28      52.0
28      56.7
28      53.7
28      54.9
28      54.9
28      54.0
28      56.5
28      56.4
