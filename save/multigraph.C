// AZ: based on code from https://root.cern.ch/doc/v610/classTMultiGraph.html

void multigraph()
{
   auto c2 = new TCanvas("c2","c2",600,400);
   TGraph *g[3];
   Double_t x[10] = {0,1,2,3,4,5,6,7,8,9};
   Double_t y[10] = {1,2,3,4,5,5,4,3,2,1};
   auto mg = new TMultiGraph();
   for (int i=0; i<3; i++) {
      g[i] = new TGraph(10, x, y);
      g[i]->SetMarkerStyle(20);
      g[i]->SetMarkerColor(i+2);
      for (int j=0; j<10; j++) y[j] = y[j]-1;
      mg->Add(g[i]);
   }
   mg->Draw("APL");
   
   //
   // Set axes titles AFTER you call mg->Draw()
   //
   mg->SetTitle("The MultiGraph Title");
   mg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
   mg->GetYaxis()->SetTitle("Coefficients");
   
   // Change the axis limits
   gPad->Modified();
   mg->GetXaxis()->SetLimits(1.5,7.5);
   mg->SetMinimum(0.);
   mg->SetMaximum(10.);
}
