void rtd()
{
   // RTD Pt1000 measurements Monterey, 370 Clay St Apt 15, Dec 23, 2014

   // T(R) data from http://www.km.kongsberg.com/ks/web/nokbg0397.nsf/AllWeb/A707D00EE0F558D6C12574E1002C2D1C/$file/tsiec751_ce.pdf?OpenElement

   Double_t R[] = {1086,  1140, 1632};
   Double_t T[] = {22.,   36.,  166.};

   TGraph* gt = new TGraph(sizeof(R)/sizeof(Double_t), R, T);
   gt->SetMarkerStyle(20);
   gt->SetFillStyle(0);
   //-- gt->SetNameTitle("gt", "T ^{o}C vs Resistance;R, Ohm;T ^{o}C");
   gt->SetNameTitle("gt", "T vs R for RTD Pt1000");
   gt->GetXaxis()->SetTitle("R, Ohm");
   gt->GetYaxis()->SetTitle("T, ^{o}C");
   gt->SetMinimum(0);

   TCanvas* can = new TCanvas;
   gt->Draw("ap");
   gt->Fit("pol1");

   TF1* f3850 = new TF1("f3850","3850e-4*(x-1000)", 1000, 2000);
   f3850->SetLineColor(4);
   f3850->SetLineStyle(2);

   TF1* f2600 = new TF1("f2600","2600e-4*(x-1000)", 1000, 2000);
   f2600->SetLineColor(8);
   f2600->SetLineStyle(2);

   f3850->Draw("same");
   f2600->Draw("same");
   TLegend* legend = can->BuildLegend();
   legend->SetFillStyle(1000);
   legend->SetFillColor(0);

   cout<< "\nThe best fit to oven data is T = 2600*(R-1000)\n" <<endl;
}
