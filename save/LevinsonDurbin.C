#include "DRS4Event.h"

// http://www.emptyloop.com/technotes/A%20tutorial%20on%20linear%20prediction%20and%20Levinson-Durbin.pdf

#include <TGraph.h>
#include <TCanvas.h>

#include <math.h>
#include <vector>

using namespace std;

// Returns in vector linear prediction coefficients calculated using Levinson Durbin
void ForwardLinearPrediction(size_t N, const Double_t* x, size_t m, Double_t* coeffs)
{
   // GET SIZE FROM INPUT VECTORS
   //-- size_t N = x.size() - 1;
   //-- size_t m = coeffs.size();

   // INITIALIZE R WITH AUTOCORRELATION COEFFICIENTS
   //-- vector<double> R( m + 1, 0.0 );
   Double_t R[100];
   for (size_t i=0; i<=m; ++i) R[i] = 0;

   for (size_t i=0; i<=m; ++i)
   {
      for (size_t j=0; j<=N-i; ++j) {
         R[i] += x[j]*x[j+i];
      }
   }

   // INITIALIZE Ak
   //-- vector<double> Ak( m + 1, 0.0 );
   Double_t Ak[100];
   for (size_t i=0; i<=m; ++i) Ak[i] = 0;
   Ak[0] = 1.0;

   // INITIALIZE Ek
   double Ek = R[0];

   // LEVINSON-DURBIN RECURSION
   for (size_t k=0; k<m; ++k) {
      // COMPUTE LAMBDA
      double lambda = 0.0;
      for (size_t j=0; j<=k; ++j) {
         lambda -= Ak[j] * R[k+1-j];
      }
      lambda /= Ek;

      // UPDATE Ak
      for (size_t n=0; n<=(k+1)/2; ++n) {
         double temp = Ak[k+1-n] + lambda*Ak[n];
         Ak[n] = Ak[n] + lambda*Ak[k+1-n];
         Ak[k+1-n] = temp;
      }

      // UPDATE Ek
      Ek *= 1.0 - lambda * lambda;
   }
   // ASSIGN COEFFICIENTS
   //-- coeffs.assign( ++Ak.begin(), Ak.end() );
   for (size_t i=0; i<m; ++i) coeffs[i] = Ak[i+1];
}
void ForwardLinearPrediction( vector<double> &coeffs, const vector<double> &x )
{
   // GET SIZE FROM INPUT VECTORS
   size_t N = x.size() - 1;
   size_t m = coeffs.size();

   // INITIALIZE R WITH AUTOCORRELATION COEFFICIENTS
   vector<double> R( m + 1, 0.0 );
   for ( size_t i = 0; i <= m; i++ )
   {
      for ( size_t j = 0; j <= N - i; j++ ) {
         R[ i ] += x[ j ] * x[ j + i ];
      }
   }

   // INITIALIZE Ak
   vector<double> Ak( m + 1, 0.0 );
   Ak[ 0 ] = 1.0;

   // INITIALIZE Ek
   double Ek = R[ 0 ];

   // LEVINSON-DURBIN RECURSION
   for ( size_t k = 0; k < m; k++ ) {
      // COMPUTE LAMBDA
      double lambda = 0.0;
      for ( size_t j = 0; j <= k; j++ ) {
         lambda -= Ak[ j ] * R[ k + 1 - j ];
      }
      lambda /= Ek;

      // UPDATE Ak
      for ( size_t n = 0; n <= ( k + 1 ) / 2; n++ ) {
         double temp = Ak[ k + 1 - n ] + lambda * Ak[ n ];
         Ak[ n ] = Ak[ n ] + lambda * Ak[ k + 1 - n ];
         Ak[ k + 1 - n ] = temp;
      }

      // UPDATE Ek
      Ek *= 1.0 - lambda * lambda;
   }
   // ASSIGN COEFFICIENTS
   coeffs.assign( ++Ak.begin(), Ak.end() );
}

// Example program using Forward Linear Prediction
int LevinsonDurbin()
{
   // CREATE DATA TO APPROXIMATE
   vector<double> original( 128, 0.0 );
   for ( size_t i=0; i<original.size(); i++) {
      original[i] = sin(i*0.01) + 0.75*sin(i*0.03) + 0.5*sin(i*0.05) + 0.25*sin(i*0.11);
   }
   // GET FORWARD LINEAR PREDICTION COEFFICIENTS
   //-- vector<double> coeffs( 4, 0.0 );
   vector<double> coeffs(16, 0);
   ForwardLinearPrediction(coeffs, original);

   // PREDICT DATA LINEARLY
   vector<double> predicted(original);
   size_t m = coeffs.size();
   for (size_t i=m; i<predicted.size(); i++) {
      predicted[i] = 0;
      for (size_t j=0; j<m; j++) {
         predicted[i] -= coeffs[j] * original[i-1-j];
      }
   }

   // CALCULATE AND DISPLAY ERROR
   double error = 0;
   for (size_t i=m; i<predicted.size(); i++) {
      printf("Index: %.2lu / Original: %.6f / Predicted: %.6f\n", i, original[i], predicted[i]);
      double delta = predicted[i] - original[i];
      error += delta*delta;
   }
   printf("Forward Linear Prediction Approximation Error: %f\n", error);

   Double_t a_i[1000];
   Double_t a_orig[1000];
   Double_t a_pred[1000];
   Int_t np = 0;

   for (size_t i=m; i<predicted.size(); i++) {
      a_i[np] = np;
      a_orig[np] = original[i];
      a_pred[np] = predicted[i];
      ++np;
   }

   TGraph* gorig = new TGraph(np,a_i,a_orig);
   gorig->SetNameTitle("gorig","original data");
   gorig->SetMarkerStyle(7);
   gorig->SetMarkerColor(1);
   gorig->SetLineColor(1);

   TGraph* gpred = new TGraph(np,a_i,a_pred);
   gpred->SetNameTitle("gpred","predinal data");
   gpred->SetMarkerStyle(7);
   gpred->SetMarkerColor(2);
   gpred->SetLineColor(2);

   new TCanvas;
   gorig->Draw("apl");
   gpred->Draw("pl");

   return 0;
}

/////////////////////// Levinson-Durbin /////////////////////////////

void moving_average(Int_t naver, Int_t np, const Double_t* y, Double_t* yaver)
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

/*
root -l mppc-72.0V-dark.dat.root
.L LevinsonDurbin.C+
linpred(8,8,128, 27,0, 420)
 */
void linpred(size_t npar=8, size_t nextra=8, size_t N=128, Int_t event=27, Int_t chan=0, Int_t offset=420, bool debug=true)
{
   TTree* tree = 0;
   DRS4Event* drs4Event = new DRS4Event(tree);
   cout<< "tree = " << tree <<endl;
   if (!tree) return;

   drs4Event->GetEntry(event);

   Float_t* x = drs4Event->chanT(chan);
   Float_t* y = drs4Event->chanV(chan);
   Double_t xuni[10000];
   Double_t yuni[10000];
   drs4Event->uniform(1024,x,y, xuni,yuni);

   Double_t input[10000];
   for (size_t i=0; i<1024; ++i) {
      input[i] = -1.*yuni[i];
   }

   if (debug) {
      TGraph* guni = new TGraph(1024,xuni,yuni);
      guni->SetNameTitle(Form("guni_evt_%d_ch_%d_npar_%lu",event,chan,npar), Form("Event %d ch %d npar %lu",event,chan,npar));
      guni->SetMarkerStyle(7);
      guni->SetMarkerColor(2);
      guni->SetLineColor(2);
      new TCanvas;
      guni->Draw("apl");
   }

   Double_t yaver[10000];
   moving_average(9, 1024,yuni, yaver);

   // CREATE DATA TO APPROXIMATE
   //--vector<double> original(Norig, 0);
   vector<double> original;
   for (size_t i=0; i<N; ++i) {
      //-- original.push_back(-1.*yuni[offset+i]);

      // input[i] = -1.*yaver[offset+i];
      //-- input[i] = -1.*yuni[offset+i];

      original.push_back(input[i+offset]);
   }

   //cout<< "original.size() = " << original.size() <<endl;

   // GET FORWARD LINEAR PREDICTION COEFFICIENTS
   vector<double> coeffs(npar, 0);
   //-- ForwardLinearPrediction(coeffs, original);

   Double_t a_coeffs[100];
   ForwardLinearPrediction(N-1, input, npar, a_coeffs);
   for (size_t i=0; i<npar; ++i) coeffs[i] = a_coeffs[i];

   // PREDICT DATA LINEARLY
   vector<double> predicted(original);
   for (size_t i=coeffs.size(); i<predicted.size(); ++i) {
      predicted[i] = 0;
      for (size_t j=0; j<coeffs.size(); ++j) {
         predicted[i] -= coeffs[j] * original[i-1-j];
      }
   }

   // CALCULATE AND DISPLAY ERROR
   double error = 0;
   for (size_t i = coeffs.size(); i<predicted.size(); ++i) {
      printf("Index: %.2lu / Original: %.6f / Predicted: %.6f\n", i, original[i], predicted[i]);
      double delta = predicted[i] - original[i];
      error += delta * delta;
   }
   printf("Forward Linear Prediction Approximation Error: %f\n", error);

   cout<< "--> look more" <<endl;

   for (size_t i=N; i<N+nextra; ++i) {
      //--predicted[i] = 0;
      predicted.push_back(0);
      //cout<< "i = " << i << " predicted.size() = " << predicted.size() <<endl;
      for (size_t j=0; j<coeffs.size(); ++j) {
         //cout<< "\tpredicted[i-1-j] = " << predicted[i-1-j] <<endl;
         //--predicted[i] -= coeffs[j] * original[i-1-j];
         predicted[i] -= coeffs[j] * predicted[i-1-j];
      }
      //-- original.push_back(-1.*yuni[offset+i]);      // to plot
      original.push_back(input[offset+i]);      // to plot
      //cout<< i << "\toriginal[i] = " << original[i] << "\t predicted[i] = " << predicted[i] <<endl;
   }

   //cout<< "predicted.size() = " << predicted.size() <<endl;

   //Double_t a_i[1000];
   Double_t a_orig[1000];
   Double_t a_pred[1000];
   Int_t np = 0;

   for (size_t i=0; i<predicted.size(); ++i) {
      //a_i[np] = np;
      //--a_orig[np] = original[i];
      a_orig[np] = original[i];
      a_pred[np] = predicted[i];
      ++np;
   }

   //TGraph* gorig = new TGraph(np,a_i,a_orig);
   TGraph* gorig = new TGraph(np,&xuni[offset],a_orig);
   gorig->SetNameTitle("gorig","original data");
   gorig->SetMarkerStyle(7);
   gorig->SetMarkerColor(1);
   gorig->SetLineColor(1);

   //TGraph* gpred = new TGraph(np,a_i,a_pred);
   TGraph* gpred = new TGraph(np,&xuni[offset],a_pred);
   gpred->SetNameTitle("gpred","predinal data");
   gpred->SetMarkerStyle(7);
   gpred->SetMarkerColor(2);
   gpred->SetLineColor(2);

   new TCanvas;
   gorig->Draw("apl");
   gpred->Draw("pl");

   return;
}
