#include <TROOT.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>

#include <iostream>

using std::cout;     using std::endl;

void dEdx()
{
   new TCanvas;
   TF1* fR = new TF1("fR", "10.*pow(x,1+[1])/([0]*(1+[1]));T_{proton}, MeV;R, mm", 0,200); fR->SetParameters(247.93,0.7637);
   fR->Draw();

   new TCanvas;
   TF1* fEx = new TF1("fEx", "pow(pow([0],1+[2])-[1]*(1+[2])*(0.1*x), 1./(1+[2]));x, mm;T_{proton}, MeV", 0,262);   fEx->SetParameters(200,247.93,0.7637);
   fEx->Draw();

   new TCanvas;
   TF1* fdEdx = new TF1("fdEdx", "[0]*pow(x,-[1]);E, MeV;dE/dx", 1,200); fdEdx->SetParameters(247.93,0.7637);
   fdEdx->Draw();

   new TCanvas;
   TF1* fdEdxx = new TF1("fdEdxx", "[1]*pow(pow([0],1+[2])-[1]*(1+[2])*(0.1*x), -[2]/(1+[2])); x, mm; dE/dx, MeV/cm", 0,262);   fdEdxx->SetParameters(200,247.93,0.7637);
   fdEdxx->Draw();
}

class FunctorEx {
public:
   bool debug;
   double a;
   double b;
   double E0;
   FunctorEx(double the_E0, bool the_debug=false): debug(the_debug), a(247.93), b(0.7637), E0(the_E0) {}
   double operator ()(double* xx, double*) {
      double& x = *xx;
      x *= 0.1;            // convert to cm
      double Ex = pow(pow(E0,b+1) -a*(1+b)*x, 1./(1.+b));
      if (debug) cout<< "Ex = " << Ex <<endl;
      return Ex > 0? Ex: 0;
   }
   double f(double x) {
      double xx[1];
      xx[0] = x;
      return (*this)(xx,0);
   }
};

void test_FunctorEx(double E0=200)
{
   FunctorEx* functorEx = new FunctorEx(E0);
   TF1* fEx = new TF1("fEx", functorEx, 0,300, 0, "FunctorEx");

   new TCanvas;
   fEx->Draw();
}

class FunctorDegrader {
public:
   //double anorm;
   bool debug;
   double a;
   double b;
   double E0;
   double t0;    // scintillator thickness, mm
   double t1;    // scintillator thickness, mm
   double t2;    // scintillator thickness, mm
   double tfoil;
   double tssd;   // total thickness of the SSD, mm
   double t_paper;
   double t_teflon_mylar;
   //FunctorDegrader(double the_E0, double anorm_0=1., bool the_debug=false): anorm(anorm_0), debug(the_debug)
   FunctorDegrader(double the_E0, bool the_debug=false): debug(the_debug)
      , a(247.93), b(0.7637), E0(the_E0)
      , t0(76), t1(76), t2(103)
   {
      tssd = 8.*0.4*2.33/1.05;         // 7.1 mm of polystyrene 
      t_paper = 0.150*0.91;            // approx. water equiv. thickness of the photo paper
      t_teflon_mylar = 0.150;          // approx. thickness of teflon+mylar
      tfoil = 1.9*13.6/1.05;           // 24.6 mm of polystyrene
   }
   double Epx(double xmm) const {
      // proton energy (in MeV) after travelling distance x (in mm)
      double x = 0.1*xmm;                                // convert in cm
      double Ex = pow(pow(E0,b+1) -a*(1+b)*x, 1./(1.+b));
      return Ex > 0? Ex: 0;
   }
   double operator ()(double* xx, double*) const {
      double& x = *xx;
      return Epx(x);
   }
   double f(double x) const {
      double xx[1];
      xx[0] = x;
      return (*this)(xx,0);
   }
   // double Signal0(double* xx, double*) const {
   //    // energy deposited in the sci0
   //    double& x_d = *xx;
   //    double t_wrap0 = tfoil+tssd + 2*t_paper + t_teflon_mylar;     // wrapping seen by proton stopped in sci0
   //    double E_d = Epx(x_d + t_wrap0);          // Energy at entrance of sci0
   //    if (E_d == 0) return 0;
   //    double Eexit = Epx(x_d+t0);
   //    if (Eexit > 0) E_d -= Eexit;
   //    return anorm*E_d;    // here proton will deposit rest of the energy into sci0
   // }
   double fSignal0(double* xx, double* par) const {
      // energy deposited in the sci0
      double x_d = *xx + par[1];
      double t_wrap0 = tfoil+tssd + 2*t_paper + t_teflon_mylar;     // wrapping seen by proton stopped in sci0
      double t_entry = x_d + t_wrap0;
      double t_exit = t_entry + t0;
      double E_d = Epx(t_entry);          // Energy at entrance of sci0
      if (E_d == 0) return 0;
      double Eexit = Epx(t_exit);
      if (Eexit > 0) E_d -= Eexit;
      return par[0]*E_d;    // here proton will deposit rest of the energy into sci0
   }
   double fSignal1(double* xx, double* par) const {
      // energy deposited in the sci1
      double& x_d = *xx;
      double t_entrance = par[1] + tfoil+tssd + x_d + 2*t_paper + 2*t_teflon_mylar + t0 + 2*t_paper + t_teflon_mylar; // wrapping seen by proton stopped in sci1
      double t_exit = t_entrance + t1;
      double Eentrance = Epx(t_entrance);
      double Eexit = Epx(t_exit);
      double Esci1 = Eentrance - Eexit;
      return par[0]*Esci1;
   }
   // double Signal1(double* xx, double*) const {
   //    // energy deposited in the sci1
   //    double& x_d = *xx;
   //    double t_entrance = tfoil+tssd + x_d + 2*t_paper + 2*t_teflon_mylar + t0 + 2*t_paper + t_teflon_mylar; // wrapping seen by proton stopped in sci1
   //    double t_exit = t_entrance + t1;
   //    double Eentrance = Epx(t_entrance);
   //    double Eexit = Epx(t_exit);
   //    double Esci1 = Eentrance - Eexit;
   //    return anorm*Esci1;
   // }
   double fSignal2(double* xx, double* par) const {
      // energy deposited in the sci2
      double& x_d = *xx;
      double t_wrap2 = par[1] + tfoil+tssd + 2*t_paper + 2*t_teflon_mylar + 2*t_paper + 2*t_teflon_mylar + 2*t_paper + t_teflon_mylar; // wrapping seen by proton entered sci1
      double Ebefore = Epx(x_d+t0+t1 + t_wrap2);          // Energy after the degrader sci0 and sci1
      return par[0]*Ebefore;    // here proton will deposit rest of the energy into sci1
   }
   // double Signal2(double* xx, double*) const {
   //    // energy deposited in the sci2
   //    double& x_d = *xx;
   //    double t_wrap2 = tfoil+tssd + 2*t_paper + 2*t_teflon_mylar + 2*t_paper + 2*t_teflon_mylar + 2*t_paper + t_teflon_mylar; // wrapping seen by proton entered sci1
   //    double Ebefore = Epx(x_d+t0+t1 + t_wrap2);          // Energy after the degrader sci0 and sci1
   //    return anorm*Ebefore;    // here proton will deposit rest of the energy into sci1
   // }
   double Bragg0(double* xx, double* par) const {
      // energy deposited in the sci0 when sci1 has no signal
      //double x_d = *xx;
      //double t_wrap0 = tfoil+tssd + 2*t_paper + t_teflon_mylar;     // wrapping seen by proton stopped in sci0
      //double t_wrap1 = t_wrap0 + t_teflon_mylar + 2*t_paper + t_teflon_mylar;     // wrapping seen by proton passed sci0
      //double E_d = Epx(x_d + t_wrap0);          // Energy at entrance of sci0
      //if (E_d == 0) return 0;
      //if (Epx(x_d+t0 + t_wrap1) > 0) return 0;  // reject events with signal in sci1
      //return anorm*E_d;    // here proton will deposit rest of the energy into sci0

      double eps = 1e-7;
      double signal1 = fSignal1(xx,par);
      if (signal1 > eps) return 0;
      return fSignal0(xx,par);
   }
   double Bragg1(double* xx, double* par) const {
      // energy deposited in the sci1 when sci1 has no signal
      //double x_d = *xx;
      //double t_wrap1 = tfoil+tssd + 2*t_paper + 2*t_teflon_mylar + 2*t_paper + t_teflon_mylar; // wrapping seen by proton stopped in sci1
      //double t_wrap2 = t_wrap1 + t_teflon_mylar + 2*t_paper + t_teflon_mylar;     // wrapping seen by proton passed sci1
      //double Ebefore = Epx(x_d+t0 + t_wrap1);          // Energy after the degrader and sci0
      //if (Ebefore == 0) return 0;
      //if (Epx(x_d+t0+t1 + t_wrap2) > 0) return 0;  // reject events with signal in sci2
      //return anorm*Ebefore;    // here proton will deposit rest of the energy into sci1

      double eps = 1e-7;
      double signal2 = fSignal2(xx,par);
      if (signal2 > eps) return 0;
      return fSignal1(xx,par);
   }
   double Bragg2(double* xx, double* par) const {
      // energy deposited in the sci2
      //double x_d = *xx;
      //double t_wrap2 = tfoil+tssd + 2*t_paper + 2*t_teflon_mylar + 2*t_paper + 2*t_teflon_mylar + 2*t_paper + t_teflon_mylar; // wrapping seen by proton entered sci1
      //double Ebefore = Epx(x_d+t0+t1 + t_wrap2);          // Energy after the degrader sci0 and sci1
      //return anorm*Ebefore;    // here proton will deposit rest of the energy into sci1
      return fSignal2(xx,par);
   }
};

FunctorDegrader* test_FunctorDegrader(double E0=200)
{
   FunctorDegrader* functorDegrader = new FunctorDegrader(E0);
   TF1* fEx = new TF1("fEx", functorDegrader, 0,300, 2, "FunctorDegrader");
   fEx->SetTitle("E_{proton}(x);x, mm;E_{proton}, MeV");

   // cout<< "fEx(150) = " << fEx(150) <<endl;        // error in C++, but OK with CINT
   // cout<< "(*fEx)(150) = " << (*fEx)(150) <<endl;  // C++ way

   new TCanvas;
   fEx->Draw();

   return functorDegrader;
}

FunctorDegrader* test_Signal0(double E0=200, double anorm=1, double offset=0)
{
   FunctorDegrader* functorDegrader = new FunctorDegrader(E0);
   TF1* fSignal0 = new TF1("fSignal0", functorDegrader, &FunctorDegrader::fSignal0, 120,300, 2, "FunctorDegrader", "fSignal0");
   fSignal0->SetTitle("Signal in sci0. Signal from sci0 in MeV;degrader thickness, mm;signal, MeV");
   fSignal0->SetLineColor(4);
   fSignal0->SetParName(0,"anorm");
   fSignal0->SetParName(1,"offset");
   fSignal0->SetParameters(anorm,offset);

   new TCanvas;
   fSignal0->Draw();

   return functorDegrader;
}

FunctorDegrader* test_Signal1(double E0=200, double anorm=1, double offset=0)
{
   FunctorDegrader* functorDegrader = new FunctorDegrader(E0);
   TF1* fSignal1 = new TF1("fSignal1", functorDegrader, &FunctorDegrader::fSignal1, 0,200, 2, "FunctorDegrader", "fSignal1");
   fSignal1->SetTitle("Signal in sci1. Signal from sci1 in MeV;degrader thickness, mm;signal, MeV");
   fSignal1->SetLineColor(4);
   fSignal1->SetParName(0,"anorm");
   fSignal1->SetParName(1,"offset");
   fSignal1->SetParameters(anorm,offset);

   new TCanvas;
   fSignal1->Draw();

   return functorDegrader;
}

TF1* get_fSignal1(double E0=200, double anorm=1, double offset=0)
{
   FunctorDegrader* functorDegrader = new FunctorDegrader(E0);
   TF1* fSignal1 = new TF1("fSignal1", functorDegrader, &FunctorDegrader::fSignal1, 0,200, 2, "FunctorDegrader", "fSignal1");
   fSignal1->SetTitle("Signal in sci1. Signal from sci1 in MeV;degrader thickness, mm;signal, MeV");
   fSignal1->SetLineColor(4);
   fSignal1->SetParName(0,"anorm");
   fSignal1->SetParName(1,"offset");
   fSignal1->SetParameters(anorm,offset);

   new TCanvas;
   fSignal1->Draw();

   return fSignal1;
}

FunctorDegrader* test_Signal2(double E0=200, double anorm=1, double offset=0)
{
   FunctorDegrader* functorDegrader = new FunctorDegrader(E0);
   TF1* fSignal2 = new TF1("fSignal2", functorDegrader, &FunctorDegrader::fSignal2, 0,120, 2, "FunctorDegrader", "fSignal2");
   fSignal2->SetTitle("Signal in sci2. Signal from sci2 in MeV;degrader thickness, mm;signal, MeV");
   fSignal2->SetLineColor(4);
   fSignal2->SetParName(0,"anorm");
   fSignal2->SetParName(1,"offset");
   fSignal2->SetParameters(anorm,offset);

   new TCanvas;
   fSignal2->Draw();

   return functorDegrader;
}

FunctorDegrader* test_Bragg0(double E0=200, double anorm=1, double offset=0)
{
   FunctorDegrader* functorDegrader = new FunctorDegrader(E0);
   TF1* fBragg0 = new TF1("fBragg0", functorDegrader, &FunctorDegrader::Bragg0, 120,300, 2, "FunctorDegrader", "Bragg0");
   fBragg0->SetTitle("Bragg in sci0. Signal from sci0 in MeV;degrader thickness, mm;signal, MeV");
   fBragg0->SetLineColor(4);
   fBragg0->SetParName(0,"anorm");
   fBragg0->SetParName(1,"offset");
   fBragg0->SetParameters(anorm,offset);

   new TCanvas;
   fBragg0->Draw();

   return functorDegrader;
}

FunctorDegrader* test_Bragg1(double E0=200, double anorm=1, double offset=0)
{
   FunctorDegrader* functorDegrader = new FunctorDegrader(E0);
   TF1* fBragg1 = new TF1("fBragg1", functorDegrader, &FunctorDegrader::Bragg1, 0,200, 2, "FunctorDegrader", "Bragg1");
   fBragg1->SetTitle("Bragg in sci1. Signal from sci1 in MeV;degrader thickness, mm;signal, MeV");
   fBragg1->SetLineColor(4);
   fBragg1->SetParName(0,"anorm");
   fBragg1->SetParName(1,"offset");
   fBragg1->SetParameters(anorm,offset);

   new TCanvas;
   fBragg1->Draw();

   return functorDegrader;
}

FunctorDegrader* test_Bragg2(double E0=200, double anorm=1, double offset=0)
{
   FunctorDegrader* functorDegrader = new FunctorDegrader(E0);
   TF1* fBragg2 = new TF1("fBragg2", functorDegrader, &FunctorDegrader::Bragg2, 0,120, 2, "FunctorDegrader", "Bragg2");
   fBragg2->SetTitle("Bragg in sci2. Signal from sci2 in MeV;degrader thickness, mm;signal, MeV");
   fBragg2->SetLineColor(4);
   fBragg2->SetParName(0,"anorm");
   fBragg2->SetParName(1,"offset");
   fBragg2->SetParameters(anorm,offset);

   new TCanvas;
   fBragg2->Draw();

   return functorDegrader;
}
