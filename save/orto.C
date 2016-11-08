#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>

#include <iostream>
using std::cout;    using std::endl;

//
// a general C function with parameter
//
Double_t orto_function(Double_t *xx, Double_t *par)
{
    Double_t x = *xx;
    Double_t sum = 0;
    for (int i=0; i<256; ++i) {
        sum += TMath::Cos(TMath::TwoPi()*par[0]*i/256) * TMath::Cos(TMath::TwoPi()*x*i/256);
    }
    return sum;
}

//
// a general C++ function object (functor)
//
class Orto_functor {
public:
    Int_t N = 256;
    Double_t fbase;
    Orto_functor(Double_t the_fbase=20.) {fbase = the_fbase;}
    Double_t operator() (Double_t *xx, Double_t *)  // the ROOT example uses double instead of Double_t in this line
    {
        Double_t x = *xx;
        Double_t sum = 0;
        for (int i=0; i<256; ++i) {
            sum += TMath::Cos(TMath::TwoPi()*fbase*i/256) * TMath::Cos(TMath::TwoPi()*x*i/256);
        }
        return sum;
    }
};

/*
Try fractional number as well as integer:
.x orto.C+(20.25)
*/
void orto(Double_t kfreq=20.)
{
    //
    // use a general C function with parameters
    //
    // TF1 *f = new TF1("f", orto_function, kfreq-5, kfreq+5, 1);
    // f->SetParameter(0, kfreq);

    //
    // use a general C++ function object (functor)
    //
    Orto_functor *orto_functor = new Orto_functor(kfreq);
    TF1 *f = new TF1("f", orto_functor, kfreq-5., kfreq+5., 0, "Orto_functor");

    f->SetTitle(Form("Ortogonality of kfreq = %0.2f with continuous k;k continuous",kfreq));
    f->SetNpx(1000);

    new TCanvas;
    f->Draw();
}
