#include <TMath.h>
#include <iostream>

using std::cout;    using std::endl;

// void unfoldResolution(Double_t Sigma12, Double_t Sigma13, Double_t Sigma23)
// {
//     Double_t sigma1 = TMath::Sqrt(0.5*(Sigma12*Sigma12 + Sigma13*Sigma13 - Sigma23*Sigma23));
//     Double_t sigma2 = TMath::Sqrt(0.5*(Sigma12*Sigma12 + Sigma23*Sigma23 - Sigma13*Sigma13));
//     Double_t sigma3 = TMath::Sqrt(0.5*(Sigma13*Sigma13 + Sigma23*Sigma23 - Sigma12*Sigma12));
// 
//     cout<< "sigma1 = " << sigma1 << " sigma2 = " << sigma2 << " sigma3 = " << sigma3 <<endl;
// }

void unfoldResolution(Double_t Sigma01, Double_t Sigma02, Double_t Sigma12, Double_t sigmaLim=16.4)
{
    Double_t arg;
    arg = 0.5*(Sigma01*Sigma01 + Sigma02*Sigma02 - Sigma12*Sigma12);
    Double_t sigma0 = arg > 0? TMath::Sqrt(arg): sigmaLim;
    arg = 0.5*(Sigma01*Sigma01 + Sigma12*Sigma12 - Sigma02*Sigma02);
    Double_t sigma1 = arg > 0? TMath::Sqrt(arg): sigmaLim;
    arg = 0.5*(Sigma02*Sigma02 + Sigma12*Sigma12 - Sigma01*Sigma01);
    Double_t sigma2 = arg > 0? TMath::Sqrt(arg): sigmaLim;

    cout<< "sigma0 = " << sigma0 << " sigma1 = " << sigma1 << " sigma2 = " << sigma2 <<endl;
}
