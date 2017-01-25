#include <TROOT.h>
#include <TMarker.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TMath.h>
#include <vector>
#include <iostream>

using std::cout;    using std::endl;

Double_t Theta(Double_t x, Double_t y, Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
    Double_t dx1 = x1 - x;
    Double_t dy1 = y1 - y;
    Double_t len1 = TMath::Sqrt(dx1*dx1 + dy1*dy1);
    Double_t dx2 = x2 - x;
    Double_t dy2 = y2 - y;
    Double_t len2 = TMath::Sqrt(dx2*dx2 + dy2*dy2);

    Double_t scalar = dx1*dx2 + dy1*dy2;
    Double_t cosTheta = scalar / len1 /len2;
    Double_t theta = TMath::ACos(cosTheta);
    return theta;
}

void rubberBand()
{
    //              0           1           2           3       4           5           6           7       8           9       10
    Double_t x[] = {-0.497967, -0.281165, 0.0779133, -0.647019, 0.0779133, -0.315041, 0.321816, 0.511518, -0.0982385, 0.477642, 0.105014};
    Double_t y[] = {-0.291667, 0.28125,   -0.302083, 0.270833,  0.28125,   0.71875,   0.5,      -0.09375, -0.708333, -0.489583, 0.697917};
    Int_t np = sizeof(x) / sizeof(Double_t);

    for (int i=0; i<np; ++i) cout<< i << "\t " << x[i] << "\t " << y[i] <<endl;

    TCanvas* can = new TCanvas;
    can->DrawFrame(-1,-1, 1,1);

    TMarker* marker = new TMarker(0, 0, 20);
    for (int i=0; i<np; ++i) marker->DrawMarker(x[i], y[i]);
    can->Modified();
    can->Update();

    // find min x as the leftmost point
    Int_t min_i = 0;
    for (int i=0; i<np; ++i) if (x[i] < x[min_i]) min_i = i;

    cout<< "x[min_i] = " << x[min_i] << " y[min_i] = " << y[min_i] <<endl;

    // Use artificial point to build one leg of the angle to find another leg
    Double_t xfake = x[min_i];
    Double_t yfake = y[min_i] - 1.;

    cout<< "Theta(x[min_i], y[min_i], xfake, yfake, -0.4, 0.6) = " << Theta(x[min_i], y[min_i], xfake, yfake, -0.4, 0.6) <<endl;

    Int_t prev_i;               // now it's artificial point (xfake, yfake)
    Int_t curr_i = min_i;
    Int_t next_i = min_i;

    std::vector<int> vertices;
    vertices.push_back(curr_i);

    Double_t theta_max = 0;
    for (int i=0; i<np; ++i)
    {
        // if (i == curr_i) continue;
        Double_t theta = Theta(x[curr_i], y[curr_i], xfake, yfake, x[i], y[i]);
        if (theta > theta_max) {
            theta_max = theta;
            next_i = i;
        }
    }
    prev_i = curr_i; 
    curr_i = next_i;
    vertices.push_back(curr_i);

    cout<< "curr_i = " << curr_i << " x[curr_i] = " << x[curr_i] << " y[curr_i] = " << y[curr_i] <<endl;

    while (next_i != min_i)
    {
        theta_max = 0;
        for (int i=0; i<np; ++i) {
            Double_t theta = Theta(x[curr_i], y[curr_i], x[prev_i], y[prev_i], x[i], y[i]);
            if (theta > theta_max) {
                theta_max = theta;
                next_i = i;
            }
        }

        prev_i = curr_i;
        curr_i = next_i;
        if (curr_i != min_i) vertices.push_back(curr_i);
    }

    // for (unsigned i=0; i<vertices.size(); ++i) cout<< i << "\t " << vertices[i] <<endl;
    for (int& i:vertices) cout<< i << "\t " << vertices[i] <<endl;

    can->DrawClone();

    TLine* line = new TLine;
    for (unsigned i=0; i<vertices.size(); ++i) {
        Int_t n1 = vertices[i];
        Int_t n2 = i < vertices.size()-1? vertices[i+1]: vertices[0];
        line->DrawLine(x[n1], y[n1], x[n2], y[n2]);
    }
}
