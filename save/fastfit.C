// Andriy Zatserklyaniy <zatserkl@gmail.com> May 22, 2017

#include <TRandom3.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <iostream>
#include <cassert>

using std::cout;    using std::endl;

void FastFitMerge(Double_t *pb, Double_t *pc, Double_t *pa)
{
    // in fact, I don't need the data points: just parameters!

    // get y at two points: x = 0 and x = 1
    Double_t x0 = 0.;
    Double_t x1 = 1.;
    Double_t yb0 = pb[0];           // + pb[1]*0
    Double_t yb1 = pb[0] + pb[1];   // pb[1] = pb[1]*1
    Double_t yc0 = pc[0];           // + pc[1]*0
    Double_t yc1 = pc[0] + pc[1];   // pc[1] = pc[1]*1

    Double_t ybc0 = 0.5*(yb0 + yc0);
    Double_t ybc1 = 0.5*(yb1 + yc1);

    pa[1] = (ybc1 - ybc0) / (x1 - x0);
    pa[0] = 0.5*(ybc0 + ybc1 - pa[1]*(x0 + x1));

    cout<< "pb[0] = " << pb[0] << " pb[1] = " << pb[1] <<endl;
    cout<< "pc[0] = " << pc[0] << " pc[1] = " << pc[1] <<endl;
    cout<< "pa[0] = " << pa[0] << " pa[1] = " << pa[1] <<endl;
    cout<<endl;

    // TODO: free memory
}

void FastFitSplitMerge(Int_t na, Double_t xa[], Double_t ya[], Double_t *pa)
{
    // Fit straight line: y = p[0] + p[1]*x

    if (na <= 2) {
        // assume that na = 2
        assert(na == 2);

        pa[1] = (ya[1] - ya[0]) / (xa[1] - xa[0]);
        pa[0] = 0.5*(ya[0] + ya[1] - pa[1]*(xa[0] + xa[1]));

        cout<< "xa[0] = " << xa[0] << " xa[1] = " << xa[1] << " ya[0] = " << ya[0] << " ya[1] = " << ya[1] << " pa[0] = " << pa[0] << " pa[1] = " << pa[1] <<endl;
        return;
    }

    static Int_t nalloc = 0;
    static Int_t nfreed = 0;

    Int_t nb = na / 2;
    Int_t nc = na - nb;
    ++nalloc;
    cout<< "allocate memory nb = " << nb << " nc = " << nc  << " nalloc = " << nalloc << " nfreed = " << nfreed <<endl;
    Double_t *xb = new Double_t[nb];
    Double_t *yb = new Double_t[nb];
    Double_t *xc = new Double_t[nc];
    Double_t *yc = new Double_t[nc];

    Double_t *pb = new Double_t[2];     // two parameters
    Double_t *pc = new Double_t[2];     // two parameters

    /// // Not the top performance because of small base for neighboring points
    ///
    /// for (int i=0; i<nb; ++i) {
    ///     xb[i] = xa[i];
    ///     yb[i] = ya[i];
    /// }
    /// for (int i=0; i<nc; ++i) {
    ///     xc[i] = xa[nb+i];
    ///     yc[i] = ya[nb+i];
    /// }

    for (int i=0, nb = nc = 0; i<na; ++i) {
        if (i % 2 == 0) {
            xb[nb] = xa[i];
            yb[nb] = ya[i];
            nb++;
        }
        else {
            xc[nc] = xa[i];
            yc[nc] = ya[i];
            nc++;
        }
    }

    FastFitSplitMerge(nb, xb, yb, pb);
    FastFitSplitMerge(nc, xc, yc, pc);
    FastFitMerge(pb, pc, pa);

    // free memory
    ++nfreed;
    cout<< "free memory nb = " << nb << " nc = " << nc << " nalloc = " << nalloc << " nfreed = " << nfreed <<endl;
    delete[] xb;
    delete[] yb;
    delete[] xc;
    delete[] yc;
    delete[] pb;
    delete[] pc;
}

void fastfit()
{
    // Linear fit

    Double_t x[100000];
    Double_t y[100000];
    Double_t ey[100000];
    Int_t np = 16;

    Double_t p[10];
                        // y = p[0] + p[1]*x
    p[0] = 0;           // inercept
    p[1] = 0.5;         // slope

    TRandom3 random;
    random.SetSeed(1);

    Double_t sigma = 1;

    for (int i=0; i<np; ++i) {
        x[i] = i;
        y[i] = p[0] + p[1]*x[i];
        y[i] += random.Gaus(y[i], sigma);
        ey[i] = sigma;
    }

    TGraphErrors *g = new TGraphErrors(np, x, y, 0, ey);
    g->SetNameTitle(
         "g",
         Form("y = p_{0} + p_{1}#upointx + #epsilon(#sigma), p_{0} = %0.2f, p_{1} = %0.2f, #sigma = %0.2f",
              p[0],p[1],sigma)
         );
    g->SetMarkerStyle(20);
    g->SetMarkerColor(9);
    g->SetLineColor(9);

    new TCanvas;
    g->Draw("ap");

    g->Fit("pol1");

    Double_t par[10];
    FastFitSplitMerge(np, x, y, par);
    for (int i=0; i<2; ++i) cout<< "par[" << i << "] = " << par[i] <<endl;

    TF1 *fline = new TF1("fline", "[0] + [1]*x", x[0], x[np-1]);
    fline->SetLineColor(8);
    fline->SetParameters(par[0], par[1]);

    fline->Draw("same");
}
