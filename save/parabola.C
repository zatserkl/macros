#include <TMath.h>
#include <TF1.h>

#include <cassert>
#include <iostream>

using std::cout;    using std::endl;

class RunningSum {
    friend std::ostream& operator << (std::ostream& os, const RunningSum& runningSum);
private:
    Int_t N;                        // capacity of the circular buffer
    Double_t* buf;                  // circular buffer for N elements
    Double_t sum;
    Int_t last;                     // position of the last inserted element
public:
    RunningSum(): N(0), buf(0), sum(0), last(-1) {
        // cout<< "default constructor RunningSum" <<endl;
    }
    RunningSum(Int_t the_N): N(the_N), sum(0), last(N-1) {
        buf = new Double_t[N];
        for (int i=0; i<N; ++i) buf[i] = 0;
    }
    RunningSum(Int_t the_N, Double_t data[]): N(the_N), sum(0), last(N-1) {
        buf = new Double_t[N];
        for (int i=0; i<N; ++i) buf[i] = data[i];
    }
    void Init(Int_t the_N) {
        N = the_N;
        assert(N > 0);              // to avoid calling more than once
        N = the_N;
        buf = new Double_t[N];
        for (int i=0; i<N; ++i) buf[i] = 0;
        last = N - 1;
    }
    void Init(Int_t the_N, Double_t data[]) {
        assert(N > 0);              // to avoid calling more than once
        N = the_N;
        buf = new Double_t[N];
        for (int i=0; i<N; ++i) buf[i] = data[i];
        last = N - 1;
    }
    ~RunningSum() {
        // cout<< "~RunningSum" <<endl;
        delete[] buf;
    }
    void Clear() {
        assert(N > 0);
        sum = 0;
        for (int i=0; i<N; ++i) buf[i] = 0;
        last = N - 1;
    }
    void Put(Double_t element)
    {
        //cout<< "RunningSum::Put: last = " << last <<endl;
        assert(N > 0);
        ++last;                     // in filled buffer this is a position of the first element
        if (last == N) last = 0;

        sum -= buf[last];           // remove value of old (first) element
        sum += element;             // add value of new element

        buf[last] = element;
    }
    RunningSum& operator << (Double_t element) {
        Put(element);
        return *this;
    }
    Double_t Sum() const {return sum;}
    Double_t Update(Double_t element) {
        Put(element);
        return Sum();
    }
};

inline std::ostream& operator << (std::ostream& os, const RunningSum& runningSum) {
    os << "N = " << runningSum.N << " sum = " << runningSum.sum << " last = " << runningSum.last;
    os << " buf:";
    for (int i=0; i<runningSum.N; ++i) os << " " << runningSum.buf[i];
    return os;
}

class ParabolaPar {
    friend std::ostream& operator << (std::ostream&, const ParabolaPar&);
public:
    Int_t N;            // the number of points in the sum
    Double_t X1;        // sum of x
    Double_t X2;        // sum of x*x
    Double_t X3;        // sum of x*x*x
    Double_t X4;        // sum of x*x*x*x
    Double_t Y;         // sum of y
    Double_t YX1;       // sum of y*x
    Double_t YX2;       // sum of y*x*x

    Double_t par[3];

    ParabolaPar() {Clear();}
    void Clear() {
        N = 0;
        X1 = 0;
        X2 = 0;
        X3 = 0;
        X4 = 0;
        Y   = 0;
        YX1 = 0;
        YX2 = 0;
    }
    void Print() const {
        cout<< "N = " << N << " X1 = " << X1 << " X2 = " << X2 << " X3 = " << X3 << " X4 = " << X4 << " Y = " << Y << " YX1 = " << YX1 << " YX2 = " << YX2 <<endl;
    }
    void Init(Int_t np, Double_t x[], Double_t y[])
    {
        Clear();
        for (int i=0; i<np; ++i) {
            ++N;
            X1 += x[i];
            X2 += x[i]*x[i];
            X3 += x[i]*x[i]*x[i];
            X4 += x[i]*x[i]*x[i]*x[i];
            Y   += y[i];
            YX1 += y[i]*x[i];
            YX2 += y[i]*x[i]*x[i];
        }
        Calc();
    }
    void Calc()
    {
        // Variables:

        // N  = Sum(i)           Y   = Sum(y)
        // X1 = Sum(x)           YX1 = Sum(y*x)
        // X2 = Sum(x*x)         YX2 = Sum(y*x*x)
        // X3 = Sum(x*x*x)    
        // X4 = Sum(x*x*x*x)

        // Determinates. NB the same minors for pair of D, D0 and pair of D1, D2.

        //      |N  X1 X2|         |X2    X3|       |X1    X2|       |X1    X2|
        // D  = |X1 X2 X3|   =   N*|        | -  X1*|        | +  X2*|        |
        //      |X2 X3 X4|         |X3    X4|       |X3    X4|       |X2    X3|
                               
        //      |Y   X1  X2|       |X2    X3|       |X1    X2|       |X1    X2|
        // D0 = |YX1 X2  X3| =   Y*|        | - YX1*|        | + YX2*|        |
        //      |YX2 X3  X4|       |X3    X4|       |X3    X4|       |X2    X3|

        //      |N   Y   X2|       |X1   YX1|      |N      Y|      |N      Y|
        // D1 = |X1  YX1 X3| =  X2*|        | - X3*|        | + X4*|        |
        //      |X2  YX2 X4|       |X2   YX2|      |X2   YX2|      |X1   YX1|

        //      |N   X1   Y|
        // D2 = |X1  X2 YX1| =
        //      |X2  X3 YX2|

        //      |N     Y X1|       |X1   YX1|      |N      Y|      |N      Y|
        //  = - |X1  YX1 X2| = -X1*|        | + X2*|        | - X3*|        |
        //      |X2  YX2 X3|       |X2   YX2|      |X2   YX2|      |X1   YX1|


        // cout<< "ParabolaPar::Calc: "; Print();

        // minors

        Double_t x2x3x3x4 = X2*X4 - X3*X3;
        Double_t x1x2x3x4 = X1*X4 - X2*X3;
        Double_t x1x2x2x3 = X1*X3 - X2*X2;

        Double_t x1yx1x2yx2 = X1*YX2 - X2*YX1;
        Double_t Nyx2yx2 = N*YX2 - Y*X2;
        Double_t Nyx1yx1 = N*YX1 - Y*X1;

        // determinants

        Double_t D  = N*x2x3x3x4 -  X1*x1x2x3x4 +  X2*x1x2x2x3;
        Double_t D0 = Y*x2x3x3x4 - YX1*x1x2x3x4 + YX2*x1x2x2x3;
        Double_t D1 =  X2*x1yx1x2yx2 - X3*Nyx2yx2 + X4*Nyx1yx1;
        Double_t D2 = -X1*x1yx1x2yx2 + X2*Nyx2yx2 - X3*Nyx1yx1;

        par[0] = D0 / D;
        par[1] = D1 / D;
        par[2] = D2 / D;
    }
};

inline std::ostream& operator << (std::ostream& os, const ParabolaPar& p) {
    os << "N = " << p.N << " X1 = " << p.X1 << " X2 = " << p.X2 << " X3 = " << p.X3 << " X4 = " << p.X4 << " Y = " << p.Y << " YX1 = " << p.YX1 << " YX2 = " << p.YX2;
    return os;
}

class RunningParabolaParPointer {
public:
    Int_t N;            // the number of points in the sum
    RunningSum* X1;     // sum of x
    RunningSum* X2;     // sum of x*x
    RunningSum* X3;     // sum of x*x*x
    RunningSum* X4;     // sum of x*x*x*x
    RunningSum* Y;      // sum of y
    RunningSum* YX1;    // sum of y*x
    RunningSum* YX2;    // sum of y*x*x

    Double_t par[3];

    RunningParabolaParPointer(Int_t the_N): N(the_N) {
        X1 = new RunningSum(N);
        X2 = new RunningSum(N);
        X3 = new RunningSum(N);
        X4 = new RunningSum(N);
        Y = new RunningSum(N);
        YX1 = new RunningSum(N);
        YX2 = new RunningSum(N);
        //Clear();
    }
    ~RunningParabolaParPointer() {
        delete X1;
        delete X2;
        delete X3;
        delete X4;
        delete Y;
        delete YX1;
        delete YX2;
    }
    void Clear() {
        // N = 0;
        // X1 = 0;
        // X2 = 0;
        // X3 = 0;
        // X4 = 0;
        // Y   = 0;
        // YX1 = 0;
        // YX2 = 0;
    }
    void Init(Int_t np, Double_t x[], Double_t y[])
    {
        Clear();
        N = 0;
        for (int i=0; i<np; ++i) {
            ++N;
            X1->Put(x[i]);
            X2->Put(x[i]*x[i]);
            X3->Put(x[i]*x[i]*x[i]);
            X4->Put(x[i]*x[i]*x[i]*x[i]);
            Y->Put(y[i]);
            YX1->Put(y[i]*x[i]);
            YX2->Put(y[i]*x[i]*x[i]);
        }
        Calc();
    }
    void Put(Double_t x, Double_t y) {
        // X1->Put(x);
        // X2->Put(x*x);
        // X3->Put(x*x*x);
        // X4->Put(x*x*x*x);
        // Y->Put(y);
        // YX1->Put(y*x);
        // YX2->Put(y*x*x);
        *X1 << x;
        *X2 << x*x;
        *X3 << x*x*x;
        *X4 << x*x*x*x;
        *Y << y;
        *YX1 << y*x;
        *YX2 << y*x*x;
        Calc();
    }
    void Calc()
    {
        // Variables:

        // N  = Sum(i)           Y   = Sum(y)
        // X1 = Sum(x)           YX1 = Sum(y*x)
        // X2 = Sum(x*x)         YX2 = Sum(y*x*x)
        // X3 = Sum(x*x*x)    
        // X4 = Sum(x*x*x*x)

        // Determinates. NB the same minors for pair of D, D0 and pair of D1, D2.

        //      |N  X1 X2|         |X2    X3|       |X1    X2|       |X1    X2|
        // D  = |X1 X2 X3|   =   N*|        | -  X1*|        | +  X2*|        |
        //      |X2 X3 X4|         |X3    X4|       |X3    X4|       |X2    X3|
                               
        //      |Y   X1  X2|       |X2    X3|       |X1    X2|       |X1    X2|
        // D0 = |YX1 X2  X3| =   Y*|        | - YX1*|        | + YX2*|        |
        //      |YX2 X3  X4|       |X3    X4|       |X3    X4|       |X2    X3|

        //      |N   Y   X2|       |X1   YX1|      |N      Y|      |N      Y|
        // D1 = |X1  YX1 X3| =  X2*|        | - X3*|        | + X4*|        |
        //      |X2  YX2 X4|       |X2   YX2|      |X2   YX2|      |X1   YX1|

        //      |N   X1   Y|
        // D2 = |X1  X2 YX1| =
        //      |X2  X3 YX2|

        //      |N     Y X1|       |X1   YX1|      |N      Y|      |N      Y|
        //  = - |X1  YX1 X2| = -X1*|        | + X2*|        | - X3*|        |
        //      |X2  YX2 X3|       |X2   YX2|      |X2   YX2|      |X1   YX1|

        Double_t x1 = X1->Sum();
        Double_t x2 = X2->Sum();
        Double_t x3 = X3->Sum();
        Double_t x4 = X4->Sum();
        Double_t y = Y->Sum();
        Double_t yx1 = YX1->Sum();
        Double_t yx2 = YX2->Sum();

        // minors

        Double_t x2x3x3x4 = x2*x4 - x3*x3;
        Double_t x1x2x3x4 = x1*x4 - x2*x3;
        Double_t x1x2x2x3 = x1*x3 - x2*x2;

        Double_t x1yx1x2yx2 = x1*yx2 - x2*yx1;
        Double_t Nyx2yx2 = N*yx2 - y*x2;
        Double_t Nyx1yx1 = N*yx1 - y*x1;

        // determinants

        Double_t D  = N*x2x3x3x4 -  x1*x1x2x3x4 +  x2*x1x2x2x3;
        Double_t D0 = y*x2x3x3x4 - yx1*x1x2x3x4 + yx2*x1x2x2x3;
        Double_t D1 =  x2*x1yx1x2yx2 - x3*Nyx2yx2 + x4*Nyx1yx1;
        Double_t D2 = -x1*x1yx1x2yx2 + x2*Nyx2yx2 - x3*Nyx1yx1;

        par[0] = D0 / D;
        par[1] = D1 / D;
        par[2] = D2 / D;
    }
};

class RunningParabolaPar {
public:
    Int_t N;            // the number of points in the sum
    RunningSum X1;      // sum of x
    RunningSum X2;      // sum of x*x
    RunningSum X3;      // sum of x*x*x
    RunningSum X4;      // sum of x*x*x*x
    RunningSum Y;       // sum of y
    RunningSum YX1;     // sum of y*x
    RunningSum YX2;     // sum of y*x*x

    Double_t par[3];

    RunningParabolaPar(Int_t the_N): N(the_N) {
        X1.Init(N);
        X2.Init(N);
        X3.Init(N);
        X4.Init(N);
        Y.Init(N);
        YX1.Init(N);
        YX2.Init(N);
    }
    RunningParabolaPar(Int_t the_N, Double_t x[], Double_t y[]): N(the_N) {
        // fills the N data points
        // cout<< "RunningParabolaPar::RunningParabolaPar: N = " << N <<endl;
        Init(N, x, y);
    }
    ~RunningParabolaPar() {
    }
    void Init(Int_t np, Double_t x[], Double_t y[])
    {
        X1.Init(N);
        X2.Init(N);
        X3.Init(N);
        X4.Init(N);
        Y.Init(N);
        YX1.Init(N);
        YX2.Init(N);
        for (int i=0; i<np; ++i) {
            X1 << x[i];
            X2 << x[i]*x[i];
            X3 << x[i]*x[i]*x[i];
            X4 << x[i]*x[i]*x[i]*x[i];
            Y << y[i];
            YX1 << y[i]*x[i];
            YX2 << y[i]*x[i]*x[i];
        }
        Calc();
    }
    void Put(Double_t x, Double_t y) {
        // X1->Put(x);
        // X2->Put(x*x);
        // X3->Put(x*x*x);
        // X4->Put(x*x*x*x);
        // Y->Put(y);
        // YX1->Put(y*x);
        // YX2->Put(y*x*x);
        X1 << x;
        X2 << x*x;
        X3 << x*x*x;
        X4 << x*x*x*x;
        Y << y;
        YX1 << y*x;
        YX2 << y*x*x;
        Calc();
    }
    void Calc()
    {
        // Variables:

        // N  = Sum(i)           Y   = Sum(y)
        // X1 = Sum(x)           YX1 = Sum(y*x)
        // X2 = Sum(x*x)         YX2 = Sum(y*x*x)
        // X3 = Sum(x*x*x)    
        // X4 = Sum(x*x*x*x)

        // Determinates. NB the same minors for pair of D, D0 and pair of D1, D2.

        //      |N  X1 X2|         |X2    X3|       |X1    X2|       |X1    X2|
        // D  = |X1 X2 X3|   =   N*|        | -  X1*|        | +  X2*|        |
        //      |X2 X3 X4|         |X3    X4|       |X3    X4|       |X2    X3|
                               
        //      |Y   X1  X2|       |X2    X3|       |X1    X2|       |X1    X2|
        // D0 = |YX1 X2  X3| =   Y*|        | - YX1*|        | + YX2*|        |
        //      |YX2 X3  X4|       |X3    X4|       |X3    X4|       |X2    X3|

        //      |N   Y   X2|       |X1   YX1|      |N      Y|      |N      Y|
        // D1 = |X1  YX1 X3| =  X2*|        | - X3*|        | + X4*|        |
        //      |X2  YX2 X4|       |X2   YX2|      |X2   YX2|      |X1   YX1|

        //      |N   X1   Y|
        // D2 = |X1  X2 YX1| =
        //      |X2  X3 YX2|

        //      |N     Y X1|       |X1   YX1|      |N      Y|      |N      Y|
        //  = - |X1  YX1 X2| = -X1*|        | + X2*|        | - X3*|        |
        //      |X2  YX2 X3|       |X2   YX2|      |X2   YX2|      |X1   YX1|

        Double_t x1 = X1.Sum();
        Double_t x2 = X2.Sum();
        Double_t x3 = X3.Sum();
        Double_t x4 = X4.Sum();
        Double_t y = Y.Sum();
        Double_t yx1 = YX1.Sum();
        Double_t yx2 = YX2.Sum();

        // minors

        Double_t x2x3x3x4 = x2*x4 - x3*x3;
        Double_t x1x2x3x4 = x1*x4 - x2*x3;
        Double_t x1x2x2x3 = x1*x3 - x2*x2;

        Double_t x1yx1x2yx2 = x1*yx2 - x2*yx1;
        Double_t Nyx2yx2 = N*yx2 - y*x2;
        Double_t Nyx1yx1 = N*yx1 - y*x1;

        // determinants

        Double_t D  = N*x2x3x3x4 -  x1*x1x2x3x4 +  x2*x1x2x2x3;
        Double_t D0 = y*x2x3x3x4 - yx1*x1x2x3x4 + yx2*x1x2x2x3;
        Double_t D1 =  x2*x1yx1x2yx2 - x3*Nyx2yx2 + x4*Nyx1yx1;
        Double_t D2 = -x1*x1yx1x2yx2 + x2*Nyx2yx2 - x3*Nyx1yx1;

        par[0] = D0 / D;
        par[1] = D1 / D;
        par[2] = D2 / D;
    }
};

void parabola_plain()
{
    Double_t xmin = -5;
    Double_t xmax = 5;
    TF1* fparabola = new TF1("fparabola", "[0] + [1]*x + [2]*x*x", xmin, xmax);
    fparabola->SetParameters(1., 2., 3.);

    Double_t x[20000];
    Double_t y[20000];

    Int_t N = 10;
    Double_t dx = (xmax - xmin) / N;

    Int_t np = 0;
    for (Double_t xcurr=xmin; xcurr<=xmax; xcurr+=dx) {
        x[np] = xcurr;
        y[np] = fparabola->Eval(xcurr);
        ++np;
    }

    ParabolaPar parabolaPar;
    parabolaPar.Clear();

    parabolaPar.N = N;
    for (int i=0; i<N; ++i) {
        parabolaPar.X1 += x[i];
        parabolaPar.X2 += x[i]*x[i];
        parabolaPar.X3 += x[i]*x[i]*x[i];
        parabolaPar.X4 += x[i]*x[i]*x[i]*x[i];
        parabolaPar.Y   += y[i];
        parabolaPar.YX1 += y[i]*x[i];
        parabolaPar.YX2 += y[i]*x[i]*x[i];
    }

    parabolaPar.Calc();
    for (int ipar=0; ipar<3; ++ipar) {
        cout<< "fparabola->GetParameter(" << ipar << ") = " << fparabola->GetParameter(ipar) << " parabolaPar.par[" << ipar << "] = " << parabolaPar.par[ipar] <<endl;
    }
}

void parabola_ParabolaPar()
{
    Double_t xmin = -5;
    Double_t xmax = 5;
    TF1* fparabola = new TF1("fparabola", "[0] + [1]*x + [2]*x*x", xmin, xmax);
    fparabola->SetParameters(1., 2., 3.);

    Double_t x[20000];
    Double_t y[20000];

    Int_t N = 9;
    Double_t dx = (xmax - xmin) / (N-1);

    Int_t np = 0;
    for (Double_t xcurr=xmin; xcurr<=xmax; xcurr+=dx) {
        x[np] = xcurr;
        y[np] = fparabola->Eval(xcurr);
        ++np;
    }

    ParabolaPar parabolaPar;
    parabolaPar.N = N;

    // for (int i=0; i<parabolaPar.N; ++i) {
    //     parabolaPar.X1 += x[i];
    //     parabolaPar.X2 += x[i]*x[i];
    //     parabolaPar.X3 += x[i]*x[i]*x[i];
    //     parabolaPar.X4 += x[i]*x[i]*x[i]*x[i];
    //     parabolaPar.Y   += y[i];
    //     parabolaPar.YX1 += y[i]*x[i];
    //     parabolaPar.YX2 += y[i]*x[i]*x[i];
    // }

    RunningSum runningX1(N);
    RunningSum runningX2(N);
    RunningSum runningX3(N);
    RunningSum runningX4(N);
    RunningSum runningY(N);
    RunningSum runningYX1(N);
    RunningSum runningYX2(N);

    // fill N elements first time

    for (int i=0; i<N; ++i) {
        runningX1 << x[i];
        runningX2 << x[i]*x[i];
        runningX3 << x[i]*x[i]*x[i];
        runningX4 << x[i]*x[i]*x[i]*x[i];
        runningY << y[i];
        runningYX1 << y[i]*x[i];
        runningYX2 << y[i]*x[i]*x[i];
    }

    parabolaPar.X1 = runningX1.Sum();
    parabolaPar.X2 = runningX2.Sum();
    parabolaPar.X3 = runningX3.Sum();
    parabolaPar.X4 = runningX4.Sum();
    parabolaPar.Y   = runningY.Sum();
    parabolaPar.YX1 = runningYX1.Sum();
    parabolaPar.YX2 = runningYX2.Sum();

    parabolaPar.Calc();
    for (int ipar=0; ipar<3; ++ipar) {
        cout<< "fparabola->GetParameter(" << ipar << ") = " << fparabola->GetParameter(ipar) << " parabolaPar.par[" << ipar << "] = " << parabolaPar.par[ipar] <<endl;
    }

    // cout<< "running sums:" <<endl;
    // cout<< runningX1 <<endl;
    // cout<< runningX2 <<endl;
    // cout<< runningX3 <<endl;
    // cout<< runningX4 <<endl;
    // cout<< runningY <<endl;
    // cout<< runningYX1 <<endl;
    // cout<< runningYX2 <<endl;
}

void parabola()
{
    Double_t xmin = -5;
    Double_t xmax = 5;
    TF1* fparabola = new TF1("fparabola", "[0] + [1]*x + [2]*x*x", xmin, xmax);
    fparabola->SetParameters(1., 2., 3.);
    for (int ipar=0; ipar<3; ++ipar) {
        cout<< "fparabola->GetParameter(" << ipar << ") = " << fparabola->GetParameter(ipar) <<endl;
    }

    Double_t x[20000];
    Double_t y[20000];

    Int_t Ndata = 15;
    //-- Int_t Ndata = 20;
    Double_t dx = (xmax - xmin) / (Ndata-1);

    for (int i=0; i<Ndata; ++i) {
        x[i] = xmin + i*dx;
        y[i] = fparabola->Eval(x[i]);
        cout<< i << "\t" << x[i] << " " << y[i] <<endl;
    }

    Int_t N = 9;

    RunningParabolaPar parabola(N, x, y);

    parabola.Calc();
    cout<< "parabola.par[0] = " << parabola.par[0] << " parabola.par[1] = " << parabola.par[1] << " parabola.par[2] = " << parabola.par[2] <<endl;

    cout<< "running parameters:" <<endl;

    for (int i=N; i<Ndata; ++i) {
        //cout<< i << "\t" << x[i] << " " << y[i] <<endl;
        parabola.Put(x[i], y[i]);
        parabola.Calc();
        cout<< "parabola.par[0] = " << parabola.par[0] << " parabola.par[1] = " << parabola.par[1] << " parabola.par[2] = " << parabola.par[2] <<endl;
    }
}
