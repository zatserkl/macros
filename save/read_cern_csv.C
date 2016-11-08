#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

using std::cout;    using std::endl;

void read_cern_csv(const char* ifname, Int_t chan1=0, Int_t chan2=0, Int_t chan3=0, Int_t chan4=0)
{
    Int_t nsamples = 0;                                 // the number of samples (data points)
    Double_t t[4][20000];                               // time array for 4 channels
    Double_t w[4][20000];                               // pulse height array for 4 channels

    Int_t nchan = 0;
    Int_t channels[4];
    if (chan1) channels[nchan++] = chan1;
    if (chan2) channels[nchan++] = chan2;
    if (chan3) channels[nchan++] = chan3;
    if (chan4) channels[nchan++] = chan4;

    if (nchan == 0) {
        cout<< "Please provide channels, e.g. read_cern_csv(\"run07.txt\", 1,4)" <<endl;
        return;
    }

    std::ifstream ifile(ifname);
    if (!ifile) {
        cout<< "Could not find file " << ifname <<endl;
        return;
    }

    TFile* ofile = new TFile(Form("%s.root",ifname), "recreate");
    TTree* tree = new TTree("wfm", "wfm tree from the CERN csv file");
    tree->SetMarkerStyle(6);
    tree->SetMarkerColor(602);
    tree->SetLineColor(602);

    std::string line;

    // read and discard line KEYSIGHT TECHNOLOGIES,DSO91204A,MY55300102,05.50.0006
    getline(ifile, line);
    // read and discard line with date
    getline(ifile, line);
    // read and discard line with some 3 numbers
    getline(ifile, line);

    while (getline(ifile, line))                                // read and discard line with the event time
    {
        for (int ich=0; ich<nchan; ++ich)
        {
            // read data for the channel
            getline(ifile, line);
            std::replace(line.begin(), line.end(), ',', ' ');   // replace all commas with spaces
            std::stringstream ss(line);
            int nsamples_old = nsamples;
            nsamples = 0;
            while (ss >> w[ich][nsamples]) nsamples++;
            if (nsamples_old > 0) {
                if (nsamples != nsamples_old) {
                    cout<< "***Error event " << tree->GetEntries() << " nsamples = " << nsamples << " while nsamples_old = " << nsamples_old <<endl;
                    break;
                }
            }

            if (tree->GetEntries() == 0 && ich == 0)            // one time action
            {
                // book the tree using the number of samples in the first channel of the first event
                for (int ichBranch=0; ichBranch<nchan; ++ichBranch) {
                    tree->Branch(Form("t%d",channels[ichBranch]), t[ichBranch], Form("t%d[%d]/D",channels[ichBranch],nsamples));
                    tree->Branch(Form("w%d",channels[ichBranch]), w[ichBranch], Form("w%d[%d]/D",channels[ichBranch],nsamples));

                    // fill the time axis
                    for (int isample=0; isample<nsamples; ++isample) t[ichBranch][isample] = isample*0.025;     // 40 GSa/s ==> 25 ps
                }

                cout<< "branches:" <<endl;
                TIter next(tree->GetListOfBranches());
                TBranch* b;
                while ((b = (TBranch*) next())) {
                    cout<< b->GetName() <<endl;
                }
            }
        }

        if (tree->GetEntries() < 10 || tree->GetEntries() % 1000 == 0) cout<< "filling event " << tree->GetEntries() <<endl;
        tree->Fill();
    }

    cout<< "Write " << tree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
    ofile->Write();

    tree->ResetBranchAddresses();
}
