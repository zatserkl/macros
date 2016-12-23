#include <TTree.h>
#include <TList.h>
#include <iostream>

using std::cout;    using std::endl;

/*
cd("sensl-3x3+31V-3x3+30V-trig2-30mV-0-2171.root")
wfm0 = wfm
cd("sensl-3x3+31V-3x3+30V-trig2-30mV-from2173.root")"
wfm1 = wfm

ofile = new TFile("sensl-3x3+31V-3x3+30V-trig2-30mV-merged.root", "create")
wfm2 = wfmMerge(wfm0, wfm1)
ofile->Write()
*/

TTree* wfmMerge(TTree* wfm0, TTree* wfm1)
{
    TList* list = new TList;
    list->Add(wfm0);
    list->Add(wfm1);
    TTree* wfm2 = TTree::MergeTrees(list);

    // NB: no need to add UserInfo()

    cout<< "wfm0->GetEntries() + wfm1->GetEntries() = wfm2->GetEntries(): " << wfm0->GetEntries() << " + " << wfm1->GetEntries() << " = " << wfm2->GetEntries() <<endl;

    return wfm2;
}
