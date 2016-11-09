// Andriy Zatserklyaniy <zatserkl@fnal.gov> Feb 6, 2016

#ifndef macro_wfmWrite_C
#define macro_wfmWrite_C

#include <TROOT.h>
#include <TEnv.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TIterator.h>

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <ctime>

/*
To run from the ROOT prompt for the filename pattern like C2data00530.txt

.L wfmWrite.C
wfmWrite("UCSC_SiPM_28V_D5_580V_W3_750V_3","data","txt", 1,2,3,0, 0,545)

To plot:
wfm->Draw("w2:t2","Entry$==0","pl")
*/

using std::cout;    using std::endl;

class TextObject {
    //
    // To write source code into TTree UserInfo Object.
    // Also use static method TextObject::print(list) to print the content of the list.
    //
public:
    TList* list;
    TextObject(TList* list_0): list(list_0) {
        // write command line
        const char* fname = gEnv->GetValue("Rint.History", "");
        // cout<< "Rint.History file is " << fname <<endl;
        std::ifstream file(fname);
        std::string line = "N/A";
        if (file) {
            std::string line_read;
            while (std::getline(file, line_read)) line = line_read;
        }
        else cout<< "\n***Warning TextObject::TextObject: could not find history file " << fname <<endl<<endl;
        // cout<< "Command line is: " << line <<endl;

        // create an object and write it into the tree
        TNamed* textObject = new TNamed("Command line", line.c_str());
        list->AddLast(textObject);
        file.close();
    }
    void writeFile(const char* fname) {
        std::ifstream file(fname);
        if (file) {
            std::stringstream ss;
            std::string line;
            while (std::getline(file, line)) ss << line << endl;

            // create an object and write it into the tree
            TNamed* textObject = new TNamed(fname, ss.str());
            list->AddLast(textObject);
            file.close();
        }
        else cout<< "\n***Warning TextObject::TextObject: could not find file " << fname <<endl<<endl;
    }
    static void print(TList* the_list)
    {
        cout<< "the number of entries in the list is " << the_list->GetEntries() <<endl;

        TIter next(the_list);
        TNamed* object;
        while ((object = (TNamed*) next())) {
            cout<< "\nobject->GetName() = " << object->GetName() <<endl;
            cout<< object->GetTitle() <<endl;
        }
    }
};

void wfmWrite(const char* dir, const char* name="Trace", const char* ext="txt", Int_t chan1=0, Int_t chan2=0, Int_t chan3=0, Int_t chan4=0, Int_t evt1=0, Int_t evt2=-1, const char* ofname="")
{
#if defined(__CINT__) && !defined(__MAKECINT__)
    cout<< "\n***Warning: This script needs to be compiled. Load it with the plus sign after the name, e.g.\n" <<endl;
    cout<< ".L " << __FILE__ << "+" <<endl;
    return;
#endif

    //  file name format string for scanf with placeholders for the dir, channel, the event number and extension
    //--old-- const char* ifname_format = "%s/C%dTrace%05d.%s";
    const char* ifname_format = "%s/C%d%s%05d.%s";

    Int_t npoints = 0;                                  // the number of samples (data points)
    Double_t t[4][20000];                               // time array for 4 channels
    Double_t w[4][20000];                               // pulse height array for 4 channels

    Int_t nchan = 0;
    Int_t channels[4];
    if (chan1) channels[nchan++] = chan1;
    if (chan2) channels[nchan++] = chan2;
    if (chan3) channels[nchan++] = chan3;
    if (chan4) channels[nchan++] = chan4;

    //
    //  open the first file to read in its header the number of data points (the value of SegmentSize)
    //

    char ifname[1000];
    sprintf(ifname, ifname_format, dir, channels[0], name, evt1, ext);  // build the filename for the first event in the first channel
    cout<< "open the first file " << ifname <<endl;
    std::ifstream ifile_first(ifname);
    if (!ifile_first) {
        cout<< "Could not find input file " << ifname <<endl;
        return;
    }

    // LECROYWP725Zi,50011,Waveform
    // Segments,1,SegmentSize,8002
    // Segment,TrigTime,TimeSinceSegment1
    // #1,07-Mar-2007 07:39:21,0                 
    // Time,Ampl

    std::string line;
    std::getline(ifile_first,line);                     // skip the line LECROYWP725Zi,50011,Waveform
    std::getline(ifile_first,line);                     // read the data line Segments,1,SegmentSize,8002 
    std::replace(line.begin(), line.end(), ',', ' ');   // replace all commas (if any) to spaces
    sscanf(line.c_str(), "Segments %*d SegmentSize %d", &npoints);
    ifile_first.close();
    cout<< "npoints = " << npoints <<endl;

    //
    //  open output ROOT file and create a tree
    //

    TFile* ofile = 0;
    if (ofname && ofname[0]) ofile = new TFile(ofname, "recreate");
    else ofile = new TFile(Form("%s.root",dir), "recreate");

    cout<< "Create output ROOT file " << ofile->GetName() <<endl;

    TTree* tree = new TTree("wfm", "wfm tree");
    tree->SetMarkerStyle(6);
    tree->SetMarkerColor(602);
    tree->SetLineColor(602);

    cout<< "Add to the " << tree->GetName() << " GetUserInfo() source " << __FILE__ <<endl;
    TextObject textObject(tree->GetUserInfo());
    textObject.writeFile(__FILE__);

    // time of the standard time_t type: the number of seconds since Jan 1, 1970 
    Int_t trigtime = 0;                             // time of the current event
    Int_t trigtime_first = 0;                       // time of the first event
    tree->Branch("trigtime", &trigtime, "trigtime/I");

    for (int ich=0; ich<nchan; ++ich) {
        tree->Branch(Form("t%d",channels[ich]), t[ich], Form("t%d[%d]/D",channels[ich],npoints));
        tree->Branch(Form("w%d",channels[ich]), w[ich], Form("w%d[%d]/D",channels[ich],npoints));
    }

    cout<< "branches:" <<endl;
    TIter next(tree->GetListOfBranches());
    TBranch* b;
    while ((b = (TBranch*) next())) {
        cout<< b->GetName() <<endl;
    }

    //
    //  loop over event files
    //

    bool stop = false;
    for (int ievent=evt1; true; ++ievent)                  // NB: an infinite loop
    //-- for (int ievent=0; ievent<10; ++ievent)
    {
        if (evt2 >= evt1 && ievent > evt2) break;
        //-- if (ievent % 100 == 0) cout<< "processing event " << ievent <<endl;
        if (false
            || ievent < 10
            || (ievent < 1000 && ievent % 100 == 0)
            || ievent % 1000 == 0
        ) cout<< "processing event " << ievent <<endl;

        stop = false;
        for (int ich=0; ich<nchan; ++ich)
        {
            sprintf(ifname, ifname_format, dir, channels[ich], name, ievent, ext);    // build the filename
            // cout<< "ifname = " << ifname <<endl;
            std::ifstream ifile(ifname);
            if (!ifile) {
                if (tree->GetEntries() == 0) cout<< "Could not find input file " << ifname <<endl;
                else cout<< "Finished: could not find the next input file " << ifname <<endl;
                stop = true;
                break;
            }

            //
            // read the event time stamp from the header
            //

            // LECROYWP725Zi,50011,Waveform
            // Segments,1,SegmentSize,8002
            // Segment,TrigTime,TimeSinceSegment1
            // #1,07-Mar-2007 07:39:21,0                 
            // Time,Ampl

            std::getline(ifile,line);    // skip the line LECROYWP725Zi,50011,Waveform
            std::getline(ifile,line);    // skip the line Segments,1,SegmentSize,8002
            std::getline(ifile,line);    // skip the line Segment,TrigTime,TimeSinceSegment1
            std::getline(ifile,line);    // read the line #1,07-Mar-2007 07:39:21,0
            if (ich == 0)                                           // use the time stamp from the first channel only
            {
                std::replace(line.begin(), line.end(), ',', ' ');   // replace all commas (if any) to spaces
                tm timeinfo;                                        // standard struct tm. See e.g. http://www.cplusplus.com/reference/ctime/tm/
                char month_str[256];
                sscanf(line.c_str(), "#%*d %2d-%3s-%4d %2d:%2d:%2d", &timeinfo.tm_mday,month_str,&timeinfo.tm_year,&timeinfo.tm_hour,&timeinfo.tm_min,&timeinfo.tm_sec);
                timeinfo.tm_mon = std::string("Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec").find(month_str) / 4;   // convert the three-letter month name into a number 0..11

                trigtime = mktime(&timeinfo);                       // time of the current event
                if (ievent == evt1) trigtime_first = trigtime;      // time of the first event
            }
            std::getline(ifile,line);    // skip the line Time,Ampl

            //
            // read the data
            //

            Int_t nlines = 0;
            for (int i=0; i<npoints; ++i)
            {
                if (std::getline(ifile,line))                           // get line from the file
                {
                    std::replace(line.begin(), line.end(), ',', ' ');   // replace all commas (if any) to spaces
                    std::istringstream ss(line);
                    ss >> t[ich][i] >> w[ich][i];                       // read the space separated data from the line
                    t[ich][i] *= 1e9;                                   // convert the time to ns
                    ++nlines;
                    //-- if (ievent == 0 && i<10) cout<< i << "\t" << t[ich][i] << " " << w[ich][i] <<endl;
                }
                else {
                    cout<< "***Warning: file " << ifname << " has " << nlines << " points instead of npoints = " << npoints <<endl;
                    break;
                }
            }
            ifile.close();
        }

        if (stop) break;
        tree->Fill();
    }

    cout<< "Write " << tree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
    int dt = trigtime - trigtime_first;
    if (dt) cout<< "Event rate = " << tree->GetEntries() << " / " << dt << " = " << tree->GetEntries()/double(dt) << " Hz" <<endl;
    ofile->Write();

    tree->ResetBranchAddresses();
}

#endif  // macro_wfmWrite_C
