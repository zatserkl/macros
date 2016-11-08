// Create a TTree for DRSOsc v.5.0.3 for the evaluation board v4 (saves in tree just one time array for all channels).
// Andriy Zatserklyaniy zatserkl@fnal.gov
// AZ: Nov 1, 2014: fixed channel time array

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <TLeaf.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>
#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cstdio>

using std::cout;      using std::endl;

struct OscTimeHeader {
   UChar_t TIME[4];                 // char word T I M E
   UChar_t BNo[2];                  // B #
   UShort_t board_number;
};

union OscTimeHeaderRecord {
   OscTimeHeader timeHeader;
   Char_t buffer[sizeof(OscTimeHeader)];
};

struct OscTime {
   UChar_t C001[4];                 // char channel number like C 0 0 1
   Float_t t[1024];
   // declare a global operator << as a friend to let it access the private data
   friend std::ostream& operator <<(std::ostream&, const OscTime&);
};
//
std::ostream& operator <<(std::ostream& os, const OscTime& oscTime)
{
    os << "C001 = ";
    for (int i=0; i<4; ++i) os << oscTime.C001[i];
    // os << '\n';
    // for (int i=0; i<1024; ++i) os << i << "\t" << oscTime.t[i] << '\n';
    return os;
}

union OscTimeRecord {
   OscTime oscTime;
   Char_t buffer[sizeof(OscTime)];
};

struct OscEventHeader {
   UChar_t EHDR[4];                 // E H D R
   UInt_t event_number;
   UShort_t year;
   UShort_t month;
   UShort_t day;
   UShort_t hour;
   UShort_t minute;
   UShort_t second;
   UShort_t millisecond;
   UShort_t reserved;
   UChar_t BNo[2];
   UShort_t board_number;
   UChar_t TNo[2];
   UShort_t trigger_cell;
};

struct OscVoltage {
   UChar_t C001[4];                 // char channel number like C 0 0 1
   UShort_t voltage[1024];
};

struct OscEventVoltage {
   OscEventHeader eventHeader;
   OscVoltage oscVoltage[4];        // reserve room for 4 channels
};

union OscEventVoltageRecord {
   OscEventVoltage oscEventVoltage;
   Char_t buffer[sizeof(OscVoltage)];
};

// class OscBin {
// private:
//    std::string ifname;     // to use input file name for messages
//    std::ifstream ifile;
//    // components of the data
//    OscTimeHeaderRecord oscTimeHeaderRecord;
//    OscTimeRecord oscTimeRecord[4];                 // 4 is the maximum number of channels
//    OscEventVoltageRecord oscEventVoltageRecord;
//    std::vector<Int_t> channels;
//    Float_t x[4][1024];
//    Float_t y[4][1024];
// public:
//    // getters
//    /// UInt_t Number() const {return oscEventVoltageRecord.eventVoltage.eventHeader.event_number;}
//    /// UInt_t Year() const {return oscEventVoltageRecord.eventVoltage.eventHeader.year;}
//    /// UInt_t Month() const {return oscEventVoltageRecord.eventVoltage.eventHeader.month;}
//    /// UInt_t Day() const {return oscEventVoltageRecord.eventVoltage.eventHeader.day;}
//    /// UInt_t Hour() const {return oscEventVoltageRecord.eventVoltage.eventHeader.hour;}
//    /// UInt_t Minute() const {return oscEventVoltageRecord.eventVoltage.eventHeader.minute;}
//    /// UInt_t Second() const {return oscEventVoltageRecord.eventVoltage.eventHeader.second;}
//    /// UInt_t Millisecond() const {return oscEventVoltageRecord.eventVoltage.eventHeader.millisecond;}
//    /// const Float_t* Time() const {return oscTimeRecord[0].time.t;}     // all times are the same for V4
//    /// const UShort_t* Voltage(Int_t ich) const {
//    ///    if (ich > 3) return 0;
//    ///    return oscEventVoltageRecord.eventVoltage.oscVoltage[ich].voltage;
//    /// }
//    /// UInt_t Nchan() const {return nchan;}
//    /// const Int_t* UsedChan() const {return usedchan;}
// 
//    // vars
//    Int_t board;
//    bool status;
//    bool operator !() const {return !status;}
//    Long64_t ifsize;
//    Int_t nchan;
//    Int_t usedchan[4];
// 
//    OscBin(): status(false) {}
//    bool Open(const char* ifname_)
//    {
//       ifname = ifname_;
//       status = true;
// 
//       // try to open input file as the oscilloscope application binary file
//       ifile.open(ifname.c_str(), std::ios::binary);
//       if (!ifile) {
//          cout<< "File not found: " << ifname <<endl;
//          status = false;
//          return false;
//       }
//       cout<< "processing file " << ifname <<endl;
// 
//       // file size
//       ifile.seekg(0, std::ios::end);
//       ifsize = ifile.tellg();
//       ifile.seekg(0);
// 
//       if (ifsize < 0) {
//          cout<< "OscBin: input file error: ifsize = " << ifsize <<endl;
//          ifile.close();
//          status = false;
//          return false;
//       }
// 
//       // read the time header
//       ifile.read(oscTimeHeaderRecord.buffer, sizeof(OscTimeHeader));
// 
//       const Char_t header_time[] = {'T', 'I', 'M', 'E'};
//       for (int ichar=0; ichar<4; ++ichar) assert(oscTimeHeaderRecord.timeHeader.TIME[ichar] == header_time[ichar]);
// 
//       board = oscTimeHeaderRecord.timeHeader.board_number;
// 
//       // read the time width data and figure out the number of channels
// 
//       const Char_t header_channel[][4] = {
//          {'C', '0', '0', '1'},
//          {'C', '0', '0', '2'},
//          {'C', '0', '0', '3'},
//          {'C', '0', '0', '4'}
//       };
//       const Char_t header_event[] = {'E', 'H', 'D', 'R'};
// 
//       for (int ichannel=0; ichannel<4; ++ichannel) {
//          // read the channel number
//          unsigned long pos = ifile.tellg();
//          ifile.read(oscTimeRecord[ichannel].buffer, sizeof(OscTime));
//          //cout<< "oscTimeRecord[ichannel].time.C001 = " << oscTimeRecord[ichannel].time.C001 <<endl;
//          cout<< "oscTimeRecord[" << ichannel << "].oscTime = " << oscTimeRecord[ichannel].oscTime <<endl;
//          bool found_channel = false;
//          for (int iheader_channel=0; iheader_channel<4; ++iheader_channel) {
//             if (std::strncmp((const char*)oscTimeRecord[ichannel].oscTime.C001, header_channel[iheader_channel],4) == 0) {
//                found_channel = true;
//                channels.push_back(iheader_channel);
//                break;
//             }
//          }
//          if (!found_channel) {
//             assert(std::strncmp((const char*)oscTimeRecord[ichannel].oscTime.C001, header_event,4) == 0);
//             // time width data finished, this is the EHDR and next data fields
//             ifile.seekg(pos, std::ios::beg);
//             break;
//          }
//       }
// 
//       cout<< "OscBin::Open: found " << channels.size() << " channels in the data" <<endl;
// 
//       assert(channels.size() > 0);
// 
//       UInt_t oscEventVoltageSize = sizeof(OscEventHeader) + channels.size()*sizeof(OscVoltage);
// 
//       //Int_t nevents = 0;
//       //for (nevents=0; ifile.read(oscEventVoltageRecord.buffer, oscEventVoltageSize); ++nevents) {}
// 
//       cout<< "read the first event" <<endl;
// 
//       ifile.read(oscEventVoltageRecord.buffer, oscEventVoltageSize);
//       
//       assert(std::strncmp((const char*)oscEventVoltageRecord.oscEventVoltage.eventHeader.EHDR, header_event, 4) == 0);
// 
//       Int_t tcell = oscEventVoltageRecord.oscEventVoltage.eventHeader.trigger_cell;
//       Float_t tsum = 0;
//       for (int i=0; i<1024; ++i) {
//          Int_t curr = (tcell + i) % 1024;
//          tsum += oscTimeRecord[0].oscTime.t[curr];
//          x[0][i] = tsum;
//          y[0][i] = (oscEventVoltageRecord.oscEventVoltage.oscVoltage[0].voltage[i] - 65535./2)/65535;
//          cout<< i << "\t curr = " << curr << "\t x[0][" << i << "] = " << x[0][i] << "\t y[0][" << i << "] = " << y[0][i] <<endl;
//       }
// 
//       TGraph* g0 = new TGraph(1024, x[0], y[0]);
//       g0->SetMarkerColor(2);
//       new TCanvas;
//       g0->Draw("ap");
// 
//       status = true;
//       return true;
//    }
// };

void ptime(Int_t entry1=0, Int_t entry2=-1, TTree* p=0)
{
   if (!p) p = (TTree*) gDirectory->Get("p");
   if (!p) {
      cout<< "Could not find tree p" <<endl;
      return;
   }

   if (entry2 < entry1) entry2 = p->GetEntries()-1;

   /*
      tm_sec	seconds after the minute	0-61*
      tm_min	minutes after the hour	0-59
      tm_hour	hours since midnight	0-23
      tm_mday	day of the month	1-31
      tm_mon	months since January	0-11
      tm_year	years since 1900	
      tm_wday	days since Sunday	0-6
      tm_yday	days since January 1	0-365
      tm_isdst	Daylight Saving Time flag	
      The Daylight Saving Time flag (tm_isdst) is greater than zero if Daylight Saving Time is in effect, zero if Daylight Saving Time is not in effect, and less than zero if the information is not available.
    * tm_sec is generally 0-59. Extra range to accommodate for leap seconds in certain systems.
    */

   struct std::tm time1;
   time1.tm_isdst = -1;       // set to "not available" otherwise mktime may corrupt the struct

   Int_t nbytes = 0;
   nbytes = p->GetEntry(entry1);
   if (nbytes <= 0) {
      cout<< "Could not load entry " << entry1 <<endl;
      return;
   }

   time1.tm_year = (Int_t) p->GetLeaf("year")->GetValue() - 1900;
   time1.tm_mon = (Int_t) p->GetLeaf("month")->GetValue() - 1;
   time1.tm_mday = (Int_t) p->GetLeaf("day")->GetValue() - 1;
   time1.tm_hour = (Int_t) p->GetLeaf("hour")->GetValue();
   time1.tm_min = (Int_t) p->GetLeaf("minute")->GetValue();
   time1.tm_sec = (Int_t) p->GetLeaf("second")->GetValue();
   Int_t millisecond1 = (Int_t) p->GetLeaf("millisecond")->GetValue();

   printf("Entry %8d date: %02d/%02d/%d %02d:%02d:%02d.%03d\n", entry1, time1.tm_mon,time1.tm_mday,time1.tm_year+1900, time1.tm_hour,time1.tm_min,time1.tm_sec,millisecond1);

   struct std::tm time2;
   time2.tm_isdst = -1;       // set to "not available" otherwise mktime may corrupt the struct

   nbytes = 0;
   nbytes = p->GetEntry(entry2);
   if (nbytes <= 0) {
      cout<< "Could not load entry " << entry2 <<endl;
      return;
   }

   time2.tm_year = (Int_t) p->GetLeaf("year")->GetValue() - 1900;
   time2.tm_mon = (Int_t) p->GetLeaf("month")->GetValue() - 1;
   time2.tm_mday = (Int_t) p->GetLeaf("day")->GetValue() - 1;
   time2.tm_hour = (Int_t) p->GetLeaf("hour")->GetValue();
   time2.tm_min = (Int_t) p->GetLeaf("minute")->GetValue();
   time2.tm_sec = (Int_t) p->GetLeaf("second")->GetValue();
   Int_t millisecond2 = (Int_t) p->GetLeaf("millisecond")->GetValue();

   printf("Entry %8d date: %02d/%02d/%d %02d:%02d:%02d.%03d\n", entry2, time2.tm_mon,time2.tm_mday,time2.tm_year+1900, time2.tm_hour,time2.tm_min,time2.tm_sec,millisecond2);

   // time difference

   time_t abs_time1 = std::mktime(&time1);   // NB: mktime may change the argument: make sure you set time1.tm_isdst = -1;
   time_t abs_time2 = std::mktime(&time2);
   Double_t dtime = (abs_time2 - abs_time1) + 0.001*(millisecond2 - millisecond1);
   //cout<< "time difference is " << dtime << " s" <<endl;
   cout<< "time difference = " << dtime << " s";
   if (entry2 > entry1) cout<< ", trigger rate = " << (entry2 - entry1) / dtime;
   cout<<endl;
}

void osctree(const char* ifname)
{
   // components of the data
   OscTimeHeaderRecord oscTimeHeaderRecord;
   OscTimeRecord oscTimeRecord[4];                 // 4 is the maximum number of channels
   OscEventVoltageRecord oscEventVoltageRecord;
   std::vector<Int_t> channels;

   // field to be saved in the tree
   Float_t x[4][1024];
   Float_t y[4][1024];
   Int_t event;
   Int_t year;
   Int_t month;
   Int_t day;
   Int_t hour;
   Int_t minute;
   Int_t second;
   Int_t millisecond;
   Int_t tc1;

   // try to open input file as the oscilloscope application binary file
   std::ifstream ifile(ifname, std::ios::binary);
   if (!ifile) {
      cout<< "File not found: " << ifname <<endl;
      return;
   }

   // read the time header
   ifile.read(oscTimeHeaderRecord.buffer, sizeof(OscTimeHeader));

   const Char_t header_time[] = {'T', 'I', 'M', 'E'};
   for (int ichar=0; ichar<4; ++ichar) assert(oscTimeHeaderRecord.timeHeader.TIME[ichar] == header_time[ichar]);

   Int_t board = oscTimeHeaderRecord.timeHeader.board_number;
   cout<< "board is " << board <<endl;

   // read the time width data and figure out the number of channels

   const Char_t header_channel[][4] = {
      {'C', '0', '0', '1'},
      {'C', '0', '0', '2'},
      {'C', '0', '0', '3'},
      {'C', '0', '0', '4'}
   };
   const Char_t header_event[] = {'E', 'H', 'D', 'R'};

   for (int ichannel=0; ichannel<4; ++ichannel) {
      // read the channel number
      unsigned long pos = ifile.tellg();
      ifile.read(oscTimeRecord[ichannel].buffer, sizeof(OscTime));
      //cout<< "oscTimeRecord[" << ichannel << "].oscTime = " << oscTimeRecord[ichannel].oscTime <<endl;
      bool found_channel = false;
      for (int iheader_channel=0; iheader_channel<4; ++iheader_channel) {
         if (std::strncmp((const char*)oscTimeRecord[ichannel].oscTime.C001, header_channel[iheader_channel],4) == 0) {
            found_channel = true;
            channels.push_back(iheader_channel);
            break;
         }
      }
      if (!found_channel) {
         assert(std::strncmp((const char*)oscTimeRecord[ichannel].oscTime.C001, header_event,4) == 0);
         // time width data finished, this is the EHDR and next data fields
         ifile.seekg(pos, std::ios::beg);
         break;
      }
   }

   cout<< "found " << channels.size() << " channels in the data" <<endl;

   if (channels.size() == 0) return;

   // create the output tree

   TFile* ofile = new TFile(Form("%s.root",ifname), "recreate");
   TTree* tree = new TTree("p", "DRS4 Oscilloscope tree for version V5");
   tree->SetMarkerStyle(6);
   tree->SetMarkerColor(46);
   tree->SetLineColor(46);

   tree->Branch("board", &board, "board/I");
   tree->Branch("event", &event, "event/I");
   tree->Branch("year", &year, "year/I");
   tree->Branch("month", &month, "month/I");
   tree->Branch("day", &day, "day/I");
   tree->Branch("hour", &hour, "hour/I");
   tree->Branch("minute", &minute, "minute/I");
   tree->Branch("second", &second, "second/I");
   tree->Branch("millisecond", &millisecond, "millisecond/I");
   tree->Branch("tc1", &tc1, "tc1/I");
   tree->Branch("t1", &x[0], "t1[1024]/F");
   // book channels in use only
   for (unsigned ich=0; ich<channels.size(); ++ich) {
      switch (channels[ich]) {
         case 0: tree->Branch("c1", &y[0], "c1[1024]/F");   break;
         case 1: tree->Branch("c2", &y[1], "c2[1024]/F");   break;
         case 2: tree->Branch("c3", &y[2], "c3[1024]/F");   break;
         case 3: tree->Branch("c4", &y[3], "c4[1024]/F");   break;
      }
   }

   UInt_t oscEventVoltageSize = sizeof(OscEventHeader) + channels.size()*sizeof(OscVoltage);

   while(ifile.read(oscEventVoltageRecord.buffer, oscEventVoltageSize))
   {
      assert(std::strncmp((const char*)oscEventVoltageRecord.oscEventVoltage.eventHeader.EHDR, header_event, 4) == 0);
      if (ifile.gcount() != oscEventVoltageSize) {
         cout<< "***Read error. Stop." <<endl;
         break;
      }

      if (tree->GetEntries() % 1000 == 0) cout<< "processing entry " << tree->GetEntries() <<endl;

      board = oscEventVoltageRecord.oscEventVoltage.eventHeader.board_number;
      event = oscEventVoltageRecord.oscEventVoltage.eventHeader.event_number;
      year = oscEventVoltageRecord.oscEventVoltage.eventHeader.year;
      month = oscEventVoltageRecord.oscEventVoltage.eventHeader.month;
      day = oscEventVoltageRecord.oscEventVoltage.eventHeader.day;
      hour = oscEventVoltageRecord.oscEventVoltage.eventHeader.hour;
      minute = oscEventVoltageRecord.oscEventVoltage.eventHeader.minute;
      second = oscEventVoltageRecord.oscEventVoltage.eventHeader.second;
      millisecond = oscEventVoltageRecord.oscEventVoltage.eventHeader.millisecond;
      tc1 = oscEventVoltageRecord.oscEventVoltage.eventHeader.trigger_cell;

      Int_t tcell = oscEventVoltageRecord.oscEventVoltage.eventHeader.trigger_cell;
      for (unsigned ichannel=0; ichannel<channels.size(); ++ichannel) {
         //--old-- Float_t tsum = 0;
         for (int i=0; i<1024; ++i) {
            Int_t curr = (tcell + i) % 1024;
            //--old-- tsum += oscTimeRecord[ichannel].oscTime.t[curr];
            //--old-- x[ichannel][i] = tsum;
            x[channels[ichannel]][i] = 0;
            for (int j=0; j<i; ++j) x[channels[ichannel]][i] += oscTimeRecord[ichannel].oscTime.t[(j+tcell)%1024];
            y[ichannel][i] = (oscEventVoltageRecord.oscEventVoltage.oscVoltage[ichannel].voltage[i]/65536. - 0.5)*1000.;    // mV
         }
      }

      tree->Fill();
   }

   // // plot the first event
   // tree->GetEntry(0);

   // TGraph* g0 = new TGraph(1024, x[0], y[0]);
   // g0->SetMarkerColor(2);
   // new TCanvas;
   // g0->Draw("ap");

   new TCanvas;
   tree->Draw(Form("c%d:t1",channels[0]+1), "Entry$<100");

   ptime();
   cout<< "Writing " << tree->GetEntries() << " entries into file " << ofile->GetName() <<endl;

   ofile->Write();
}
