// Andriy Zatserklyaniy, zatserkl@fnal.gov

// DRS4 v4

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TObjArray.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>
#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <cstring>

using std::cout;      using std::endl;

struct OscTime {
   Char_t header[4];
   UInt_t number;
   UShort_t year;
   UShort_t month;
   UShort_t day;
   UShort_t hour;
   UShort_t minute;
   UShort_t second;
   UShort_t millisecond;
   UShort_t reserved;
   Float_t b1_t[1024];
};

struct OscChannel {
   Char_t header[4];
   UShort_t voltage[1024];
};

struct OscEvent {
   OscTime time;
   OscChannel channel[4];
};

union OscRecord {
   OscEvent event;
   Char_t buffer[sizeof(OscEvent)];
   OscRecord() {for (unsigned i=0; i<sizeof(OscEvent); ++i) buffer[i] = 0;}
};

class OscBin {
private:
   std::string ifname;     // to use input file name for messages
   std::ifstream ifile;
   OscRecord oscRecord;
public:
   // getters
   UInt_t Number() const {return oscRecord.event.time.number;}
   UInt_t Year() const {return oscRecord.event.time.year;}
   UInt_t Day() const {return oscRecord.event.time.day;}
   UInt_t Hour() const {return oscRecord.event.time.hour;}
   UInt_t Minute() const {return oscRecord.event.time.minute;}
   UInt_t Second() const {return oscRecord.event.time.second;}
   UInt_t Millisecond() const {return oscRecord.event.time.millisecond;}
   const Float_t* Time() const {return oscRecord.event.time.b1_t;}
   const UShort_t* Voltage(Int_t ich) const {
      if (ich > 3) return 0;
      return oscRecord.event.channel[ich].voltage;
   }
   UInt_t Nchan() const {return nchan;}
   const Int_t* UsedChan() const {return usedchan;}

   // vars
   bool status;
   bool operator !() const {return !status;}
   Long64_t ifsize;
   Int_t nchan;
   Int_t usedchan[4];

   OscBin(): status(false) {}
   bool Open(const char* ifname_)
   {
      // use like:
      //
      // OscBin oscBin;
      // oscBin.Open("osc.bin");
      // if (!oscBin) {
      //    cout<< "Input file does not seem to be the DRS4 oscilloscope application binary file" <<endl;
      //    return;
      // }

      ifname = ifname_;
      status = true;

      // try to open input file as the oscilloscope application binary file
      ifile.open(ifname.c_str(), std::ios::binary);
      if (!ifile) {
         cout<< "File not found: " << ifname <<endl;
         status = false;
         return false;
      }
      cout<< "processing file " << ifname <<endl;

      // file size
      ifile.seekg(0, std::ios::end);
      ifsize = ifile.tellg();
      ifile.seekg(0);

      if (ifsize < 0) {
         cout<< "OscBin: input file error: ifsize = " << ifsize <<endl;
         ifile.close();
         status = false;
         return false;
      }

      // read the first event to figure out the number of channels

      const Char_t header_event[] = {'E', 'H', 'D', 'R'};
      const Char_t header_channel[][4] = {
         {'C', '0', '0', '1'},
         {'C', '0', '0', '2'},
         {'C', '0', '0', '3'},
         {'C', '0', '0', '4'}
      };

      Long64_t size = sizeof(OscTime) + 4*sizeof(OscChannel);
      while (ifsize < size) size -= sizeof(OscChannel); // min file size: sizeof(OscTime) + 1*sizeof(OscChannel)

      if (ifsize < Long64_t(sizeof(OscTime) + sizeof(OscChannel))) {
         cout<< "OscBin: input file " << ifname << " is too short to contain even single DRS4 osclloscope event" <<endl;
         ifile.close();
         status = false;
         return false;
      }

      ifile.read(oscRecord.buffer, size);

      // check the event header
      if (std::strncmp(oscRecord.event.time.header, header_event, 4) != 0) {
         status = false;
         cout<< "OscBin: No event header found, the input file " << ifname << " does not seem to be the DRS4 osclloscope application binary file" <<endl;
         return false;
      }

      // read channels' headers to figure out how many them and their numbers
      for (nchan=0; nchan<4; ++nchan)
      {
         Int_t channel = 0;
         for (int ich=0; ich<4; ++ich) {
            if (std::strncmp(oscRecord.event.channel[nchan].header, header_channel[ich], 4)==0) {
               channel = ich+1;
               usedchan[ich] = channel;
               // cout<< "OscBin::Open: channel = " << channel << " usedchan[" << ich << "] = " << usedchan[ich] <<endl;
            }
         }
         if (channel == 0) {
            // this is not a channel header
            break;
         }
      }
      if (nchan == 0) {
         status = false;
         cout<< "OscBin: No channels found, the input file " << ifname << " does not seem to be the DRS4 osclloscope application binary file" <<endl;
         return false;
      }

      // cout<< "OscBin: oscilloscope event contains " << nchan << " channels" <<endl;

      // return to the beginning of the file
      ifile.seekg(0);

      status = true;
      return true;
   }
   bool SetEventPointer(Int_t record)
   {
      // set file get pointer to read record record (counting from 0)
      if (!status) return false;
      Long64_t record_size = sizeof(OscEvent);
      if (record_size*(record+1) <= ifsize) ifile.seekg(record*record_size);
      else {
         cout<< "File " << ifname << " has no event " << record <<endl;
         return false;
      }
      // cout<< "OscBin::SetEventPointer: ifile.tellg() = " << ifile.tellg() <<endl;
      return true;
   }
   Long64_t GetEventPointer()
   {
      // set file get pointer to read record record (counting from 0)
      if (!status) return false;
      Long64_t pos = ifile.tellg();
      // cout<< "OscBin::GetEventPointer: ifile.tellg() = " << pos <<endl;
      return pos;
   }
   bool ReadEvent()
   {
      // reads current event

      if (!status) return false;
      Long64_t record_size = sizeof(OscTime) + nchan*sizeof(OscChannel);
      Long64_t curr = ifile.tellg();
      if (curr < 0) {
         cout<< "OscBin::ReadNext: read error: curr = " << curr <<endl;
         return false;
      }
      if (ifsize < curr+record_size) return false;
      else {
         ifile.read(oscRecord.buffer, record_size);
         // check data consistency
         const Char_t header_event[] = {'E', 'H', 'D', 'R'};
         const Char_t header_channel[][4] = {
            {'C', '0', '0', '1'},
            {'C', '0', '0', '2'},
            {'C', '0', '0', '3'},
            {'C', '0', '0', '4'}
         };
         // event header
         if (std::strncmp(oscRecord.event.time.header, header_event, 4) != 0) {
            cout<< "OscBin::ReadNext: corrupted event header" <<endl;
            return false;
         }
         // channel headers
         for (int ich=0; ich<nchan; ++ich) {
            if (std::strncmp(oscRecord.event.channel[ich].header, header_channel[ich], 4) != 0) {
               cout<< "OscBin::ReadNext: corrupted channel header for ich = " << ich <<endl;
               return false;
            }
         }
      }
      return true;
   }
};

/*
TODO: 
1) exect EHDR and C001,C002, etc.
2) read complete 4-bytes word and union with two 2-bytes words
3) exceptions for readout error
4) easy way to set brach addresses to read the tree (MakeClass may work)
   -- how to do with different number of channels?
*/

struct PulseBuffer
{
   Int_t event;
   Int_t year;
   Int_t month;
   Int_t day;
   Int_t hour;
   Int_t minute;
   Int_t second;
   Int_t millisecond;
   Int_t tc1;              // not in use, for compatibility with pulse tree only
   Float_t b1_t[1024];
   Float_t b1_c1[1024];
   Float_t b1_c2[1024];
   Float_t b1_c3[1024];
   Float_t b1_c4[1024];
   Int_t tc2;              // not in use, for compatibility with pulse tree only
   Float_t b2_t[1024];
   Float_t b2_c1[1024];
   Float_t b2_c2[1024];
   Float_t b2_c3[1024];
   Float_t b2_c4[1024];
   Int_t usedchan[4];
   TBranch* b_event;
   TBranch* b_year;
   TBranch* b_month;
   TBranch* b_day;
   TBranch* b_hour;
   TBranch* b_minute;
   TBranch* b_second;
   TBranch* b_millisecond;
   TBranch* b_tc1;              // not in use, for compatibility with pulse tree only
   TBranch* b_b1_t;
   TBranch* b_b1_c1;
   TBranch* b_b1_c2;
   TBranch* b_b1_c3;
   TBranch* b_b1_c4;
   TBranch* b_tc2;              // not in use, for compatibility with pulse tree only
   TBranch* b_b2_t;
   TBranch* b_b2_c1;
   TBranch* b_b2_c2;
   TBranch* b_b2_c3;
   TBranch* b_b2_c4;
   TBranch* b_usedchan;
   std::map<TBranch*, void*> baddr;
   void clear() {
      event = 0;
      year = 0;
      month = 0;
      day = 0;
      hour = 0;
      minute = 0;
      second = 0;
      millisecond = 0;
      tc1 = 0;
      for (int ich=0; ich<4; ++ich) usedchan[ich] = 0;
      for (int ipoint=0; ipoint<1024; ++ipoint) {
         b1_t[ipoint] = 0;
         b1_c1[ipoint] = 0;
         b1_c2[ipoint] = 0;
         b1_c3[ipoint] = 0;
         b1_c4[ipoint] = 0;
         b2_t[ipoint] = 0;
         b2_c1[ipoint] = 0;
         b2_c2[ipoint] = 0;
         b2_c3[ipoint] = 0;
         b2_c4[ipoint] = 0;
      }
   }
   ~PulseBuffer() {
   }
   PulseBuffer():
      b_event(0)
      , b_year(0)
      , b_month(0)
      , b_day(0)
      , b_hour(0)
      , b_minute(0)
      , b_second(0)
      , b_millisecond(0)
      , b_tc1(0)              // not in use, for compatibility with pulse tree only
      , b_b1_t(0)
      , b_b1_c1(0)
      , b_b1_c2(0)
      , b_b1_c3(0)
      , b_b1_c4(0)
      , b_tc2(0)              // not in use, for compatibility with pulse tree only
      , b_b2_t(0)
      , b_b2_c1(0)
      , b_b2_c2(0)
      , b_b2_c3(0)
      , b_b2_c4(0)
      , b_usedchan(0)
   {
      clear();
   }
   PulseBuffer(const PulseBuffer& buf) {
      operator =(buf);
   }
   PulseBuffer& operator =(const PulseBuffer& buf) {
      if (&buf == this) return *this;
      event = buf.event;
      year = buf.year;
      month = buf.month;
      day = buf.day;
      hour = buf.hour;
      minute = buf.minute;
      second = buf.second;
      millisecond = buf.millisecond;
      for (int ich=0; ich<4; ++ich) usedchan[ich] = buf.usedchan[ich];
      for (int ipoint=0; ipoint<1024; ++ipoint) {
         b1_t[ipoint] = buf.b1_t[ipoint];
         b1_c1[ipoint] = buf.b1_c1[ipoint];
         b1_c2[ipoint] = buf.b1_c2[ipoint];
         b1_c3[ipoint] = buf.b1_c3[ipoint];
         b1_c4[ipoint] = buf.b1_c4[ipoint];
         b2_t[ipoint] = buf.b2_t[ipoint];
         b2_c1[ipoint] = buf.b2_c1[ipoint];
         b2_c2[ipoint] = buf.b2_c2[ipoint];
         b2_c3[ipoint] = buf.b2_c3[ipoint];
         b2_c4[ipoint] = buf.b2_c4[ipoint];
      }
      b_event = buf.b_event;
      b_year = buf.b_year;
      b_month = buf.b_month;
      b_day = buf.b_day;
      b_hour = buf.b_hour;
      b_minute = buf.b_minute;
      b_second = buf.b_second;
      b_millisecond = buf.b_millisecond;
      b_tc1 = buf.b_tc1;
      b_b1_t = buf.b_b1_t;
      b_b1_c1 = buf.b_b1_c1;
      b_b1_c2 = buf.b_b1_c2;
      b_b1_c3 = buf.b_b1_c3;
      b_b1_c4 = buf.b_b1_c4;
      b_tc2 = buf.b_tc2;
      b_b2_t = buf.b_b2_t;
      b_b2_c1 = buf.b_b2_c1;
      b_b2_c2 = buf.b_b2_c2;
      b_b2_c3 = buf.b_b2_c3;
      b_b2_c4 = buf.b_b2_c4;
      b_usedchan = buf.b_usedchan;
      return *this;
   }
   void UsedChan() const {
      cout<< "DRS4 channels: ";
      for (int ich=0; ich<4; ++ich) {
         if (usedchan[ich]) cout<< usedchan[ich] << "  ";
      }
      cout<<endl;
   }
   void Book(TTree* tree)
   {
      // books branches of the new tree
      b_event = tree->Branch("event", &event, "event/I");
      b_year = tree->Branch("year", &year, "year/I");
      b_month = tree->Branch("month", &month, "month/I");
      b_day = tree->Branch("day", &day, "day/I");
      b_hour = tree->Branch("hour", &hour, "hour/I");
      b_minute = tree->Branch("minute", &minute, "minute/I");
      b_second = tree->Branch("second", &second, "second/I");
      b_millisecond = tree->Branch("millisecond", &millisecond, "millisecond/I");
      b_tc1 = tree->Branch("tc1", &tc1, "tc1/I");
      b_b1_t = tree->Branch("b1_t", &b1_t, "b1_t[1024]/F");
      // channel in use: from 1 to 4 (its number), not in use: 0
      b_usedchan = tree->Branch("usedchan", &usedchan, "usedchan[4]/I");
      // book used channels only
      if (usedchan[0] > 0) b_b1_c1 = tree->Branch("b1_c1", &b1_c1, "b1_c1[1024]/F");
      if (usedchan[1] > 0) b_b1_c2 = tree->Branch("b1_c2", &b1_c2, "b1_c2[1024]/F");
      if (usedchan[2] > 0) b_b1_c3 = tree->Branch("b1_c3", &b1_c3, "b1_c3[1024]/F");
      if (usedchan[3] > 0) b_b1_c4 = tree->Branch("b1_c4", &b1_c4, "b1_c4[1024]/F");
   }
   void Connect(TTree* tree)
   {
      clear();
      //SaveBranchAddresses(tree);
      // connects tree buffers with variables to use for event-by-event analysis
      // rest of the event
      cout<< "&event = " << &event <<endl;
      if ((b_event = tree->GetBranch("event"))) b_event->SetAddress(&event);
      if ((b_year = tree->GetBranch("year"))) b_year->SetAddress(&year);
      if ((b_month = tree->GetBranch("month"))) b_month->SetAddress(&month);
      if ((b_day = tree->GetBranch("day"))) b_day->SetAddress(&day);
      if ((b_hour = tree->GetBranch("hour"))) b_hour->SetAddress(&hour);
      if ((b_minute = tree->GetBranch("minute"))) b_minute->SetAddress(&minute);
      if ((b_second = tree->GetBranch("second"))) b_second->SetAddress(&second);
      if ((b_millisecond = tree->GetBranch("millisecond"))) b_millisecond->SetAddress(&millisecond);
      if ((b_tc1 = tree->GetBranch("tc1"))) b_tc1->SetAddress(&tc1);
      if ((b_b1_t = tree->GetBranch("b1_t"))) b_b1_t->SetAddress(&b1_t);
      // channel in use: from 1 to 4 (its number), not in use: 0
      if ((b_usedchan = tree->GetBranch("usedchan"))) b_usedchan->SetAddress(&usedchan);
      // connect used channels only
      if ((b_b1_c1 = tree->GetBranch("b1_c1"))) b_b1_c1->SetAddress(&b1_c1);
      if ((b_b1_c2 = tree->GetBranch("b1_c2"))) b_b1_c2->SetAddress(&b1_c2);
      if ((b_b1_c3 = tree->GetBranch("b1_c3"))) b_b1_c3->SetAddress(&b1_c3);
      if ((b_b1_c4 = tree->GetBranch("b1_c4"))) b_b1_c4->SetAddress(&b1_c4);
      if ((b_tc1 = tree->GetBranch("tc1"))) b_tc1->SetAddress(&tc1);
      if ((b_b2_t = tree->GetBranch("b2_t"))) b_b2_t->SetAddress(&b2_t);
      // channel in use: from 1 to 4 (its number), not in use: 0
      if ((b_usedchan = tree->GetBranch("usedchan"))) b_usedchan->SetAddress(&usedchan);
      // connect used channels only
      if ((b_b2_c1 = tree->GetBranch("b2_c1"))) b_b2_c1->SetAddress(&b2_c1);
      if ((b_b2_c2 = tree->GetBranch("b2_c2"))) b_b2_c2->SetAddress(&b2_c2);
      if ((b_b2_c3 = tree->GetBranch("b2_c3"))) b_b2_c3->SetAddress(&b2_c3);
      if ((b_b2_c4 = tree->GetBranch("b2_c4"))) b_b2_c4->SetAddress(&b2_c4);
   }
   void Disconnect(TTree* tree) {
      tree->ResetBranchAddresses();
   }
   TGraph* plot(Int_t ch=0) const      // NB: init to non-existant channel
   {
      if (ch == 0) {
         cout<< "Usage: plot(channel)" <<endl;
         UsedChan();
         return 0;
      }

      const Int_t color[] = {2, 4, 6, 8};

      const Float_t* v = 0;
      if (ch == 1 and usedchan[0]) v = b1_c1;
      if (ch == 2 and usedchan[1]) v = b1_c2;
      if (ch == 3 and usedchan[2]) v = b1_c3;
      if (ch == 4 and usedchan[3]) v = b1_c4;
      if (!v) {
         cout<< "Channel not found: " << ch <<endl;
         UsedChan();
         return 0;
      }
      TGraph* gr = new TGraph(1024,b1_t,v);
      gr->SetNameTitle(Form("gr_evt_%d_chan_%d",event,ch), Form("gr_evt_%d_chan_%d;time, ns;V",event,ch));
      gr->SetMarkerStyle(6);
      gr->SetMarkerColor(color[ch]);
      gr->Draw("ap");
      return gr;
   }
};

class Osc {
public:
   static PulseBuffer buf;
   std::string ifname;
   TTree* tree;
   OscBin oscBin;                // oscilloscope binary file reader

   Osc(TTree* tree_): tree(tree_)
   {
      TFile* ifile = tree->GetCurrentFile();
      if (ifile) ifname = ifile->GetName();

      buf.Connect(tree);                 // connect tree buffers with static fields of class Osc
      tree->GetEntry(0);

      cout<< "DRS4 channels: ";
      for (int ich=0; ich<4; ++ich) {
         if (buf.usedchan[ich]) cout<< buf.usedchan[ich] << "  ";
      }
      cout<<endl;
   }

   Osc(const char* ifname_): ifname(ifname_), tree(0)
   {
      if (oscBin.Open(ifname.c_str())) {
         cout<< "Successfully opened DRS4 oscilloscope binary file " << ifname.c_str() <<endl;
         return;
      }

      // try to open as a root file
      TFile* ifile = TFile::Open(ifname.c_str());
      if (!ifile) {
         cout<< "Not a root file: " << ifname <<endl;
         return;
      }
      tree = (TTree*) ifile->Get("pulse");
      if (!tree) {
         cout<< "There is no DRS4 oscilloscope tree in the file " << ifname <<endl;
         return;
      }
      cout<< "Found DRS4 oscilloscope tree \"pulse\" in the ROOT file " << ifname <<endl;
      cout<< "The number of entries is " << tree->GetEntries() <<endl;

      buf.Connect(tree);                 // connect tree buffers with static fields of class Osc
      tree->GetEntry(0);

      cout<< "DRS4 channels: ";
      for (int ich=0; ich<4; ++ich) {
         if (buf.usedchan[ich]) cout<< buf.usedchan[ich] << "  ";
      }
      cout<<endl;
   }

   void Convert(Int_t entry1=0, Int_t entry2=0)
   {
      if (tree) return;
      if (!oscBin) return;

      cout<< "DRS4 channels: ";
      for (int ich=0; ich<4; ++ich) {
         buf.usedchan[ich] = oscBin.UsedChan()[ich];
         if (buf.usedchan[ich]) cout<< buf.usedchan[ich] << "  ";
      }
      cout<<endl;

      TFile* ofile = TFile::Open(Form("%s.root",ifname.c_str()), "recreate");

      tree = new TTree("pulse", "pulse tree for DRS4 oscilloscope data");
      tree->SetMarkerStyle(6);
      tree->SetMarkerColor(46);
      tree->SetLineColor(46);

      buf.Book(tree);                             // book tree

      Int_t ientry = entry1;
      oscBin.SetEventPointer(ientry);

      while (oscBin.ReadEvent())
      {
         if (tree->GetEntries() % 100 == 0) cout<< "processing ientry " << ientry << " oscBin.Number = " << oscBin.Number() <<endl;

         for (int ich=0; ich<4; ++ich) {
            buf.usedchan[ich] = oscBin.UsedChan()[ich];
            // cout<< ich << "\b1_t " << usedchan[ich] <<endl;
         }
         buf.event = oscBin.Number();
         buf.year = oscBin.Year();
         buf.day = oscBin.Day();
         buf.hour = oscBin.Hour();
         buf.minute = oscBin.Minute();
         buf.second = oscBin.Second();
         buf.millisecond = oscBin.Millisecond();
         for (int ipoint=0; ipoint<1024; ++ipoint) buf.b1_t[ipoint] = oscBin.Time()[ipoint];
         for (unsigned ich=0; ich<oscBin.Nchan(); ++ich) {
            Float_t* v = 0;
            Int_t channel = oscBin.UsedChan()[ich];
            if (channel > 0) {
               if (channel == 1) v = buf.b1_c1;
               if (channel == 2) v = buf.b1_c2;
               if (channel == 3) v = buf.b1_c3;
               if (channel == 4) v = buf.b1_c4;
            }
            if (!v) {
               cout<< "something wrong: v == 0" <<endl;
               return;
            }
            for (int ipoint=0; ipoint<1024; ++ipoint) {
               v[ipoint] = (Float_t(oscBin.Voltage(ich)[ipoint]) - 32767.5) / 65536.;
            }
         }

         tree->Fill();
         if (entry2 > 0 && ientry >= entry2) break;
         ++ientry;
      }

      cout<< "write " << tree->GetEntries() << " entries into file " << ofile->GetName() <<endl;
      ofile->Write();

      new TCanvas;
      tree->Draw("b1_c1:b1_t","Entry$<100");
   }
};

void showbra(TTree* tree) {
   for (int i=0; i<tree->GetListOfBranches()->GetEntries(); ++i) {
      TBranch* branch = (TBranch*) tree->GetListOfBranches()->At(i);
      void* address = branch->GetAddress();
      cout<< "branch " << branch->GetName() << " address: " << address <<endl;
   }
}

// define static members
PulseBuffer Osc::buf;

TTree* pulsetree(TTree* tree)
{
   Osc* osc = new Osc(tree);
   return osc->tree;
}

TTree* pulsetree(const char* ifname, Int_t entry1=0, Int_t entry2=0)
{
   Osc* osc = new Osc(ifname);
   osc->Convert(entry1,entry2);
   return osc->tree;
}
