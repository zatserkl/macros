/*
 
   Name:           drstree.C
   Created by:     Stefan Ritt <stefan.ritt@psi.ch>
   Date:           July 30th, 2014
   Modified by:    Andriy Zatserklyaniy <zatserkl@fnal.gov> Nov 7, 2014: Create branches for active channels only
 
   Purpose:        Example program under ROOT to read a binary data file written 
                   by the DRSOsc program. Decode time and voltages from waveforms 
                   and display them as a graph. Put values into a ROOT Tree for 
                   further analysis.
 
                   To run it, do:
 
                   - Crate a file test.dat via the "Save" button in DRSOsc
                   - start ROOT
                   root [0] .L drstree.C+
                   root [1] drstree("test.dat");
 
*/
 

#include <cstring>
#include <cstdio>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "Getline.h"

typedef struct {
  char           time_header[4];
  char           bn[2];
  unsigned short board_serial_number;
} THEADER;

typedef struct {
  char           event_header[4];
  unsigned int   event_serial_number;
  unsigned short year;
  unsigned short month;
  unsigned short day;
  unsigned short hour;
  unsigned short minute;
  unsigned short second;
  unsigned short millisecond;
  unsigned short range;
  char           bs[2];
  unsigned short board_serial_number;
  char           tc[2];
  unsigned short trigger_cell;
} EHEADER;

/*-----------------------------------------------------------------------------*/

void drstree(const char *filename, int nevents_max=1000000) {
   THEADER th;
   EHEADER eh;
   char hdr[4];
   unsigned short voltage[1024];
   double waveform[4][1024], time[4][1024];
   float bin_width[4][1024];
   char rootfile[1024];
   int i, j, ch, n, chn_index;
   double t1, t2, dt;

   // open the binary waveform file
   FILE *f = fopen(Form("%s", filename), "r");
   if (f == NULL) {
      printf("Cannot find file \'%s\'\n", filename);
      return;
   }

   //open the root file
   strcpy(rootfile, filename);
   strcat(rootfile, ".root");
   TFile *outfile = new TFile(rootfile, "RECREATE");
   
   // define the drs tree
   TTree *drs = new TTree("drs",filename);
   drs->SetMarkerColor(602);
   drs->SetLineColor(602);
   drs->SetFillColor(602);
   drs->SetFillStyle(3001);

   // read time header
   fread(&th, sizeof(th), 1, f);
   printf("Found data for board #%d\n", th.board_serial_number);

   std::vector<int> channels;

   // read time bin widths
   memset(bin_width, sizeof(bin_width), 0);
   for (ch=0 ; ch<5 ; ch++) {
      fread(hdr, sizeof(hdr), 1, f);
      if (hdr[0] != 'C') {
         // event header found
         fseek(f, -4, SEEK_CUR);
         break;      
      }
      i = hdr[3] - '0' - 1;
      printf("Found timing calibration for channel #%d\n", i+1);
      fread(&bin_width[i][0], sizeof(float), 1024, f);

      channels.push_back(i+1);
   }

   // event header fields
   drs->Branch("event", &eh.event_serial_number, "event/i");
   drs->Branch("year", &eh.year, "year/s");
   drs->Branch("month", &eh.month, "month/s");
   drs->Branch("day", &eh.day, "day/s");
   drs->Branch("hour", &eh.hour, "hour/s");
   drs->Branch("minute", &eh.minute, "minute/s");
   drs->Branch("second", &eh.second, "second/s");
   drs->Branch("millisecond", &eh.millisecond, "millisecond/s");
   drs->Branch("range", &eh.range, "range/s");
   drs->Branch("board", &eh.board_serial_number, "board/s");
   drs->Branch("tcell", &eh.trigger_cell, "tcell/s");
   //
   for (unsigned ichannel=0; ichannel<channels.size(); ++ichannel) {
      int channel = channels[ichannel];
      drs->Branch(Form("t%d",channel), time[channel-1], Form("t%d[1024]/D",channel));  
      drs->Branch(Form("w%d",channel), waveform[channel-1], Form("w%d[1024]/D",channel));  
   }

   // loop over all events in data file
   for (n=0; n<nevents_max ; n++) {
      // read event header
      i = (int)fread(&eh, sizeof(eh), 1, f);
      if (i < 1)
         break;
         
      if (eh.event_serial_number < 1000 && eh.event_serial_number%100 == 0) printf("Found event #%d\n", eh.event_serial_number);
      if (eh.event_serial_number%1000 == 0) printf("Found event #%d\n", eh.event_serial_number);

      // reach channel data
      for (ch=0 ; ch<5 ; ch++) {
         i = (int)fread(hdr, sizeof(hdr), 1, f);
         if (i < 1)
            break;
         if (hdr[0] != 'C') {
            // event header found
            fseek(f, -4, SEEK_CUR);
            break;      
         }
         chn_index = hdr[3] - '0' - 1;
         fread(voltage, sizeof(short), 1024, f);
         
         for (i=0 ; i<1024 ; i++) {
            // convert data to volts
            waveform[chn_index][i] = (voltage[i] / 65536. + eh.range/1000.0 - 0.5);
            
            // calculate time for this cell
            for (j=0,time[chn_index][i]=0 ; j<i ; j++)
              time[chn_index][i] += bin_width[chn_index][(j+eh.trigger_cell) % 1024];            
         }
      }
    
      // align cell #0 of all channels
      t1 = time[0][(1024-eh.trigger_cell) % 1024];
      for (ch=1 ; ch<4 ; ch++) {
         t2 = time[ch][(1024-eh.trigger_cell) % 1024];
         dt = t1 - t2;
         for (i=0 ; i<1024 ; i++)
            time[ch][i] += dt;
      }

      // fill root tree
      drs->Fill();
   }

   // print number of events
   printf("Wrote %d events into the output file %s\n", n, outfile->GetName());
   
   // save and close root file
   drs->Write();
   drs->ResetBranchAddresses();
   //outfile->Close();
}
