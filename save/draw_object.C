// bertrand.bellenot@cern.ch roottalk@cern.ch Dec 2, 2014
// Here is a simple example macro showing how to read and draw an object
// by using its TKey from the TFile’s list of TKeys

#include "Riostream.h"
#include "TFile.h"
#include "TList.h"
#include "TKey.h"

void draw_object(const char *file_name, const char *obj_name)
{
   // first open the file
   TFile *file = TFile::Open(file_name);
   if (!file || file->IsZombie()) {
      std::cout << "Cannot open " << file_name << "! Aborting..." << std::endl;
      return;
   }
   // get the list of keys
   TList *list = (TList *)file->GetListOfKeys();
   if (!list) {
      std::cout << "Cannot get the list of TKeys! Aborting..." << std::endl;
      return;
   }
   // try to find the proper key by its object name
   TKey *key = (TKey *)list->FindObject(obj_name);
   if (!key) {
      std::cout << "Cannot find a TKey named" << obj_name << "! Aborting..." << std::endl;
      return;
   }
   // finally read the object itself
   TObject *obj = ((TKey *)key)->ReadObj();
   if (!obj) {
      std::cout << "Cannot read the object named " << obj_name << "! Aborting..." << std::endl;
      return;
   }
   // and draw it
   obj->Draw();
}
