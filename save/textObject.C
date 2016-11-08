#include <TROOT.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TList.h>
#include <TNamed.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

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
        cout<< "Rint.History file is " << fname <<endl;
        std::ifstream file(fname);
        if (!file) {
            cout<< "***Error TextObject::TextObject: could not find history file " << fname <<endl;
            exit(1);
        }

        std::string line_read;
        std::string line;
        while (std::getline(file, line_read)) line = line_read;
        cout<< "Command line is: " << line <<endl;

        // create an object and write it into the tree
        TNamed* textObject = new TNamed("Command line", line.c_str());
        list->AddLast(textObject);
        file.close();
    }
    void writeFile(const char* fname) {
        std::ifstream file(fname);
        if (!file) {
            cout<< "***Error TextObject::TextObject: could not find file " << fname <<endl;
            exit(1);
        }

        std::stringstream ss;
        std::string line;
        while (std::getline(file, line)) ss << line << endl;

        // create an object and write it into the tree
        TNamed* textObject = new TNamed(fname, ss.str());
        list->AddLast(textObject);
        file.close();
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

TList* textObject()
{
    TList* list = new TList;
    TextObject textObject(list);

    const char* fname = __FILE__;
    textObject.writeFile(fname);

    cout<< "\nlist content\n" <<endl;

    textObject.print(list);

    return list;
}
