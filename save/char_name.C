#include <iostream>

using std::cout;    using std::endl;

void char_name_print(const char* name="")
{
    if (name && *name) cout<< name <<endl;
    else if (!name) cout<< "name is zero" <<endl;
    else cout<< "name is empty" <<endl;
}

void char_name()
{
    char_name_print(0);
    char_name_print();
    char_name_print("some name");
}
