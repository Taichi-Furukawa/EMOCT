#include <iostream>
#include "libs/table.h"

using namespace std;
using namespace arch;

int main() {
    cout << "Hello, World!" << endl;
    table<float> ta(4,4);
    ta.fill(1.0f);

    cout<<ta<<endl;





    return 0;
}