#include <iostream>
#include "libs/table.h"
#include "GACT.h"

using namespace std;
using namespace arch;

int main() {
    cout << "Hello, World!" << endl;
    GACT ga(1000,1000,1.0,0.5);
    ga.Evolution();





    return 0;
}