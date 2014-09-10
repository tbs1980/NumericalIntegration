#include "kronrodLaurieGautschi.h"
#include "kronrodPiessensBurkardt.h"

#include <iostream>
#include <iomanip>

using namespace Kronrod;
using namespace Eigen;
using namespace std;

int main()
{
    Array<double,Dynamic,2> ans;
    int n = 10;
    ans = multiPrecisionKronrod(n);

    cout << "Laurie/Gautschi:" << endl;
    cout << setprecision(15);
    cout << ans << endl;
    cout << endl << endl;

    cout << "PiessensBurkardt" << endl;

    Array<double, Dynamic, 1> xGK;
    Array<double, Dynamic, 1> wGK;
    Array<double, Dynamic, 1> wG;

    kronrod(n, xGK,  wGK, wG);

    cout << xGK << endl << endl;
    cout << wGK << endl << endl;
    cout << wG << endl << endl;
    return 0;
}
