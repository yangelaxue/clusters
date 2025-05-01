#include "TurbGen.h"

using namespace std;

/* Global variables */
enum {X,Y,Z};
double NDIM = 3.;
int N[3] = {64, 64, 64};
double L[3] = {1., 1., 1.};
double KMIN = 2.; // minimum wavenumber in units of 2pi / L[X]
double KMAX = 20.; // maximum wavenumber in units of 2pi / L[X]
double KMID = 1e38; // middle wavenumber in units of 2pi / L[X]
int spect_form = 2; // 0: band, 1: parabolic, 2: power-law
double power_law_exp = -5/3;    // power-law exponent -2: Burgers, -5/3: Kolmogorov, 1.5: would Kazantsev
double power_law_exp_2 = -5/3;
double angles_exp = 1.; // if power-law spectrum, spectral sampling of modes.
                        // number of modes chosen  increases as K^angles_exp
double sol_weight = .5; // 1.0: solenoidal driving, 0.0: compressive driving, 0.5: natural mixture
int random_seed = 10000; // random seed for this instance

int seed = -random_seed;
// seed = -seed;






int main() {
    TurbGen tg = TurbGen();
    // tg.test();
    // tg.verbose = 2;

    std::string parameter_file = "TurbGen.par";
    tg.init_driving(parameter_file, 0.);
    // printf("Print this!\n");

    // for (int i=0;i<10;i++) {
    //     double out = tg.ran1(&seed);
    //     double out2 = tg.ran2(&seed);
    //     std::cout<<out<<" "<<out2<<std::endl;
    // }

    return 0;
};