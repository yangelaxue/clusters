/*
.h files include all the functions and variables that we would like to use.
*/

// Libraries regarding io functions
#include <iostream>
#include <iomanip>      // manip = manipulate // allows for std::stringstream and << operators to work.
#include <fstream>      // file stream // Deals with file streams (regarding evolfile)
// #include <sstream>      // string stream
// Libraries with useful functions
#include <iterator> // (regarding evolfile)
// #include <algorithm>
// Easy use of variable types
#include <vector>
// #include <string>
// C-style libraries
#include <cstring>      // Manipuate strings in C-style
// #include <cstdio>       // Manipulate standard io in C-style
// #include <cstdlib>      // Includes C functions
// #include <cstdarg>
#include <cfloat>
#include <cmath>


namespace NameSpaceTurbGen {
    // Defines the maximum number of modes allowed to be generated.
    static const int tgd_max_nmodes = 100000;
}


class TurbGen
{
    // Define all the variables which will be used internally within the class.
    private:
        enum {X,Y,Z};   // X, Y and Z will have specific labels assigned to them, By default, it is 0, 1 and 2 respectively
                        // This makes tracking indices easier
        std::vector<double> mode[3], aka[3], akb[3], OUphases, ampl;    // mode arrays[3], phases[3] and amplitudes
                                                                        // [3] indicates that variables will store 3 vectors
        int nmodes; // Number of modes to be generated
        double ndim; // ndim = 1, 1.5, 2, 2.5 or 3
        int ncmp; // ncmp = 1 for ndim=1, ncmp = 2 for ndim=2, ncmp = 3 for ndim=1.5,2.5 or 3
        int random_seed, seed; // random_seed is the original starting seed; seed gets updated by calls to random number generator
        int spect_form; // 0: band, 1: parabola, 2: power law
        int nsteps_per_t_turb; // Number of driving steps per turnover time
        int step; // internal OU step number
        double L[3]; //domain physical length L[dim] with dim = X, Y or Z
        double velocity; // velocity dispersion
        double t_decay; // turbulent turnover time, equivalent to T in exp(-t/T) (autocorrelation time-scale)
        double dt; // time-step for OU update and for generating driving patterns
        double power; // driving energy injection rate (~ v^3/L)
        double OUvar; // OUvar = sqrt(power/t_decay)
        double kmin, kmax; // minimum and maximun wave modes to drive or generate
        double kmid; // partitions the power-law spectrum into [kmin,kmid] and (kmid,kmax]
        double sol_weight, sol_weight_norm; // coefficient for decomposing turbulence and driving into solenoidal and compressive modes
        double power_law_exp, power_law_exp_2; // for a power-law spectrum. power_law_exp_2 can be used for (kmid,kmax]
        double angles_exp; // for a power-law spectrum. Determines how sparsely to sample k-space (angles_exp=2. for dense sampling)
        double ampl_factor[3]; // scale amplitude by this factor (default 1., 1., 1.)
        int ampl_auto_adjust; // switch (0,1) to turn off/on automatic amplitude adjustment
        // Variables not directly used in constructing turbulence:
        int verbose; // Shell output detail. 0: not output, 1: standard output, 2: detailed output
        std::string ClassSignature;
        std::string parameter_file;
        std::string evolfile;
        // MPI related variables:
        int PE; // MPI task for printf purposes, if provided.

    /* Class related functions! */
    // These are equivalent to the __init__ functions in Python (I think)
    public: TurbGen(void) {
        Constructor(0);
    };
    public: TurbGen(const int PE) {
        Constructor(PE);
    };
    private: void Constructor(const int PE) {
        ClassSignature = "TurbGen: ";
        this->PE = PE;
        verbose = 2; // default verbose value
        evolfile = "TurbGen.dat";
    };
    // This deletes everything class related (deallocates memory).
    // This is called a destructure, preceeded by a ~.
    public: ~TurbGen() { // This is run automatically at the end of the script.
        if (verbose>1) std::cout<<ClassSignature<<"destrutor called."<<std::endl;
    };


    /* UNIMPLEMENTED FUNCTIONS
    



    public: void set_verbose(const int verbose) {
        // Sets the level of detail to print to interface while running this programme.
        this->verbose = verbose;    // in this case, this->verbose refers to the to variable verbose on this class,
                                    // as opposed to the input variable verbose
    };

    private: std::string FuncSig(const std::string func_name) {
        // create and return function signature
        return func_name+": ";
    };
    
    */

    /* get functions */
    // these variables are private otherwise
    public: 
        double get_turnover_time(void) {
            return t_decay;
        };
        int get_nsteps_per_t_turb(void) {
            return nsteps_per_t_turb;
        };
        int get_number_of_compoents(void) {
            return ncmp;
        };
        std::vector<std::vector<double>> get_modes(void) {  // returns a 2-D vector
                                                                    // typically, std::vectors have dynamic sizes.
            std::vector<std::vector<double>> ret;               // standing for return
            ret.resize((int)ndim);                              // (int) always rounds down
            for (int d=0; d<(int)ndim; d++) ret[d] = mode[d];   // ret is (int)ndim dimensions, mode is 3 dimensions
            return ret;
        };
        std::vector<double> get_amplitudes(void) {
            return ampl;
        };
    /* set functions */
    private: void set_number_of_components(void) {
        if ((ndim!=1) && (ndim!=1.5) && (ndim!=2) && (ndim!=2.5) && (ndim!=3)) {
            // TODO: Raise Error!!! 
        }
        if (ndim==1) ncmp = 1;
        else if (ndim==2) ncmp = 2;
        else if (ndim==1.5 || ndim==2.5 || ndim==3) ncmp=3;
        else  {
            printf("ERROR: number of dimensions must be 1, 1.5, 2, 2.5 or 3.\n");
            exit(-1); // exit(-1) terminates the programme, with -1 letting us know that an error was encountered.
        }        
    }
    private: void set_solenoidal_weight_normalisation(void) {
        // sqrt and pow come from the <cmath> library
        sol_weight_norm = sqrt(3./ncmp)*sqrt(3.)*1./sqrt(1.-2.*sol_weight+ncmp*pow(sol_weight,2.));
    };

    /**********************************/
    /* Initialise Turbulence Driving! */
    /**********************************/
    public: int init_driving(std::string parameter_file) {
        return init_driving(parameter_file, 0.0); // call with time = 0.0
    };
    public: int init_driving(std::string parameter_file, const double time) {
        /*
        Initialise Turbulence Driving.
        Does not drive turbulence itself.
        */

        // See if parameter file exists before reading from its values.
        this->parameter_file = parameter_file;
        FILE * fp = fopen(parameter_file.c_str(), "r");
        if (fp==NULL) {
            printf("Parameter file does not exist.\n");
            exit(-1);
        }
        else {
            fclose(fp);
        }

        // Read from parameter file and define variables.
        double KMIN, KMAX, KDRIV;
        std::vector<double> ret;
        ret = read_from_parameter_file("ndim", "d"); ndim = ret[0];
        set_number_of_components();
        // Deal with multiple values for a given parameter.
        std::string ncmp_str; // variable ncmp as string type
        std::stringstream dummystream; dummystream << ncmp; dummystream >> ncmp_str; dummystream.clear();
        // read remaining parameters from file
        ret = read_from_parameter_file("L", ncmp_str+"d");
        for (unsigned int d = 0; d < ncmp; d++) L[d] = ret[d];
        ret = read_from_parameter_file("velocity", "d"); velocity = ret[0]; // Target turbulent velocity dispersion
        ret = read_from_parameter_file("KDRIV", "d"); KDRIV = ret[0]; // driving wavenumber in 2pi/L; sets t_decay below
        ret = read_from_parameter_file("KMIN", "d"); KMIN = ret[0]; // min wavenumber in 2pi/L
        ret = read_from_parameter_file("KMAX", "d"); KMAX = ret[0]; // max wavenumber in 2pi/L
        ret = read_from_parameter_file("sol_weight", "d"); sol_weight = ret[0]; // solenoildal weight
        ret = read_from_parameter_file("spect_form", "i"); spect_form = (int)ret[0]; // spectral form
        ret = read_from_parameter_file("power_law_exp", "d"); power_law_exp = ret[0]; // power-law exponent (if spect_form=2)
        power_law_exp_2 = power_law_exp; // driving does not support a 2nd PL section (yet)
        ret = read_from_parameter_file("angles_exp", "d"); angles_exp = ret[0]; // angles sampling exponent (if spect_form=2)
        ret = read_from_parameter_file("ampl_factor", ncmp_str+"d"); // adjust driving amplitude by factor in x[y[z]] (to adjust to target velocity)
        for (unsigned int d = 0; d < ncmp; d++) ampl_factor[d] = ret[d];
        ret = read_from_parameter_file("ampl_auto_adjust", "i"); ampl_auto_adjust = (int)ret[0]; // automatic amplitude adjustment switch
        ret = read_from_parameter_file("random_seed", "i"); random_seed = (int)ret[0]; // random seed
        ret = read_from_parameter_file("nsteps_per_t_turb", "i"); nsteps_per_t_turb = (int)ret[0]; // number of pattern updates per t_decay


        // DEFINE QUANTITIES!

        // Characteristic k-scales
        kmin = (KMIN-DBL_EPSILON) * 2*M_PI/L[X];
        kmax = (KMAX+DBL_EPSILON) * 2*M_PI/L[X];
        kmid = kmax; // Driving does not support a second power-law.
        // Characteristic time scales
        t_decay = L[X] / KDRIV / velocity;
        dt = t_decay / nsteps_per_t_turb;
        // Initiate driving
        step = -1; // For the OU process
        seed = random_seed;

        // Characteristic velocity field properties
        const double ampl_coeff = .15; // This yields good user-to-target velocity dispersion.
        power = pow(ampl_coeff*velocity, 3.) / L[X];
        OUvar = sqrt(power/t_decay);

        // Some checks regarding auto-adjusting amplitudes
        if ((ampl_auto_adjust==1) && time>0.) {
            if (!read_ampl_factor_from_evol_file(time)) return -1;
        }

        // This is such that user inputs of ampl_factor are more user friendly.
        for (unsigned int d=0; d<ncmp; d++) {
            ampl_factor[d] = pow(ampl_factor[d], 1.5);
        }

        // Initialse a single turbulent realisation with these given input parameters.
        set_solenoidal_weight_normalisation();
        init_modes();
        OU_noise_init();
        get_decomposition_coeffs();

        if (PE == 0) {
            std::ofstream outfilestream(evolfile.c_str(), std::ios::app);
            outfilestream   << std::setw(24) << "#01_time" << std::setw(24) << "#02_time_in_t_turb"
                            << std::setw(24) << "#03_ampl_factor_x" << std::setw(24) << "#04_ampl_factor_y" << std::setw(24) << "#05_ampl_factor_z"
                            << std::setw(24) << "#06_v_turb_prev_x" << std::setw(24) << "#07_v_turb_prev_y" << std::setw(24) << "#08_v_turb_prev_z"
                            << std::endl; // This creates column labels, ready for columns to be filled.
            outfilestream.close();
            outfilestream.clear();
        }

        printf("Turbulence Driving is Initiated!\n");

        return 0;
    }
    
    /*********************/
    /* Drive Turbulence! */
    /*********************/

    public: bool check_for_update(const double time) {
        /* Do not automatically adjust power amplitude */
        double v_turb[3]; for (int d=0; d<3; d++) v_turb[d] = -1.;
        return check_for_update(time, v_turb);
    }

    public: bool check_for_update(const double time, const double v_turb[]) {
        int step_requested = floor(time/dt); // dt is defined as t_decay/nsteps_per_t_turb in init_driving.
        if (step_requested<=step) {
            printf("Requested step is less than or equal to current step."); // step=-1 upon init_driving.
            return false;
        }

        // Update OUphase.
        for (int is=step; is<step_requested; is++) {
            OU_noise_update();
            step++;
        }

        /*
        Adjust driving amplitude to read user-defined target velocity dispersion v_turb
        */

        get_decomposition_coeffs(); // Calculate new solenoidal and compressive coefficients from OUphases.
        double time_gen = step*dt;
        if (PE==0) write_to_evol_file(time, ampl_factor, v_turb);

        return true;
    }

    /***************/
    /* Return grid */
    /***************/

    void get_turb_vector_unigrid(const double pos_beg[], const double pos_end[], const int n[], float *return_grid[]) {

        // Compute grid increments dx, dy and dz.
        double del[3] = {1., 1., 1.};
        for (int d=0; d<(int)ndim; d++) if (n[d]>1) del[d] = (pos_end[d]-pos_beg[d]) / (n[d]-1);

        //stratch variables
        double a[3];
        double real, imag;
        // loop over all points in return grid
        for (int k=0; k<n[Z]; k++) {
            for (int j=0; j<n[Y]; j++) {
                for (int i=0; i<n[X]; i++) {
                    // clear
                    a[X] = 0.; a[Y] = 0.; a[Z] = 0.;

                    // Loop over all k-modes
                    for (int m=0; m<nmodes; m++) {
                        // Calculate real and imag component of exp(i(\vec{k} \cdot \vec{x})) = cos(\vec{k} \cdot \vec{x}) + isin(\vec{k} \cdot \vec{x})
                        real = cos(mode[X][m]*(pos_beg[X]+i*del[X]) + mode[Y][m]*(pos_beg[Y]+j*del[Y]) + mode[Z][m]*(pos_beg[Z]+k*del[Z]));
                        imag = sin(mode[X][m]*(pos_beg[X]+i*del[X]) + mode[Y][m]*(pos_beg[Y]+j*del[Y]) + mode[Z][m]*(pos_beg[Z]+k*del[Z]));

                        a[X] += ampl[m] * (aka[X][m]*real - akb[X][m]*imag);
                        if (ncmp > 1) a[Y] += ampl[m] * (aka[Y][m]*real - akb[Y][m]*imag);
                        if (ncmp > 2) a[Z] += ampl[m] * (aka[Z][m]*real - akb[Z][m]*imag);
                    }

                    

                    
                }
            }
        }

    }
    
    /*****************************************/
    /* Initialise Single Turbulent Instance! */
    /*****************************************/
    public: int init_single_realisation(
        const double ndim, const double L[3], const double KMIN, const double KMAX,
        const int spect_form, const double power_law_exp, const double angles_exp,
        const double sol_weight, const int random_seed) {
        /*
        Initialise modes for a single power-law model.
        */ 
       return init_single_realisation(
        ndim, L, KMIN, KMAX, KMAX,
        spect_form, power_law_exp, power_law_exp, angles_exp,
        sol_weight, random_seed);
    }
    public: int init_single_realisation(
        const double ndim, const double L[3], const double KMIN, const double KMID, const double KMAX,
        const int spect_form, const double power_law_exp, const double power_law_exp_2, const double angles_exp,
        const double sol_weight, const int random_seed) {
        /*
        Initialise modes for a single turbulent instance, as use for an initial condition.
        The actual field is not created here.
        */ 
        
        // Define class variables
        this->ndim = ndim;
        this->L[X] = L[X];
        this->L[Y] = L[Y];
        this->L[Z] = L[Z];
        this->kmin = (KMIN-DBL_EPSILON) * 2*M_PI/L[X];
        this->kmax = (KMAX+DBL_EPSILON) * 2*M_PI/L[X];
        this->kmid = (KMID+DBL_EPSILON) * 2*M_PI/L[X];
        this->spect_form = spect_form;
        this->power_law_exp = power_law_exp;
        this->power_law_exp_2 = power_law_exp_2;
        this->angles_exp = angles_exp;
        this->sol_weight = sol_weight;
        this->random_seed = random_seed;
        seed = this->random_seed; // Sets the dynamic seed as the initial random seed.
        OUvar = 1.;
        for (unsigned int d=0; d<3; d++) ampl_factor[d] = 1.; // Ampl factor is applied after the field is already created. Default to 1.
        // Set derived parameters
        set_number_of_components(); // defines ncmp
        set_solenoidal_weight_normalisation(); // defines sol_weight_norm
        // Initialise modes
        init_modes();
        // Initialise random phases
        OU_noise_init();
        // Calculate real and coefficients.
        get_decomposition_coeffs(); // This is based on how compressive or solenoidal the turbulence is.

        printf("Success!\n");
        
        return 0;      
    };

    /***********************************************/
    /* Helpers for Tnitialising a Turbulent Field! */
    /***********************************************/

    private: int init_modes(void) {
        /*
        Initiate modes for a single realisation.
        */

        int ikmin[3], ikmax[3], ik[3], tot_nmodes_full_sampling;
        double k[3], ka, kc, amplitude, parab_prefact;

        // For use in a power-law spectrum
        int iang, nang; // iang is a counter, nang is number to modes (angles) per k-shell
        double rand, phi, theta;

        // Define charactertic kc. Means more for a parabolic spectra.
        kc = (spect_form=1) ? 0.5*(kmin+kmax) : kmin;

        // For use for a paraboloid spectrum
        parab_prefact = -4./pow(kmax-kmin,2.);

        // Calculate the total number of modes to be generated in the case of full sampling (angles_exp=2.)
        ikmin[X] = 0;
        ikmin[Y] = 0;
        ikmin[Z] = 0;
        ikmax[X] = 256;
        ikmax[Y] = ((int)ndim>1) ? 256 : 0;
        ikmax[Z] = ((int)ndim>2) ? 256 : 0;
        
        nmodes = 0;
        for (ik[X]=-ikmax[X]; ik[X]<=ikmax[X]; ik[X]++) {
            k[X] = ik[X] * 2*M_PI / L[X];
            for (ik[Y]=-ikmax[Y]; ik[Y]<=ikmax[Y]; ik[Y]++) {
                k[Y] = ik[Y] * 2*M_PI / L[Y];
                for (ik[Z]=-ikmax[Z]; ik[Z]<=ikmax[Z]; ik[Z]++) {
                    k[Z] = ik[Z] * 2*M_PI / L[Z];
                    ka = sqrt(k[X]*k[X] + k[Y]*k[Y] + k[Z]*k[Z]);
                    if ((ka>=kmin) && (ka<=kmax)) nmodes++;
                }
            }
        }
        tot_nmodes_full_sampling = nmodes;
        nmodes = 0; // nmodes is a counter. Reset this counter.

        if (spect_form!=2) { // band and parabolic spectra always uses the full number of modes.
            if (tot_nmodes_full_sampling>NameSpaceTurbGen::tgd_max_nmodes) {
                printf("Too many modes");
                exit(-1);
            }
        }

        /* Create modes for a band or parabolic spectrum! */
        if (spect_form!=2) {
            for (ik[X]=-ikmax[X]; ik[X]<=ikmax[X]; ik[X]++) {
                k[X] = ik[X] * 2*M_PI / L[X];
                for (ik[Y]=-ikmax[Y]; ik[Y]<=ikmax[Y]; ik[Y]++) {
                    k[Y] = ik[Y] * 2*M_PI / L[Y];
                    for (ik[Z]=-ikmax[Z]; ik[Z]<=ikmax[Z]; ik[Z]++) {
                        k[Z] = ik[Z] * 2*M_PI / L[Z];
                        ka = sqrt(k[X]*k[X] + k[Y]*k[Y] + k[Z]*k[Z]);

                        if ((ka>=kmin) && (ka<=kmax)) {
                            if (spect_form==0) amplitude = 1.;
                            if (spect_form==1) amplitude = fabs(parab_prefact*pow(ka-kc,2.)+1); // fabs is a float absolute function

                            amplitude = sqrt(amplitude) * pow(kc/ka,((int)ndim-1.)/2.);
                            
                            // Append newlt calculated modes and amplitudes to the end of respective arrays.
                            ampl.push_back(amplitude); // pusk_back places new element at the end of the array. Much like append in Python.
                            mode[X].push_back(k[X]);
                            if ((int)ndim>1) mode[Y].push_back(k[Y]);
                            if ((int)ndim>2) mode[Z].push_back(k[Z]);

                            nmodes++;

                            // update if verbose
                            // if (nmodes % 1000 == 0) if (verbose) std::cout<< "... " << nmodes << " modes created out of " << tot_nmodes_full_sampling << " ..." <<std::endl;
                        }
                    }
                }
            }
        }

        /* Create modes for a power-law spectrum! */
        if (spect_form==2) {
            printf("Power-law spectrum used");
            
            // Set initial random seed from the dynamic seed.
            int seed_init = -seed;
            double rand = ran2(&seed_init); // The first random needs to be generated to begin the sequence,
                                            // this number is not used.

            // Loop over all k-shells in [KMIN,KMAX]
            ikmin[0] = std::max(1, (int)round(kmin*L[0]/2/M_PI));
            ikmax[0] = (int)round(kmax*L[0]/2/M_PI);
            for (ik[0]=ikmin[0]; ik[0]<=ikmax[0]; ik[0]++) {
                
                // Loop over all modes to be generated per k-shell.
                nang = pow(2.,(int)ndim)*ceil(pow((double)ik[0],angles_exp));
                for (iang=0; iang<=nang; iang++){
                    
                    /* These determined at time of running */
                    // Define phi
                    phi = 2*M_PI * ran2(&seed);
                    if (ndim==1) phi = (phi<M_PI) ? 0 : M_PI;
                    // Define theta
                    theta = (ndim>2) ? acos(1.-2.*ran2(&seed)) : M_PI/2.;

                    // Randomises the length of the k-mode within the k-shell.
                    rand = ik[0] + ran2(&seed) - .5;
                    
                    // Define k[X], k[Y] and k[Z].
                    k[X] = 2*M_PI/L[X] * round(rand*sin(theta)*cos(phi));
                    k[Y] = (ndim>1) ? 2*M_PI/L[Y] * round(rand*sin(theta)*sin(phi)) : 0.;
                    k[Z] = (ndim>2) ? 2*M_PI/L[Z] * round(rand*cos(phi)) : 0.;
                    // Rounding ensures that field is periodic.

                    // Calculate the amplitude of this k-mode.
                    ka = sqrt(k[X]*k[X] + k[Y]*k[Y] + k[Z]*k[Z]);

                    if ((ka>=kmin) && (ka<=kmax)) {
                        if (nmodes>NameSpaceTurbGen::tgd_max_nmodes){
                            printf("Too many modes. Exting...\n");
                            exit(-1);
                        }

                        // Calculate the power-spectrum amplitude.
                        amplitude = (ka<kmid) ? pow(ka/kc,power_law_exp) : pow(kmid/kmin,power_law_exp)*pow(ka/kmid,power_law_exp_2);
                        amplitude = sqrt( 2.*sqrt(3.) * amplitude * pow((double)ik[0],(int)ndim-1) / (double)(nang) ) * pow(kc/ka,((int)ndim-1)/2.);

                        nmodes++; // Increase counter

                        // Add elements to the end of vectors.
                        ampl.push_back(amplitude);
                        mode[X].push_back(k[X]);
                        if ((int)ndim > 1) mode[Y].push_back(k[Y]);
                        if ((int)ndim > 2) mode[Z].push_back(k[Z]);

                        // Print updates
                        if (nmodes % 1000 == 0) std::cout<< "... " << nmodes << " modes generated..." <<std::endl;
                    }
                }
            }
        }
        return 0;
    };

    private: void OU_noise_init(void) {
        /*
        Initialise a pseudo-random sequence for the OU process
        */
        OUphases.resize(nmodes*ncmp*2); // For each component of each mode, generate a random coefficient for their real and imag components.
        for (int m=0; m<nmodes; m++) {
            for (int d=0; d<ncmp; d++) {
                for (int ir=0; ir<2; ir++) {
                    OUphases[2*ncmp*m+2*d+ir] = OUvar * get_normal_ran(); // This is sequential. OUvar is set different by different methods.
                }
            }
        }
    };

    private: void OU_noise_update(void) {
        /*
        The OU process for driving turbulence is

        x_(n+1) = f*x_n + sigma*sqrt(1-f^2)*z_n

        where a = dv/dt and x is the Fourier transform of a, f=e^(-dt/t_decay),
        z_n is a random Gaussian number and sigma is the desired variance of the OU sequence.
        */
        const double damping_factor = exp(-dt/t_decay);
        for (int m=0;m<nmodes;m++) {
            for (int d=0;d<ncmp;d++) {
                for (int ir=0;ir<2;ir++) {
                    OUphases[2*ncmp*m+2*d+ir] = OUphases[2*ncmp*m+2*d+ir]*damping_factor
                                                + OUvar*sqrt(1.-pow(damping_factor,2.))*get_normal_ran();
                }
            }
        }
    }

    private: void get_decomposition_coeffs(void) {
        /*
        Get the coefficents for the real and imaginary components of the forcing field.
        This is based on how compressive or solenoidal the turbulent mode is.
        This does not take into consideration of the amplitude of forcing of each mode.
        */

        // These are the real and imaginary coefficients for each dimension of each mode.
        for (int d=0; d<3; d++) {
            aka[d].resize(nmodes);
            akb[d].resize(nmodes);
        }

        double ka, kb, kk, diva, divb, curla, curlb;
        for (int m=0; m<nmodes; m++) {
            // Loop over all modes.
            ka = 0.;
            kb = 0.;
            kk = 0.;
            for (int d=0; d<ncmp; d++){
                // Loop over all components of each mode.
                double this_mode = mode[std::min(d,(int)ndim-1)][m]; // This deals with fractional dimensions.
                kk = kk + this_mode*this_mode;                  // get k^2
                ka = ka + this_mode*OUphases[2*ncmp*m+2*d+1];   // get k\cdot<random phase>
                kb = kb + this_mode*OUphases[2*ncmp*m+2*d+0];   // get k\cdot<another random phase>
            }
            for (int d=0; d<ncmp; d++) {
                // Loop over all components of each mode.
                double this_mode = mode[std::min(d,(int)ndim-1)][m];
                diva = this_mode*ka/kk;
                divb = this_mode*kb/kk;
                curla = OUphases[2*ncmp*m+2*d+0] - divb;
                curlb = OUphases[2*ncmp*m+2*d+1] - diva;
                aka[d][m] = sol_weight*curla + (1.-sol_weight)*divb;
                akb[d][m] = sol_weight*curlb + (1.-sol_weight)*diva;
            }
        }
    };











    /********************/
    /* HELPER FUNCTIONS */
    /********************/

    /*
    Random Number Generators
    */
    private: double get_normal_ran(void) {
        /*
        Returns a uniform random variable with unit variance.
        */
        double r1 = ran1(&seed);
        double r2 = ran1(&seed);
        double g1 = sqrt(2.*log(1./r1)) * cos(2.*M_PI*r2);
        return g1;
    }
    private: double ran1(int *idum) {
        /*
        Uniform random number generator.
        Will generate a sequence of random numbers with the same seed,
        unless the seed changes, in which case the algorithm restarts.
        */ 

        static const int IA=16807, IM=2147483647, IQ=127773, IR=2836;
        static const double AM=1./IM, RNMX=1.-1.2e-7;
        if (*idum<=0) *idum = std::max(-*idum, 1);
        int k = *idum/IQ;
        *idum = IA*(*idum-k*IQ)-IR*k;
        if (*idum<0) *idum = *idum+IM;
        int iy = *idum;
        double ret = std::min(AM*iy, RNMX);
        return ret;
    };
    private: double ran2(int *idum) {
        /*
        More robust uniform random number generator.
        */
        static const int IM1=2147483563, IM2=2147483399, IMM1=IM1-1, IA1=40014, IA2=40692, IQ1=53668,
                         IQ2=52774, IR1=12211, IR2=3791, NTAB=32, NDIV=1+IMM1/NTAB;
        static const double AM=1./IM1, RNMX=1.-1.2e-7;
        static int idum2 = 123456789;
        static int iy = 0;
        static int iv[NTAB];
        int j, k;
        if (*idum <= 0) {
            *idum = std::max(-*idum, 1);
            idum2 = *idum;
            for (j = NTAB+7; j >= 0; j--) {
                k = *idum/IQ1;
                *idum = IA1*(*idum-k*IQ1)-k*IR1;
                if (*idum < 0) *idum += IM1;
                if (j < NTAB) iv[j] = *idum;
            }
            iy = iv[0];
        }
        k = *idum/IQ1;
        *idum = IA1*(*idum-k*IQ1)-k*IR1;
        if (*idum < 0) *idum += IM1;
        k = idum2/IQ2;
        idum2 = IA2*(idum2-k*IQ2)-k*IR2;
        if (idum2 < 0) idum2 += IM2;
        j = iy/NDIV;
        iy = iv[j]-idum2;
        iv[j] = *idum;
        if (iy < 1) iy += IMM1;
        double ret = std::min(AM*iy, RNMX);
        return ret;
    };

    /*
    File Readers
    */
    private: std::vector<double> read_from_parameter_file(std::string var_name,std::string type) {

        std::vector<double> ret;
        FILE * fp; // *fp is a pointer to a file.
        char * line = NULL; // *line is a pointer to a char, char is similar (but different) to a string, storing a single character.
        size_t len = 0; // This stores a positive integer size value.
        ssize_t read; // This stores an integer size value.

        // Check if parameter_file exists before proceeding.
        fp = fopen(parameter_file.c_str(), "r");
        if (fp==NULL) {
            printf("Supplied parameter_file does not exist.\n");
            exit(-1);
        }

        // Search for variable value.
        bool found = false; // check if variable var_name is found.
        // Loop over all lines until variable name is found. If found==true, break.
        while ((read=getline(&line,&len,fp)) != -1) { // this defines the line variable.
            if (strncmp(line, var_name.c_str(), strlen(var_name.c_str())) == 0) {   // strlen returns the length of a string, 
                                                                                    // strncmp compares strings of a given length.
                // The variable has been found!
                char *substr1 = strstr(line, "="); // extract everything after and including "=" from the line variable.
                char *substr2 = strstr(substr1, "!"); // deal with comment '! ...'
                char *substr3 = strstr(substr1, "#"); // deal with comment '# ...'
                int end_index = strlen(substr1);
                if ((substr2==NULL) && (substr3==NULL)) {
                    end_index -= std::max(strlen(substr2),strlen(substr3));
                }
                else {
                    if (substr2 != NULL) end_index -= strlen(substr2);
                    if (substr3 != NULL) end_index -= strlen(substr3);
                }
                
                char dest[100];
                memset(dest, '\0', sizeof(dest));
                strncpy(dest, substr1+1, end_index-1); // Copy string from substring1 into destination of given length.
                if ((type=="i") || (type=="1i")) {
                    ret.resize(1); ret[0] = double(atoi(dest)); // atoi converts string to integer.
                }
                if ((type=="d") || (type=="1d")) {
                    ret.resize(1); ret[0] = double(atof(dest)); // atof converts string to float.
                }
                if ((type=="2d") || (type=="3d")) { // The integer preceeding data type tell's us the length of the parameter value.
                    int ncmp;
                    std::stringstream dummystream;
                    dummystream << type[0]; dummystream >> ncmp; dummystream.clear(); // This equates ncmp with the first int in type.
                    ret.resize(ncmp);
                    std::string cmpstr = dest;
                    int cnt = 0; // Counts number of ',' in string
                    for (size_t offset = cmpstr.find(","); offset != std::string::npos; offset = cmpstr.find(",", offset+1)) {
                        ++cnt; // This increments the value of cnt.
                    }
                    if ((cnt != 0) && (cnt != ncmp-1)) { // error check
                        printf("Number of values does not match description.\n");
                        exit(-1);
                    }
                    for (int d=0; d<ncmp; d++) {
                        std::string substr = cmpstr; // copy into substr in case no commas provided
                        if (cnt == ncmp-1) { // if user has specified ncmp components; else copy the same value into each component
                            substr = cmpstr.substr(0, cmpstr.find(",")); // extract substring before first ','
                            cmpstr = cmpstr.substr(cmpstr.find(",")+1, cmpstr.size()-cmpstr.find(",")); // reduce cmpstr to remainder for next loop iter
                        }
                        dummystream << substr; dummystream >> ret[d]; dummystream.clear(); // copy into ret
                    }
                }
                found = true;
            }
            if (found) break;
        }
        fclose(fp);
        if (line) free(line); // de-allocates memory to this variable.
        if (!found) {
            std::cout<<var_name<<std::endl;
            printf("Variable has not been found\n");
            exit(-1);
        }
        return ret;
    }

    public: bool read_ampl_factor_from_evol_file(const double time) { // TODO: COMMENT
        return true;
        // std::ifstream infilestream(evolfile.c_str());
        // if (infilestream.is_open()) {
        //     // read each line from file into a vector
        //     std::vector<std::string> lines;
        //     std::string line;
        //     while (std::getline(infilestream, line)) lines.push_back(line);
        //     infilestream.close();
        //     bool initial_copy_done = false;
        //     // loop backwards through the vector to find the right time entry
        //     for (int i = lines.size() - 1; i >= 0; i--) {
        //         // skip header line(s)
        //         if (lines[i].find("time") != std::string::npos) continue;
        //         // split line
        //         std::istringstream buffer(lines[i]);
        //         std::vector<std::string> line_split((std::istream_iterator<std::string>(buffer)), std::istream_iterator<std::string>());
        //         double time_in_file = double(atof(line_split[0].c_str())); // time
        //         double ampl_factor_in_file[3];
        //         for (unsigned int d = 0; d < ncmp; d++) ampl_factor_in_file[d] = double(atof(line_split[d+2].c_str()));
        //         if (!initial_copy_done) { // copy the last valid ampl_factor from the evolution file
        //             for (unsigned int d = 0; d < ncmp; d++) ampl_factor[d] = ampl_factor_in_file[d];
        //             initial_copy_done = true;
        //         }
        //         if (time_in_file >= time-1.1*dt) { // search backwards and update ampl_factor with last valid time
        //             for (unsigned int d = 0; d < ncmp; d++) ampl_factor[d] = ampl_factor_in_file[d];
        //         }
        //         else break;
        //     }
        // }
        // else {
        //     return false;
        // }
        // return true;
    }; // read_ampl_factor_from_evol_file

    public: bool write_to_evol_file(const double time, const double ampl_factor[], const double v_turb[]) { // TODO: Comment!
        std::ofstream outfilestream(evolfile.c_str(), std::ios::app);
        outfilestream.precision(16);
        outfilestream   << std::scientific << std::setw(24) << time << std::setw(24) << time/t_decay
                        << std::setw(24) << pow(ampl_factor[X],1.0/1.5) << std::setw(24) << pow(ampl_factor[Y],1.0/1.5) << std::setw(24) << pow(ampl_factor[Z],1.0/1.5)
                        << std::setw(24) << v_turb[X] << std::setw(24) << v_turb[Y] << std::setw(24) << v_turb[Z]
                        << std::endl;
        outfilestream.close();
        outfilestream.clear();
        return true;
    }; // write_to_evol_file
        
        
        // & here points to the address of the variables.
            

};