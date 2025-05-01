"""
My version of TurbGen originally written in CPP by Christoph Federrath.
Includes
    - the generation of a single instance of a turbulent field,
    - simulating evolution in a turbulent field using an Ornstein-Uhlenbeck process.

TODO Finish commenting on this
TODO Remove redundant stuff
TODO Tidy
Data: Mar 2025
Author: Angela Xue
"""

#%% Imports

from math import pi as PI
import numpy as np
from numpy.random import random

#%% Other functions

#%% Start Class

class TurbGen:

    """
    Class which generates a turbulent field given a parameter file.
    Assumes uniformly spaced grids, wtih all sides of the grid having equal length and resolution.
    Allows for 1, 2, or 3 dimensions.

    Global Parameters
    -----------------
    RANDOM_SEED : int
        Seed for randomly generated field.
    NSTEPS_PER_T_TURB : int
        Number of time steps per t_decay
    VELOCITY : float
        Target turbulent velocity dispersion

    K_MIN : float
        Minimum wave number in turbulent field in units of 2PI/L.
    K_MAX : float
        Maximum wave number in turbulent field in units of 2PI/L.
    K_DRIV : float
        Driving turbulent mode.
    K_MID : float TODO Obsolete?
        Middle wave number in turbulent field in units of 2PI/L,
        for if we want to split up the spectrum into [K_MIN,K_MID] and [K_MID,K_MAX].

    POWER_LAW_EXP :float
        Power law in spectrum.
        Burgers: -2, Kolmogorov: -5/2, Kazantsev: 3/2.
    POWER_LAW_EXP_2 : float
        Second powerlaw is [K_MID,K_MAX] is used. TODO Obsolete?
    ANGLES_EXP : float
        Spectral sampling of angles.
        Number of modes (angles) in k-shell surface increases as k^angles_exp.
        2.: for full sampling, 0.: for healpix-type sampling

    SOL_WEIGHT : float
        Solenoidal weight.
        1.: solenoidal driving, 0.: compressive driving, .5: natural mixture

    FILE_NAME = "turbulence_output.h5" TODO

    METHODS
    -------

    init_modes
        Creates lists of kx, ky and kz modes and their corresponding power-laws.
    
    _get_decomposition_coeffs
        Before generating the turbulent field, calculate the real and imaginary coefficients
        using the list of modes and amplitubes found in init_modes
    
    get_turb_vector_unigrid
        Creates a turbulent field on a uniform grid.

    """

    def __init__(self,parameter_file:str):
        """
        Initiate class with a parameter file which holds all global variables.
        """
        
        with open(parameter_file,'r') as f:
            lines = f.readlines()

            for line in lines:
                if line.startswith("RANDOM_SEED "):
                    self.RANDOM_SEED = int(line.split('=')[1].split('#')[0])
                
                # These concern turbulence driving.
                elif line.startswith("NSTEPS_PER_T_TURB "):
                    self.NSTEPS_PER_T_TURB = int(line.split('=')[1].split('#')[0])
                elif line.startswith("VELOCITY "):
                    self.VELOCITY = float(line.split('=')[1].split('#')[0])
                elif line.startswith("AMPL_FACTOR "):
                    AMPL_FACTOR = float(line.split('=')[1].split('#')[0])
                elif line.startswith("AMPL_AUTO_ADJUST "):
                    self.AMPL_AUTO_ADJUST = float(line.split('=')[1].split('#')[0])
                elif line.startswith("K_DRIV "):
                    self.K_DRIV = float(line.split('=')[1].split('#')[0])

                # These concern turbulence generation.
                elif line.startswith("NDIM "):
                    self.NDIM = int(line.split('=')[1].split('#')[0])
                elif line.startswith("N "):
                    N = int(line.split('=')[1].split('#')[0])
                elif line.startswith("L "):
                    L = float(line.split('=')[1].split('#')[0])
                elif line.startswith("K_MIN "):
                    self.K_MIN = float(line.split('=')[1].split('#')[0])
                elif line.startswith("K_MAX "):
                    self.K_MAX = float(line.split('=')[1].split('#')[0])
                elif line.startswith("K_MID "):
                    self.K_MID = float(line.split('=')[1].split('#')[0])
                elif line.startswith("POWER_LAW_EXP "):
                    POWER_LAW_EXP = line.split('=')[1].split('#')[0]
                    if '/' in POWER_LAW_EXP:
                        self.POWER_LAW_EXP = float(POWER_LAW_EXP.split('/')[0])/float(POWER_LAW_EXP.split('/')[1])
                    else:
                        self.POWER_LAW_EXP = float(POWER_LAW_EXP)
                elif line.startswith("POWER_LAW_EXP_2 "):
                    POWER_LAW_EXP_2 = line.split('=')[1].split('#')[0]
                    if '/' in POWER_LAW_EXP_2:
                        self.POWER_LAW_EXP_2 = float(POWER_LAW_EXP_2.split('/')[0])/float(POWER_LAW_EXP_2.split('/')[1])
                    else:
                        self.POWER_LAW_EXP_2 = float(POWER_LAW_EXP_2)
                elif line.startswith("ANGLES_EXP "):
                    self.ANGLES_EXP = float(line.split('=')[1].split('#')[0])
                elif line.startswith("SOL_WEIGHT "):
                    self.SOL_WEIGHT = float(line.split('=')[1].split('#')[0])
                else: pass

        if self.NDIM==1:
            self.N = (N,1,1)
            self.L = (L,1,1)
        elif self.NDIM==2:
            self.N = (N,N,1)
            self.L = (L,L,1)
        elif self.NDIM==3:
            self.N = (N,N,N)
            self.L = (L,L,L)
        else:
            raise ValueError("Dimensions other than 1, 2 or 3 are not supported.")
        self.delt = tuple(_L/_N for _L,_N in zip(self.L,self.N))
        self.AMPL_FACTOR = (AMPL_FACTOR**1.5,)*self.NDIM # amplitude is actually ~velocity^1.5

        # set_number_of_components()
        self.ncmp = self.NDIM

        self._set_solenoidal_weight_normalisation()
        self._set_pos_beg_end()

    def _set_pos_beg_end(self,):
        #TODO Maybe change for MPI
        # Define center of cells at each end of the grid.
        self.pos_beg = [delt/2 for delt in self.delt]
        self.pos_end = [L-delt/2 for L,delt in zip(self.L,self.delt)]
    
    def _set_solenoidal_weight_normalisation(self):
        """As described."""
        self.sol_weight_norm = np.sqrt(3/self.ncmp) * np.sqrt(3)
        self.sol_weight_norm /= np.sqrt(1 - 2*self.SOL_WEIGHT + self.ncmp*self.SOL_WEIGHT**2)

    def _set_turnover_time(self,):
        """
        The characteristic time scale for the driving turbulent mode is ~ L/v
        where L is the injection scale, here determined by K_DRIV.
        """

        self.t_decay = self.L[0] / self.K_DRIV / self.VELOCITY

    def _set_dt(self,):
        self.dt = self.t_decay/self.NSTEPS_PER_T_TURB

    def init_modes(self,k_min,k_max,k_mid):
        """
        Initiate k-modes that make up the turbulent field.

        Loop over every k-mode
        - calculate number of modes to generate for this k-shell
        - uniformly randomise the angles for generated k-modes
        - calculate kx, ky, kz for each mode in this k-shell
        - append k-modes between k_min and k_max to a list of all k-
        - for each k-mode in the k-shell, calculate the corresponding power-law amplitude.
        - append these amplitudes to a list of amplitudes.
        """

        kc = k_min # Characteristic k for scaling the amplitude below.
        
        self.nmodes = 0 # reset
        
        # Define empty arrays.
        self.ampl_arr = np.empty((0,1)) # Amplitude
        self.mode_arr = np.empty((0,3)) # k-vectors

        # Reset random seed.
        np.random.seed(self.RANDOM_SEED)

        # Loop over all k sizes.
        for ik in range(int(self.K_MAX+1)):
            nang = int(2**self.NDIM * np.ceil(ik**self.ANGLES_EXP)) # Number of angles to calculate per k-shell

            # Uniform sample of phi
            phi = 2*PI * random(nang)
            if self.NDIM==1:
                phi[phi<PI] = 0
                phi[phi!=0] = PI

            # 'Uniform' sample of theta
            theta = PI/2 * np.ones(nang) if self.NDIM<=2 else np.arccos(1.0 - 2.0*random(nang))

            # k vector
            rand = ik + random(nang) - .5
            kx = 2*PI/self.L[0] * np.round(rand*np.sin(theta)*np.cos(phi))
            ky = 2*PI/self.L[0] * np.round(rand*np.sin(theta)*np.sin(phi)) if self.NDIM>1 else np.zeros(nang)
            kz = 2*PI/self.L[0] * np.round(rand*(np.cos(theta))) if self.NDIM>2 else np.zeros(nang)

            ka = (kx**2 + ky**2 + kz**2)**.5

            # Add wave modes to list.
            where = np.where((ka>=k_min)*(ka<=k_max))
            # print(where)
            self.nmodes += len(where[0])

            amplitude = np.zeros(len(where[0]))
            
            amplitude[ka[where]<k_mid] = (ka[where][ka[where]<k_mid]/kc)**self.POWER_LAW_EXP
            amplitude[ka[where]>=k_mid] = (k_mid/k_min)**self.POWER_LAW_EXP * (ka[where][ka[where]>=k_mid]/k_mid)**self.POWER_LAW_EXP_2
            
            amplitude = np.sqrt(amplitude * ik**(self.NDIM-1) * 4.0*np.sqrt(3.0) / nang) * (kc/ka[where])**((self.NDIM-1)/2)
                
            self.ampl_arr = np.concatenate([self.ampl_arr,amplitude[:,np.newaxis]])

            modes = np.concatenate([kx[where], ky[where], kz[where]],axis=-1).reshape((3,-1)).T
            # print(ka[where])
            # print(kx[where],ky[where])
            # print(where)
            # print(modes)
            # print(np.sum([_kx**2 for _kx in modes],axis=-1)**.5)
            # print()
            self.mode_arr = np.concatenate([self.mode_arr,modes])
        
        # if self.NDIM==2:
        #     self.mode_arr = self.mode_arr[:,:2]
        # elif self.NDIM==1:
        #     self.mode_arr = self.mode_arr[:,0]

        self.ampl_arr = self.ampl_arr.flatten()

    def _OU_noise_init(self,OUvar):
        """
        Generate random Ornstein-Uhlenbeck variances.
        The shape of OUphases allows for the calculation of nmodes x NDIM x 2
        """

        OUphases = OUvar * np.random.normal(0,1,(self.nmodes,3,2))
        if self.NDIM<=1:
            OUphases[:,1][:] = 0
        if self.NDIM<=2:
            OUphases[:,2][:] = 0
        self.OUphases= OUphases


    def _get_decomposition_coeffs(self,):
        """
        The wave modes are decomposed into real and imaginary components.
        Additionally, there will be a mixture of compressive (\\Delta x v = 0) modes
        and non-compressive or solenoidal modes (\\Delta\\cdot v = 0) modes.

        This calculates the coefficients of the real and imaginary components of the turbulent field.
        """

        kk = np.sum([kx**2 for kx in self.mode_arr],axis=-1) # k^2 for each k-mode.
        ka = np.sum([kx*self.OUphases[:,d,1] for d,kx in enumerate(self.mode_arr.T)],axis=0) # k x OUrand for each phase
        kb = np.sum([kx*self.OUphases[:,d,0] for d,kx in enumerate(self.mode_arr.T)],axis=0) # k x OUrand for each phase
        # All these have shape==(<self.nmodes>,)
        
        diva = np.array([kx*ka/kk for kx in self.mode_arr.T])
        divb = np.array([kx*kb/kk for kx in self.mode_arr.T])
        curla = np.array([self.OUphases[:,d,0] - divb[d] for d in range(3)])
        curlb = np.array([self.OUphases[:,d,1] - divb[d] for d in range(3)])
        # All these have shape==(self.ncmp,self.nmodes)

        self.aka = np.array([self.SOL_WEIGHT*curla[d] + (1-self.SOL_WEIGHT)*divb[d] for d in range(3)]).T
        self.akb = np.array([self.SOL_WEIGHT*curlb[d] + (1-self.SOL_WEIGHT)*diva[d] for d in range(3)]).T
        # print(self.aka[:,2])
        # print(self.akb[:,2])
        # All these have shape==(self.nmodes,self.ncmp)

    def init_single_realisation(self,):
        """
        This creates a single realisation or instant of a turbulent field.
        """

        # Convert ks into physical units.
        k_min = self.K_MIN * 2*PI/self.L[0]
        k_max = self.K_MAX * 2*PI/self.L[0]
        k_mid = self.K_MID * 2*PI/self.L[0]

        # Initiate modes.
        self.init_modes(k_min,k_max,k_mid)
        OUvar = 1. # Ornstein-Uhlenbeck variance
        self._OU_noise_init(OUvar)
        self._get_decomposition_coeffs()

        # # set_solenoidal_weight_normalisation()
        # sol_weight_norm = self._set_solenoidal_weight_normalisation()

    def get_turb_vector_unigrid(self):
        """
        Generate turbulent field on a uniform grid.
        """

        # Precompute amplitude including normalisation factors for each mode.
        ampl = 2. * self.sol_weight_norm

        # Predefine trigonometric values for more efficient use of time.
        # These will have shape==(self.nmodes,self.N[d]) for d=0,1,2
        sinxi = np.sin(np.tile(self.mode_arr[:,0],(self.N[0],1)).T * (self.pos_beg[0]+np.tile(np.arange(self.N[0]),(self.nmodes,1))*self.delt[0]))
        cosxi = np.cos(np.tile(self.mode_arr[:,0],(self.N[0],1)).T * (self.pos_beg[0]+np.tile(np.arange(self.N[0]),(self.nmodes,1))*self.delt[0]))
        try:
            sinyj = np.sin(np.tile(self.mode_arr[:,1],(self.N[1],1)).T * (self.pos_beg[1]+np.tile(np.arange(self.N[1]),(self.nmodes,1))*self.delt[1]))
            cosyj = np.cos(np.tile(self.mode_arr[:,1],(self.N[1],1)).T * (self.pos_beg[1]+np.tile(np.arange(self.N[1]),(self.nmodes,1))*self.delt[1]))
        except:
            sinyj = np.zeros((len(self.nmodes,self.N[1])))
            cosyj = np.ones((len(self.nmodes,self.N[1])))
        try:
            sinzk = np.sin(np.tile(self.mode_arr[:,2],(self.N[2],1)).T * (self.pos_beg[2]+np.tile(np.arange(self.N[2]),(self.nmodes,1))*self.delt[2]))
            coszk = np.cos(np.tile(self.mode_arr[:,2],(self.N[2],1)).T * (self.pos_beg[2]+np.tile(np.arange(self.N[2]),(self.nmodes,1))*self.delt[2]))
        except:
            sinzk = np.zeros((len(self.nmodes,self.N[2])))
            coszk = np.ones((len(self.nmodes,self.N[2])))

        # print(sinxi.shape, sinyj.shape, sinzk.shape)
        # print(cosxi.shape, cosyj.shape, coszk.shape)
        
        # # Generate grid of turbulent velocities.
        # # Define real and imaginary components of
        # # e^{i \vec{k} \cdot \vec{x}} = cos(kx*x + ky*y + kz*z) + i sin(kx*x + ky*y + kz*z)
        sinxi = np.tile(sinxi[:,:,np.newaxis,np.newaxis],(1,1,self.N[1],self.N[2]))
        cosxi = np.tile(cosxi[:,:,np.newaxis,np.newaxis],(1,1,self.N[1],self.N[2]))
        sinyj = np.tile(sinyj[:,np.newaxis,:,np.newaxis],(1,self.N[0],1,self.N[2]))
        cosyj = np.tile(cosyj[:,np.newaxis,:,np.newaxis],(1,self.N[0],1,self.N[2]))
        sinzk = np.tile(sinzk[:,np.newaxis,np.newaxis,:],(1,self.N[0],self.N[1],1))
        coszk = np.tile(coszk[:,np.newaxis,np.newaxis,:],(1,self.N[0],self.N[1],1))

        real = (cosxi*cosyj - sinxi*sinyj) * coszk - (sinxi*cosyj + cosxi*sinyj) * sinzk
        imag = (cosxi*sinzk + sinyj*coszk) * cosxi + (cosyj*coszk - sinyj*sinzk) * sinxi

        vx = ampl * (np.transpose(np.tile(self.aka[:,0],(*self.N,1)),(3,0,1,2))*real - np.transpose(np.tile(self.akb[:,0],(*self.N,1)),(3,0,1,2))*imag).sum(axis=0)
        vy = ampl * (np.transpose(np.tile(self.aka[:,1],(*self.N,1)),(3,0,1,2))*real - np.transpose(np.tile(self.akb[:,1],(*self.N,1)),(3,0,1,2))*imag).sum(axis=0)
        vz = ampl * (np.transpose(np.tile(self.aka[:,2],(*self.N,1)),(3,0,1,2))*real - np.transpose(np.tile(self.akb[:,2],(*self.N,1)),(3,0,1,2))*imag).sum(axis=0)

        self.v = np.array([vx,vy,vz])

    def init_driving(self,):

        k_min = self.K_MIN * 2*PI/self.L[0]
        k_max = self.K_MAX * 2*PI/self.L[0]
        k_mid = k_max # Driving does not support different power exponents

        self._set_turnover_time() # Set t_decay
        self._set_dt() # Set dt
        self.step = -1 # Initiate step at -1

        ampl_coeff = 0.15   # This is the default amplitude coefficient ('ampl_coeff') that often
                            # (based on Mach ~ 1, naturally-mixed driving) leads to a good match of the
                            # user-to-target velocity dispersion.

        energy = (ampl_coeff*self.VELOCITY)**3. / self.L[0] # Energy input rate => driving amplitude ~ sqrt(energy/t_decay).
                                                            # Note that energy input rate ~ velocity^3 / L_box.
        self.OUvar = np.sqrt(energy/self.t_decay) # Set Ornstein-Uhlenbeck (OU) variance
        
        # Initiate modes.
        self.init_modes(k_min,k_max,k_mid) # Defines self.mode_arr and self.ampl_arr given self.RANDOM_SEED
        self._OU_noise_init(self.OUvar) # Define OUphases scaled by OUvar
        self._get_decomposition_coeffs() # Defines coefficients ready to define a turbulent grid.

    def _OU_noise_update(self,OUvar):
        """
        Updates the OU factor by the sequence
            x_n+1 = f x_n + sigma * sqrt(1-f**2) z_n
        where f = exp(-dt/ts) and z_n is a Gaussian random
        variable drawn from a Gaussian with unit variance, and
        sigma is the desired variance of the OU sequence.
        """

        damping_factor = np.exp(-self.dt/self.t_decay)
        OUphases = self.OUphases
        self.OUphases = OUphases * damping_factor + (1-damping_factor**2)**.5 * OUvar * np.random.normal(0.,1.,(self.nmodes,self.ncmp,2))

    def check_for_update(self,time,return_steps=[-1]):
        """
        Check to update the field.

        Parameters
        ----------
            time : float
                The simulation time.
        """

        # Automatic amplitude adjustment.
        # if self.AMPL_AUTO_ADJUST==1 and v_turb[0]>0:
        #     v_turb_for_ampl_adjust = v_turb

        #     if time>.1*self.t_decay:
        #         ampl_factor_new = self.AMPL_FACTOR * (self.VELOCITY/v_turb_for_ampl_adjust/self.ncmp)**1.5
        #         ampl_adjust_timescale = (ampl_factor_new/self.AMPL_FACTOR)**.5 * self.t_decay
        #         damping_factor = np.exp(-dt/ampl_adjust_timescale)
        #         ampl_factor_adjusted = self.AMPL_FACTOR * damping_factor + (1.0-damping_factor) * ampl_factor_new
        # print(return_steps)

        v_t = []

        step_requested = int(np.floor(time / self.dt))
        for step in range(self.step,step_requested):
            self.step += 1
            self._OU_noise_update(self.OUvar)
            if self.step in return_steps:
                self._get_decomposition_coeffs()
                self.get_turb_vector_unigrid(self.pos_beg,self.pos_end)
                v_t.append(self.v.copy())
        self._get_decomposition_coeffs()
        self.get_turb_vector_unigrid(self.pos_beg,self.pos_end)
        v_t.append(self.v.copy())

        return np.array(v_t)
        

        # print(ampl_factor_new,ampl_adjust_timescale,damping_factor,ampl_factor_adjusted)

    def evolve_turb_grid(self,t_end,return_steps=[-1]):
        """
        Initiate a turbulent field, then evolve it.
        """

        self.init_driving()
        v_t = self.check_for_update(t_end,return_steps)

        return v_t


def main():

    import os, struct
    
    # Set up file and directories

    save_dir = "./turbulent field"

    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    parameter_file = os.path.join(save_dir,"parameters.txt")

    # Generate turbulent field

    tg = TurbGen(parameter_file)

    tg.init_single_realisation()
    pos_beg, pos_end = tg._set_pos_beg_end()
    tg.get_turb_vector_unigrid()

    # Shift and scale velocity field

    mean = np.mean(tg.v,axis=(1,2,3))
    mean2 = np.mean(tg.v**2,axis=(1,2,3))

    std = (mean2-mean*2)**.5
    v = np.array([(_v-_mean)/_std for _v,_mean,_std in zip(tg.v,mean,std)])

    # Save turbulent velocities and associated grid0.out file for PLUTO

    f_names = ["vx10.dbl", "vx20.dbl", "vx30.dbl"]
    n_points = np.prod(tg.N)

    for i,f_name in enumerate(f_names):
        with open(os.path.join(save_dir, f_name),'wb') as f_o:   
            f_o.write(struct.pack('<'+'d'*n_points,*((v[i]).flatten())))

    with open(os.path.join(save_dir, "grid0.out"),'w') as f_out:
        f_out.write("# GEOMETRY:   CARTESIAN\n")
        for d in range(tg.NDIM):
            f_out.write(f"{tg.N[d]}\n")
            for i in range(tg.N[d]):
                # xL = -tg.L[d]/2 + (tg.L[d]*(i))/(tg.N[d])
                # xR = -tg.L[d]/2 + (tg.L[d]*(i+1))/(tg.N[d])
                xL = -tg.L[d]/2 + (tg.L[d]*(i-.5))/(tg.N[d]-1)
                xR = -tg.L[d]/2 + (tg.L[d]*(i+.5))/(tg.N[d]-1)
                if i==0 or i==tg.N[d]-1:
                    print(xL,xR)
                f_out.write("{} {:.6f} {:.6f}\n".format(i+1,xL,xR))


if __name__=="__main__":
    main()




        
    
# %%
