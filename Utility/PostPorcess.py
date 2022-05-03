""" Histogramming and post processing methods """
import numpy as np
import logging
import scipy.io as sio 
from shutil import copyfile as cp
import matplotlib.pyplot as plt
logger = logging.getLogger(__name__)
import itertools
from Utility.constants import FixedPointNumber as FPN, constants
import Utility.constants
#I know they might not
#warnings.filterwarnings(action='ignore', category=ConvergenceWarning)

class PostProcess:
    algorithm = property(lambda self: self.__alg, lambda self, val: self.__set_alg(val))
    sim_res = property(lambda self: self.__sim_res, lambda self, val: self.__set_sim_res(val))

    def __init__(self, sim_res, spad_per_tdc, Tmod, NN=10000, algorithm = 'Histogram') -> None:
        ''' Post Processin class 
        default , forms histogram according to the time& res given by the TDC

        Args:
            algorithm: What algorithm to use, defauls is histogramming
            sim_res: resulution the simulation is performed.
            spad_per_tdc: How many SPADs is assigned to TDC cluster
            TMod: Time between each illumination
        
        '''
        logger.info('Initiating PostProc instance')
        
        self.__set_alg(algorithm)
        self.__set_sim_res(sim_res)
        self.__spad_per_tdc = spad_per_tdc
        self.__tmod = Tmod
        self.__NN = NN
        self.n = np.zeros(spad_per_tdc) #times calculation is called
        #allocate memory
        
        if(algorithm == 'Histogram'):
            self.__hist_mem = np.zeros([spad_per_tdc, int(Tmod//sim_res)], dtype=int)
        elif(algorithm == 'All'):
            self.__hist_mem = np.zeros([spad_per_tdc, int(Tmod//sim_res)], dtype=int)
            self.__mem = np.zeros([spad_per_tdc], dtype=int) #int?
            self.__mem_back = np.zeros([spad_per_tdc], dtype=int) #int?
            self.convergence_track = []
            for _ in range(spad_per_tdc):
                self.convergence_track.append([])

            self._lambda = [1/4] * self.__spad_per_tdc
            self.dyn_lambda = 2000
            self.lambda_min = 1/1024
            self.arr_ctr = [0] * self.__spad_per_tdc
        else:
            self.__mem = np.zeros([spad_per_tdc], dtype=int) #int?
            self.__mem_back = np.zeros([spad_per_tdc], dtype=int)
            self.convergence_track = []
            for _ in range(spad_per_tdc):
                self.convergence_track.append([])
    def __set_alg(self, alg):
        self.__alg = alg

    def __set_sim_res(self, sim_res):
        self.__sim_res = sim_res

    def calculate_histogram_parameters(self,arrival_time,time_res):
        time_res_adjust_factor = int(time_res//self.__sim_res)
        arrival_time_adjusted = time_res_adjust_factor*int(arrival_time//time_res_adjust_factor)
        return time_res_adjust_factor, arrival_time_adjusted

    def update_histogram(self,arrival_pixel,arrival_time,time_res):
        time_res_adjust_factor, arrival_time_adjusted = self.calculate_histogram_parameters(arrival_time, time_res)
        self.__hist_mem[arrival_pixel,arrival_time_adjusted:arrival_time_adjusted+time_res_adjust_factor] += 1

    def get_histogram(self):
        return self.__hist_mem

    def iir_filter_config(self,iir_lambda):
        self._lambda = iir_lambda
        # memory:

    def iir_filter(self,arrival_pixel, arrival_time, time_res):
        # y[k] = (1-lambda)y[k-1] + lambda x[k]
        # lambda range: 2^0, 2^-1, 2^-2, ... (can be realized with right shift on hw)
        y_k = (1 - self._lambda[arrival_pixel]) * self.__mem[arrival_pixel] + self._lambda[arrival_pixel] * arrival_time
        self.arr_ctr[arrival_pixel] += 1
        if self.arr_ctr[arrival_pixel] > self.dyn_lambda and self.dyn_lambda != 0 and self._lambda[arrival_pixel] > self.lambda_min :
            self._lambda[arrival_pixel] /= 2
            self.arr_ctr[arrival_pixel] = 0
        #update
        self.__mem[arrival_pixel] = y_k
        self.convergence_track[arrival_pixel].append(y_k)

    ## EM
    
    def run_alg(self,arrival_pixel_out, arrival_time_out,time_res=40e-12):
        if self.__alg == 'Histogram':
            self.update_histogram(arrival_pixel_out, arrival_time_out,time_res)
        elif self.__alg == 'IIR':
            self.iir_filter(arrival_pixel_out, arrival_time_out, time_res)
        elif self.__alg == 'All':
            self.update_histogram(arrival_pixel_out, arrival_time_out,time_res)
            self.iir_filter(arrival_pixel_out, arrival_time_out, time_res)

    def save_output(self, tof_list, time_steps, outstr, isipy, NN, config_file, nbins=100, secondary_offset=10e-9, figformat='png',save=True,trial=0):
        if self.algorithm == 'Histogram':
            matdict = {}
            for nums in range(self.__spad_per_tdc):
                matdict['pixel{}_all_ar_tof{}'.format(
                    nums, int(tof_list[nums]*1e9))] = self.__hist_mem[nums, :]
                fig = plt.figure()
                plt.bar(time_steps,self.__hist_mem[nums, :],width=self.sim_res)
                plt.savefig('./{}/all_ar_tof_{}.{}}'.format(outstr,nums,figformat))
                if(not isipy):
                    plt.close(fig)
            # TODO: undefined because of earlier comment
            sio.savemat('./{}/N{}_strt{}_stp{}.mat'
                        .format(outstr, NN, int(tof_list[0]*1e9), int(1e9*tof_list[-1])), matdict)
        elif self.algorithm == 'IIR':
            pass

        elif self.algorithm == 'All':
            #For histogram 
            matdict = {}
            for nums in range(self.__spad_per_tdc):
                matdict['pixel{}_all_ar_tof{}'.format(
                    nums, int(tof_list[nums]*1e9))] = self.__hist_mem[nums]
                matdict[f'pixel{nums}_iir_convergance'] = self.convergence_track[nums]
                if(save):
                    time_line = time_steps #np.linspace(0,self.__tmod,num=nbins)
                    x_line = time_line # * 3e8 / 2 #lightSpeed
                    time_to_dist =  3e8 /2
                    figsize = (3,2)
                    fig, ax = plt.subplots(1,figsize=(6,2))
                    #plt.tight_layout()
                    ax.bar(x_line*time_to_dist, self.__hist_mem[nums],width=self.sim_res*time_to_dist)
                    ax.set_xlabel('Distance (m)',fontsize=12)
                    ax.set_ylabel('Counts',fontsize=12)

                    plt.suptitle('Histogram', fontsize=14)
                    plt.tight_layout()
                    plt.savefig('./{}/all_ar_tof{}_{}-{}+{}_{}_{}.'.format(outstr,nums,int(tof_list[0]*1e9),int(tof_list[-1]*1e9),figsize[0],figsize[1],trial) + figformat)

                    fig2 = plt.figure(figsize=figsize)
                    plt.plot(np.linspace(0,len(self.convergence_track[nums]),len(self.convergence_track[nums])), np.array(self.convergence_track[nums])*self.sim_res*3e8/2)
                    plt.axhline(tof_list[nums] *3e8/2 ,linestyle='--',color='red')
                    plt.axhline((tof_list[nums] + secondary_offset) *3e8/2 ,linestyle='--',color='red')
                    plt.ylabel('Distance (m)', fontsize=12)
                    plt.suptitle('IIR filter convergene in time', fontsize=14)
                    plt.xlabel('Sample', fontsize=12)
                    plt.tight_layout()
                    plt.savefig('./{}/convergeIIR_{}_{}_{}_{}.{}'.format(outstr,nums,figsize[0],figsize[1],trial,figformat))

                    if(True):
                        plt.close(fig)
                        plt.close(fig2)
            sio.savemat('./{}/N{}_strt{}_stp{}_{}.mat'
                        .format(outstr, NN, int(tof_list[0]*1e9), int(1e9*tof_list[-1]), trial), matdict)
            #plot convergence

        cp(config_file,f'./{outstr}/configs.yaml')       
