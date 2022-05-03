""" Partial histogram algorithm file"""

# python libraries
import cProfile
import logging
import multiprocessing as mp

# third-party libraries
import numpy as np
import tqdm

import matplotlib.pyplot as plt

from Utility.sipm import spad_wise_profiler

# first-party libraries
try: 
    import Utility.arrival_checker as avc
    import Utility.event_generator as ev
    from Utility.sipm import spad_wise as spad_wise
    from Utility.tdc import TDC
    from Utility.PostPorcess import PostProcess

except ModuleNotFoundError: # means you're probably trying to run partial_hist test code on its own
    import os
    import sys
    sys.path.append(os.path.abspath('..'))  # add directory above to path
    import Utility.arrival_checker as avc
    import Utility.event_generator as ev
    from Utility.tdc import TDC
    from Utility.PostPorcess import PostProcess


logger = logging.getLogger(__name__)
# logger.setLevel(level=logging.DEBUG)

class PHR:

    # read/write properties
    period = property(lambda self: self.__period, lambda self, val: self.__set_hist_params(val, None, None))
    n_bits = property(lambda self: self.__n_bits, lambda self, val: self.__set_hist_params(None, val, None))
    n_steps = property(lambda self: self.__n_steps, lambda self, val: self.__set_hist_params(None, None, val))
    # read-only properties 
    n_bit_per_hist = property(lambda self: self.__n_bit_per_hist)
    tdc_res = property(lambda self: self.__tdc_res)


    def __init__(self, period, n_bits=10, n_steps=3):
        ''' Partial Histogram Class
        Defaults run setup identical to Zhang et. al. 2019, where the algorithm is taken from

        Args:
            period (float): laser period in seconds, typically order ~ns
            n_bits (int): number of bits in the TDC
            n_steps (int): number of partial histogram steps
        '''
        logger.info('Initiating PHR instance')

        self.__period = period # total laser period
        self.__set_hist_params(self.__period, n_bits, n_steps)

        #self.__phr_res = [7, 4, 0]  # smallest bit for each phr window

    def __set_hist_params(self, period, n_bits, n_steps):
        ''' Set the histogram parameters, in particular:
            - set the class parameters, if provided
            - set the number of bits per histogram iteration

        Args:
            period (float): laser period in seconds, typically order ~ns
            n_bits (int): number of bits in the TDC
            n_steps (int): number of partial histogram steps
        '''
        if period is not None: self.__period = period
        if n_bits is not None: self.__n_bits = n_bits
        if n_steps is not None: self.__n_steps = n_steps

        # set the minimum TDC resolution
        self.__tdc_res = self.__period/(2**self.__n_bits) 

        # calculate the bits per step and remainder
        base = self.__n_bits//self.__n_steps
        remainder = self.__n_bits % self.__n_steps 
    
        # set the bits per step in an array
        self.__n_bit_per_hist = [base]*self.__n_steps # set the base number of bits per step
        self.__n_bit_per_hist[-1] += remainder  # add the remainder to the last step

def run_algorithm(PHR, N, const, lambda_bg, lambda_sum, sipm_arr, cd_arr, first_only=True): # need to add more inputs
    ''' Partial histogram algorithm from Zhang et. al. 2019
        First a histogram of bits [9:7] is constructed, then [6:4], then [3:0]
    Args
        N (int): number of iterations to perform
        const (constants.constants): object full of constants from main sim -> should replace
        lambda_sum (np.array): arary of photon arrival rates
        sipm_arr (sipm.SiPM): array of silicon photomultipler objects
        cd_arr (coincidence_detector.CD): array of coincidence detectors
    '''
    profiler = False
    iterationProfiler = True
    ## PRE-ALLOCATED ARRAYS
    # copied definitions from montesim.py temporarily
    
    _all_cd_t = [[]] * const.n_spad_per_tdc 
    _all_cd_s = [[]] * const.n_spad_per_tdc
    arr_t = [[]] * const.n_spad_per_tdc 
    arr_s = [[]] * const.n_spad_per_tdc
   
    #Create TDC Objects (defines how many TDC and what the deadline is)
    tdc_a = TDC(const.n_tdc_per_gp, const.n_spad_per_tdc, const, first_only=first_only)

    # Instantiate Post Process module (After TDC)
    postProcModule = PostProcess(const.sim_res, const.n_spad_per_tdc, const.T, N, const.postAlg) #It instantiates memory as well
    ## CONSTANTS FOR THE SIMULATION
    num_iter_per_phr = N//PHR.n_steps # number of iterations to perform for each PHR
    pool = mp.Pool(processes=const.n_spad_per_tdc) # set up the number of processes in the pool

    ## VARIABLES FOR THE SIMULATION
    bin_offset = np.zeros(const.n_spad_per_tdc) # init bin offset to zero for each SiPM
    bin_size = int(const.T//const.sim_res) # init bin size to whole window
    lambda_sum_part = lambda_sum # init to whole array for first iteration #!should always be whole

    #For gating the arrivals if we want. 
    min_keep_index = np.zeros(const.n_spad_per_tdc)
    max_keep_index = np.zeros(const.n_spad_per_tdc) + bin_size

    #performance analysis
    if(iterationProfiler):
        p = cProfile.Profile()
        p.enable()
    for phr_iter in range(PHR.n_steps):

        # compute smallest tdc resolution in this window
        # i.e. with 10 total bits, if we've already done iterations of 3 and 3 bits, the smallest bit is
        #   10 - sum([3, 3]) = 4
        # hence the smallest tdc resolution is min(tdc_res) * (2**min_bit_res)
        bit_res = PHR.n_bits - sum(PHR.n_bit_per_hist[0:phr_iter+1])
        tdc_window_res = PHR.tdc_res * (2**bit_res)
        logger.debug(f'bit_per_hist: {PHR.n_bit_per_hist}')
        logger.debug(f'bit_res: {bit_res}')
        logger.debug(f'tdc_res: {PHR.tdc_res}')
        logger.debug(f'TDC Set to res: {tdc_window_res}')
        # compute the number of bins in this window
        phr_n_bins = 2**(PHR.n_bit_per_hist[phr_iter])
        logger.debug(f'starting with {phr_n_bins} total bins')
        # divide the last bin size by the current number of bins
        bin_size = int(bin_size//phr_n_bins)

        logger.debug(f'Partial histogram with TDC resolution {tdc_window_res*1e12} ps')



        # perform num_iter_per_phr to form each partial histogram
        for n in tqdm.trange(num_iter_per_phr, desc='partial-hist-iter-{}'.format(phr_iter), ascii=True):
            if const.mt==2:
                # assemble the arrays into a work package
                if profiler:
                    output = [pool.apply_async(spad_wise_profiler,
                                                args=(shared_macro_pix, const, lambda_bg, lambda_sum_part[shared_macro_pix], sipm_arr[shared_macro_pix], cd_arr[shared_macro_pix],first_only,)) 
                            for shared_macro_pix in range(const.n_spad_per_tdc)]
                    result = [p.get() for p in output]  # Collect processes' output
                else:
                    output = [pool.apply_async(spad_wise,
                                                args=(shared_macro_pix, const, lambda_bg, lambda_sum_part[shared_macro_pix], sipm_arr[shared_macro_pix], cd_arr[shared_macro_pix],first_only,)) 
                            for shared_macro_pix in range(const.n_spad_per_tdc)]
                    result = [p.get() for p in output]  # Collect processes' output

                for i in range(const.n_spad_per_tdc):  # store the outputs
                    arr_t[i] = result[i][0]
                    arr_s[i] = result[i][1]
                    sipm_arr[i].spad_avail = result[i][2] #since we have created a new process, they worked on a copy of these classes and have not updated the original. So we update the original
                    sipm_arr[i].spad_free_times = result[i][3]
            else:
                #Normal single sipm at a time (No multi-processing)
                for shared_macro_pix in range(const.n_spad_per_tdc):
                    arr_t[shared_macro_pix], arr_s[shared_macro_pix], _, _ = spad_wise(shared_macro_pix, const, lambda_bg, lambda_sum_part[shared_macro_pix], sipm_arr[shared_macro_pix], cd_arr[shared_macro_pix],first_only)

            #IF we are gating events let's delete them
            for i in range(const.n_spad_per_tdc):
                # masking the events falling outside bin of interest
                _all_cd_t[i] = [x for x in arr_t[i] if x >= min_keep_index[i] and x < max_keep_index[i]]
                _all_cd_s[i] = arr_s[i] #!do we need them?
            
            # check for TDC dead time (serial)

            ####
            if n == 0:#!why?
                continue

            #update tdc block deadtime
            tdc_a.tdc_dead_checker(_all_cd_t, const, postProcModule,time_res=tdc_window_res)

            # update TDC deadtime
            tdc_a.freeTimeUpdate_NextRount(const.T, const)

        if(iterationProfiler):
            p.disable()
            p.dump_stats(f'profile_it_{const.device}_{const.n_spads_per_pix}.prof')

        # define the collector array
        collector = np.zeros([const.n_spad_per_tdc, phr_n_bins])
        target_bin = np.zeros(const.n_spad_per_tdc)
        SiPM_hists_a = postProcModule.get_histogram()

        # partial histogram  technique (Zooming over maximum bin)
        # for each sipm, find the maximum bin
        for sipm in range(const.n_spad_per_tdc):
            for hist_bin in range(phr_n_bins):
                # calculate the max and min bin index for each
                min_bin_index = int(bin_offset[sipm] + hist_bin*bin_size)
                max_bin_index = int(bin_offset[sipm] + (hist_bin+1)*bin_size)
                collector[sipm,hist_bin] = np.sum(SiPM_hists_a[sipm, min_bin_index:max_bin_index])
            # the target bin is the max collector bin
            target_bin[sipm] = np.argmax(collector[sipm])
            logger.debug("target_bin for sipm-{} is: {}".format(sipm, target_bin[sipm]))
            # compute the points to keep EVENTS #! not in lambda_sum
            min_keep_index[sipm] = int(bin_offset[sipm] + target_bin[sipm]*bin_size)
            max_keep_index[sipm] = int(bin_offset[sipm] + (target_bin[sipm] + 1)*bin_size)

            # set the selected bin in lambda_sum_part
            #!In the approach used by Chao Zhang et. al. they just ignore incoming. 
            #! this does not result in zeroing lambda sum (The probability of having events does not change)
            ## lambda_sum_part[sipm, min_keep_index:max_keep_index] = lambda_sum[sipm, min_keep_index:max_keep_index] 
            # store the new bin offset for next iteration
            bin_offset[sipm] = bin_offset[sipm] + target_bin[sipm]*bin_size
            logger.debug(f'zooming on {min_keep_index[sipm]}:{max_keep_index[sipm]}, offset={bin_offset}')
    
    time_steps = np.linspace(0, const.T-const.sim_res,
                             num=int(const.T//const.sim_res)) 
    #let's save the histogram
    postProcModule.save_output(const.t_tof, time_steps, const.outstr, 0, N, const.systemConfig, secondary_offset = const.secondary_step,save = const.save_figures ,trial = const.trial_number)


    return SiPM_hists_a #, SiPM_hists_f, SiPM_hists_fncd

#used foe debugging this module
if __name__ == "__main__":
    pass
    # example code, could be added to top-level file
    #phr = PHR(100e-9, n_bits=10) # supply PHR with the TDC resolution
    #SiPM_hists_a, SiPM_hists_f, SiPM_hists_fncd = phr.run_algorithm(NN, const, lambda_sum_org, SiPM_array, CD_array)
