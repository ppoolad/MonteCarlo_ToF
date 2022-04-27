import cProfile
import numpy as np 
##import Utility.constants as const
import Utility.event_generator as ev 
import Utility.arrival_checker as avc
import logging
import cProfile
logger = logging.getLogger(__name__)
logger.setLevel(level=logging.DEBUG)
from numba import jit, cuda

#SIPM Class defines system's pixel structure. 
class SiPM:
    """
        Creating the SiPM object. This class defines the SiPM's stucture in our system

        ---------------
        Args:
            nspad: How many SPADs create a SIPM?
            const: Config file to read constant parameters (SPAD pitch, deadtime, etc.)
    """
    def __init__(self,nspad,const):
        print("SiPM Created")
        self.nspad = nspad
        self.free_time = np.zeros(self.nspad)
        self.logger = logging.getLogger('SPAD_logger')
        self.logger.setLevel(logging.DEBUG)
        self.logger.info('Initiating SiPM instance')
        self.miss_counts = []
        self.spad_recharge = const.spad_recharge
        self.spad_avail = [1 for _ in range(nspad)]
        self.spad_free_times = [-1 for _ in range (nspad)]

    """
    This function checks which of the incoming photons get detected by the SPADs. 
    ---------------------------
    args: 
        time_arr: A sorted list of times photons have hit the sensor's surface
        spad_arr: SPAD index for the time array above (e.g. the photon hit at t=2ps was on SPAD number 2)
        const: configurations
    """
    def all_arrival_valid_checker_list(self, time_arr, spad_arr, const):
        valid_arr = []
        valid_spad = []
        found_first = [0 for _ in range(self.nspad)]
        curr_step = 0
        curr_spad = 0
        last_step = 0
        timeline_length = int(const.T//const.sim_res)
        for i in range(len(time_arr)):
            curr_step = time_arr[i] #Read Time
            curr_spad = spad_arr[i] #Which SPAD?

                #vvvv this part might be unnecessary (test it) -- v If the spad get's recharged befor current event/
            if (self.spad_free_times[curr_spad] != -1) and (self.spad_free_times[curr_spad] <= curr_step): 
                if(self.spad_recharge): # SPAD gets re-charged
                    self.spad_avail[curr_spad] = 1 #Means it is availalbe 
                    self.spad_free_times[curr_spad] = -1
                elif((found_first[curr_spad] == False)): #SPAD does not get recharged but it's the first one!
                    self.spad_avail[curr_spad] = 1 #? Why we are making it available?
                    self.spad_free_times[curr_spad] = -1      
                else: #SPAD Should not get recharged and it's not the firts one : Discard!
                    if(sum(found_first) == self.nspad ): #all SPADs found their first
                        break
                    else:
                        continue              

            if self.spad_avail[curr_spad] >= 1: #SPAD is AVAILABLE
                # use one of the SPADs
                self.spad_avail[curr_spad] = 0 #Got it! Not available anymore :D
                found_first[curr_spad] = True  #In case there is no recharge shceme.
                valid_arr.append(curr_step)    #Add it to validated list.#! Potentially performance limiting
                valid_spad.append(curr_spad)
                # set dead for some time
                self.spad_free_times[curr_spad] = curr_step + ev.dd_time_gen(const.spad_dead, const.spad_sigma, const) #It will be dead until it gets recharged again. 
            last_step = curr_step

        for i in range(len(self.spad_free_times)):
            # the zero events will be automatically cleared on the next cycle
            self.spad_free_times[i] = max(-1, self.spad_free_times[i] - timeline_length) #Check if spad will remain dead for the next iteration!

            if self.spad_free_times[i] == -1: #It won't :)
                self.spad_avail[i] = 1 #Make it available
            else:
                self.spad_avail[i] = 0 #It's still dead        
        
        valid_arr = np.array(valid_arr) #! Copy! Potentially performance limiting
        valid_spad = np.array(valid_spad)

        return valid_arr, valid_spad, self.spad_avail, self.spad_free_times

#Run spad wise but wrap in around a profiler
#Useful to profile it separately since it is executed inside a process pool
def spad_wise_profiler(shared_macro_pix, const, lambda_bg, lambda_sum, SiPM_array, CD_array, all_events=True):

    p = cProfile.Profile()
    p.enable()
    _all_cd_t, _all_cd_s, nxt_spad_avail, nxt_spad_free_times = spad_wise(shared_macro_pix, const, lambda_bg, lambda_sum, SiPM_array, CD_array, all_events)
    p.disable()
    p.dump_stats(f'profile{shared_macro_pix}_{const.device}_{const.n_spads_per_pix}.prof')
    return _all_cd_t, _all_cd_s, nxt_spad_avail, nxt_spad_free_times

def spad_wise(shared_macro_pix, const, lambda_bg, lambda_sum, SiPM_array, CD_array, all_events=True):
    '''Processing function for each SPAD.
            An event list is created, checked against dead times, coincidence are detected'''
    logger.debug(
        '-------------MACRO PIXEL {}---------------'.format(shared_macro_pix))

    #if torch gpu is selected
    if const.device == 'torchGPU':
        #print("torch")
        event_ar_t, event_ar_spad = ev.pytorch_generator(lambda_bg, lambda_sum, const.n_spads_per_pix, const.sim_res, const.tof_echo, const.echo_split)
    
    #if cpu is selected
    else:
        event_ar_t, event_ar_spad = ev.rnd_gen_vect_dotile(lambda_bg, lambda_sum, const.n_spads_per_pix, const.sim_res, const.tof_echo, const.echo_split)

    #Validate against deadtime
    all_ar_valid_t, all_ar_valid_spad, nxt_spad_avail, nxt_spad_free_times = SiPM_array.all_arrival_valid_checker_list(
        event_ar_t, event_ar_spad, const)

    #Coincidence detection
    _all_cd_t, _all_cd_s = CD_array.coincidence_checker(all_ar_valid_t, all_ar_valid_spad, all_events=all_events)

    return _all_cd_t, _all_cd_s, nxt_spad_avail, nxt_spad_free_times