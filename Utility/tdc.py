import numpy as np 
##import Utility.constants as const
import Utility.event_generator as ev 
import logging

class TDC:
    def __init__(self,ntdc,npix,const,first_only=False):
        self.logger = logging.getLogger('TDC_logger')
        self.logger.info('Initiating TDC instance')
        self.ntdc = ntdc
        self.availableTDC = ntdc
        self.occupiedTDC = 0
        self.npix = npix
        self.freeTime = np.zeros(ntdc)
        self.TDCState = np.zeros(ntdc,dtype=bool) #TDC busy state 0 free, 1 busy
        self.TDCMiss = np.zeros([(int(const.T//const.sim_res))])
        self.first_only = first_only

    def tdc_dead_checker(self,sipm_time, const, PProc,time_res=40e-12):
        self.logger.debug('unflatter sipm_time: {}'.format(sipm_time))
        if type(sipm_time) == list:
            _flat_list = []
            _flat_index = []
            for idx,sublist in enumerate(iter(sipm_time)):
                self.logger.debug(sublist)
                for item in sublist:
                    _flat_list.append(item)
                    _flat_index.append(idx)
            flat_list=np.array(_flat_list)
            flat_index=np.array(_flat_index)
            self.logger.debug('{}:{}'.format(_flat_list,_flat_index))
            sipm_time_sorted = np.sort(flat_list)
            sipm_time_idx = flat_index[flat_list.argsort()]
            self.logger.debug('sipm_time: {},{}'.format(sipm_time_sorted,sipm_time_idx))
            #np.array([item for sublist in sipm_time for item in sublist])
        else:
            flat_list = sipm_time
            sipm_time_sorted = sorted(flat_list)
            sipm_time_idx = flat_list.argsort()
            self.logger.debug('sipm_time: {},{}'.format(sipm_time_sorted,sipm_time_idx))
        for idx in range(len(flat_list)):
            arrival_time = sipm_time_sorted[idx]
            arrival_pixel = sipm_time_idx[idx]
            self.logger.debug('%d from %d',arrival_time,arrival_pixel)
            if arrival_time != -1: #There is an arrival
                self.update_tdc_avail(arrival_time)
                if self.availableTDC > 0: #There is TDC available
                    # print(self.TDCState)
                    # print(self.availableTDC)
                    # print('npwhere: ',np.where(self.TDCState == False))
                    which_to_allocate = np.where(self.TDCState == False)[0][0] #Alocate first available
                    self.logger.debug('Allocating to: %d', which_to_allocate)
                    self.freeTime[which_to_allocate] = arrival_time + const.tdc_dd_dist
                    self.TDCState[which_to_allocate] = True #occupy it
                    self.availableTDC -=  1
                    self.occupiedTDC += 1
                    # Dynamic resolution
                    #arrival_times_out.append(arrival_time)
                    #arrival_pixel_out.append(arrival_pixel)
                    PProc.run_alg(arrival_pixel, arrival_time, time_res)
                    #time_res_adjust_factor = int(time_res//const.sim_res)
                    #self.logger.debug('Adjust factor for time res:{}'.format(time_res_adjust_factor))
                    #arrival_time_adjusted = time_res_adjust_factor*int(arrival_time//time_res_adjust_factor)
                    #self.logger.debug('adjusted time for TDC res:{}'.format(arrival_time_adjusted))
                    #self.logger.debug('Range to be incr{}:{}'.format(arrival_time_adjusted,arrival_time_adjusted+time_res_adjust_factor))
                    #tdc_out_times[arrival_pixel,arrival_time_adjusted:arrival_time_adjusted+time_res_adjust_factor] += 1
                    #print('pix: {}, t:{}, {}'.format(arrival_pixel,arrival_time,tdc_out_times[arrival_pixel,arrival_time]))
                else:
                    #self.TDCMiss = np.append(self.TDCMiss,arrival_time)
                    self.TDCMiss[arrival_time] += 1
                    self.logger.debug('All busy, missed!')
                    #tdc_out_times[arrival_pixel] = np.append(tdc_out_times[arrival_pixel],-1) #Missed
            # else:
            #     #tdc_out_times[arrival_pixel] = np.append(tdc_out_times[arrival_pixel],-1) #Missed
            #     pass
        #return arrival_pixel_out,arrival_times_out

    def update_tdc_avail(self,current_time):
        for tdcs in range(self.ntdc):
            my_freetime = self.freeTime[tdcs]
            my_status = self.TDCState[tdcs]
            self.logger.debug('now: %d free: %d ',current_time, my_freetime)
            if (current_time > my_freetime) and my_status == True: #Am I even Busy ?
                self.logger.debug('Freed: %d',tdcs)
                self.availableTDC += 1
                self.occupiedTDC -= 1
                self.TDCState[tdcs] = False #Free it
                self.freeTime[tdcs] = 0.0 #Reset it

    def freeTimeUpdate_NextRount(self,Tmod,const):
        for tdcs in range(self.ntdc):
            if self.TDCState[tdcs] == True:
                if self.freeTime[tdcs] >=  int(Tmod//const.sim_res):
                    self.logger.debug("Async TDC busy. %d ",self.freeTime[tdcs])
                    self.freeTime[tdcs] = self.freeTime[tdcs] - int(Tmod//const.sim_res) #carrying tdc state to next round
                    self.logger.debug("%d",self.freeTime[tdcs])
                else:
                    self.logger.debug("Finished round. %d to: ",self.freeTime[tdcs])
                    self.freeTime[tdcs] = 0.0 #carrying tdc state to next round
                    self.logger.debug("%d",self.freeTime[tdcs])
                    self.TDCState[tdcs] = False
                    self.availableTDC += 1
                    self.occupiedTDC -= 1
