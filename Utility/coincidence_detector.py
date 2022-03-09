import numpy as np
import logging

#logger = logging.getLogger(__name__)

class CD:

    cd_factor = property(lambda self: self._cd_factor, lambda self, val: setattr(self, '_cd_factor', val))
    cd_n_spad = property(lambda self: self._cd_n_spad, lambda self, val: setattr(self, '_cd_n_spad', val))
    cd_window = property(lambda self: self._cd_window, lambda self, val: setattr(self, '_cd_window', val)) #CD Window is actually window in sec/sim_res ?

    def __init__(self, cd_factor, cd_window,cd_n_spad):
        ''' Coincidence detector class
        Args:
            cd_factor (int): threshold for coincidence detection
            cd_window (cd_window): window in simulation timesteps over which to check for valid coincident events
        '''
        self.logger = logging.getLogger('CD_logger')
        self.logger.setLevel(logging.DEBUG)

        self.logger.info('Initiating CD instance')
        self._cd_factor = cd_factor
        self._cd_window = cd_window
        self._cd_n_spad = cd_n_spad

    def coincidence_checker(self, arrival_times, arrival_spads, all_events=True):
        ''' Check for valid coincidence in the provided array of times and SPADs
        I.e., if a list of events for SPADs 1-4 are provided, then these SPADs
        are equivalently fed into a coincidence detector.

        Implementation is based on parallel AND gates, with cd_window equal to pulse length. See Beer, 2018.

        Args:
            arrival_times(numpy.ndarray): array of photon arrival times at the sensor in simulation timesteps
            arrival_spads(numpy.ndarray): array of SPAD numbers corresponding to the above arrival times
                                          i.e. if arrival_time[0] = 35 and arrival_spad[0] = 3 then SPAD 3 had an arrival at 35 timesteps
            all_events(bool): Whether to return all events or only the first event

        Returns:
            (list, list): arrival_times and arrival_spads with only the events meeting coincidence detection
        '''
        if self._cd_factor == 1:
            # coincidence detection not enabled
            if all_events:
                # return everything
                #print('bypass CD')
                return arrival_times, arrival_spads
            else:
                return [arrival_times[0]], [arrival_spads[0]]

        arrival_cd = [] # output list for arrival times
        spad_cd = [] # output list for spads 

        spad_w = [] # list of spads currently in the window, could be made faster by omitting
        time_w = [] # list of times currently in the window

        spad_in_list = [0 for _ in range(self.cd_n_spad)]
        len_list = 0

        for ind, time in enumerate(arrival_times):
            # self.logger.debug('Read {} from {}'.format(time,ind))
            spad = arrival_spads[ind]   # get the current SPAD index

            # while (len(time_w) > 0) and (time > time_w[0] + self._cd_window):   # iterate until we're current
            while len_list > 0 and (time > time_w[0] + self._cd_window):
                # self.logger.debug('Window passed, delete first {}:{}'.format(time_w,spad_w))
                spad_in_list[spad_w[0]] = 0  # remove this spad from occupancy list
                del time_w[0] # delete the first time
                del spad_w[0] # delete the first spad
#                 time_w.clear() # Window past, start over
#                 spad_w.clear()
                len_list -= 1  # reduce the length of the list
            
            # if spad in spad_w:
            if spad_in_list[spad]:
                # spad has triggered an event in this window already, ignore it
                # need to do this after deleting stale times to ensure we don't miss any
                # self.logger.debug('mine, ignore it')
                continue

            # self.logger.debug('Add it to the window list')
            time_w.append(time) # add current time to list
            spad_w.append(spad) # add current SPAD to list
            len_list += 1

            # if len(spad_w) >= self._cd_factor: # CHECK THIS STEP, CAN LEAD TO VERY CLOSE ONES
            if len_list >= self._cd_factor:
                # self.logger.debug('FOUND @ {}. first in {}'.format(time_w[-1],time_w[0]))
                #found = True
                # if len(spad_w) > self._cd_factor:
                #     continue #ignore if widnow met but there are more events
                    #for now we ignore it, we can use them to increase our confidence though
                arrival_cd.append(time) # returns last time
                spad_cd.append(spad) # returns last SPAD
            
            if all_events == False and len(arrival_cd) > 0:
                # exit out early because we only want the first event
                return arrival_cd, spad_cd
                
        return arrival_cd, spad_cd
