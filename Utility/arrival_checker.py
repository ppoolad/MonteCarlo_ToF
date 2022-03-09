#Arrival checker

import numpy as np 
##import Utility.constants as const
import Utility.event_generator as ev 
import logging

module_logger = logging.getLogger('arrival')
def first_arrival_valid_checker(event_ar_t,event_ar_spad):
    #for spads in range(const.n_spads_per_pix):
        #which_are_mine = np.where(event_ar_spad == spads)[0]
    events_of_0 = np.where(event_ar_spad == 0)[0]
    if len(events_of_0) >= 1:
        events_of_0 = events_of_0[0]
        zeroth_spad_event = event_ar_t[events_of_0]
    else:
        zeroth_spad_event = -1
    _,uniques_idx = np.unique(event_ar_spad,return_index=True)
    uniques_idx = np.sort(uniques_idx)
    event_far_spad = event_ar_spad[uniques_idx]
    event_far_t = event_ar_t[uniques_idx] #get only first one
    return event_far_t, event_far_spad, zeroth_spad_event


def all_arrival_valid_checker(event_ar_t,event_ar_spad,const):
    valid_event_t = np.array([],'int')
    valid_event_spad = np.array([],'int')
    for spads in range(const.n_spads_per_pix): #iterate over each spad
        which_are_mine = np.where(event_ar_spad == spads)[0]
        if(len(which_are_mine) > 0): # Do I have any ?
            my_arrivals = event_ar_t[which_are_mine]
            #print(spads,my_arrivals)
            #valid_times = np.array([my_arrivals[0]]) #first always valid
            #valid_event_t = np.append(valid_event_t,valid_times)
            #time_unitl_last_valid = 0.0
            if(len(my_arrivals) > 1): #I have more than 1 arrival
                finish = False
                idx = 1
                while not(finish):
                    my_time = my_arrivals[idx-1] #I am the reference
                    other_times = my_arrivals[idx:] - my_time #distance of others with me
                    not_valids = np.where(other_times < ev.dd_time_gen(const.spad_dead, const.spad_sigma,const))[0]+idx #Who is too close?
                    my_arrivals = np.delete(my_arrivals, not_valids) #delete those from the list
                    idx += 1 #next
                    if(idx >= len(my_arrivals)):
                        finish = True #are we done ?
            valid_event_t = np.append(valid_event_t,my_arrivals)
            valid_event_spad = np.append(valid_event_spad, [(spads)]*len(my_arrivals))
    
    p = valid_event_t.argsort()
    valid_event_t.sort()
    q = valid_event_spad[p]
    return valid_event_t,q

def coincidence_detector(CD_factor, CD_Window, arrival_times, arrival_spads):
    #We have a list of event, sorted by time and spad number. we want to check in CDf happens in CDw. 
    if CD_factor == 1:
        found = True
        return arrival_times, arrival_spads, found
    else:
        found = False
        cd_time = -1
        cd_spads = -1
        for idx,times in enumerate(arrival_times):
            who_am_i = arrival_spads[idx]
            #print('i am: ', who_am_i, idx)
            find_other_arrivals = np.where(arrival_spads == who_am_i)[0]
            #print('mu others: ',find_other_arrivals)
            exclude_my_others = np.delete(arrival_times,find_other_arrivals[1:])[idx:]
            exclude_my_others_sp = np.delete(arrival_spads,find_other_arrivals[1:])[idx:]
            #print(arrival_times,arrival_spads,'\n',exclude_my_others)
            distance_to_others = exclude_my_others - times
            how_many_in_window = np.where(distance_to_others <= CD_Window)[0]
            #print(how_many_in_window)
            if(len(how_many_in_window) >= CD_factor):
                # print(how_many_in_window)
                #print('Found coincidence!',idx)
                # print(arrival_times)
                cd_time = arrival_times[idx+CD_factor-1]
                cd_spads = exclude_my_others_sp[how_many_in_window]
                found = True
                break
        return cd_time,cd_spads,found

