#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:07:43 2020

@author: Pooya Poolad
A poisson event generator.
"""
try:
    import numpy as np
    import numpy.random as rnd
    #import matplotlib.pyplot as plt
    from numba import jit,cuda
    from numba.cuda.random import create_xoroshiro128p_states, xoroshiro128p_uniform_float32

except ModuleNotFoundError: # means you're probably trying to run partial_hist test code on its own
    import os
    import sys
    sys.path.append(os.path.abspath('..'))  # add directory above to path
    import numpy as np
    import numpy.random as rnd
    #import matplotlib.pyplot as plt
    from numba import jit,cuda
    from numba.cuda.random import create_xoroshiro128p_states, xoroshiro128p_uniform_float32
from pdb import set_trace; 

import torch
import cupy as cp
# #@torch.no_grad()
# @torch.jit.script

def pytorch_generator(lambda_bg, lambda_sum, n_spad, sim_res, tof_echo, echo_split):
    #defualts to gpu 1// you might want to change this
    torch.cuda.set_device('cuda:0')
    torch.no_grad()
    cuda1 = torch.device('cuda:0')

    #copy to gpu
    signal_lambdax = torch.tensor(lambda_sum[0],device=cuda1,dtype=float)
    
    #calculate the summation on gpu
    signal_lambda = (signal_lambdax + lambda_bg) * sim_res
    
    #secondary reflections
    signal_lambda_secx = torch.tensor(lambda_sum[1], device=cuda1,dtype=float)
    signal_lambda_sec = (signal_lambda_secx + lambda_bg) * sim_res
    n_spad_per_sipm = n_spad
    if tof_echo:
        lambda_sum_org = torch.zeros([ len(signal_lambda), n_spad_per_sipm], device=cuda1)
        split_spad = np.rint(np.array(echo_split) * n_spad_per_sipm).astype(int)
        #TODO: Only 2 Echo supported here for now #
        lambda_sum_org[:,0:split_spad[0]] = torch.tile(signal_lambda.reshape(signal_lambda.size()[0],1),[1,split_spad[0]])
        lambda_sum_org[:,split_spad[0]:split_spad[0]+split_spad[1]] = torch.tile(signal_lambda_sec.reshape(signal_lambda_sec.size()[0],1),[1,split_spad[1]])

    else:
        lambda_sum_org = torch.tile((lambda_bg + signal_lambda) * sim_res,[n_spad_per_sipm,1], device=cuda1)

    #generate random numbers    
    rnd_list = torch.rand(lambda_sum_org.shape, device=cuda1) #torch.cuda.FloatTensor(lambda_sum_org.shape).uniform_()
    
    #get the event list
    ev_list = torch.where(rnd_list < lambda_sum_org)
    times, spads = ev_list[0].cpu().numpy(), ev_list[1].cpu().numpy()
    return (times, spads)

def rnd_gen_vect_random(p_tile):
    p_tile = p_tile
    rnd_list = rnd.random(p_tile.shape)
    return np.where(rnd_list < p_tile)

def rnd_gen_vect_dotile(lambda_bg, lambda_sum, n_spad, sim_res, tof_echo, echo_split):
    p_tile = create_lambda_t(lambda_bg, lambda_sum[0], n_spad, sim_res, tof_echo, lambda_sum[1], echo_split)

    return rnd_gen_vect_random(p_tile)

def create_lambda_t(lambda_bg, signal_lambda, n_spad_per_sipm, sim_res, tof_echo=0, signal_lambda_sec = None, echo_split= [0.8, 0.2]):

    if tof_echo:
        #lambda_sum_org_untiled = np.zeros([len(time_steps)])
        lambda_sum_org = np.zeros([ len(signal_lambda), n_spad_per_sipm])
        split_spad = np.rint(np.array(echo_split) * n_spad_per_sipm).astype(int)
        start_chunk = 0
        for echo,chunks in enumerate(split_spad):
            #first chunk is primary
            if echo == 0:
                signal_gamma = (signal_lambda + lambda_bg) * sim_res
            else:
                signal_gamma = ( signal_lambda_sec[echo] + lambda_bg) * sim_res
        
        #Tile the lambda because we have n_spad_per_sipm independent spads (Some have different distro)
            lambda_sum_org[:,start_chunk:start_chunk+chunks] = np.tile(signal_gamma.reshape(signal_gamma.size,1),[1,chunks])
            start_chunk += chunks
    else:
        #Tile the lambda because we have n_spad_per_sipm i.i.d. spad
        lambda_sum_org = np.tile((lambda_bg + signal_lambda) * sim_res,[n_spad_per_sipm,1])
    
    return lambda_sum_org


def rnd_gen_vect_pretiled(lambda_sum, n_spad, sim_res):
    #use this if you have tiled lambdas in the main. 
    return rnd_gen_vect_random(lambda_sum)

def dd_time_gen(dd_avg,std_dev,const,num_out=False):
    num = rnd.normal(loc=dd_avg,scale=std_dev)
    indexed = num//const.sim_res
    if(num_out):
        return int(indexed),num
    else:
        return int(indexed)

#For testing this module independently
if __name__ == "__main__":
    pass
    #Test program
    #launching cuda core
    # n_spad = 4
    # tof_echo = 1
    # echo_split = [0.8, 0.2]
    # lambda_sum = [1e12]
    # lambda_sig = np.zeros((600,1))
    # lambda_sig[100:120] +=  4e12
    # lambda_sum.append(lambda_sig)
    # lambda_sig_sec = np.zeros((600,1))
    # lambda_sig_sec[200:220] +=  5e11
    # lambda_sum.append(lambda_sig_sec)
    # sim_res = 10e-12
    # cuda_generator_launcher(lambda_sum, n_spad, sim_res, tof_echo, echo_split)

