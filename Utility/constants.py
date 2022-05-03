#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 17:41:20 2020

@author: pooladpo

pyhisical constants and parameters
"""
from os.path import join, split, splitext
import numpy as np
import scipy.constants
import Utility.config_util as confutil
class constants:
    def __init__(self, ttof, _T,configPath='./configs/montesim_test_config_complete.yaml'):
        # Read configurations from YAML file
        config_requires = ["Lens", "Laser", "Sensor","System"]
        myArgs  = confutil.parse_yaml_config(configPath, required=config_requires)
        try:
            self.device = myArgs.System['Device']
        except:
            self.device = 'cpu'

        #should it save png figures?
        try:
            self.save_figures = myArgs.System['save_figures']
        except:
            self.save_figures = True
        # Setup configurations
        self.c = scipy.constants.c  #Speed of light
        self.T = _T #Total time
        self.t_tof = ttof
        self.dist = self.t_tof*self.c/2
        self.outstr = 'output' #default value
        self.mu = self.t_tof
        self.laser_pulse_length = myArgs.Laser['Laser_pulse_length'] #s FWHM
        self.laser_sigma = self.laser_pulse_length / (2*np.sqrt(2*np.log(2)))
        self.reflectivity = myArgs.System['object_reflectivity']
        self.mt = myArgs.System['mt']
        
        self.rpt_freq = 1/self.T
        self.fps = myArgs.System['frames_per_second'] #Hz
        self.N = int(self.rpt_freq // self.fps)
        ##How many MonteCarlo Trial?
        
        self.trial_number = 0
        try:
            self.nTrial = myArgs.System['n_montecarlo_trial']
        except:
            self.nTrial = 1

        self.t = 0.0 #?
        self.sim_res = myArgs.System['sim_res']
        self.tdc_sigma = myArgs.System['tdc_sigma']
        self.spad_sigma = myArgs.Sensor['Spad_sigma']
        try: 
            self.SPAD_DCR  = myArgs.Sensor['DCR']
        except:
            self.SPAD_DCR  = 30000 #30khz count per second

        self.etha_pde = myArgs.Sensor['Sensor_pdp']
        self.wvlngth  = myArgs.Laser['WaveLength']
        
        #self.e_photon = scipy.constants.h * scipy.constants.c / self.wvlngth #hc/lambda #2.091e-19
        self.sigma = np.sqrt(self.laser_sigma**2 + self.tdc_sigma**2 + self.spad_sigma**2)
        self.bg_lambda = myArgs.Solar['constant_rate']
        self.opticConfig = configPath #? Probably won't need
        try:
            self.weather_factor = myArgs.Solar['weather_factor']
        except:
            self.weather_factor = 1
        self.systemConfig = configPath
        # self.n_spads_per_pix = 1 # useless, is modified in update_CD
        self.n_spad_per_tdc = 1 # useless, is updated by the top level code to be equal to the number of TOF steps
        self.n_tdc_per_gp  = myArgs.System['n_tdc_per_cluster']
        self.tdc_dead = myArgs.System['tdc_dead_time']
        self.spad_dead = myArgs.Sensor['SPAD_deadtime']
        self.spad_recharge = myArgs.System['SPAD_recharge']
        self.spad_dd_dist = self.spad_dead//self.sim_res
        self.tdc_dd_dist = self.tdc_dead//self.sim_res
        self.CD_Window = myArgs.System['CD_Window'] #CD Stuff will be updated by update_CD method
        self.CD_Factor = myArgs.System['coincidence']
        self.overlapping_cd = myArgs.System['coincidence_auto_reset']
        self.cd_wdist = self.CD_Window//self.sim_res
        self.first_only = myArgs.System['first_only']
        self.alg = myArgs.System['alg']
        self.postAlg = myArgs.System['Post_Process_Alg']
        self.nhistogram_bins = myArgs.System['nbins']

        try:
            self.tf_order = myArgs.System['Tap_order']
        except:
            self.tf_order = 'inc'
        if(self.alg == 'partial-hist'):
            self.phr_iter = myArgs.System['phr_iter']
        else:
            self.phr_iter = 1

        #Multiple reflection support
        try:
            self.tof_echo = myArgs.System['ToF_Echo']
        except:
            self.tof_echo = 0

        #Distance of multiple reflections
        try:
            self.secondary_step = myArgs.System['secondary_step']
        except:
            self.secondary_step = 10e-9

        try:
            self.echo_split = myArgs.System['echo_split']
        except:
            self.echo_split = [0.8, 0.2]

        try:
            self.em_kernels = myArgs.System['EMKernels']
        except:
            self.em_kernels = 3

        try:
            self.buffer_size = myArgs.System['EMBatch']
        except:
            self.buffer_size = 200
        try:
            self.max_iter_cache = myArgs.System['EMIter']
        except:
            self.max_iter_cache = 100

        try:
        #Gateing
            self.gate_start = myArgs.System['gate_start']
            self.gate_end = myArgs.System['gate_end']
        except:
            self.start_gate = 0
            self.end_gate = self.T
        #########

    def update_CD(self,CD_W=2e-9,CD_Factor=[2,4]):
        self.CD_Window = CD_W
        self.CD_Factor = CD_Factor[0]
        self.n_spads_per_pix = int(CD_Factor[1])
        self.cd_wdist = self.CD_Window//self.sim_res

    def set_background(self,lambda_bg):
        self.bg_lambda = lambda_bg

    def set_laser(self,pmax,fwhm):
        self.laser_pulse_length = fwhm #s FWHM
        self.laser_sigma = self.laser_pulse_length / (2*np.sqrt(2*np.log(2)))
        self.laser_max_p = pmax #Watt
        self.laser_output_E = pmax* self.laser_sigma * np.sqrt(2*np.pi)
        self.laser_E = self.reflectivity * (self.laser_output_E /4 / np.pi / self.dist**2 / (np.tan(self.fov * np.pi/180)** 2)) * self.pix_area ##self.reflectivity * (self.P_30cm * (30e-2)**2 / (2*self.dist**2) * self.laser_pulse_length)

    def print_constants(self):
        attr_list = [a for a in dir(self) if not a.startswith('__') and not callable(getattr(self, a))]
        print("Constants configuration:")
        for attr in attr_list:
            print(f"\t{attr}: {getattr(self,attr)}")



class FixedPointNumber:
    ''' Fixed point Class
    This class handles fixed point implementation and conversions
    Args:
        fp: Floating point number
        n: number of bits 
    
    
    '''
    real = property(lambda self: self._real, lambda self, val: self.set_value(val))
    
    def __init__(self, fp, n, m=1) -> None:

        
        self.n_fraction_bits = n
        self.m_integer_bits = m
        self.real  = fp

    def get_quantized(self) -> float:
        fixed_point = (self.real/(2**-self.n_fraction_bits)).astype('int')
        quantized = fixed_point * 2**-self.n_fraction_bits
        return quantized
    def get_fp_integer(self) -> int:
        fixed_point = (self.real/(2**-self.n_fraction_bits)).astype('int')
        return fixed_point

    def set_value(self,real) -> None:
        self._real  = real
        #self.overflow = 0 if real < (2 ** self.m_integer_bits) else 1

    quantized = property(get_quantized)
    fixedpoint_int = property(get_fp_integer)