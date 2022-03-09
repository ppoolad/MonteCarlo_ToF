#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 18:08:06 2020

@author: pooladpo
Montecarlo Simulator
"""
import argparse
import numpy as np
import logging
from Utility.constants import constants
from Utility.sipm import SiPM
from Utility.coincidence_detector import CD
from Utility.partial_hist import PHR
from Utility.partial_hist import run_algorithm as phr_run
from Utility.optics import Optic_system
import cProfile
from pstats import Stats, SortKey
import Utility.config_util as confutil


#from IPython.core.debugger import set_trace #for debug
logging.getLogger('matplotlib').setLevel(logging.ERROR)
logging.getLogger('Utility.partial_hist').setLevel(logging.DEBUG)
logger = logging.getLogger(__name__)
logger.setLevel(level=logging.DEBUG)

def simulate(const, NN, nbins, mt, sv, outstr, alg, pos=0, windowMode=0, windowTaps=10,config_file='./configs/montesim_test_config_complete.yaml'):

    '''Main body of the simulator'''
    T = const.T  # Time Range
    time_steps = np.linspace(0, const.T-const.sim_res,
                             num=int(const.T//const.sim_res))  # Time Vector

    mu = const.t_tof  # ToF Vectors

    # allocate for signal powers (laser and the background)
    signal_lambda = np.zeros([len(mu), len(time_steps)])
    
    #Optic system object
    optics = Optic_system(config_file)

    #How much photons coming from the background
    lambda_bg = optics.calculate_lambda_bg()
    #How much photons (max) coming from the signal?
    lambda_sig_max = optics.calculate_lambda_sig(const.dist,const.reflectivity)

    #Also, we might have secondary targets
    lambda_sig_max_secondary = [[]]*const.tof_echo


    #How much photons coming from  the secondary targets? (Targets behind our target)
    for echo in range(const.tof_echo):
        #calculate maximum peak reflected by the scondary reflections
        lambda_sig_max_secondary[echo] = optics.calculate_lambda_sig(const.dist + ((echo+1)*const.secondary_step)*const.c/2,const.reflectivity)

    # Now lets create the \lambda(t) envelope for the signals.
    lambda_sum_org = np.zeros([len(mu), const.tof_echo+1, len(time_steps)]) #Let's out them all here

    for idx, tofs in enumerate(mu):
        #primary reflection
        lambda_sum_org[idx,0,:] = optics.create_gamma_t_envelope(lambda_sig_max[idx], time_steps, tofs, optics._Laser.laser_sigma) #Create the gaussian laser temporal response
        for echo in range(const.tof_echo):
            #secondary ones
            lambda_sum_org[idx,echo+1,:] = optics.create_gamma_t_envelope(lambda_sig_max_secondary[echo][idx], time_steps, tofs + (echo+1)*const.secondary_step, optics._Laser.laser_sigma) #for the reflections if any

        #now if we store it somewhere we won't need to copy it in event generator
        #UNCOMMENT HERE IF YOUR WANT TO SEE THE LAMBDA
        # fig, ax= plt.subplots(1, figsize=(6,4))
        # ax.plot(time_steps,lambda_sum_org[0,0,:] + lambda_bg)
        # ax.set_ylim(2e8)
        # ax.set_xlabel('t (s)',fontsize=12)
        # ax.set_ylabel('$\lambda$ (ph/s)',fontsize=12)

        # #ax.legend(['Reconstruction','Histogram (Normalized)'])
        # plt.suptitle('Photon arrival rate')
        # plt.tight_layout()
        # plt.savefig('./{}/lambda_0.png'.format(outstr))
    
    #Create SIPM Objects
    SiPM_array = []
    
    for _ in range(const.n_spad_per_tdc):
        SiPM_array.append(SiPM(const.n_spads_per_pix, const))

    #Coincidence Detector Module Object
    #Each SIPM need a CD module
    CD_array = [CD(const.CD_Factor,const.cd_wdist,const.n_spads_per_pix)
                for _ in range(const.n_spad_per_tdc)]

    if NN == 0:
        NN = const.N
        print('going for full sim')

    #How many bits is the TDC code?
    tdc_nbits = int(np.log2(nbins))

    #partial histogramming post processing technique.
    #This unit can to simple histogram acquisition + partial histogram technique proposed by Zhang et. al. (JSSC 2019)
    phr = PHR(const.T, n_bits=tdc_nbits, n_steps=const.phr_iter)
    #run algorithm
    #it also saves data
    SiPM_hists_a = phr_run(phr,NN, const, lambda_bg, lambda_sum_org, SiPM_array, CD_array, first_only=const.first_only)


#TODO: Matbe find a better way of passing all the parameters? why not parsing args here?
def main(alg, n_steps=5, start=10, step_res=1e-8, Tmod=50e-9, N=1000, DBG='', nbins=100, mt=0, sv=1, outstr='outputs', cd=[2, 4], lambda_bg=100, WinMode=0, WinTaps=10,config_file='./configs/montesim_test_config_complete.yaml'):
    cpu_count = 39  # Let's keep 1 core for background tasks.
    print('Using {} CPU cores'.format(cpu_count))
    tof_steps = np.arange(start=start*1e-9, step=step_res,
                            stop=start*1e-9+n_steps*step_res-10e-12)
    myconst =  constants(tof_steps, Tmod,config_file)
    myconst.n_spad_per_tdc = n_steps
    myconst.update_CD(CD_Factor=cd)  # this also updates the 
    myconst.set_background(lambda_bg)
    myconst.outstr = outstr
    print('to be tested: ', tof_steps)
    myconst.print_constants()
    simulate(myconst, N, nbins, mt, sv, outstr, alg,
                windowMode=WinMode, windowTaps=WinTaps,config_file=config_file)


if __name__ == '__main__':
    # set up the argument parser
    parser = argparse.ArgumentParser(description='Process some integers.')


    subparser = parser.add_subparsers(dest="subparser_name")

    # config parser
    parser_cfg = subparser.add_parser('config', help='Launch the simulation with a config file')
    parser_cfg.add_argument("-f",
                        "--filename",
                        dest="filename",
                        required=True,
                        help="Path to a yaml config file for simulation")
    parser_cfg.add_argument("-o",
                        "--output_dir",
                        dest="outfolder",
                        required=True,
                        help="Path to output directory")
    parser_cfg.add_argument("-p",
                        "--profiling",
                        dest="do_profiling",
                        action="store_true",
                        help="profile the code")

    parser_cfg.set_defaults(do_profiling=False)
    # required parameters in the config file, should match cmdline args
    config_requires = ["System"]
    # set up the command line arguments #no more command line, use YAML File

    # parse the arguments
    args = parser.parse_args()
    outfolder = ''
    if args.subparser_name == None:
        parser.print_help()
    elif args.subparser_name == 'config':
        # overwriting args could be a mistake, but we'll try it
        fileName= args.filename
        outfolder= args.outfolder
        do_profiling = args.do_profiling
        sys_args = confutil.parse_yaml_config(args.filename, required=config_requires)
        args = sys_args
    else:
        raise("Commanline config is not supported anymore, use the YAML file")

    argvs = vars(args)
    argv=argvs['System']
    #getting the default output folder if not defined
    if outfolder != '' :
        argv['outfolder'] = outfolder

    # set up debugging
    DBG = argv['dbg_log']
    if DBG != '':
        logging.basicConfig(level=logging.DEBUG, filename=DBG, filemode='w+',
                            format='%(asctime)s - %(name)s - %(levelname)s -\t %(message)s')
    logger.info('Hi!') # hello

    # convert fraction string into regular integer    
    cd_prep =   argv['coincidence'].split('/')
    cd_factor = [int(i) for i in cd_prep]
    print("simulating %d steps of %gs until %gs for %d iteration" %
          (argv['steps'], argv['step_res'], argv['Tmod'], argv['NIter']))
    
    # call the main function
    #If we want to measure performance of the mother-process
    if do_profiling:
        with cProfile.Profile() as pr:
            main(alg=argv['alg'], n_steps=argv['steps'], start=argv['start_tof'], step_res=argv['step_res'], Tmod=argv['Tmod'], N=argv['NIter'],
                DBG=argv['dbg_log'], nbins=argv['nbins'], mt=argv['mt'], sv=argv['sv'], outstr=argv['outfolder'], cd=cd_factor,
                lambda_bg=sys_args.Solar['constant_rate'], WinMode=argv['Windowing'], WinTaps=argv['WindowTaps'],config_file=fileName)

        with open('profiling_stats.txt', 'w') as stream:
            stats = Stats(pr, stream=stream)
            stats.strip_dirs()
            stats.sort_stats('time')
            stats.dump_stats('.prof_stats')
            stats.print_stats()        

    else:
        main(alg=argv['alg'], n_steps=argv['steps'], start=argv['start_tof'], step_res=argv['step_res'], Tmod=argv['Tmod'], N=argv['NIter'],
            DBG=argv['dbg_log'], nbins=argv['nbins'], mt=argv['mt'], sv=argv['sv'], outstr=argv['outfolder'], cd=cd_factor,
            lambda_bg=sys_args.Solar['constant_rate'], WinMode=argv['Windowing'], WinTaps=argv['WindowTaps'],config_file=fileName)
