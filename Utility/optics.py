'''
This class provides means to calculate environment conditions
and photon arrival rates based on the lens properties and SPAD specs.


The inputs to the simulator are:
TBD

Depends on
    - Field of View
    - F Number
    - Solar irradiance
    - ..
Author: Pooya Poolad
Created: 2021-03-17

'''
#%%
import numpy as np
import logging
import argparse
import Utility.config_util as confutil
import Utility.logging_util as logutil # just for easy logger formatting
import pandas as pd
import matplotlib.pyplot as plt
import scipy.constants as const
logger = logutil.setup_logging(__name__, debug_level=logutil.INFO)
#%%
class SolarModel:
    def __init__(self,config):
        self.constant_rate = config['constant_rate']
        if self.constant_rate <= 0:
            try:
                path = config['path_to_irradiance']
            except:
                path = './Utility/Data/astmg173.xlsx'

            try: 
                self.weather_factor = config['weather_factor']
            except:
                self.weather_factor = 1

            spectral_irradiance = pd.read_excel(path,'SMARTS2',dtype=np.float64,header=1, usecols='A:D', nrows=1500)
            keys = spectral_irradiance.keys()
            self.direct_irradiance = spectral_irradiance[keys[3]] # W*m-2*nm-1
            self.tilted_irradiance = spectral_irradiance[keys[2]] # W*m-2*nm-1
            if config['sun_angle'] == 'direct':
                self.target_irradiance = self.direct_irradiance
            else:
                self.target_irradiance = self.tilted_irradiance
            self.wavelengths = spectral_irradiance[keys[0]] #nm
        
    def plotDirect(self):
        plt.figure()
        plt.plot(self.wavelengths, self.direct_irradiance)
        plt.title('Direct Solar Spectral Irradiance')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Irradiance (W/m2/nm)')
        
    def plotIndirect(self):
        plt.figure() 
        plt.plot(self.wavelengths, self.tilted_irradiance)
        plt.title('Indirect Solar Spectral Irradiance')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Irradiance (W/m2/nm)')



class Laser:
    def __init__(self,config):
        self.P_max = config['Laser_Pmax']
        self.WaveLength = config['WaveLength']
        self.Beam_Div_h = config['Beam_Div_h']
        self.Beam_Div_v = config['Beam_Div_v']
        self.Beam_shape = config['Beam_Shape']
        self.Laser_pulse_length = config['Laser_pulse_length']
        self.laser_sigma = self.Laser_pulse_length / (2*np.sqrt(2*np.log(2)))

    def calculatePeakPowerAtTarget(self,Distance_vector):
        #self.laser_output_E = self.P_max* self.laser_sigma * np.sqrt(2*np.pi)
        #self.laser_E = (self.laser_output_E / np.pi / Distance_vector**2 / (np.tan(self.Beam_Div * np.pi/180)** 2)) ##self.reflectivity * (self.P_30cm * (30e-2)**2 / (2*self.dist**2) * self.laser_pulse_length)
        if self.Beam_shape == 'square':
            #Etx = Etx/(4.d2.tan(fov_h)tan(fov_v))
            Coverage_area = 4 * Distance_vector**2 * np.tan(self.Beam_Div_v/2 * np.pi/180 ) * np.tan(self.Beam_Div_h/2 * np.pi/180)
           # M_preScene = self.P_max / (4 * Distance_vector**2 * np.tan(self.Beam_Div_v/2 * np.pi/180 ) * np.tan(self.Beam_Div_h/2 * np.pi/180)) #Radiant Exposure
        elif self.Beam_shape == 'circle':
            Coverage_area = (np.pi * Distance_vector**2 * np.tan(self.Beam_Div_v/2 * np.pi/180 ) * np.tan(self.Beam_Div_h/2 * np.pi/180))
            #M_preScene = self.P_max / (np.pi * Distance_vector**2 * np.tan(self.Beam_Div_v/2 * np.pi/180 ) * np.tan(self.Beam_Div_h/2 * np.pi/180)) #Radiant Exposure
        else:
            logger.error("Not a valid Shape")
            return -1
        M_preScene = self.P_max / Coverage_area
        #print(Coverage_area)
        #print(M_preScene)
        return M_preScene

    def calculatePeakPowerReflected(self,reflectivity,Distance_vector):
        M_scene = reflectivity * self.calculatePeakPowerAtTarget(Distance_vector)
        return M_scene
        
class Sensor:
    def __init__(self,config):
            self.FF = config['Sensor_FF']
            self.Area = config['Pixel_A']
            self.PDP  = config['Sensor_pdp']


class Lens:
    def __init__(self,config):
        self.Diameter = config['Lens_D']
        self.Focal = config['Lens_f']
        self.Lens_transmittance = config['Lens_t']
        self.Filter_transmittance = config['filter_t']
        self.Filter_fwhm = config['filter_fwhm']
        self.Filter_fc = config['filter_fc']

    def LambertianFactor(self,Distance_vector):
        sin2alpha = self.Diameter**2 / (self.Diameter**2 + 4 * Distance_vector**2)
        return sin2alpha

    def LensOutput(self,Distance_vector):
        lens_projection = self.LambertianFactor(Distance_vector) * self.Lens_transmittance * self.Filter_transmittance
        return lens_projection

    def LensOutputFiltered(self, input_wavelength, input_spirradiance):
        dlambda = np.diff(input_wavelength) #because it's in nm
        e_photon = const.c * const.Planck / (input_wavelength*1e-9)

        filter_reponse = self.Filter_transmittance * np.exp(-(input_wavelength - self.Filter_fc*1e9)**2/(2*(self.Filter_fwhm*1e9)**2/8/np.log(2)))  #gaussian filter with a fc and fwhm
        phpspnm_dlambda = input_spirradiance
        
        e_f_dlambda =  filter_reponse[1:] * dlambda #irradiance is w/m2/nm
        print(f'Solar w/m2: {np.sum(e_f_dlambda * phpspnm_dlambda)}')
        #plt.plot(input_wavelength,filter_reponse*np.max(phpspnm_dlambda))
        #plt.plot(input_wavelength,input_spirradiance/e_photon)
        self.irradiance = np.sum( phpspnm_dlambda[1:]/e_photon * e_f_dlambda) #ph/s/m2
        return self.irradiance


class Optic_system:
    def __init__(self,configFile):
        # required parameters in the config file, should match cmdline args
        config_requires = ["Lens", "Laser", "Sensor"]
        myArgs  = confutil.parse_yaml_config(configFile, required=config_requires)
        self._Laser   =   Laser(myArgs.Laser)
        self._Sensor  =   Sensor(myArgs.Sensor)
        self._Lens    =   Lens(myArgs.Lens)
        self._solar   =   SolarModel(myArgs.Solar)
        self.Plank = const.Planck
        self.C = const.c
        

    def calculateLaserRefPower(self,distance_vector,reflectivity):
        #Calculate Power Envelope
        ref_Pmax = self._Laser.calculatePeakPowerReflected(reflectivity,distance_vector)
        lens_factor = self._Lens.LensOutput(distance_vector)
        P_pix_sig = ref_Pmax * lens_factor * self._Sensor.FF * self._Sensor.PDP * self._Sensor.Area / (self._Lens.Focal**2)
        return P_pix_sig

    '''

    '''

    def calculateEnvRefPower(self):
        P_sol = self._Lens.LensOutputFiltered(self._solar.wavelengths,self._solar.target_irradiance) * self._Lens.Lens_transmittance
        P_pix_env = P_sol * self._Sensor.FF * self._Sensor.Area * self._Sensor.PDP * self._Lens.Diameter**2 / (4* self._Lens.Focal**2) #this is an approx refer to diaries for complete 
        return P_pix_env

    def calculate_lambda_bg(self):
        if(self._solar.constant_rate > 0):
            self.lambda_bg = self._solar.constant_rate
        else:
            self.lambda_bg = self.calculateEnvRefPower() / self._solar.weather_factor
        #print(self.lambda_bg)
        print(f'BG Ph/s : {self.lambda_bg}')
        return self.lambda_bg

    def calculate_lambda_sig(self,distance_vector,reflectivity=0.75):
        E_ph = self.Plank * self.C / self._Laser.WaveLength
        self.lambda_sig_max = self.calculateLaserRefPower(distance_vector,reflectivity) / E_ph
        
        return self.lambda_sig_max
    # def calculate_lambda_sig_envelope(self,tof_vector,reflectivity = 0.75):
    #     distance_vector = tof_vector * self.C / 2
    #     lambda_max = self.calculate_lambda_sig(distance_vector,reflectivity)
    #     lambda_sig = lambda_max * np.exp()
    @staticmethod
    def create_gamma_t_envelope(lambda_sig_max, time_steps, tofs, sigma):
        signal_lambda = lambda_sig_max * \
                np.exp(- (time_steps.T - tofs)**2 / (2 * sigma**2))
        return signal_lambda

    
    def plot_sweep_lambdas_distance(self,d_init,d_end,d_step):
        d_vect = np.linspace(d_init,d_end, num=(d_end - d_init)//d_step)
        self.calculate_lambda_bg()
        self.calculate_lambda_sig(d_vect)
        print(f'Lambda bg is {self.lambda_bg}')
        plt.figure()
        plt.plot(d_vect,self.lambda_sig_max)
        horiz_line = np.array([self.lambda_bg for i in range(len(d_vect))])
        plt.semilogy(d_vect, horiz_line)
        plt.figure()
        plt.semilogy(d_vect,self.lambda_sig_max/horiz_line)

'''
    In case module is runned standalone
    Input:
        yaml configuration file including all the constants
        and path to solar spectral irradiance file
'''
#%%
if __name__ == "__main__":
    # config parser
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers(dest="subparser_name")

    parser_cfg = subparser.add_parser('config', help='Launch the simulation with a config file')
    parser_cfg.add_argument("-f",
                        "--filename",
                        dest="configFile",
                        required=True,
                        help="Path to a yaml config file for simulation")
    # parse arguments
    try:
        args = parser_cfg.parse_args()
    except:
        args = argparse.Namespace()
        args.configFile = './configs/ISCAS_Configs/config2.yaml'
    optic = Optic_system(args.configFile)
    optic.plot_sweep_lambdas_distance(1,100,1)

# %%
