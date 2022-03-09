# MonteCarlo numeric simulator for ToF SPAD arrays
This is a numeric simulator to simulate behavior of SPAD arrays used for ToF measurements. 

- [MonteCarlo numeric simulator for ToF SPAD arrays](#montecarlo-numeric-simulator-for-tof-spad-arrays)
  - [Setup](#setup)
  - [Methodology](#methodology)
  - [Algorithms](#algorithms)
  - [Features](#features)
    - [Commands](#commands)
  - [Usage](#usage)
    - [Prequisites](#prequisites)
    - [Run](#run)
      - [Config File Config](#config-file-config)
      - [Command Line Config](#command-line-config)



## Setup
This setup consists of Silicon Photon Multipliers (SiPMs) each could contain $n$ SPAD in close proximity of each other. The assumption is SPADs in each SiPM or Macro-Pixel all receive photons with the same distributions (i.e. ToF is the same for them). Each $m$ Macro-Pixel then could share $t$ Time to Digital Converter (TDC). Unlike SPADs inside a SiPM, each SiPM can have different ToF. Finally, each run of this simulator in fact simulates $M$ Macro-Pixel each having $N$ SPADs sharing $T$ TDCs. 

Here an ASCII schematic of the circuit can be seen.

    | SPAD1 | SPAD2|               |
    |:-----:|:----:|      SiPM2    |
    | SPAD3 | SPADn|               |
    |:------------:|:-------------:|
    |      ...     |     SiPM3     |

## Methodology
1. Photon rate (lambda) for all the arrays will be calculated before hand as: lambda_sum(t) = lambda_bg + lambda_sig(t)
2. Having the photon rate, events can be generated using bernoulli trials assuming this is a poisson process and approximate it as bernoulli. Each SPAD (n*m) will get independent events.
3. Events are passed to SPAD object so ones that violate SPAD dead time are eliminated.
4. Output is passed to Coincidence detector unit to check for n event in w time window.
5. Events after Coincidence detector will be passed to TDC objects. They will be distributed in a round robin fashion and missed events will be eliminated.
6. Stages 1-5 will reiterate x times, x is a function of frame rate and illumination period. 
7. A histogram will be generated for each Macro-Pixel. User can save as .mat file or plot it.

## Algorithms
Two different high-level algorithmic approaches are currently implemented. The main `montesim.py` file can call both algorithms by passing the `--alg` parameter on the command line. Passing `time-focusing` or `partial-hist` runs either the time focusing algorithm or partial histogram algorithm, respectively. Future updates will work on increasing the number of parameters that can be changed from the command line.

## Features
- MultiProcessing is used for SiPMs (each SiPM is independent and will be calculated in parallel)
- Utility/constants.py has the physical constants and user can config them.
### Commands
   Please refer to .YAML File for the list of available options
## Usage
  `python montesim.py config -f path/to/config/file.yaml -o /path/to/save/outputs`
### Prequisites 
- python3.8
- numpy
- multiprocessing
- logging 
- scipy 
- matplotlib
- tdqm
### Run

Various simulator constants can be changed by modifying `Utility/constants.py`. Outside of these constants, the simulator [commands](#commands) can be passed to the simulator in two different ways. Arguments can be passed in via the command line, or through a `.yaml` configuration file. To select one of these two options, the simulator needs to be called with either `config` or `cmdline`. 

#### Config File Config

 python montesim.py config -f <filename>
 
 Would call the simulator and specify that a config file should be used. When a config file is used, the `-f` or `--filename` parameter should be passed to the simulator, rather than the commands specified in the [Commands](#commands) section. In this case, the argument following `-f` should be a path to the `.yaml` config file containing all the same arguments as in [Commands](#commands). 

#### Command Line Config

`python3.7 montesim.py cmdline --NIter 1000 --mt 2 --sv 2`

Would call the simulator, specify that command line configuration commands should be used, and parse the relevant config options. If an argument is *not* specified, then the simulator will assume some default parameters. Be aware of these defaults when running your simulation. 
    
