#!/usr/bin/env python

##########################################################################
##              adapted from Charlie Payne's goOak.py
##              things I should consistently change are:
##                  * wrkdir, MNU, smax, omega_norm_max
##                  file paramater list syle, ARGS['e3max']
##                  * (A,Z) in...
##                  valence_space, custom_valence_space
##                  * hw in [..., e in [..., MNU['Decay'] in [...
##                  2bme, 3bme, LECs
##              -Antoine Belley
##########################################################################

import os
import argparse
import sys
from time import time,sleep
from datetime import datetime
from interaction import *


#Paths were to find the files from the IMSRG and fro the job submission. You should change this to match your directory
# if str(os.environ['HOSTNAME']) == 'oak.arc.ubc.ca':
#   imasms = '/global/home/belley/Scripts_decay/'                           # ' ' ' ' ' ' nuqsub.py script lives
# elif str(os.environ['HOSTNAME'])[7:] == 'cedar.computecanada.ca':
#   imasms = '/home/belleya/projects/def-holt/belleya/Scripts_decay/'       # ' ' ' ' ' ' nuqsub.py script lives
# else:
#   print('This cluster is not known')
#   print('Add the paths to execute_M2nu.py in this cluster')
#   print('Exiting...')
#   exit(1)
# sys.path.append(imasms)
# from cluster_class import cedar, oak


parser = argparse.ArgumentParser()
parser.add_argument('ZI',     help = 'Atomic (proton) number of the initial nucleus (I)', type=int)
parser.add_argument('A',      help = 'Mass number of the decay', type=int)
parser.add_argument('BB',     help = '"OS", "HF", "3N", "NN"')
parser.add_argument('int',    help = 'The desired interraction to use. Can be \'BARE\', \'magic\', \'N3LO\', \'N4LO\', \'N2LOsat\', \'PWA\' or \'Delta\'')
parser.add_argument('emax',   help = 'Maximum energy to limit valence space.', type = int)
parser.add_argument('e3max',  help = 'Size of the 3-body model space', type = int)
parser.add_argument('hw',     help = 'Frequency for the harmonic oscillator basis.', type = float)
parser.add_argument('Decay',  help = 'The decay that will be used for the operators')
parser.add_argument('-t', '--time', action = 'store', help = 'Wall time for the submission on cedar')
parser.add_argument('-x', '--extra', action = 'store', help = 'Extra string to add to the name of the output directory')
parser.add_argument('-s', '--scratch', action = 'store', help = 'Name of a directory to write Omega operators so they don\'t need to be stored in memory. If not given, they\'ll just be stored in memory')
parser.add_argument('-r', '--reference', action = 'store', help = 'Name of the reference state if not default')
args = parser.parse_args()


"""USER SPECIFIC INFO'S! DON'T FORGET TO CHANGE THOSE"""
#Path to the executables, make sure to change this one to match your path
# exe = '/Users/antoinebelley/bin/imsrg++'
exe = '/Users/antoinebelley/Documents/TRIUMF/imsrg/src/imsrg++'
#Don't forget to change this. I don't want emails about your calculations...
mail_address = 'antoine.belley@mail.mcgill.ca'
#Directory that contain all the result directories
wrkdir = '/Users/antoinebelley/Documents/TRIUMF/results/'
# wrkdir = '/Users/antoinebelley/Documents/TRIUMF/Benchmark_TBME/'

#Lists to find the isotope
#Some files use with capital letters and some with lower case
#It's a bit of mess...
ELEM = ['blank', 'h', 'he',
        'li', 'be', 'b',  'c',  'n',  'o', 'f',  'ne',
        'na', 'mg', 'al', 'si', 'p',  's', 'cl', 'ar',
        'k',  'ca', 'sc', 'ti', 'v',  'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr',
        'rb', 'sr', 'y',  'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn', 'sb', 'te', 'i',  'xe',
        'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm']
ELEM2 = ['n','H','He','Li','Be','B','C','N',
       'O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K',
       'Ca','Sc','Ti','V','Cr','Mn','Fe','Co',  'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y',
       'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',  'Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
       'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb']# ,'Bi','Po','At','Rn','Fr','Ra','Ac','Th','U','Np','Pu']
nucI = f'{ELEM[args.ZI]}{args.A}'

os.chdir(wrkdir)
#MNU is a directory containing the specific info for the MNU decays
MNU = {}
#Name of the diretory containing the results
outname = f'output_{nucI}_{args.int}'
if args.extra:
  outname = f'{outname}_{args.extra}'
outname =f'{outname}/'
if not os.path.exists(outname): os.mkdir(outname)
MNU['outname'] = outname
MNU['dirout']  = f'{wrkdir}{outname}'
MNU['Reduced'] = 'R'
MNU['Ec']      = '7.72' #7.72 for Ca48 but value shouldn't affect the results that much
MNU['SRC']     = 'none' #Can be change to 'AV18', 'CD-Bonn', 'Miller-Spencer'
MNU['int']     = f'{args.int}'
MNU['BB']      = f'{args.BB}'
MNU['gap'] = 'off'
#Verifiy that we have the BARE int for the OS basis
if MNU['BB']=='OS' and MNU['int'] != 'BARE':
  print('OS needs BARE interaction. Switching to BARE...')
  MNU['int'] = 'BARE'
MNU['Decay'] = f'{args.Decay}'

ARGS = {}

ARGS['input_op_fmt'] = 'miyagi'
ARGS['A'] = f'{args.A}'
if not args.reference:
  ARGS['reference'] = f'{ELEM2[args.ZI]}{args.A}'
else:
  ARGS['reference'] = f'{args.reference}'
# Set the space (sp)
if args.A == 3 and args.ZI == 1:
 ARGS['valence_space'] = 's-shell'
 MNU['Ec'] = 2.74
if args.A == 6 and args.ZI == 2:
 ARGS['valence_space'] = 'p-shell'
 MNU['Ec'] = 2.74
elif args.A == 8 and args.ZI == 2:
  ARGS['valence_space'] = 'p-shell' # AKA: p
  MNU['Ec'] = 3.17
if args.A == 9 and args.ZI == 4:
 ARGS['valence_space'] = 'p-shell'
 MNU['Ec'] = 2.74
if args.A == 14 and args.ZI == 7:
 ARGS['valence_space'] = 'p-shell'
 MNU['Ec'] = 2.74
elif args.A == 10 and (args.ZI == 2 or args.ZI == 4 or args.ZI == 6):
  ARGS['valence_space'] = 'p-shell' # AKA: p
  MNU['Ec'] = 3.54
elif args.A == 14 and args.ZI == 6:
  ARGS['valence_space'] = 'psd5-shell' # AKA: p
  ARGS['custom_valence_space'] = "He4,p0p3,n0p3,p0p1,n0p1,p1s1,n1s1,p0d5,n0d5"
  MNU['Ec'] = 4.19
elif args.A == 14 and args.ZI == 8:
  ARGS['valence_space'] = 'p-shell' # AKA: p
  MNU['Ec'] = 4.19
elif args.A == 22 and args.ZI == 8:
  ARGS['valence_space'] = 'sd-shell' # AKA: sd
  MNU['Ec'] = 5.25
elif args.A<40 and args.ZI<20:
  ARGS['valence_space'] = 'sd-shell' # AKA: sd
  MNU['Ec'] = 1.12*args.A**(1/2)
elif args.A == 40 and args.ZI == 20:
  # ARGS['valence_space'] = 'fp-shell' # AKA: fp
  ARGS['valence_space'] = 'sd-shell'  # AKA: fp
  #ARGS['valence_space'] = 'Ca40' # for Jiangming
elif args.A >= 42 and (args.ZI == 20 or args.ZI == 22 or args.ZI == 24):
  ARGS['valence_space'] = 'fp-shell' # AKA: fp
  MNU['Ec'] = 7.72
  # MNU['Ec'] = 0.143
  #ARGS['valence_space'] = 'Ca48' # for Jiangming
elif args.A == 76 and args.ZI == 32:
  ARGS['valence_space'] = 'pf5g9' # this is just a label when custom_valence_space is set
  ARGS['custom_valence_space'] = 'Ni56,p0f5,n0f5,p1p3,n1p3,p1p1,n1p1,p0g9,n0g9' # AKA: Ni56 core with jj44pn
  MNU['Ec'] = 9.41
elif args.A == 82 and args.ZI == 34:
  ARGS['valence_space'] = 'pf5g9' # this is just a label when custom_valence_space is set
  ARGS['custom_valence_space'] = 'Ni56,p0f5,n0f5,p1p3,n1p3,p1p1,n1p1,p0g9,n0g9' # AKA: Ni56 core with jj44pn
  MNU['Ec'] = 9.41
elif args.A == 130 and args.ZI == 52:
  ARGS['valence_space'] = 'sdg'
  ARGS['custom_valence_space'] = 'Sn100,p0g7,n0g7,p1d3,n1d3,p1d5,n1d5,p2s1,n2s1,p0h11,n0h11'
  MNU['Ec'] = 12.76
elif args.A == 136 and args.ZI == 54:
  ARGS['valence_space'] = 'sdg'
  ARGS['custom_valence_space'] = 'Sn100,p0g7,n0g7,p1d3,n1d3,p1d5,n1d5,p2s1,n2s1,p0h11,n0h11'
  MNU['Ec'] = 13.06
else:
  ARGS['valence_space'] = '0hw-shell'



#Print verifictation to screen
print('outname   = ',MNU['outname'])
print('dirout    = ',MNU['dirout'])
print('Reduced   = ',MNU['Reduced'])
print('Ec        = ',MNU['Ec'])
print('SRC       = ',MNU['SRC'])
print('BB        = ',MNU['BB'])
print('int       = ',MNU['int'])
print('gap?      = ',MNU['gap'])
print('Decay     = ',MNU['Decay'])


#ARGS is a directory containing the values to pass as argument for the imsrg exe

if MNU['int'] == 'BARE' or MNU['BB'] == 'HF':
  ARGS['smax'] = '0'
else:
  ARGS['smax'] = '500'
ARGS['dsmax'] = '0.3'
# ARGS['eta_criterion'] = 1e-6
### Norm of Omega at which we split off and start a new transformation
ARGS['omega_norm_max'] = '0.5'
### Solution method
ARGS['method'] = 'magnus'
if MNU['BB'] == '3N':
  ARGS['write_HO_ops'] = 'False'
  ARGS['write_HF_ops'] = 'False'
#Name of a directory to write Omega operators so they don't need to be stored in memory. If not given, they'll just be stored in memory.
#if args.emax == 14:
#  ARGS['scratch'] = '/home/belleya/scratch/'
#elif args.scratch:
#  ARGS['scratch'] = f'{args.scratch}'
if args.e3max>=26 or args.int == 'Sample34': 
  ARGS["no2b_precision"] = "half"

run_multiple = False
#Set the parameters for the interraction and decay
if MNU['int'] == 'BARE':
  print('running BARE-style...')
  ARGS['basis'] = 'oscillator'
  BARE.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
  if MNU['BB'] != 'OS':
    print('for MNU[int] = BARE, choose MNU[BB] = OS')
    print('exiting...')
    exit()
  print('running with OS basis (no IMSRG evolution)...')
elif MNU['int'] == 'mag16':
  print('running magic interraction...')
  mag16.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
elif MNU['int'] == 'magic':
  print('running magic interraction...')
  if args.hw ==12:
    if args.e3max>=26:
      magic12high.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
    else:
      magic12.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
  elif args.hw ==16:
    if args.e3max>=26:
      magic16high.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
    else:
      magic16.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
elif MNU['int'] == 'magicwocoul':
  magic16woCoulomb.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
elif MNU['int'] == 'magiciso':
  magic16iso.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
elif MNU['int'] == 'EM2':
  print('running EM2.0/2.0 interraction...')
  EM2.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
elif MNU['int'] == 'N2LOsat':
  print('running N2LO Sat int...')
  N2LO_SAT.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
elif MNU['int'] == 'LO_EMN':
  print('running LO...')
  if args.hw == 15:
    LO_EM500.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
elif MNU['int'] == 'N2LO_EMN':
  print('running N2LO...')
  if args.hw == 15:
    N2LO_EM500.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
elif MNU['int'] == 'NLO_EMN':
  print('running N2LO...')
  if args.hw == 15:
    NLO_EM500.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
elif MNU['int'] == 'N3LO_EMN':
  print('running N3LO...')
  if args.hw == 15:
    N3LO_EMN500.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
elif MNU['int'] == 'N3LO':
  print('running N3LO...')
  if args.hw == 16:
      N3LO_EM500_LNL.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
  elif args.hw == 20:
      N3LO_EM500_LNL20.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
  elif args.hw == 15:
    N3LO_EM500_LNL15.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
elif MNU['int'] == 'N4LO':
  print('running N4LO...')
  N4LO_EM500_LNL.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
elif MNU['int'] == 'PWA':
  print('running PWA int...')
  PWA.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
elif MNU['int'] == 'Delta':
  print('running PWA int...')
  if args.e3max >= 26:
    DELTA_NNLO_G0high.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
  else:
    DELTA_NNLO_GO.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'])
elif MNU['int'] == 'Sample34':
  print('running 34 sample int...')
  run_multiple = True
else:
  print('this *int* has not been set up yet!')
  print('exiting...')
  exit()

ARGS['BetaCM'] = '5'
ARGS['hwBetaCM'] = f'{args.hw}'
ARGS['input_op_fmt'] = 'miyagi'
ARGS['hw'] = f'{args.hw}'
ARGS['emax'] = f'{args.emax}'
ARGS['e3max'] = f'{args.e3max}'
ARGS['denominator_delta'] =  10
#ARGS["3bme_type"] = "no2b"
ARGS['freeze_occupations'] =  "false"

# #Print the parameters choosen
# print('smax      = ',ARGS['smax'])
# print('o.norm.m. = ',ARGS['omega_norm_max'])
# print('method    = ',ARGS['method'])
# print('e3max     = ',ARGS['e3max'])
# print('Z ========= ',args.ZI)
# print ('A ========= ',args.A)
# print('ref       = ',ARGS['reference'])
# print('v.sp.     = ',ARGS['valence_space'])
# print('=======================')
# print(f'hw        = {args.hw}')
# print('=======================')
# print('-----------------------')
# print(f'e         = {args.emax}')
# print('-----------------------')
# print('LECs      = ',ARGS['LECs'])


# Make an estimate of how much time to request and the memory needed. Only used for slurm at the moment.
if   args.emax <  5 : 
  time_request = '00:10:00'
  vmem         = 125
elif args.emax <=  8 :
  time_request = '03:00:00'
  vmem         = 187
elif args.emax <= 10 : 
  time_request = '06:00:00'
  vmem         = 187
elif args.emax <= 12 :
  time_request = '24:00:00'
  vmem         = 187
else:
  time_request = '72:00:00'
  vmem         = 250
if args.time:
  time_request = args.time

if run_multiple==True:
  #Submit jobs for the 34 sample interactions.
  for sample in [2007, 2460, 3448, 3373, 1141, 3319, 3098, 1429, 3895, 1172, 90, 723, 2125, 2245, 750, 1469, 1177, 493, 500, 3621, 606, 3480, 3260, 3813, 3105, 2411, 3350, 774, 2618, 4117, 922, 1173, 1802, 3472]:
  #for sample in [1141,3813,4117]:
    Sample34.set_args(ARGS, MNU['BB'], args.Decay, args.hw, MNU['SRC'], MNU['Ec'],sample)
    print('smax      = ',ARGS['smax'])
    print('o.norm.m. = ',ARGS['omega_norm_max'])
    print('method    = ',ARGS['method'])
    print('e3max     = ',ARGS['e3max'])
    print('Z ========= ',args.ZI)
    print ('A ========= ',args.A)
    print('ref       = ',ARGS['reference'])
    print('v.sp.     = ',ARGS['valence_space'])
    print('=======================')
    print(f'hw        = {args.hw}')
    print('=======================')
    print('-----------------------')
    print(f'e         = {args.emax}')
    print('-----------------------')
    # print('LECs      = ',ARGS['LECs'])
    # print('-----------------------')
    print('Sample    = ',sample)
        #Jobname and logname for the job submission
    jobname  = f'{ARGS["valence_space"]}_{MNU["int"]}_{sample}_{MNU["BB"]}_{ARGS["LECs"]}_{ARGS["method"]}_{ARGS["reference"]}_e{ARGS["emax"]}_E{ARGS["e3max"]}_s{ARGS["smax"]}_hw{ARGS["hw"]}_A{ARGS["A"]}_{MNU["Decay"]}_'
    ### Some optional parameters that we probably want in the output name if we're using them
    if 'lmax3' in ARGS:  jobname  += 'l%d'%(ARGS['lmax3']) + '_'
    if 'eta_criterion' in ARGS: jobname += 'eta%s'%(ARGS['eta_criterion']) + '_'
    if 'core_generator' in ARGS: jobname += '' + ARGS['core_generator'] + '_'
    if 'gapE' in ARGS: jobname += 'gap' + ARGS['gapE'] + '_'
    if 'BetaCM' in ARGS: jobname += 'BCM' + ARGS['BetaCM'] + '_'
    ARGS['flowfile'] = MNU['outname'] + 'BCH_' + jobname + 'flow_' + datetime.fromtimestamp(time()).strftime('%y%m%d%H%M.dat')
    ARGS['intfile']  = MNU['outname'] + jobname
    logname = f'imsrg_log/{jobname}'
    cmd = ' '.join([exe] + ['%s=%s'%(x,ARGS[x]) for x in ARGS])

    ### Make a directory for the log files, if it doesn't already exist
    if not os.path.exists('imsrg_log'): os.mkdir('imsrg_log')

    #Submit the job if we're running in batch mode, otherwise just run in the current shell
    if str(os.environ['HOSTNAME']) == oak.host:
      oak.submit_job(cmd,jobname,mail_address)
    elif str(os.environ['HOSTNAME'])[7:] == cedar.host: 
      cedar.submit_job(cmd,jobname, logname, mail_address, vmem = vmem, time = time_request)
    else:
      print('Cluster is not known. Please create an instance for it.')
      print('Exiting...')
      exit(1)
    sleep(0.1)

else:
  print('smax      = ',ARGS['smax'])
  print('o.norm.m. = ',ARGS['omega_norm_max'])
  print('method    = ',ARGS['method'])
  print('e3max     = ',ARGS['e3max'])
  print('Z ========= ',args.ZI)
  print ('A ========= ',args.A)
  print('ref       = ',ARGS['reference'])
  print('v.sp.     = ',ARGS['valence_space'])
  print('=======================')
  print(f'hw        = {args.hw}')
  print('=======================')
  print('-----------------------')
  print(f'e         = {args.emax}')
  print('-----------------------')
  # print('LECs      = ',ARGS['LECs'])
  #Jobname and logname for the job submission
  if len(MNU["Decay"])>100: MNU["Decay"] = "r_depedance_GT_DGT"
  jobname  = f'{ARGS["valence_space"]}_{MNU["int"]}_{MNU["BB"]}_{ARGS["method"]}_{ARGS["reference"]}_e{ARGS["emax"]}_E{ARGS["e3max"]}_s{ARGS["smax"]}_hw{ARGS["hw"]}_A{ARGS["A"]}_{MNU["Decay"]}_'
  if args.extra:
    jobname  += f'{args.extra}_'
  ### Some optional parameters that we probably want in the output name if we're using them
  if 'lmax3' in ARGS:  jobname  += 'l%d'%(ARGS['lmax3']) + '_'
  if 'eta_criterion' in ARGS: jobname += 'eta%s'%(ARGS['eta_criterion']) + '_'
  if 'core_generator' in ARGS: jobname += '' + ARGS['core_generator'] + '_'
  if 'gapE' in ARGS: jobname += 'gap' + ARGS['gapE'] + '_'
  if 'BetaCM' in ARGS: jobname += 'BCM' + ARGS['BetaCM'] + '_'
  ARGS['flowfile'] = MNU['outname'] + 'BCH_' + jobname + 'flow_' + datetime.fromtimestamp(time()).strftime('%y%m%d%H%M.dat')
  ARGS['intfile']  = MNU['outname'] + jobname
  logname = f'imsrg_log/{jobname}'
  cmd = ' '.join([exe] + ['%s=%s'%(x,ARGS[x]) for x in ARGS])
  print(cmd)

  ### Make a directory for the log files, if it doesn't already exist
  if not os.path.exists('imsrg_log'): os.mkdir('imsrg_log')

  #Submit the job if we're running in batch mode, otherwise just run in the current shell
  # if str(os.environ['HOSTNAME']) == oak.host:
  #   oak.submit_job(cmd,jobname,mail_address)
  # elif str(os.environ['HOSTNAME'])[7:] == cedar.host: 
  #   cedar.submit_job(cmd,jobname, logname, mail_address, vmem = vmem, time = time_request)
  # else:
  #   print('Cluster is not known. Please create an instance for it.')
  #   print('Exiting...')
  #   exit(1)
  # sleep(0.1)

  os.system(cmd)

