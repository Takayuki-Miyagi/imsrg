 ###File contianing the definition of the intercation class in the instances
###that have been implemented so far. 
###Author: Antoine Belley


class Interaction():
  """Class that describes what files should be used for a particular interaction"""
  def __init__(self, file2e1max, file3e1max, file2bme, file3bme, LECs, Threebme_type, file2c = None, fileIf = None,  reg_cutoff=500, reg_power=3):
    self.file2e1max = file2e1max
    self.file3e1max = file3e1max
    self.file2bme   = file2bme
    self.file3bme   = file3bme
    self.LECs       = LECs
    self.Threebme_type  = Threebme_type
    self.file2c     = file2c     #File for the 2bc
    self.fileIf     = fileIf     #File for the induced force
    self.reg_power = reg_power
    self.reg_cutoff = reg_cutoff

  def set_args(self, Dict, BB, Decay, hw, SRC , Ec, sample=None):
    """Set the required argument for the imsrg for the specific interaction and decay"""
    Dict['file2e1max'] = self.file2e1max
    Dict['file3e1max'] = self.file3e1max
    if BB != 'OS':
      if sample == None:
        Dict['2bme']       = self.file2bme%(hw)
      else:
        Dict['2bme']       = self.file2bme%(sample,hw)
    else:
      Dict['2bme'] = 'none'
    Dict['LECs']  = self.LECs
    if BB == 'OS':
      Dict['3bme'] = 'none'
    if BB == '2N':
      Dict['3bme'] = 'none'
    if BB == '3N' or BB == 'HF':
      if sample == None:
        Dict['3bme'] = self.file3bme%(hw)
      else:
        Dict['3bme'] = self.file3bme%(sample,hw)
      Dict['3bme_type'] = self.Threebme_type
    if Decay == '2':
      Dict['Operators'] = 'GamowTeller'
    elif Decay == '2c':
      if self.file2c != None:
        Dict['OperatorsFromFile'] = self.file2c
      else:
        print('This interaction has no file for the 2bc')
        print('Exiting...')
        exit(1)
    elif Decay == '2IF':
      if self.fileIf != None:
        Dict['OperatorsFromFile'] = self.fileIf
      else:
        print('This interaction has no file for the induced force')
        print('Exiting...')
        exit(1)
    elif Decay == 'GT' or Decay == 'F' or Decay == 'T' or Decay == 'C':
      Dict['Operators'] = f'M0nu_{Decay}_{Ec}_{SRC}'
    elif Decay == 'M0nu':
      Dict['Operators'] = f'M0nu_GT_{Ec}_{SRC},M0nu_F_{Ec}_{SRC},M0nu_T_{Ec}_{SRC},M0nu_C_{self.reg_cutoff}_{self.reg_power}'
    elif Decay == 'M0nuHeavy':
      Dict['Operators'] =  f'M0nuHeavy_GT_{SRC}_AA,M0nuHeavy_GT_{SRC}_AP,M0nuHeavy_GT_{SRC}_PP,M0nuHeavy_F_{SRC}_VV,M0nuHeavy_T_{SRC}_AP,M0nuHeavy_T_{SRC}_PP'
    elif Decay == 'Test_op_from_file':
      # Dict['OperatorsFromFile'] = "0vbbFVV_SR^0_2_0_2^/Users/antoinebelley/Documents/TRIUMF/NuHamil-public/0vbbFermiVV_Ec_7.72_Range_SR-N2LO_TwBME-HO_NN-only_N3LO_EM500_bare_hw16_emax4_e2max8.me2j.gz"
      Dict['OperatorsFromFile'] = "0vbbN2LO^0_2_0_2^/Users/antoinebelley/Documents/TRIUMF/NuHamil-public/0vbbGamowTellerPP_Ec_7.72_Range_SR_-N2LO-Local3-500_TwBME-HO_NN-only_N3LO_EM500_srg2.0_hw15_emax14_e2max28.me2j.gz"
    elif Decay == "None":
      pass
    else:
      Dict['Operators'] = Decay
    
###########################################################################################
# me2j_dir = "/home/belleya/projects/def-holt/shared/me2j"
# me3j_dir = '/home/belleya/projects/def-holt/shared/me3j'
me2j_dir = "/Users/antoinebelley/Documents/TRIUMF/Interactions"
me3j_dir = '/Users/antoinebelley/Documents/TRIUMF/Interactions'
#Instances of the interraction implemented so far
BARE    = Interaction(file2e1max = '0 file2e2max=0 file2lmax=0', file3e1max = '0 file3e2max=0 file3e3max=0', 
                      file2bme = 'none', file3bme = 'none', LECs = 'none', Threebme_type='none')

mag16   = Interaction(file2e1max = '14 file2e2max=28 file2lmax=14', file3e1max = '14 file3e2max=28 file3e3max=16',
                      file2bme = '/home/belleya/projects/def-holt/shared/exch/ME_share/vnn_hw%d.00_kvnn10_lambda1.80_mesh_kmax_7.0_100_pc_R15.00_N15.dat_to_me2j.gz', 
                      file3bme = '/home/belleya/projects/def-holt/shared/exch/ME_share/jsTNF_Nmax_16_J12max_8_hbarOmega_%d.00_Fit_cutoff_2.00_nexp_4_c1_1.00_c3_1.00_c4_1.00_cD_1.00_cE_1.00_2pi_0.00_2pi1pi_0.00_2picont_0.00_rings_0.00_J3max_9_new_E3_16_e_14_ant_EM1.8_2.0.h5_to_me3j.gz',
                      LECs     = 'EM1.8_2.0', Threebme_type='full')

magic12   = Interaction(file2e1max = '16 file2e2max=32', file3e1max = '16 file3e2max=32 file3e3max=24',
                     file2bme = '/home/belleya/projects/def-holt/shared/me2j/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw%d_emax16_e2max32.me2j.gz', 
                     file3bme = '/home/belleya/projects/def-holt/shared/me3j/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw%d_ms16_32_24.stream.bin', 
                     LECs     = 'EM1.8_2.0', Threebme_type = 'no2b') #For hw=12

magic12high   = Interaction(file2e1max = '16 file2e2max=32', file3e1max = '16 file3e2max=32 file3e3max=28',
                     file2bme = '/home/belleya/projects/def-holt/shared/me2j/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw%d_emax16_e2max32.me2j.gz', 
                     file3bme = '/home/belleya/projects/def-holt/shared/me3j/NO2B_half_ThBME_EM1.8_2.0_3NFJmax15_IS_hw%d_ms16_32_28.stream.bin', 
                     LECs     = 'EM1.8_2.0', Threebme_type = 'no2b') #For hw=12, e3max>= 26

magic16   = Interaction(file2e1max = '18 file2e2max=36 file2lmax=18', file3e1max = '18 file3e2max=36 file3e3max=24',
                      file2bme = f'{me2j_dir}/TwBME-HO_NN-only_N3LO_EM500_srg1.80_hw%d_emax18_e2max36.me2j.gz', 
                      file3bme = f'{me3j_dir}/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw%d_ms18_36_24.stream.bin', 
                      LECs     = 'EM1.8_2.0', Threebme_type = 'no2b') #for hw=16

magic16woCoulomb   = Interaction(file2e1max = '14 file2e2max=28 file2lmax=14', file3e1max = '18 file3e2max=32 file3e3max=24',
                      file2bme = f'{me2j_dir}/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw%d_emax14_e2max28_wocoul.me2j.gz', 
                      file3bme = f'{me3j_dir}/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw%d_ms18_36_24.stream.bin', 
                      LECs     = 'EM1.8_2.0', Threebme_type = 'no2b') #for hw=16

magic16iso   = Interaction(file2e1max = '14 file2e2max=28 file2lmax=14', file3e1max = '18 file3e2max=32 file3e3max=24',
                      file2bme = f'{me2j_dir}/TwBME-HO_NN-only_N3LO_EM500_iso_srg1.8_hw%d_emax14_e2max28_wocoul.me2j.gz', 
                      file3bme = f'{me3j_dir}/NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw%d_ms18_36_24.stream.bin', 
                      LECs     = 'EM1.8_2.0', Threebme_type = 'no2b') #for hw=16

magic16high   = Interaction(file2e1max = '18 file2e2max=36 file2lmax=18', file3e1max = '16 file3e2max=32 file3e3max=28',
                      file2bme = '/home/belleya/projects/def-holt/shared/me2j/TwBME-HO_NN-only_N3LO_EM500_srg1.80_hw%d_emax18_e2max36.me2j.gz', 
                      file3bme = '/home/belleya/projects/def-holt/shared/me3j/NO2B_half_ThBME_EM1.8_2.0_3NFJmax15_IS_hw%d_ms16_32_28.stream.bin', 
                      LECs     = 'EM1.8_2.0', Threebme_type = 'no2b') #for hw=16

EM2   = Interaction(file2e1max = '16 file2e2max=32', file3e1max = '16 file3e2max=32 file3e3max=24',
                      file2bme = '/project/6006601/shared/me2j/TwBME-HO_NN-only_N3LO_EM500_srg2.0_hw%d_emax16_e2max32.me2j.gz', 
                      file3bme = '/project/6006601/shared/me3j/NO2B_ThBME_EM2.0_2.0_3NFJmax15_IS_hw%d_ms16_32_24.stream.bin', 
                      LECs     = 'EM2.0_2.0', Threebme_type = 'no2b')

N2LO_SAT = Interaction(file2e1max = '16 file2e2max=32 file2lmax=16', file3e1max= '18 file3e2max=36 file3e3max=24',
                       file2bme = '/home/belleya/projects/def-holt/shared/me2j/TwBME-HO_NN-only_N2LO_sat_bare_hw%d_emax16_e2max32.me2j.gz',
                       file3bme = '/home/belleya/projects/def-holt/shared/me3j/NO2B_ThBME_N2LOsat_3NFJmax15_IS_hw%d_ms18_36_24.stream.bin',
                       LECs = 'N2LOsat', Threebme_type = 'no2b')

PWA      = Interaction(file2e1max = '16 file2e2max=32', file3e1max= '18 file3e2max=36 file3e3max=24',
                       file2bme = '/home/belleya/projects/def-holt/shared/me2j/TwBME-HO_NN-only_N3LO_EM500_srg2.0_hw%d_emax16_e2max32.me2j.gz',
                       file3bme = '/home/belleya/projects/def-holt/shared/me3j/NO2B_ThBME_N2LOsat_3NFJmax15_IS_hw%d_ms18_36_24.stream.bin',
                       LECs = 'PWA', Threebme_type = 'no2b')

LO_EM500 = Interaction(file2e1max = '16 file2e2max=32 file2lmax=16', file3e1max= '16 file3e2max=32 file3e3max=22',
                             file2bme = '/Users/antoinebelley/Documents/TRIUMF/Interactions/TwBME-HO_NN-only_LO_EMN500_srg2.0_hw%d_emax16_e2max32.me2j.gz',
                             file3bme = '/Users/antoinebelley/Documents/TRIUMF/Interactions/NO2B_ThBME_srg2.0_ramp46-9-44-15-42_LO_EMN500_IS_hw%dfrom30_ms16_32_22.stream.bin',
                             LECs = 'LO_EM500', Threebme_type = 'no2b')

NLO_EM500 = Interaction(file2e1max = '16 file2e2max=32 file2lmax=16', file3e1max= '16 file3e2max=32 file3e3max=22',
                             file2bme = '/Users/antoinebelley/Documents/TRIUMF/Interactions/TwBME-HO_NN-only_NLO_EMN500_srg2.0_hw%d_emax16_e2max32.me2j.gz',
                             file3bme = '/Users/antoinebelley/Documents/TRIUMF/Interactions/NO2B_ThBME_srg2.0_ramp46-9-44-15-42_NLO_EMN500_IS_hw%dfrom30_ms16_32_22.stream.bin',
                             LECs = 'NLO_EM500', Threebme_type = 'no2b')

N2LO_EM500 = Interaction(file2e1max = '16 file2e2max=32 file2lmax=16', file3e1max= '16 file3e2max=32 file3e3max=22',
                             file2bme = '/Users/antoinebelley/Documents/TRIUMF/Interactions/TwBME-HO_NN-only_N2LO_EMN500_srg2.0_hw%d_emax16_e2max32.me2j.gz',
                             file3bme = '/Users/antoinebelley/Documents/TRIUMF/Interactions/NO2B_ThBME_srg2.0_ramp46-9-44-15-42_N2LO_EMN500_3NFJmax15_c1_-0.74_c3_-3.61_c4_2.44_cD_-0.36_cE_-0.189_LNL2_650_500_IS_hw%dfrom30_ms16_32_22.stream.bin',
                             LECs = 'N2LO_EM500', Threebme_type = 'no2b')

N3LO_EMN500 = Interaction(file2e1max = '16 file2e2max=32 file2lmax=16', file3e1max= '16 file3e2max=32 file3e3max=22',
                             file2bme = '/Users/antoinebelley/Documents/TRIUMF/Interactions/TwBME-HO_NN-only_N3LO_EMN500_srg2.0_hw%d_emax16_e2max32.me2j.gz',
                             file3bme = '/Users/antoinebelley/Documents/TRIUMF/Interactions/NO2B_ThBME_srg2.0_ramp46-9-44-15-42_N3LO_EMN500_3NFJmax15_c1_-1.2_c3_-4.43_c4_2.67_cD_-1.33_cE_-0.413_LNL2_650_500_IS_hw%dfrom30_ms16_32_22.stream.bin',
                             LECs = 'N3LO_EM500', Threebme_type = 'no2b')

N3LO_EM500_LNL = Interaction(file2e1max = '16 file2e2max=32 file2lmax=16', file3e1max= '16 file3e2max=32 file3e3max=22',
                             file2bme = '/home/belleya/projects/def-holt/shared/me2j/TwBME-HO_NN-only_N3LO_EM500_srg2.0_hw%d_emax16_e2max32.me2j.gz',
                             file3bme = '/home/belleya/projects/def-holt/shared/me3j/NO2B_ThBME_srg2.0_ramp46-9-44-15-42_N3LO_EM500_3NFJmax15_c1_-0.81_c3_-3.2_c4_5.4_cD_0.7_cE_-0.06_LNL2_650_500_IS_hw%dfrom30_ms16_32_22.stream.bin',
                             LECs = 'N3LO_EM500_LNL', Threebme_type = 'no2b')

N3LO_EM500_LNL15 = Interaction(file2e1max = '18 file2e2max=36 file2lmax=18', file3e1max= '16 file3e2max=32 file3e3max=28',
                             file2bme = '/home/belleya/projects/def-holt/shared/me2j/TwBME-HO_NN-only_N3LO_EM500_srg2.00_hw%d_emax18_e2max36.me2j.gz',
                             file3bme = '/home/belleya/projects/def-holt/shared/me3j/NO2B_half_ThBME_srg2.0_ramp48_N3LO_EM500_JJmax13_c1_-0.81_c3_-3.2_c4_5.4_cD_0.7_cE_-0.06_LNL2_650_500_IS_hw%dfrom30_ms18_36_28.stream.bin',
                             LECs = 'N3LO_EM500_LNL', Threebme_type = 'half')


N3LO_EM500_LNL20 = Interaction(file2e1max = '14 file2e2max=28 file2lmax=14', file3e1max= '14 file3e2max=28 file3e3max=14',
                             file2bme = '/project/6006601/shared/me2j/TwBME-HO_NN-only_N3LO_EM500_srg2.0_hw%d_emax14_e2max28.me2j.gz',
                             file3bme = '/project/6006601/shared/me3j/ThBME_srg2.0_ramp40-5-36-7-32-9-28-11-24_N3LO_EM500_c1_-0.81_c3_-3.2_c4_5.4_cD_0.7_cE_-0.06_LNL2_650_500_IS_hw%d_ms14_28_14.me3j.gz',
                             LECs = 'N3LO_EM500_LNL', Threebme_type = 'full')

N4LO_EM500_LNL = Interaction(file2e1max = '16 file2e2max=32 file2lmax=16', file3e1max= '16 file3e2max=32 file3e3max=22',
                             file2bme = '/home/belleya/projects/def-holt/shared/me2j/TwBME-HO_NN-only_N4LO_EMN500_srg2.0_hw%d_emax16_e2max32.me2j.gz',
                             file3bme = '/home/belleya/projects/def-holt/shared/me3j/NO2B_ThBME_srg2.0_ramp46-9-44-15-42_N4LO_EMN500_3NFJmax15_c1_-0.73_c3_-3.38_c4_1.69_cD_-1.8_cE_-0.31_LNL2_650_500_IS_hw%dfrom30_ms16_32_22.stream.bin',
                             LECs = 'N4LO_EM500_LNL', Threebme_type = 'no2b')

DELTA_NNLO_GO = Interaction(file2e1max = '18 file2e2max=36, file2lmax=18', file3e1max= '16 file3e2max=32 file3e3max=24',
                             file2bme = '/home/belleya/projects/def-holt/shared/me2j/TwBME-HO_NN-only_DNNLOgo_bare_hw%d_emax18_e2max36.me2j.gz',
                             file3bme = '/home/belleya/projects/def-holt/shared/me3j/NO2B_ThBME_DNNLOgo_3NFJmax15_IS_hw%d_ms16_32_24.stream.bin',
                             LECs = 'N2LOgo', Threebme_type = 'no2b', reg_power=4, reg_cutoff=394)


DELTA_NNLO_G0high = Interaction(file2e1max = '18 file2e2max=36, file2lmax=18', file3e1max= '16 file3e2max=32 file3e3max=28',
                             file2bme = '/home/belleya/projects/def-holt/shared/me2j/TwBME-HO_NN-only_DNNLOgo_bare_hw%d_emax18_e2max36.me2j.gz',
                             file3bme = '/home/belleya/projects/def-holt/shared/me3j/NO2B_half_ThBME_DNNLOgo_3NFJmax15_IS_hw%d_ms16_32_28.stream.bin',
                             LECs = 'N2LOgo', Threebme_type = 'no2b', reg_power=4, reg_cutoff=394)

Sample34       = Interaction(file2e1max = '16 file2e2max=32, file2lmax=16', file3e1max= '16 file3e2max=32 file3e3max=28',
                             file2bme = '/home/belleya/projects/def-holt/shared/me2j/TwBME-HO_NN-only_DeltaGO394_sample%d_bare_hw%d_emax16_e2max32.me2j.gz',
                             file3bme = '/home/belleya/projects/def-holt/shared/me3j/delta_samples/NO2B_half_ThBME_3NFJmax15_DN2LOGO394_sample%d_IS_hw%d_ms16_32_28.stream.bin',
                             LECs = 'N2LOgo', Threebme_type = 'no2b', reg_power=4, reg_cutoff=394)



