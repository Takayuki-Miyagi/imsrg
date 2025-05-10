import os
import numpy as np

ELEM = ['n','H','He','Li','Be','B','C','N',
       'O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K',
       'Ca','Sc','Ti','V','Cr','Mn','Fe','Co',  'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y',
       'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',  'Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
       'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb']


# ZA = [[52,130],[54,136]]
# ZA = [[2,6],[2,8],[4,10],[6,14],[8,22],[20,48]]
# ZA=[[4,10],[6,14]]
# ZA = [[24,50]]
ZA = [[8,22]]
# ZA = [[8,22],[20,48]]
# ZA = [[2,6]]
BB  = ['3N']
interactions = ['magic']
emax = [4]
e3max = [16]
hw = [16]
rs = np.concatenate([np.linspace(0,1,20),np.linspace(1,7,20)])
for az in ZA:
  op = ""
  EC = 1.12*az[1]**(1/2)
  for r in rs:
      op += f'DGTR_{r:.2f},M0nuR_{EC:.2f}_{r:.2f},DGTRLocal_{r:.2f},'#
  for bb in BB:
      for inter in interactions:
          for e in emax:
              for e3 in e3max:
                  for freq in hw:
                      cmd =  f"python3 goImsrgMnu.py {az[0]} {az[1]} {bb} {inter} {e} {e3} {freq} {op}"
                      # print(cmd)
                      os.system(cmd)