#!/usr/bin/env python3

import pyIMSRG

ms = pyIMSRG.ModelSpace(1,'He6','He6')
ut = pyIMSRG.UnitTest(ms)
pyIMSRG.Commutator.SetUseIMSRG3(True)
passed = ut.TestCommutators()

print('passed? ',passed)

exit(not passed)
