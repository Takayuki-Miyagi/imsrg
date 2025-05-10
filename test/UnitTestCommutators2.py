#!/usr/bin/env python3

import pyIMSRG

ms = pyIMSRG.ModelSpace(2,'He6','He6')
ut = pyIMSRG.UnitTest(ms)
pyIMSRG.Commutator.SetUseIMSRG3(True)
pyIMSRG.Commutator.SetUseIMSRG3N7(True)
passed = ut.TestCommutators()

print('passed? ',passed)

exit(not passed)
