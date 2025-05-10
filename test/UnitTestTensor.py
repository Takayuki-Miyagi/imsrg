#!/usr/bin/env python3

import pyIMSRG

ms = pyIMSRG.ModelSpace(1,'He6','He6')
ut = pyIMSRG.UnitTest(ms)
passed = ut.TestCommutators_Tensor()

print('passed? ',passed)

exit(not passed)
