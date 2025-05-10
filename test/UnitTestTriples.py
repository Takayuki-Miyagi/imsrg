#!/usr/bin/env python3

import pyIMSRG

ms = pyIMSRG.ModelSpace(3,'C14','C14')
ut = pyIMSRG.UnitTest(ms)
passed = ut.TestPerturbativeTriples()

print('passed? ',passed)
exit(not passed)
