#!/usr/bin/env python3

import pyIMSRG

ms = pyIMSRG.ModelSpace(2,'He6','He6')
ms.PreCalculateSixJ()
ut = pyIMSRG.UnitTest(ms)

print('scalar 3b')
Rando = ut.RandomOp(ms,0,0,0,3,+1)
passed = ut.TestNormalOrdering(Rando)

print('lambda=1 2b')
Rando = ut.RandomOp(ms,1,0,0,2,+1)
passed &= ut.TestNormalOrdering(Rando)

print('parity-changing 2b')
Rando = ut.RandomOp(ms,0,0,1,2,+1)
passed &= ut.TestNormalOrdering(Rando)

print('charge-changing 2b')
Rando = ut.RandomOp(ms,0,1,0,2,+1)
passed &= ut.TestNormalOrdering(Rando)

print('passed? ',passed)

exit(not passed)
