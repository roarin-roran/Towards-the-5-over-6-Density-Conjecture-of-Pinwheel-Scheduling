import sys
import solver_foresight


inputGardenStrings = sys.argv[1].split(',')
inputGardenInts = []
for i in inputGardenStrings:
    inputGardenInts.append(int(i))


oursolver_foresight = solver_foresight.solver_foresight(inputGardenInts, False, True, False)
oursolver_foresight.solve()