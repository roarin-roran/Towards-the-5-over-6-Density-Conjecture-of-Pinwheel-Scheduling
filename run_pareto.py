import sys
import pareto
import time 


desired_surface_size = int(sys.argv[1])


startTime = time.time()
sampleP = pareto.pareto(desired_surface_size)
sampleP.run()

solveTimeCost = time.time() - startTime
print("time taken:", solveTimeCost)
