# graph
import solver_graph
# naive
import solver_naive
# opt
import solver_opt
# foresight
import solver_foresight

from fractions import Fraction
import random
import math
import sys

class pinwheel_sovler_tester:
	def __init__(self):
		random.seed(77) 

	def generatePWSP(self):
		budget = Fraction(1,1)
		pWSP = []

		while True:
			newTask = math.floor(1/random.random())
			newTaskDensity = Fraction(1, newTask)

			if newTaskDensity <= budget:
				pWSP.append(newTask)
				budget = budget - newTaskDensity
			else:
				break

		if budget > 0:
			pWSP.append(math.ceil(1/budget))
			budget = budget - Fraction(1, math.ceil(1/budget))
		
		if budget < 0:
			print("generation is broken,", pWSP, "is overly dense!")

		pWSP.sort()

		return pWSP

	def generateManyPWSPs(self, desiredNumberOfPSPs):
		manyPSPs = []

		while len(manyPSPs) < desiredNumberOfPSPs:
			newPSP = self.generatePWSP()

			if newPSP not in manyPSPs:
				manyPSPs.append(newPSP)

		return manyPSPs

	def indefinitelySolvePSPs(self, testingGraph, testingNaive, twoToTheK):
		if testingGraph:
			f = open("round_one.csv", "w")
			f.write("problem, graph, naive, opt, foresight, \n")
		elif testingNaive:
			f = open("round_two.csv", "w")
			f.write("problem, naive, opt, foresight, \n")
		else:
			f = open("round_three.csv", "w")
			f.write("problem, opt, foresight, \n")

		solvedPWSPs = []

		while(True):
			newPWSP = self.generatePWSP()

			if twoToTheK:
				cap = 2**(len(newPWSP) - 1)
				
				for i in range(len(newPWSP)):
					if newPWSP[i] > cap:
						newPWSP[i] = cap



			if newPWSP not in solvedPWSPs:
				solvedPWSPs.append(newPWSP)

				print(newPWSP)
				problem = str(newPWSP)
				problem = problem.replace(",", ";")

				f.write(problem + ",")

				if testingGraph:
					print("graph")
					graphSolver = solver_graph.solver_graph(newPWSP)
					graphSolver.solve()
					print(graphSolver.solveTimeCost)
					f.write(str(graphSolver.solveTimeCost) + ",")

				if testingNaive or testingGraph:
					print("naive")
					naiveSolver = solver_naive.solver_naive(newPWSP, False, False)
					naiveSolver.solve()
					print(naiveSolver.solveTimeCost)
					f.write(str(naiveSolver.solveTimeCost) + ",")

				print("opt")
				optSolver = solver_opt.PinwheelSolver(newPWSP, False, False)
				optSolver.solve()
				print(optSolver.solveTimeCost)
				f.write(str(optSolver.solveTimeCost) + ",")


				print("foresight")
				foresightSolver = solver_foresight.solver_foresight(newPWSP, False, False, False)
				foresightSolver.solve()
				print(foresightSolver.solveTimeCost)
				f.write(str(foresightSolver.solveTimeCost) + ",")

				print()
				f.write("\n")

	def round_one(self):
		self.indefinitelySolvePSPs(testingGraph = True, testingNaive = True, twoToTheK = True)

	def round_two(self):
		self.indefinitelySolvePSPs(testingGraph = False, testingNaive = True, twoToTheK = True)

	def round_three(self):
		self.indefinitelySolvePSPs(testingGraph = False, testingNaive = False, twoToTheK = True)




ourPWST = pinwheel_sovler_tester()

our_input = int(sys.argv[1])

if our_input == 1:
	ourPWST.round_one()
elif our_input == 2:
	ourPWST.round_two()
elif our_input == 3:
	ourPWST.round_three()
else:
	print("invalid input - please use 1, 2 or 3 to choose between the three rounds")

