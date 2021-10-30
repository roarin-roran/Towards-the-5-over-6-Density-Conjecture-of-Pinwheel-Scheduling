from fractions import Fraction
from collections import OrderedDict
from math import floor, ceil
from operator import itemgetter
import multiprocessing as mp
import time
import solver_foresight
import solver_opt
import solver_naive

# constants
# solvers
NAIVE     = 0
OPTIMISED = 1
FORESIGHT = 2
# approximation methods
NO_APPROX   = 0
P_ONE       = 1
P_FIVE      = 2
P_K_MINUS_1 = 3

CACHE_SIZE = False

class pareto_w_density:
	class AllParetoSurfaces:
		class ParetoSurface:
			# initiate ParetoSurface
			def __init__(self, numberOfTasks):
				self.numberOfTasks = numberOfTasks

				# set the pareto surface up so that elements can be compared to it in reverse order of symmertry
				self.symmetries = OrderedDict()
				for i in reversed(range(numberOfTasks)):
					self.symmetries[i+1] = {}

			def getNumberOfProblems(self):
				return len(self.getAllProblems())


		# AllParetoSurfaces
		def __init__(self, lengthLimit, allSolutionsUsed):
			self.lengthLimit = lengthLimit
			self.allSolutionsUsed = allSolutionsUsed

			# set up each lengths array as it's own pareto surface of the given length
			self.lengths = OrderedDict()
			for i in reversed(range(lengthLimit)):
				self.lengths[i+1] = self.ParetoSurface(i+1)

		# add a node to the appropriate pareto surface.
		def push(self, pSNode):
			paretoSurface = self.lengths[len(pSNode.problem)]

			symmetryGroup = paretoSurface.symmetries[pareto_w_density.getSymmetry(pSNode.problem)]

			# follow/build the node path to the location of the added node
			currentNode = symmetryGroup
			for i in range(len(pSNode.problem) - 1):
				if pSNode.problem[i] not in currentNode:
					currentNode[pSNode.problem[i]] = {}

				currentNode = currentNode[pSNode.problem[i]]

			currentNode[pSNode.problem[-1]] = pSNode

		# remove a node from the pareto surface. if it doesn't exist - throws errors
		def pop(self, problem):

			paretoSurface = self.lengths[len(problem)]

			symmetryGroup = paretoSurface.symmetries[pareto_w_density.getSymmetry(problem)]

			allPathNodes = []

			# find all nodes on the path to the deleted node
			currentNode = symmetryGroup
			for i in range(len(problem)):
				allPathNodes.append(currentNode)

				try:
					currentNode = currentNode[problem[i]]
				except:
					print("there was an error when accessing element", i, "of", problem)
					print("this is mostly likely caused by the problem not existing in the following list:")

					self.printAllPSPs(len(problem))

					print("exiting")

					exit(1)

			# find and remove the thread that only contains the problem to be deleted
			for i in range(len(allPathNodes)):

				if len(allPathNodes[-1-i]) > 1:
					allPathNodes[-1-i].pop(problem[-1-i])

					return currentNode

			# if it hasn't returned, then this is the only problem in its symmetry group!
			symmetryGroup.pop(problem[0])

		# compares a solved PSP with the existing pareto surface, 
		# if that PSP dominates any points, they're removed
		# if that PSP is undominated, it is added
		def compareSolved(self, solvedPSP, soln):
			symmetry = pareto_w_density.getSymmetry(solvedPSP)

			paretoSurface = self.lengths[len(solvedPSP)]

			# trash collection - all dominated problems, waiting to be removed
			self.allDominatedProblems = []

			for i in paretoSurface.symmetries:
				symmetryGroup = paretoSurface.symmetries[i]
				
				if i > symmetry:
					dominator = self.getFirstDominator(solvedPSP, symmetryGroup)
		
					if dominator:
						break
				elif i == symmetry:
					self.getAllDominated(solvedPSP, solvedPSP, symmetryGroup)
			
					dominator = self.getFirstDominator(solvedPSP, symmetryGroup)
					
					if dominator:
						break
					
				else:
					self.getAllDominated(solvedPSP, solvedPSP, symmetryGroup)


			for i in self.allDominatedProblems:
				self.pop(i)


			if dominator:
				outputSolnNo = dominator.solutionNumber
			else:
				pSNode = pareto_w_density.PSNode(solvedPSP, soln, self.allSolutionsUsed)
				self.push(pSNode)

				outputSolnNo = pSNode.solutionNumber

			return outputSolnNo


		# returns the most symmetric dominator, if there is one
		def compareUnsolved(self, newPSP):
			dominator = False

			paretoSurface = self.lengths[len(newPSP)]

			for i in paretoSurface.symmetries:
				symmetryGroup = paretoSurface.symmetries[i]

				dominator = self.getFirstDominator(newPSP, symmetryGroup)
		
				if dominator:
					break
			
			return dominator

		# recursively find and return the first dominating element, if there is one - else return false.
		def getFirstDominator(self, solvedPSPFragment, subTree):
			foundDominator = False

			# if there are remaining elements
			if len(solvedPSPFragment) > 0:
				# copy to prevent leaks
				solvedPSPFragment = solvedPSPFragment.copy()
				elementOfComparison = solvedPSPFragment.pop(0)

				# recurse on all compatible subtrees of the inputted subtree
				for key in subTree:
					if key <= elementOfComparison:
						foundDominator = self.getFirstDominator(solvedPSPFragment, subTree[key])

					# if a dominator is found in this step, return it
					if foundDominator:
						break
				# getting to the end means we found a dominator! - is this true??
			else:
				# if there are no elements left, the subtree is just the PSNode at the leaf
				return subTree
			
			return foundDominator
				
		# recursively finds all problems
		def getAllDominated(self, solvedPSP, solvedPSPFragment, subTree):
			# if there are remaining elements
			if len(solvedPSPFragment) > 0:
				# copy to prevent leaks
				solvedPSPFragment = solvedPSPFragment.copy()
				elementOfComparison = solvedPSPFragment.pop(0)

				# recurse on all compatible subtrees of the inputted subtree
				for key in subTree:
					if key >= elementOfComparison:
						self.getAllDominated(solvedPSP, solvedPSPFragment, subTree[key])

			else:
				# problems will dominate themselves without this line
				if subTree.problem == solvedPSP:
					pass 
				else:
					self.allDominatedProblems.append(subTree.problem)

		def printAllPSPs(self, length):
			allPSPs = []

			paretoSurface = self.lengths[length]
			for i in paretoSurface.symmetries:
				self.recursivelyAccessAllPSPs(paretoSurface.symmetries[i], length, allPSPs)

			print()
			j = 1
			for i in allPSPs:
				print(j, ":", i[0])
				j = j+1

		def recursivelyAccessAllPSPs(self, pTNode, length, outputList):
			if type(pTNode) == dict:
				for i in pTNode:
					self.recursivelyAccessAllPSPs(pTNode[i], length - 1, outputList)
			else:
				outputList.append([pTNode.problem, pTNode.solution])

	# node on the pareto surface
	class PSNode:
		def __init__(self, problem, solution, allSolutionsUsed):
			self.problem = problem
			self.solution = solution

			solnNovel = True
			for i in allSolutionsUsed:
				if problem == allSolutionsUsed[i].problem:
					solnNovel = False 
					print("non-unique PSNode generated!")
					print(problem)
					exit(0)
					break

			self.solutionNumber = len(allSolutionsUsed)
			allSolutionsUsed[self.solutionNumber] = self

			#input()

	# node in the pareto trie - the certificate, used slightly differently from the ParetoSurface
	class CertificateTrieNode:
		def __init__(self, problemGiven, problemSolved, parent = False):
			self.problemGiven  = problemGiven
			self.problemSolved = problemSolved

			self.children = []

			if parent:
				parent.children.append(self)

		def printCTNode(self):
			print([self.problemGiven, self.problemSolved], end = '')
			print(" (", end = '')

			j = False
			for i in self.children:
				if j:
					print(", ", end = '')
				else:
					j = True
				print(i.problemGiven, end = '')

			print(")", end = '')
			print()

	
	# * * * * * * * ** * * * * * * *
	# * * * * * * * ** * * * * * * *
	# * * * * * INITIATION * * * * *
	# * * * * * * * ** * * * * * * *
	# * * * * * * * ** * * * * * * *


	# initiate pareto_w_density
	def __init__(self, densityLimit, lengthLimit, solver, approxMethod, keepCerts, outputFilename):
		print("initiating for lengthLimit", lengthLimit)

		self.densityLimit = densityLimit
		self.lengthLimit = lengthLimit
		self.solver = solver
		self.approxMethod = approxMethod
		self.keepCerts = keepCerts
		self.outputFilename = outputFilename

		# a list of all solutions ever used by the solver, maintained as part of the certificate.
		self.allSolutionsUsed = {}

		# initiate Pareto surface
		self.pSurfaces = self.AllParetoSurfaces(lengthLimit, self.allSolutionsUsed)

		# initialise cap
		self.maxTaskSeperation = 2**(lengthLimit-1)

		self.timeSpentSolving = 0

		self.numberOfConsideredProblems = 0
		self.numberOfUpdates = 0

		# initiate the Pareto surface with an approximation
		self.unfoldedSurface = []

		if approxMethod == NO_APPROX:
			pass
		elif approxMethod == P_ONE:
			self.approximatePateroSurface_P_One()
		elif approxMethod == P_FIVE:
			self.approximatePateroSurface_P_Five()
		elif approxMethod == P_K_MINUS_1:
			self.approximatePateroSurface_P_K_Minus_One()

		# findme

		if CACHE_SIZE:
			# initialise the cache
			# first, find a suitable first value.
			for i in self.allSolutionsUsed:
				if len(self.allSolutionsUsed[i].problem) == self.lengthLimit:
					self.cache = self.allSolutionsUsed[i]
					break
			self.cacheSymmetry = pareto_w_density.getSymmetry(self.cache.problem)
			self.numRepeats = 0
			self.maxNumRepeats = 0


	# * * * * * * * * * * * * * * * * * * * * *
	# * * * * * SURFACE APPROXIMATION * * * * *
	# * * * * * * * * * * * * * * * * * * * * *

	# preprocessing step, recursively generating the 5/6 surface for all smaller problems
	# down to k=5, for which the full pareto surface is provided.
	def approximatePateroSurface_P_K_Minus_One(self):
		ParetoSurface = []
		#self.generateUnfoldedSurface()

		# if the full pareto surface is unavailable, recursively approximate it.
		if self.lengthLimit > 5:
			# find the 5/6 surface for k-1
			nextSmallestSurfaceSolver = pareto_w_density(self.densityLimit, 
														 self.lengthLimit - 1, 
														 self.solver, 
														 self.approxMethod, 
														 self.keepCerts, 
														 self.outputFilename)
			nextSmallestSurfaceSolver.solve()

			# consider all elements in the next smallest surface, 
			for length in nextSmallestSurfaceSolver.pSurfaces.lengths:
				allPSPs = []
				paretoSurface = nextSmallestSurfaceSolver.pSurfaces.lengths[length]
				for i in paretoSurface.symmetries:
					nextSmallestSurfaceSolver.pSurfaces.recursivelyAccessAllPSPs(paretoSurface.symmetries[i], length, allPSPs)

				for i in allPSPs:

					# add the element to the relevent pareto surface
					self.compareSolved(i[0], i[1])

					# unfold the children of this element
					allChildren = self.unfoldAllChildren(i[0])

					# generate the unfolded children of this element
					for k in allChildren:
						self.recursivelyGenerateUnfoldedSurface(k, i[0], i[1])

			# compare all elements of this unfolded surface with the pareto surface
			for i in self.unfoldedSurface:
				self.compareSolved(i[0], i[1])
		

		# if the full pareto surface is known, use it.
		elif self.lengthLimit == 5:
			self.approximatePateroSurface_P_Five()

	# light version of the preprocessing, maintained mostly for comparisons. faster for large k,
	# but probably much slower when rapid searching is added. but simple and reliable.
	def approximatePateroSurface_P_Five(self):
		self.approximatePateroSurface_P_One()

		ParetoSurface = []

		ParetoSurface.append([[3,4,5,8],		[0,1,3,0,2,1,0,2]])
		ParetoSurface.append([[3,5,5,5],		[0,1,2,0,3,1,0,2,3]])

		ParetoSurface.append([[2,6,8,10,16],	[0,1,0,2,0,4,0,1,0,3,0,2,0,1,0,3]])
		ParetoSurface.append([[2,6,10,10,10],	[0,1,0,2,0,3,0,1,0,4,0,2,0,1,0,3,0,4]])
		ParetoSurface.append([[3,4,5,14,14],	[0,1,2,0,3,1,0,2,0,1,4,0,2,1]])
		ParetoSurface.append([[3,4,6,10,16],	[0,1,2,0,3,1,0,2,0,1,3,0,2,1,0,4]])
		ParetoSurface.append([[3,4,6,11,11],	[0,1,2,0,4,1,0,2,1,0,3]])
		ParetoSurface.append([[3,4,8,8,8],		[0,1,3,0,4,1,0,2]])
		ParetoSurface.append([[3,5,5,9,9],		[0,1,4,0,2,1,0,3,2]])
		ParetoSurface.append([[3,5,6,7,12],		[0,1,3,0,2,1,0,3,1,0,2,4]])
		ParetoSurface.append([[3,5,7,7,9],		[0,1,2,0,3,1,0,4,2,0,1,3,0,2,1,0,4,3]])
		ParetoSurface.append([[3,5,7,8,8],		[0,1,2,0,3,1,0,4,0,2,1,0,3,4]])
		ParetoSurface.append([[4,4,5,7,12],		[0,1,2,3,0,1,4,2,0,1,3,2]])
		ParetoSurface.append([[4,5,5,6,10],		[0,1,2,4,0,3,1,2,0,3]])
		ParetoSurface.append([[4,5,5,7,7],		[0,1,4,2,0,3,1,0,2,4,1,0,3,2]])

		if self.lengthLimit > 5:
			for i in ParetoSurface:
				allChildren = self.unfoldAllChildren(i[0])

				for j in allChildren:
					self.recursivelyGenerateUnfoldedSurface(j, i[0], i[1])

			for i in self.unfoldedSurface:
				self.compareSolved(i[0], i[1])
		
		for i in ParetoSurface:
			self.compareSolved(i[0], i[1])


	# minimal version of preprocessing - unfolds [0]
	def approximatePateroSurface_P_One(self):
		# unfold [1]
		self.generateUnfoldedSurface()

		# compare all elements to the pareto surface
		for i in self.unfoldedSurface:
			self.compareSolved(i[0], i[1])

	# * * * * * * * * * * * * * * * * * * *
	# * * * * * SURFACE UNFOLDING * * * * *
	# * * * * * * * * * * * * * * * * * * *

	# recursively generates all PWS instances generatable by unfolding
	def generateUnfoldedSurface(self):
		print("generating unfolded surface (this may take a few minutes)")

		self.unfoldedSurface = []

		self.unfoldedSurface.append([[1], [0]])
		
		allIntegerInputs = self.unfoldAllChildren([1])
		for i in allIntegerInputs:
			self.recursivelyGenerateUnfoldedSurface(i)

		self.unfoldedSurface = sorted(self.unfoldedSurface, key=itemgetter(0))

		for i in range(len(self.unfoldedSurface)):
			print("unfoldedSurface member", i, ": ", self.unfoldedSurface[i])

		for i in self.unfoldedSurface:
			self.compareSolved(i[0], i[1])

	# recursively function for the unfolding surface
	def recursivelyGenerateUnfoldedSurface(self, prefix, parent = [1], parentSolution = [0]):
		# if the prefix is unexamined, proceed to solve it, record it and find it's children
		if self.prefixNotInOutputList(prefix):
			solution = pareto_w_density.unfoldSolution(parent, prefix, parentSolution)

			self.unfoldedSurface.append([prefix, solution])

			# if the prefix is not full length, generate its descendents
			if len(prefix) != self.lengthLimit:
				
				# generate all children of this prefix				
				allChildren = self.unfoldAllChildren(prefix)

				# for all distinct elements in the prefix
				for newPrefix in allChildren:

					# if the new prefix has already been examined, don't recurse
					if self.prefixNotInOutputList(newPrefix):
						prefix 		= prefix.copy()
						newPrefix 	= newPrefix.copy()
						solution 	= solution.copy()

						#self.recursivelyGenerateUnfoldedSurface(newPrefix, solution)
						self.recursivelyGenerateUnfoldedSurface(newPrefix, prefix, solution)

	# unfolds a given prefix, generating all of its immediate unfolded descendents
	def unfoldAllChildren(self, prefix):
		allChildren = []

		# for all distinct elements in the prefix
		for i in range(len(prefix)):

			# proceed to unfolding only if this is the first time you've seen this character
			proceedToUnfoldingStage = True
			if i != 0:
				if prefix[i-1] == prefix[i]:
					proceedToUnfoldingStage = False

			if proceedToUnfoldingStage:
				# unfold the ith character
				oldElement = prefix[i]

				# replace the element with each possible unfolding
				for j in range(self.lengthLimit - len(prefix)):
					
					# generate the new problem
					newPrefix = prefix.copy()
					newPrefix[i] = (j+2)*oldElement
					# add the placeholders for the new element
					for k in range(j+1):
						newPrefix.append((j+2)*oldElement)

					# sort the problem
					newPrefix.sort()

					allChildren.append(newPrefix)

		return allChildren

	# returns true iff the given input is not in the outputList
	def prefixNotInOutputList(self, prefix):
		prefixNotInOutputList = True
		for i in self.unfoldedSurface:
			if i[0] == prefix:
				prefixNotInOutputList = False

		return prefixNotInOutputList



	# * * * * * * * * *  * * * * * * * *
	# * * * * * * * * *  * * * * * * * *
	# * * * * * PARETO SURFACE * * * * *
	# * * * * * * * * *  * * * * * * * *
	# * * * * * * * * *  * * * * * * * *


	def addToParetoSurfaces(self, pSNode):
		return self.pSurfaces.push(pSNode)

	def popFromParetoSurface(self, problem):
		return self.pSurfaces.pop(problem)

	def printAllPSPs(self, length):
		return self.pSurfaces.printAllPSPs(length)

	def compareSolved(self, solvedPSP, soln):
		return self.pSurfaces.compareSolved(solvedPSP, soln)

	def compareUnsolved(self, newPSP):
		self.numberOfConsideredProblems += 1

		return self.pSurfaces.compareUnsolved(newPSP)

	# * * * * * * * * * * * * * *
	# * * * * * UTILITY * * * * *
	# * * * * * * * * * * * * * *

	# returns the number of self similar terminating characters of a pinwheel scheduling problem
	@staticmethod
	def getSymmetry(pWSProblem):
		symmetry = 0 
		for i in reversed(pWSProblem):
			if i == pWSProblem[-1]:
				symmetry += 1
			else:
				break

		return symmetry

	# returns the density of an inputted pinwheel scheduling instance.
	@staticmethod
	def getDensity(partialSolution):
		density = Fraction(0,1)

		for i in partialSolution:
			density += Fraction(1,i)

		return density

	@staticmethod
	def isFirstDominant(firstProblem, secondProblem):
		if len(firstProblem) != len(secondProblem):
			print("invalid comparison!")

		for i in range(len(firstProblem)):
			if firstProblem[i] > secondProblem[i]:
				return False

		return True

	# * * * * * * * * * * * * * *
	# * * * * * * * * * * * * * *
	# * * * * * SOLVING * * * * *
	# * * * * * * * * * * * * * *
	# * * * * * * * * * * * * * *

	# returns the smallest denominator that can be inverted and added to a prefix
	# without raising the total density above the limit
	def getNextTaskSeperation(self, prefix):
		if prefix[-1] == self.maxTaskSeperation:
			nextTaskSeperation = self.maxTaskSeperation
		elif self.densityLimit - self.getDensity(prefix) < 0.00000001:
			nextTaskSeperation = self.maxTaskSeperation
		else:
			spaceRemaining = Fraction(1, self.densityLimit - self.getDensity(prefix))

			nextTaskSeperation = ceil(spaceRemaining)
			
			# the nextTaskSeperation cannot be less than the last task of the prefix
			if nextTaskSeperation < prefix[-1]:
				nextTaskSeperation = prefix[-1]

			# if this is not the last task to be added to the prefix, the density needs to be
			# less than 5/6, not less than or equal to
			if len(prefix) < self.lengthLimit - 1:
				
				newDensity = Fraction(1, nextTaskSeperation) + self.getDensity(prefix)
				if newDensity == self.densityLimit:
					nextTaskSeperation += 1

			if nextTaskSeperation > self.maxTaskSeperation:
				nextTaskSeperation = self.maxTaskSeperation

		return nextTaskSeperation

	# given any node, finds and solves all descendent leaf nodes
	def recursivelySolveParetoTree(self, prefix, parent):
		nextTaskSeperation = self.getNextTaskSeperation(prefix)
		numLevelsToClimb = 0

		# explore the children of this node until an appropriately symmetric solvable one is found.
		while numLevelsToClimb == 0:
			currentPrefix = prefix.copy()
			currentPrefix.append(nextTaskSeperation)

			# if our prefix is not full length, proceed recursively
			if len(prefix) != self.lengthLimit:
				# if the modified prefix is exactly full length, pass the parent back in, this'll end up getting
				# recorded by the if clause. this node will have no children
				if len(currentPrefix) == self.lengthLimit:
					numLevelsToClimb = self.recursivelySolveParetoTree(currentPrefix, parent) - 1
				# if the prefix is still not full length, keep going
				else:
					pTNode = self.CertificateTrieNode(currentPrefix, False, parent)

					numLevelsToClimb = self.recursivelySolveParetoTree(currentPrefix, pTNode) - 1
				nextTaskSeperation += 1
			# if our prefix is full length, solve it
			else:	
				# first, check for domination - a pre-existing solution to the problem.
				# dominatingNode = self.compareUnsolved(prefix)
				
				dominatingNode = self.findSolution(prefix)




				if dominatingNode:
					numLevelsToClimb = max(self.getSymmetry(prefix), self.getSymmetry(dominatingNode.problem))
					print("( k=",self.lengthLimit,")  \tdominated :	\t", prefix, "\tby:\t", dominatingNode.problem)

					pTNode = self.CertificateTrieNode(prefix, dominatingNode.solutionNumber, parent)

				else:
					print("solve:		\t\t\t", prefix)
					pTNode = self.parralelFoldAndSolve(prefix, parent)

					# if the prefix is unsolvable, foldAndSolve will exit
					print("solved: \t\t\t\t", prefix, "\tby:\t", self.lastPinwheelInstanceSolved)
					
					numLevelsToClimb = max(self.getSymmetry(prefix), self.getSymmetry(self.lastPinwheelInstanceSolved))

				# constant to fight the removal of 1 each time				
				numLevelsToClimb += 1

		return numLevelsToClimb

	def findSolution(self, prefix):
		# findme
		# print("have ", self.numRepeats, " max ", self.maxNumRepeats)
		
		if CACHE_SIZE:
			self.maxNumRepeats = CACHE_SIZE * self.lengthLimit

			# if we're not timed out
			if self.numRepeats < self.maxNumRepeats:
				if self.isFirstDominant(self.cache.problem, prefix):
					# print("the solution we have works!", self.cache.problem, " dominates ", prefix)
					
					self.numRepeats += 1
				else:
					#print("the solution we have doesn't work", self.cache.problem, " doesn't dominate ", prefix)
					
					newProblem = self.compareUnsolved(prefix)

					self.numRepeats = 0
					if newProblem:
						self.cache = newProblem
						self.cacheSymmetry = pareto_w_density.getSymmetry(self.cache.problem)
						
						#self.maxNumRepeats -= 1
					else:
						#self.maxNumRepeats = 0
						return False
			# if we are timed out
			else:
				# refresh the cache
				newProblem = self.compareUnsolved(prefix)

				# if we have a dominant
				if newProblem:
					self.cache = newProblem
					newSymmetry = pareto_w_density.getSymmetry(self.cache.problem)
				
					#if newSymmetry > self.cacheSymmetry:
					#	self.maxNumRepeats -= 1
					#else:
					#	self.maxNumRepeats += 1

					self.cacheSymmetry = newSymmetry
				else:
					# if not, a new cache will be provided
					self.numRepeats = 0
					#self.maxNumRepeats = 0
					return False

			# print("final output", self.cache.problem)
			# input()

			return self.cache
		else:
			return self.compareUnsolved(prefix)

	# solves thw whole pareto trie by solving the subtrees of all input nodes
	def solve(self):
		# start the clock
		startTime = time.time()

		root = self.CertificateTrieNode([], False)

		# for all possible first options.
		# note that [1] is excluded on density grounds - this is done here
		# to simplify the certificate checker
		for i in range(self.lengthLimit - 1):
			bough = self.CertificateTrieNode([i+2], False, root)

			self.recursivelySolveParetoTree([i+2], bough)

		self.timeTaken = time.time() - startTime

		print()
		print("complete!")
		print("time cost:           ", self.timeTaken)
		searchPercentage = (100 * (self.timeTaken - self.timeSpentSolving)/self.timeTaken)
		print("search percentage:   ", round(searchPercentage), "%")
		print("problems considered: ", self.numberOfConsideredProblems)
		print("problems solved:     ", self.numberOfUpdates)

		print()

		f = open(self.outputFilename, "a")
		f.write(str(self.lengthLimit))
		f.write(",")
		f.write(str(self.timeTaken))
		f.write(",")
		f.write(str(self.timeSpentSolving))
		f.write(",")
		f.write(str(self.numberOfUpdates))
		#f.write(",")
		#f.write(str(self.pSurfaces.lengths[self.lengthLimit].getNumberOfProblems()))
		f.write("\n")
		f.close()

		return [root, self.allSolutionsUsed]



	# * * * * * * *  * * * * * * *
	# * * * * * FOLDIING * * * * *
	# * * * * * * *  * * * * * * *

	# returns an array of all potentially solvable foldings of an input, including the input itself
	@staticmethod
	def generateFoldedProblems(pWSProblem):
		# stores the various foldings of the problem
		# initialise with the origional problem
		foldedProblems = [[pWSProblem, 1, pWSProblem[-1]]]

		# determine how many poles should be combined - start from 1 so at least 2 should be folded
		# initialise
		foldingNumber = 1
		foldedProblem = pWSProblem.copy()

		# round and record all possibly solvable foldings
		while True:
			# the element that all later elements will be rounded to
			foldedElement = foldedProblem[-1-foldingNumber]

			# round all elements
			for i in range(foldingNumber + 1):
				foldedProblem[-1- i] = foldedElement

			# don't try to fold the first element - causes overrun
			# we also should not see anything that's solvable by RR and isn't
			# in the folded problems array we initialise with
			if foldedElement == pWSProblem[0]:
				break
			else:
				# count the length of the self-similar suffix
				while foldedProblem[-1-foldingNumber] == foldedElement:
					foldingNumber += 1

			# record
			foldedElement = floor(foldedProblem[-1]/foldingNumber)

			# avoid dividing by zero and also necesarily unsolvable instances
			if foldedElement > 1:
				foldedProblems.append([foldedProblem.copy(), foldingNumber, foldedElement])
			else:
				break

		# fold the problems:
		for i in foldedProblems:
			# remove the poles to be folded
			foldedProblem = i[0][:(len(pWSProblem) - i[1])]
			
			# add the folded pole
			foldedProblem.append(i[2])
			foldedProblem.sort()
			i[0] = foldedProblem.copy()

		# remove any overly dense problems
		for i in range(len(foldedProblems)):
			if pareto_w_density.getDensity(foldedProblems[i][0]) > 1:
				foldedProblems = foldedProblems[:i]
				break

		return foldedProblems

	# folds, solves and records a PWS problem, using the most symmetric possible solution.
	# as the origional input is included, if no input is solvable there's a hypothesis refutation.
	# legacy function used chiefly for bug fixing.
	def foldAndSolveSingleThread(self, pWSProblem, parent):
		foldedProblemDescriptions = self.generateFoldedProblems(pWSProblem)

		# attempt problems in order of solution symmetry
		for i in reversed(foldedProblemDescriptions):
			unfoldedSolutionOutput = self.solveFoldedProblemFromDescription(i)
			
			if unfoldedSolutionOutput[0] == True:

				newSolnNumber = self.compareSolved(unfoldedSolutionOutput[1], unfoldedSolutionOutput[2])

				self.lastPinwheelInstanceSolved = unfoldedSolutionOutput[1]

				pTNode = self.CertificateTrieNode(foldedProblemDescriptions[0][0], newSolnNumber, parent)

				self.timeSpentSolving += unfoldedSolutionOutput[3]

				if unfoldedSolutionOutput[4]:
					self.numberOfUpdates += 1

				break

		if unfoldedSolutionOutput[0] == True:
			return pTNode
		else:
			# if we can't solve any folded version, try and solve the actual input
			print("input", pWSProblem, "is unsolvable! some hypothesis has been refuted")
			exit(0)

	# solves a folded problem from a listable input, used for parralelisation.
	def solveFoldedProblemFromDescription(self, foldedProblemDescription):
		foldedProblem = foldedProblemDescription[0]
		foldingNumber = foldedProblemDescription[1]
		foldedElement = foldedProblemDescription[2]

		unfoldedAndSolvedProblem = []
		unfoldedSolution = []
		outputSolutionNumber = False

		density = self.getDensity(foldedProblem).numerator/self.getDensity(foldedProblem).denominator
		print("working on:", foldedProblemDescription, "which has density", density)

		if self.solver == NAIVE:
			pinwheelSolver = solver_naive.solver_naive(foldedProblem.copy(), False, False)
		elif self.solver == OPTIMISED:
			pinwheelSolver = solver_opt.PinwheelSolver(foldedProblem.copy(), False, False)
		elif self.solver == FORESIGHT:
			pinwheelSolver = solver_foresight.solver_foresight(foldedProblem.copy(), False, False, False)



		if pinwheelSolver.solve() == 1:

			foldedAndSolvedProblem = pinwheelSolver.getSolnQuality()
			foldedSolution		   = pinwheelSolver.getSolution()

			outputSolutionNumber = self.compareSolved(foldedAndSolvedProblem, foldedSolution)

			unfoldedAndSolvedProblem = self.unfoldSolvedProblem(foldedProblem, foldedAndSolvedProblem, foldingNumber, foldedElement)
			unfoldedSolution = self.unfoldSolution(foldedAndSolvedProblem, unfoldedAndSolvedProblem, foldedSolution)

			solvability = True

		else:
			print("missing feature: failure surface")
			solvability = False
				
		return [solvability, unfoldedAndSolvedProblem, unfoldedSolution, pinwheelSolver.solveTimeCost, outputSolutionNumber]

	# finds a direct or dominating solution to a given problem, in parralel
	def parralelFoldAndSolve(self, pWSProblem, parent):
		foldedProblemDescriptions = self.generateFoldedProblems(pWSProblem)

		print()
		print("solve, in parralel:")

		p = mp.Pool(len(foldedProblemDescriptions))

		for unfoldedSolutionOutput in p.imap_unordered(self.solveFoldedProblemFromDescription, foldedProblemDescriptions, chunksize=1):
			if unfoldedSolutionOutput[0] == True:
				p.terminate()
				print("terminating")
				print()

				newSolnNumber = self.compareSolved(unfoldedSolutionOutput[1], unfoldedSolutionOutput[2])

				self.lastPinwheelInstanceSolved = unfoldedSolutionOutput[1]

				pTNode = self.CertificateTrieNode(foldedProblemDescriptions[0][0], newSolnNumber, parent)

				self.timeSpentSolving += unfoldedSolutionOutput[3]

				if unfoldedSolutionOutput[4]:
					self.numberOfUpdates += 1

				break

		# just to be safe
		p.terminate()

		return pTNode

	# unfold the problem, based on the hardest problem solved by a schedule, 
	# and what was done to generate it
	@staticmethod
	def unfoldSolvedProblem(foldedProblem, foldedAndSolvedProblem, foldingNumber, foldedElement):
		foldedProblem = foldedProblem.copy()
		foldedAndSolvedProblem = foldedAndSolvedProblem.copy()

		# remove the folded element
		indexToUnfold = foldedProblem.index(foldedElement)

		valueToUnfold = foldedAndSolvedProblem.pop(indexToUnfold)

		# add its unfoldings
		unfoldedAndSolvedProblem = foldedAndSolvedProblem.copy()
		for i in range(foldingNumber):
			unfoldedAndSolvedProblem.append(valueToUnfold*foldingNumber)
		unfoldedAndSolvedProblem.sort()

		return unfoldedAndSolvedProblem
		
	# unfold the solution to a pinwheel scheduling problem, using a previous unfolding of the problem
	@staticmethod	
	def unfoldSolution(foldedAndSolvedProblem, unfoldedProblem, foldedSolution):
		if len(foldedAndSolvedProblem) == len(unfoldedProblem):
			return(foldedSolution)
		else:
			# find the unfolded task 
			for i in range(len(foldedAndSolvedProblem)):
				if foldedAndSolvedProblem[i] != unfoldedProblem[i]:
					unfoldLocation = i
					break

			# map the elements in the folded problem onto elements in the unfolded problem
			mapping = []
			j = 0
			for i in range(len(foldedAndSolvedProblem)):
				if i == unfoldLocation:
					mappingAddition = []

					blockingProblem = foldedAndSolvedProblem.copy()
					for k in range(len(unfoldedProblem)):
						if unfoldedProblem[k] in blockingProblem:
							blockingProblem.remove(unfoldedProblem[k])
						else:
							mappingAddition.append(k)
					mapping.append(mappingAddition)
				else:
					while j < len(unfoldedProblem) - 1:
						if foldedAndSolvedProblem[i] == unfoldedProblem[j]:
							break
						else:
							j = j + 1

					mapping.append(j)
					j = j + 1

			# construct the solution repeat unit
			unfoldedSolnSection = foldedSolution.copy()

			for i in range(len(unfoldedSolnSection)):
				if unfoldedSolnSection[i] == unfoldLocation:
					unfoldedSolnSection[i] = '*'
				else:
					unfoldedSolnSection[i] = mapping[unfoldedSolnSection[i]]

			# construct the solution blank
			unfoldedSolution = []
			for i in mapping[unfoldLocation]:
				unfoldedSolution = unfoldedSolution + unfoldedSolnSection


			# cycle through the solution replacing the wildcards with the appropriate values
			j = 0
			for i in range(len(unfoldedSolution)):
				if unfoldedSolution[i] == '*':
					unfoldedSolution[i] = mapping[unfoldLocation][j] 
					j = j + 1
					if j >= len(mapping[unfoldLocation]):
						j = 0

			return unfoldedSolution



