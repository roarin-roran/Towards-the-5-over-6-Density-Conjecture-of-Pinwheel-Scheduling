from fractions import Fraction

# checks the certificate of a pareto with density instance.
class certificate_checker:
	def __init__(self, certificate, densityLimit, lengthLimit):
		self.certificateTrie  = certificate[0]
		self.certificateSolns = certificate[1]
		self.densityLimit = densityLimit
		self.lengthLimit  = lengthLimit

	def check(self):
		print("checking certificate")

		print("checking solutions")
		self.checkAllSolutions()
		print("all solutions ok")

		print("checking trie")
		self.checkBranchNode(self.certificateTrie)
		print("trie ok")

		# all test failures cause the program to exit
		print("certificate successfully tested")
	
	def checkAllSolutions(self):
		for i in self.certificateSolns:
			self.checkSolution(self.certificateSolns[i])

	def checkSolution(self, currentNode):
		# initialise
		solnQuality = [0]*len(currentNode.problem)
		gapSinceCut = [0]*len(currentNode.problem)

		# double soln to get wrap around
		for i in currentNode.solution*2:
			# cut
			gapSinceCut[i] = 0

			for j in range(len(gapSinceCut)):
				# grow
				gapSinceCut[j] += 1
				# increment max
				if gapSinceCut[j] > solnQuality[j]:
					solnQuality[j] = gapSinceCut[j]

		solnQualityInsufficient = False
		for i in range(len(solnQuality)):

			if solnQuality[i] > currentNode.problem[i]:
				solnQualityInsufficient = True
				break

		if solnQualityInsufficient:
			print("issue with the solution to ", currentNode.problem, "it only has quality", solnQuality)
			print("solution: ", currentNode.solution)
			exit(0)

	def recursivelycheck(self, currentNode):
		for i in currentNode.children:
			# if the child is a leaf, check it as a leaf
			if len(i.problemGiven) == self.lengthLimit:
				self.checkLeafNode(i)
			# else, check it as a branch
			else:
				self.checkBranchNode(i)

	def checkLeafNode(self, currentNode):
		#print("checking leaf node",  currentNode.problemSolved, currentNode.solution)
		self.checkDomination(currentNode)

	def checkBranchNode(self, currentNode):
		#print("checking branch node", currentNode.problemGiven)

		# first child is correct wrt density
		#print("checking first child of", currentNode.problemGiven)			
		self.checkFirstChild(currentNode.children[0])
		# continuous, with no gaps
		self.checkContinuityOfChildren(currentNode)
		# last char is RR complete (either problem given or problem sovled.)
		self.checkRRCompleteness(currentNode)

		self.recursivelycheck(currentNode)

	def checkDomination(self, currentNode):
		#print("checking domination", currentNode.problemGiven, currentNode.problemSolved)
		# solved problem is leq input problem

		problemSolved = self.certificateSolns[currentNode.problemSolved].problem

		for i in range(len(currentNode.problemGiven)):
			if currentNode.problemGiven[i] < problemSolved[i]:
				print(i)
				print("domination issue!", currentNode.problemGiven, "is not solved by", problemSolved)
				exit(0)

	# uses density and escalation to ensure that there should be no earlier child
	def checkFirstChild(self, firstChild):
		#print("checking first child: ", firstChild.problemGiven)


		# to be correct, this node must have acceptable density and it's hypothetical
		# replacement must not be permissible

		if self.hasAcceptableDensity(firstChild.problemGiven):

			# if this node has unit length, we can make a hypothetical child
			if len(firstChild.problemGiven) == 1:
				hypotheticalChild = True
			# if this unit has two selfsimilar tasks at the end, we can't reduce the last one
			# and this node is fine as a first child
			elif firstChild.problemGiven[-1] == firstChild.problemGiven[-2]:			
				hypotheticalChild = False
				firstChildOK = True
			# otherwise, we can make a hypothetical child
			else:
				hypotheticalChild = True

			# if a hypothetical child can sensibly exist, make one.
			if hypotheticalChild:
				# tighten the last task of the hypothetical child
				hypotheticalChild = firstChild.problemGiven.copy()
				hypotheticalChild[-1] = hypotheticalChild[-1] - 1

				# result must be invalid, or else we whould have the hypothetical child
				if self.hasAcceptableDensity(hypotheticalChild):
					firstChildOK = False
				else:
					firstChildOK = True

		# exit if the test is failed, silent pass
		if not firstChildOK:
			print("the node", firstChild.problemGiven, "should have a predecessor")
			exit(0)

	@staticmethod
	def checkContinuityOfChildren(currentNode):
		firstChildAddition = currentNode.children[0].problemGiven[-1]

		for i in range(len(currentNode.children) - 1):
			if currentNode.children[i+1].problemGiven[-1] != firstChildAddition + i + 1:
				print("noncontinuity in the child array of", currentNode.problemGiven)
				exit(0)

	def checkRRCompleteness(self, currentNode):
		rHLDescentant = self.recursivelyFindRightHandLeafDescendant(currentNode)

		symmetryOfPG = 0
		for i in reversed(rHLDescentant.problemGiven):
			if i == rHLDescentant.problemGiven[-1]:
				symmetryOfPG = symmetryOfPG + 1

		symmetryOfPS = 0

		try:
			problemSolved = self.certificateSolns[rHLDescentant.problemSolved].problem
		except Exception as e:
			print("issue with the dominant of", currentNode.problemGiven)
			print("problem solved (", currentNode.problemSolved, ") does not exist, causing:")
			print(repr(e))

			exit(0)


		for i in reversed(problemSolved):
			if i == problemSolved[-1]:
				symmetryOfPS = symmetryOfPS + 1

		symmetryOfCurrentNode = max(symmetryOfPS, symmetryOfPG)

		if len(currentNode.problemGiven) + symmetryOfCurrentNode < self.lengthLimit:
			print("issue with the round robin completeness of", currentNode.problemGiven, end='')
			print(", which has rr descendants of", rHLDescentant.problemGiven, "and", rHLDescentant.problemSolved)
			exit(0)

	def recursivelyFindRightHandLeafDescendant(self, currentNode):
		if len(currentNode.problemGiven) == self.lengthLimit:
			return currentNode
		else:
			return self.recursivelyFindRightHandLeafDescendant(currentNode.children[-1])

	def hasAcceptableDensity(self, solutionFragment):
		density = certificate_checker.getDensity(solutionFragment)
		roundingNumber = 2**(self.lengthLimit - 1)

		if density < self.densityLimit:
			return True
		elif density == self.densityLimit:
			if len(solutionFragment) == self.lengthLimit:
				return True
			else:
				return False
		elif solutionFragment[-1] == roundingNumber:
			solutionFragmentCopy = solutionFragment.copy()
			while solutionFragmentCopy[-1] == roundingNumber:
				solutionFragmentCopy.pop()
			
			#print(solutionFragment, solutionFragmentCopy)
			#input()

			return self.hasAcceptableDensity(solutionFragmentCopy)
		else:
			return False






			#yay! infinite recursion!

	@staticmethod
	def getDensity(solutionFragment):
		density = Fraction(0,1)

		for i in solutionFragment:
			density += Fraction(1,i)

		return density

	@staticmethod
	def recursivelyPrint(currentNode):
		tabbing = len(currentNode.problemGiven)*"\t"

		print(tabbing, currentNode.problemGiven)
		if len(currentNode.children) > 0:
			print(tabbing, "(")
			for i in currentNode.children:
				certificate_checker.recursivelyPrint(i)
			print(tabbing, ")")
	