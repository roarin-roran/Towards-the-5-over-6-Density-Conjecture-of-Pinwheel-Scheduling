import time 


# an optimised solver for pinwheel systems with wildcards - 
# doesn't use minimum solution length, because this seems to slow things down for some reason
class solver_wildcard_opt:
	class moveStatePair:
		def __init__(self, move, state, bBSArray):
			self.move = move
			self.state = state
			self.bBSArray = bBSArray


	def __init__(self, knownPoles):
		self.knownPoles = knownPoles

		# solutions with a wildcard
		fullSolns = []
		#solutions without a wildcard
		staticSolns = []

		# start the exploration:
		self.partialSolution = []

		# * * * first day * * *

		# this is the state after day 1
		firstState = [1]*len(self.knownPoles)
		# the move made on day 1 - this is a free choice 
		# and we chose a wildcard
		firstMove = "*"

		initialBBSArray = self.constructInitialBBSArray()

		initialBBSArray = []

		self.partialSolution.append(self.moveStatePair("*", firstState, initialBBSArray))

		# * * * exploration variables * * *
		self.moveDownNext = True
		self.nextDownMove = 0

		# record all solutions that don't work for 
		self.tightSolutions = []

		self.solvable = 0

		self.setNaiveMinSolnLength()


	def constructInitialBBSArray(self):
		bBSArray = []

		for i in range(len(self.knownPoles) - 1):
			if self.knownPoles[i+1] == self.knownPoles[i]:
				bBSArray.append(i+1)
				
		# print enough stuff to check the bBSArray is correct
		if len(bBSArray) > 0 and False:
			print(self.knownPoles)
			print(bBSArray)
			input()

		return bBSArray

	# sets the minimum solution length - no point in testing any solution
	# shorter than this
	def setNaiveMinSolnLength(self):
		# first check the density, because any pinwheel case that would
		# fail on density grounds has an infinite minimum solution length
		naiveDensity = 0
		for i in self.knownPoles:
			naiveDensity += 1/i

		if naiveDensity > 1:
			self.solvable = -1
		else:
			# the length of the round robin solution.
			self.minSolnLength = len(self.knownPoles)
			minNumber = [1]*self.minSolnLength


			actedThisItteration = True
			while actedThisItteration:
				actedThisItteration = False
				for i in range(len(minNumber)):
					while (minNumber[i]/self.minSolnLength) < \
							(1/self.knownPoles[i]):
						minNumber[i] += 1
						self.minSolnLength += 1
						actedThisItteration = True





			#self.minSolnLength = 0




			print("minSolnLength = ", self.minSolnLength)

	def printPartialSolution(self):
		for i in self.partialSolution:
			print("move:\t", i.move, "\tstate:\t", i.state)

	# makes the next down move, tests it and determines what to do next
	def down(self):
		# test the next down move
		if self.nextDownMove != "*":
			if isinstance(self.nextDownMove, int):
				if self.nextDownMove < 0 or self.nextDownMove >= len(self.knownPoles):
					print("y'move's bad guv: bad int", self.nextDownMove)
					exit(0)
			else:
				print("y'move's bad guv: float or character")
				exit(0)

		# having passed input testing, increment all poles
		currentState = []
		for i in self.partialSolution[-1].state:
			currentState.append(i+1)

		# test for repetition
		if self.nextDownMove == self.partialSolution[-1].move:
			#print("repetition detected!", self.nextDownMove)
			# if the repeated move is a *, replace it with a 0
			if self.nextDownMove == "*":
				self.nextDownMove = 0
			# else, if it can be incremented increment it
			elif self.nextDownMove < (len(self.knownPoles) - 1):
				self.nextDownMove += 1
			# else ban it
			else:
				self.moveDownNext = False

		# ban moves which are blocked by symmetry
		if self.moveDownNext:
			if self.nextDownMove in self.partialSolution[-1].bBSArray:
				self.moveDownNext = False

		if self.moveDownNext:
			# having incremented all poles, test it for failure
			if self.testForFailure():
				self.moveDownNext = False
			else:
				currentBBSArray = self.partialSolution[-1].bBSArray
				
				# having tested the state for failure, cut the input
				if self.nextDownMove != "*":
					currentState[self.nextDownMove] = 0

					# if the move that's just been played frees up a new move,
					# remove that move from the BBSArray
					if self.nextDownMove + 1 in currentBBSArray:
						currentBBSArray.remove(self.nextDownMove + 1)

				# test for success - if this finds a wildcarded
				# solution we'll exit, otherwise we'll just backtrack.
				if self.testForSuccess(currentState):
					self.moveDownNext = False
				# if no solution is found, append the state and keep moving.
				else:
					self.partialSolution.append(self.moveStatePair(self.nextDownMove, currentState, currentBBSArray))
					self.nextDownMove = "*"
			#print(currentState)

	# tests for immediate failure of constraint violations.
	def testForFailure(self):
		# look at every pole
		for i in range(len(self.knownPoles)):
			# if it's too high, this pole causes a failure
			if self.partialSolution[-1].state[i] >= self.knownPoles[i]:
				return True

		# if we get here, all is good.
		return False

	# test for success
	def testForSuccess(self, currentState):
		for i in range(len(self.partialSolution)):

			if self.partialSolution[i].state == currentState:
				# this is in a cyclic state
				solution = self.partialSolution[(i-len(self.partialSolution)):]
				
				# solution is a list of moveStatePairs, test the moves for a wildcard
				wildcardInSoln = False
				for j in solution:
					if j.move == "*":
						self.solutionWithWildCards = solution
						self.solvable = 1
						wildcardInSoln = True
						break

				return True
			i+=1
		# no cycles found - not a cyclic state
		return False

	def testForSuccessWithFullMinimumSolutionLength(self, currentState):
		
		# there's no reason to test a solution that's too short to possibly
		# be complete
		if len(self.partialSolution) >= self.minSolnLength:
			#i = self.minSolnLength - 1

			#maxValue = len(self.partialSolution) - self.minSolnLength 
			#while i < maxValue:
			for i in range(len(self.partialSolution) - 2*self.minSolnLength + 1):

				if self.partialSolution[i].state == currentState:
					# this is in a cyclic state
					solution = self.partialSolution[(i-len(self.partialSolution)):]
					
					# solution is a list of moveStatePairs, test the moves for a wildcard
					wildcardInSoln = False
					for j in solution:
						if j.move == "*":
							self.solutionWithWildCards = solution
							self.solvable = 1
							wildcardInSoln = True
							break

					return True

				i += 1
		# no cycles found - not a cyclic state
		return False

	def up(self):
		if len(self.partialSolution) > 1:
			poppedMove = self.partialSolution.pop(-1).move
		else:
			#print("this is unsolvable!")
			self.solvable = -1
			return

		# the wildcard is attempted before the numbers, so if
		# a wildcard is unplayable, try some numbers
		if poppedMove == "*":
			self.nextDownMove = 0
			self.moveDownNext = True
		# if one number doesn't work, try the next one
		elif poppedMove < (len(self.knownPoles) - 1):
			self.nextDownMove = poppedMove + 1
			self.moveDownNext = True
		# if all numbers have been tried, we have to keep climbing
		elif poppedMove >= (len(self.knownPoles) - 1):
			self.moveDownNext = False
		# if it's not a number or a *, it's invalid
		else:
			print("invalid input to up!", poppedMove)
			exit(0)

		#print("popped move: ", poppedMove)

	def solve(self):
		while self.solvable == 0:
			#self.printPartialSolution()

			if self.moveDownNext:
				#print("down")
				self.down()
			else:
				#print("up")
				self.up()


		if self.solvable == 1:
			print("input of", self.knownPoles, "is solvable, with a solution:")
			justMoves = []
			for i in self.solutionWithWildCards:
				justMoves.append(i.move)
			self.solutionWithWildCards = justMoves

			print(justMoves)
			return True
		elif self.solvable == -1:
			print("input of", self.knownPoles, "is tight.")
			return False

	def findMaxWildcardGap(self):
		# check that the function is not called out of sequence.
		if not self.solvable == 1:
			print("no solution, but findMaxWildcardGap was called")
			return
		
		# run throught the solution twice to include interactions cut by the start of the solution
		doubleSolution = self.solutionWithWildCards + self.solutionWithWildCards

		# use a start of 1 to count for testing before cutting
		timeSinceWildcard = 1
		maxTimeSinceWildcard = 0
		for i in doubleSolution:
			# if this is the wildcard, reset the time.
			if i == "*":
				timeSinceWildcard = 1
			# otherwise, increment it.
			else:
				timeSinceWildcard += 1

			if timeSinceWildcard > maxTimeSinceWildcard:
				maxTimeSinceWildcard = timeSinceWildcard

		return maxTimeSinceWildcard

#startTime = time.time()

#ourPWS = solver_wildcard_opt([3,3,6,12,24,48,96,200])
#ourPWS.solve()
#print(ourPWS.findMaxWildcardGap())

#solveTimeCost = time.time() - startTime
#print("time taken:", solveTimeCost)

