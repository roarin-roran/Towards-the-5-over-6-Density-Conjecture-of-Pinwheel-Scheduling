

class solver_wildcard_naive:
	class moveStatePair:
		def __init__(self, move, state):
			self.move = move
			self.state = state

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

		self.partialSolution.append(self.moveStatePair("*", firstState))

		# * * * exploration variables * * *
		self.moveDownNext = True
		self.nextDownMove = 0

		# record all solutions that don't work for 
		self.tightSolutions = []

		self.solvable = 0

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

		# having incremented all poles, test it for failure
		if self.testForFailure():
			self.moveDownNext = False
		else:
			# having tested the state for failure, cut the input
			if self.nextDownMove != "*":
				currentState[self.nextDownMove] = 0

			# test for success - if this finds a wildcarded
			# solution we'll exit, otherwise we'll just backtrack.
			if self.testForSuccess(currentState):
				self.moveDownNext = False
			# if no solution is found, append the state and keep moving.
			else:
				self.partialSolution.append(self.moveStatePair(self.nextDownMove, currentState))
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

				#if not wildcardInSoln:
					#self.tightSolutions.append(solution)
						



				#print("current state: ", currentState)
				#for j in solution:
				#	print(j.move)

				return True
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
		i = 1
		while self.solvable == 0:
			#print()
			#print("day", i)
			i += 1
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


ourPWS = solver_wildcard_naive([3,3,9,9,12])
ourPWS.solve()
print(ourPWS.findMaxWildcardGap())

