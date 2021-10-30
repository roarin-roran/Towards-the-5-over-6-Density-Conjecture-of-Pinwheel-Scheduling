import networkx as nx
import time

class solver_graph:
	def __init__(self, pWSP):
		self.pWSP = tuple(pWSP)

		self.DG = nx.DiGraph()

		self.allStatesVisited = []

		initialState = tuple([0] * len(pWSP))

		#self.recursivelyExploreAllChildren(initialState)

		# graph:
		#j = 1
		#for i in list(self.DG.nodes):
		#	print("state", j, ":\t", i)
		#	j += 1


	# if this state is unseen, explore all of its children, else do nothing.
	def recursivelyExploreAllChildren(self, state):
		if state not in self.allStatesVisited:
			# don't revisit this state
			self.allStatesVisited.append(state)

			# need a list to modify and a tuple to put into the graph
			stateList = list(state)

			# grow
			for i in range(len(stateList)):
				stateList[i] += 1

			# failure testing
			for i in range(len(stateList)):
				if stateList[i] > self.pWSP[i]:
					return

			# branch and recursively continue, creating edges.
			for i in range(len(stateList)):
				newState = stateList.copy()
				newState[i] = 0
				newState = tuple(newState)

				self.DG.add_edge(state, newState)
				self.recursivelyExploreAllChildren(newState)

	def nonRecursivelyExploreAllChildren(self):
		allUnexploredStates = []
		allExploredStates = []

		initialState = [0] * len(self.pWSP)
		allUnexploredStates.append(initialState)

		while len(allUnexploredStates) > 0:
			inputState = allUnexploredStates[0].copy()
			currentState = inputState.copy()
			allUnexploredStates.pop(0)

			if tuple(currentState) not in allExploredStates:
				# don't revisit the state
				allExploredStates.append(tuple(currentState))

				# grow
				for i in range(len(currentState)):
					currentState[i] += 1

				skip = False
				# failure testing
				for i in range(len(currentState)):
					if currentState[i] > self.pWSP[i]:
						#print("jaammamamamamaam")
						skip = True

				if not skip:
					# branch and psuedo-recursively continue, creating edges.
					for i in range(len(currentState)):
						newState = currentState.copy()
						newState[i] = 0

						self.DG.add_edge(tuple(inputState), tuple(newState))
						allUnexploredStates.append(newState)




	def solve(self):
		startTime = time.time()

		self.nonRecursivelyExploreAllChildren()

		if len(max(nx.algorithms.components.strongly_connected_components(self.DG), key=len)) > 1:
			print(self.pWSP, " is strongly connected! - solvable")
		else:
			print(self.pWSP, " is not strongly connected - unsolvable") 
		
		self.solveTimeCost = time.time() - startTime		


#ourPGS = solver_graph([2,2])
#ourPGS.solve()
#print(ourPGS.solveTimeCost)