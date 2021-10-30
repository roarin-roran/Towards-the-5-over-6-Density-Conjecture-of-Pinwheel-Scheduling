# constructs the pareto and failure surface for a given number of tasks, using naive methods.




import pareto_surface_manager
import solver_foresight
import solver_wildcard_opt as PWS
import time 




# https://www.yworks.com/yed-live/?file=https://gist.githubusercontent.com/roarin-roran/46f0376f1bb57562306fbd12b6fb19e8/raw/0dde64c5c52e26c87847321d212af42bb812a8b3/trie%20exploration%20graph%201-4
class pareto:
	class node:
		def __init__(self, recurrenceVector, parent):
			self.recurrenceVector = recurrenceVector
			self.parent = parent

			self.satisfiable = None
			self.soln = None
			self.looseness = None
			self.looseSoln = None
			self.complete = None
			self.resolved = False

			self.children = []

		# returns the satisfiability of the recurrence vector of this node
		# if this is unknown, make a pinwheel solver and find it.
		def getSatisfiability(self):
			# if sat is not known
			if self.satisfiable == None:
				# make a pinwheel solver for this vector
				pinwheelSolver = solver_foresight.solver_foresight(self.recurrenceVector, False, False, False)

				# solve it
				PWSReturn = pinwheelSolver.solve()

				# if it's solvable, record the solution and that it has one
				if PWSReturn == 1:
					self.satisfiable = True
					self.soln = pinwheelSolver.soln
				# if it's unsolvable, record that
				elif PWSReturn == -1:
					self.satisfiable = False
					print("I'm unsolvable")
				# catch illegal returns
				else:
					print("pinwheelSolver returned a forbidden value:", PWSReturn)
					exit(1)

			return self.satisfiable

		# returns a solution if one exists or false if no solution can exist
		def getSoln(self):
			# if the solution has not yet been sought, find it
			if self.satisfiable == None:
				self.getSatisfiability()

			# return the solution if it exists, otherwise false
			if self.satisfiable:
				return self.soln
			elif self.satisfiable == False:
				return False
			# catch an illegal value
			else:
				print("illegal value of self.satisfiable after running getSatisfiability:", self.satisfiable)

		# returns true if the node is complete (not subject to further exploration because all solutions up to the
		# target number of tasks has been found).
		def getCompleteness(self):
			# to be complete, a node must be satisfiable
			if self.getSatisfiability():
				# if it's satisfiable and full length, a node is complete
				if len(self.recurrenceVector) == self.parent.numberOfTasks:
					self.complete = True
				# if it's less than full length, it can still be complete if it has exactly one child
				# and that child is complete
				elif len(self.recurrenceVector) < self.parent.numberOfTasks:
					# if it has exactly one child, the point may be complete
					if len(self.children) == 1:
						self.complete = self.children[0].getSatisfiability()
					# if it doesn't have exactly one child it's not complete
					else:
						self.complete = False

				# if a node is over length, someat is broken
				else:
					print("error: longer recurrenceVector than numberOfTasks")
					exit(1)
			# if the node is not satisfiable, it's not complete
			else:
				self.complete = False

			return self.complete

		# returns true iff a node is loose: solvable with wildcards
		def getLooseness(self):
			# if looseness is undefined, find it.
			if self.looseness == None:
				if self.satisfiable == None:
					self.getSatisfiability()

				# if this input is known to be unsolvable, it's tight
				if self.satisfiable == False:
					self.looseness = False
				elif self.satisfiable == True:
					pinwheelSolver = PWS.solver_wildcard_opt(self.recurrenceVector)

					self.looseness = pinwheelSolver.solve()
					
					if self.looseness:
						self.looseSoln = pinwheelSolver.solutionWithWildCards

				# catch illegal value.
				else:
					print("illegal value of self.satisfiable after running getSatisfiability:", self.satisfiable)

			return self.looseness

		# returns the loose solution if one can exist, else returns false
		def getLooseSoln(self):
			# if looseness is undefined, find it
			if self.looseness == None:
				self.getLooseness()

			if self.looseness:
				return self.looseSoln
			else:
				return False

		# creates the full subtree of this node, returns true iff this node is compete, else false
		def resolve(self):
			print("beginning of node:", self.recurrenceVector)

			# can only call this once
			if self.resolved:
				print("error, attempted to resolve the same node multiple times")
				exit(1)
			else:
				self.resolved = True

			# test satisfiability
			if not self.getSatisfiability():
				print("I'm unsatisfiable: ", self.getSatisfiability())

				newFSP = self.parent.SE.FSPoint(self.recurrenceVector)
				self.parent.SE.compareWithSurface(newFSP, False)
				return False

			# test length - if it's full length, it's complete in its own right
			if len(self.recurrenceVector) == self.parent.numberOfTasks:
				newPSP = self.parent.SE.PSPoint(self.getSoln(), self.recurrenceVector)
				self.parent.SE.compareWithSurface(newPSP, True)

				self.complete = True
				return True
			# if it's not full length its children will determine its completeness
			elif len(self.recurrenceVector) < self.parent.numberOfTasks:
				# if this node can have children, have children until one is complete
				if self.getLooseness():
					# node hasn't been resolved before - should have no children
					anyCompleteChildren = False
					newChar = self.recurrenceVector[-1]
					while not anyCompleteChildren:
						print("itterating loop: ", self.recurrenceVector, "+ *")
						#input()

						newRV = self.recurrenceVector.copy()
						newRV.append(newChar)
						newChar += 1

						print("I'm creating a child: ", newRV)

						self.children.append(self.parent.node(newRV, self.parent))

						# resolve the child, stopping if the child is complete
						anyCompleteChildren = self.children[-1].resolve()

						if anyCompleteChildren:
							print("done!!")
						
					if len(self.children) == 1:
						self.complete = True
						return True
					else:
						return False

				# if this node is barren, compare with the failure surface and return false
				else:
					newFSP = self.parent.SE.FSPoint(self.recurrenceVector)
					self.parent.SE.compareWithSurface(newFSP, False)
					return False

				
			# if a node is over length, someat is broken
			else:
				print("error: longer recurrenceVector than numberOfTasks")
				exit(1)


	def __init__(self, numberOfTasks):
		self.numberOfTasks = numberOfTasks

		self.SE = pareto_surface_manager.pareto_surface_manager(numberOfTasks)

	@staticmethod
	def isRecurrenceVectorChildComplete(self, recurrenceVector):
		if len(recurrenceVector) == self.numberOfTasks:
			return True
		elif len(recurrenceVector) < self.numberOfTasks:
			return False
		else:
			print("bad input to isRecurrenceVectorChildComplete - recurrenceVector too long")
			exit(1)

	# a series explorer designed to be easy to parralelise
	# creates all the depth 1 nodes and explores from them
	def explore(self):
		for i in range(self.numberOfTasks):
			newNode = self.node([i + 1], self)
			newNode.resolve()

	# tests all functions 
	@staticmethod
	def test():
		ourTE = pareto(4)

		# test get satisfiability, get soln, get looseness and get loose soln

		ourNode = ourTE.node([2,2,2], ourTE)
		if ourNode.getSatisfiability() != False:
			print("problem with nodes, failed test 1")
		if ourNode.getSoln() != False:
			print("problem with nodes, failed test 2")
		if ourNode.getLooseSoln() != False:
			print("problem with nodes, failed test 3")
		if ourNode.getLooseness() != False:
			print("problem with nodes, failed test 4")


		ourNode = ourTE.node([2,3], ourTE)
		if ourNode.getLooseSoln() != False:
			print("problem with nodes, failed test 5")
		if ourNode.getSatisfiability() != True:
			print("problem with nodes, failed test 6")
		if ourNode.getSoln() != [1,0]:
			print("problem with nodes, failed test 7")
		if ourNode.getLooseness() != False:
			print("problem with nodes, failed test 8")

		ourNode = ourTE.node([2,4,8], ourTE)
		if ourNode.getLooseness() != True:
			print("problem with nodes, failed test 9")
		if ourNode.getSoln() != [1, 0, 2, 0, 1, 0]:
			print("problem with nodes, failed test 10")
		if ourNode.getSatisfiability() != True:
			print("problem with nodes, failed test 11")
		if ourNode.getLooseSoln() != [1, 0, 2, 0, 1, 0, '*', 0]:
			print("problem with nodes, failed test 12")

		# test completeness
		if ourNode.getCompleteness() != False:
			print("problem with nodes, failed test 13")

		ourNode = ourTE.node([2,4,8,8], ourTE)
		if ourNode.getCompleteness() != True:
			print("problem with nodes, failed test 14")

		print("missing tests: completeness")

	def run(self):
		print("runnning for", self.numberOfTasks, "tasks")

		self.explore()
		print("Parato surface: ")
		self.SE.printParetoSurface()
		print()
		print("Failure surface: ")
		self.SE.printFailureSurface()

		print("all done boss")

#pareto.run(4)