# *********************************************
# * * * * * * * * * * * * * * * * * * * * * * *
# * * * * * * * * * * NAIVE * * * * * * * * * *
# * * * * * * * * * * * * * * * * * * * * * * *
# *********************************************

# * * * * * * * * * * DISCRETE INSTANCE SOLVING * * * * * * * * * *
# the solver_naive class is initiated with a single
# instance of pinwheel scheduling and explores it,
# recording a sollution if there is one

from collections import deque
import numpy as np
import sys
import time


class solver_naive:
    "An automated solution finder for one instance of pinwheel scheduling, which \
    uses and produces either a solution or (in testing mode) a witness table \
    (with testing mode disabled, this is generated but not saved)"

    class MoveStatePair:
        "a node on the exploration tree, consisting of the last move played and the \
        current state of the system"
        def __init__(self, m, s):
            self.move = m
            self.state = s

        def getMove(self):
            return self.move

        def getState(self):
            return self.state

    def __init__(self, wheels, testingMode, reccomendedPrints):
        self.wheels = wheels
        # wheels should be stored in descending order
        self.wheels.sort()
        self.wheels.reverse()

        self.testingMode = testingMode
        self.reccomendedPrints = reccomendedPrints

        # pole heights are all initially 0,
        # but over the course of the first day all grow to height 1
        gaps = []
        for i in range(len(wheels)):
            gaps.append(1)

        # the gardener cuts the first pole on the first day:
        move = 0
        gaps[0] = 0

        firstMSP = self.MoveStatePair(move, gaps)
        self.nodeStack = deque()
        self.nodeStack.append(firstMSP)

        # needed for exploration
        self.moveDownNext = True
        self.nextDownMove = 1
        # -1 for no, 1 for yes and 0 for we don't know yet
        self.solvable = 0
        self.soln = []

        # a list of prefixes which are unsolvable
        if testingMode:
            self.witnessList = []

        if testingMode:
            print("initiated successfully:", wheels)

        # cost, either by counting nodes or by counting comparisons in
        # testForSuccess. testForSuccess will dominate here, because
        # it grows as n^2 per node instead of the constant m per node
        self.nodeCost = 0
        self.failureTestCost = 0
        self.successTestCost = 0

    # **********************************
    # * * * * * * * * ** * * * * * * * *
    # * * * * * CORE FUNCTIONS * * * * *
    # * * * * * * * * ** * * * * * * * *
    # **********************************

    # popping (up), pushing (down) and testing for success and failure

    # manages the exploration of the tree with a loop and switch statement.
    # as each move is informed by the last, a recursive process was initially
    # used here - this should work better
    def solve(self):
        self.startTime = time.time()

        # while we don't know whether this pinwheel instance is solvable or not
        while(self.solvable == 0):

            # if we want to go down and there's somewhere to go, go down
            if self.moveDownNext and self.nextDownMove < len(self.wheels):
                self.down(self.nextDownMove)
            # otherwise, go up
            else:
                self.up()

        self.failureTestCost = len(self.wheels)*self.nodeCost

        self.solveTimeCost = time.time() - self.startTime

        if self.solvable == 1:
	        self.sortOutput()

        return self.solvable

    # sorts the output in descending order, including a redefinition of the problem and relabelling the solution
    def sortOutput(self):
        unsortedProblem = self.getSolnQuality()

        sortedProblem = unsortedProblem.copy()
        sortedProblem.sort()
        #sortedProblem.reverse()

        # map the unsorted problem onto the sorted problem, use this to fix the solution.
        mapping = []
        sortedProblemCopy = sortedProblem.copy()
        for i in unsortedProblem:
            for j in range(len(sortedProblemCopy)):
                if i == sortedProblemCopy[j]:
                    mapping.append(j)
                    sortedProblemCopy[j] = -1
                    break

        # use the mapping to relabel the solution
        for i in range(len(self.soln)):
            self.soln[i] = mapping[self.soln[i]]

        self.wheels = sortedProblem.copy()

    # find the hardest PWS problem that is solved by this solution.
    def getSolnQuality(self):
        if self.solvable != 1:
            print("unsolved!")
            exit(1)
        else:
            self.numPoles = len(self.wheels)

            # initialise
            self.solnQuality = [0]*self.numPoles
            gapSinceCut = [0]*self.numPoles

            # double soln to get wrap around
            for i in self.soln*2:
                # cut
                gapSinceCut[i] = 0

                for j in range(len(gapSinceCut)):
                    # grow
                    gapSinceCut[j] += 1
                    # increment max
                    if gapSinceCut[j] > self.solnQuality[j]:
                        self.solnQuality[j] = gapSinceCut[j]

            return self.solnQuality

    def getSolution(self):
        if self.solvable != 1:
            print("invalid input - unsolved or unsolvable system", self.solvable, self.dBGMaxSeperations)
            exit(0)
        else:
            return self.soln 

    # moves down the tree, adding a new node if testing is passed.
    def down(self, inputMove):
        # if testing mode is on, test that the input letter is valid before
        # proceeding.
        if self.testingMode:
            self.testInputLetter(inputMove)

        gaps = self.nodeStack[-1].getState()
        # increment all hole seperations
        gaps = [i+1 for i in gaps]


        # if one gap gets too big, this prefix is a bust -
        # abandon the last addition
        if self.testForFail():
            self.moveDownNext = False
        # if all gaps are small enough
        else:
            if(self.testingMode):
                print("down", inputMove, self.getAllMoves())

            # cut the inputted pole
            gaps[inputMove] = 0

            # add the inputted move and its consequences to the stack
            self.nodeStack.append(self.MoveStatePair(inputMove, gaps))

            # see whether this is a complete solution
            if(self.testForSuccess()):
                if self.testingMode or self.reccomendedPrints:
                    print("success! solution found:", self.soln)
            # if it's not, keep going down, starting with the first
            # letter again
            else:
                self.moveDownNext = True
                self.nextDownMove = 0


        self.nodeCost += 1

    # returns true if any gap is too large, false otherwise
    def testForFail(self):
        testFailed = False
        j = 0
        # for all gaps in the current state
        for i in self.nodeStack[-1].getState():
            # check it's not larger than its maximum size.
            testFailed = testFailed or (i >= self.wheels[j])
            j += 1

        return testFailed

    # tests whether the solution prefix produced so far is complete.
    # if it is, record it and use self.solvable to communicate completion.
    # if not, nothing happens.
    def testForSuccess(self):

        # tests 
        i = 0
        while i < len(self.nodeStack) - 1:
            #print(i, self.nodeStack[i].getState())

            if np.array_equal(self.nodeStack[i].getState(),
                              self.nodeStack[-1].getState()):
                self.recordSolln(i)
                self.solvable = 1
            i += 1

            
            self.successTestCost += len(self.wheels)

        return self.solvable

    # removes a day from the stack, informing
    def up(self):
        if(self.testingMode):
            print("up")
            self.witnessList.append(self.getAllMoves())
        if len(self.nodeStack) == 1:
            if self.testingMode or self.reccomendedPrints:
                print("no solution!")
            self.solvable = -1
        else:
            lastNode = self.nodeStack.pop()
            lastMove = lastNode.getMove()

            self.nextDownMove = lastMove + 1
            self.moveDownNext = True

    # ***************************************
    # * * * * * * * * * * * * * * * * * * * *
    # * * * * * ANCILIARY FUNCTIONS * * * * *
    # * * * * * * * * * * * * * * * * * * * *
    # ***************************************

    # assorted support functions

    def recordSolln(self, repetitionStart):
        soln = []
        i = repetitionStart
        while i < len(self.nodeStack)-1:
            soln.append(self.nodeStack[i].getMove())
            i += 1

        # note that this result is also tested by bgtp solver in normal
        # operation, which verifies the performance of everything below
        # it including this - this test is therefore testing mode only
        if self.testingMode:
            if solver_naive.testPinwheelSoln(self.wheels, soln, True):
                self.soln = soln
            else:
                print("no solution recorded,", soln,
                      "does not solve", self.wheels)
        else:
            self.soln = soln

    def getAllMoves(self):
        allMoves = []
        for i in self.nodeStack:
            allMoves.append(i.getMove())
        return allMoves

    def printAllMoves(self):
        for i in self.nodeStack:
            print(i.getMove())

    def printAllMoveStatePairs(self):
        for i in self.nodeStack:
            print(i.getMove(), i.getState())

    def getSoln(self):
        if self.solvable == 0:
            print("bug warning! tried to get nonexisting solution - unsolved")
            self.solve()
        elif self.solvable == -1:
            print("bug warning! tried to get nonexisting solution \
                  - unsolvable")
            return []
        return self.soln


    # ***************************
    # * * * * * * * * * * * * * *
    # * * * * * TESTING * * * * *
    # * * * * * * * * * * * * * *
    # ***************************

    # testing and verification functions

    # a silent pass test for the validity of the input letter.
    def testInputLetter(self, inputMove):
        if (not isinstance(inputMove, int)) \
                or (inputMove < 0) \
                or (inputMove >= len(self.wheels)):

            print("invalid input to down",
                  not isinstance(inputMove, int),
                  (inputMove < 0),
                  (inputMove >= len(self.wheels)))
            print(inputMove, self.wheels)
            exit(1)

    # as implied by reportSuccess, this will be a silent pass noisy fail
    # test by default
    @staticmethod
    def testPinwheelSoln(wheels, solution, reportSuccess):
        # the time since this wheel was last spiked
        gaps = [0]*len(wheels)

        passedAllTests = True

        # a valid solution must include all elements from wheels.
        for i in range(len(wheels)):
            if i not in solution:
                print("wheel", i, "not adressed!")
                passedAllTests = False

        if passedAllTests:
            # a valid solution must include only elements from wheels.
            for i in solution:
                if i not in range(len(wheels)):
                    print("there is no wheel", i, "!")
                    passedAllTests = False

        if passedAllTests:
            # a valid solution must not have any gap larger than the
            # urgency given in wheels
            gaps = [0] * len(wheels)
            for i in solution * 2:
                gaps = [i+1 for i in gaps]
                if gaps[i] > wheels[i]:
                    print("wheel", i, "has too large a gap")
                    passedAllTests = False
                gaps[i] = 0

        if passedAllTests and reportSuccess:
            print("soln tested and passed!")

        return passedAllTests

    @staticmethod
    def runTestCases():
        print(" - - - - - testing pinwheel solver- - - - - ")
        print("should succeed:")
        oursolver_naive = solver_naive([6, 4, 2], True, True)
        oursolver_naive.solve()
        print()

        print("should fail:")
        oursolver_naive = solver_naive([4, 3, 2], True, True)
        oursolver_naive.solve()
        # print(oursolver_naive.witnessList)
        print()

        print("should succeed:")
        oursolver_naive = solver_naive([1], True, True)
        oursolver_naive.solve()
        print()

        print("should succeed:")
        oursolver_naive = solver_naive([5, 10, 5, 10, 3], True, True)
        oursolver_naive.solve()
        print()

#oursolver_naive = solver_naive([2,2], False, True)
#oursolver_naive.runTestCases()