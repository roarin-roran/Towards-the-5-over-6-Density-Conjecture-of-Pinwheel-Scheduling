# *******************************************
# * * * * * * * * * * * * * * * * * * * * * *
# * * * * * * * * * * OPT * * * * * * * * * *
# * * * * * * * * * * * * * * * * * * * * * *
# *******************************************



from collections import deque
import numpy as np
import sys
import time

# **********************************************
# * * * * * * * * * * * ** * * * * * * * * * * *
# * * * * * * * * * * SOLVER * * * * * * * * * *
# * * * * * * * * * * * ** * * * * * * * * * * *
# **********************************************

class Solver:
    "a parent class for exact, greedy and cutfastest solvers: if it needs to \
    perform success testing, it extends this class - allowing multiple classes \
    to be modified simultaniously"

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

    # tests whether the solution prefix produced so far is complete.
    # if it is, record it and use self.solvable to communicate completion.
    # if not, nothing happens.
    def testForSuccess(self):
        # there's no reason to test a solution that's too short to possibly
        # be complete
        if len(self.nodeStack) >= self.minSolnLength:
            # the shortest way to cut all poles is a partial round robin
            i = self.minSolnLength - 1
            while i < (len(self.nodeStack) - self.minSolnLength):
                if np.array_equal(self.nodeStack[i].getState(),
                                  self.nodeStack[-1].getState()):
                    self.recordSolln(i)
                    self.solvable = 1
                i += 1

        self.successTestCost += self.numPoles

        return self.solvable



    # * * * * * * * * * * * * * * * * * * * *
    # * * * * * ANCILIARY FUNCTIONS * * * * *
    # * * * * * * * * * * * * * * * * * * * *

    def recordSolln(self, repetitionStart):
        soln = []
        i = repetitionStart
        while i < len(self.nodeStack)-1:
            soln.append(self.nodeStack[i].getMove())
            i += 1

        if len(soln) < self.minSolnLength:
            print("WARNING: self reported solution length is wrong")
            sys.exit()
        else:
            self.soln = soln

    def getAllMoves(self):
        allMoves = []
        for i in self.nodeStack:
            allMoves.append(i.getMove())
        return allMoves

    def printAllMoves(self):
        for i in self.nodeStack:
            print(i.getMove(), end ='')
        print()

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


# ******************************************************
# * * * * * * * * * * * * * ** * * * * * * * * * * * * *
# * * * * * * * * * * DBGT EXPLORING * * * * * * * * * *
# * * * * * * * * * * * * * ** * * * * * * * * * * * * *
# ******************************************************

class solver_opt(Solver):
    "functions shared by the core PinwheelSolver class, as well as the less \
    important SchedulableDensityMeasurerer class, which both explore discretized \
    bamboo gardens. the core methods of this class will be popping and pushing \
    from the move stack and testing for constraint failure - though optimisations \
    and some small methods will also appear here."
    
    # * * * * * *  * * * * * *
    # * * * * * CORE * * * * *
    # * * * * * *  * * * * * *

    # SDM builds on the init function used by pinwheel
    def sharedInitComponents(self, wheels, testingMode, reccomendedPrints):
        self.wheels = wheels
        # wheels should be stored in descending order
        self.wheels.sort()
        self.wheels.reverse()

        self.numPoles = len(self.wheels)

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

        # setup the node stack
        firstMSP = self.MoveStatePair(move, gaps)
        self.nodeStack = deque()
        self.nodeStack.append(firstMSP)

        # needed for exploration
        self.moveDownNext = True
        self.nextDownMove = 1
        # -1 for no, 1 for yes and 0 for we don't know yet
        self.solvable = 0
        self.soln = []

        # an array of the elements of the wheels array which are identical
        # to an earlier element. these should not be played before that element
        # to avoid testing symetric situations more than once
        self.blockedBySymmetry = []
        self.anyBlockedBySymmetry = False
        self.constructSymetricBlockArray()

        # a list of prefixes which are unsolvable
        if testingMode:
            self.witnessList = []

        if testingMode:
            print("initiated successfully:", wheels)

        # cost, by counting comparisons in success or failure testing.
        self.successTestCost = 0
        self.failureTestCost = 0

        # set the minimum solution length, which is used to avoid unnecesary
        # success testing
        self.setNaiveMinSolnLength()

        self.nodeCost = 0

    # moves down the tree, adding a new node if testing is passed.
    def down(self, inputMove, defaultMove):
        self.gaps = self.nodeStack[-1].getState()
        # increment all hole seperations
        self.gaps = [i+1 for i in self.gaps]

        # if one gap gets too big, this prefix is a bust -
        # abandon the last addition
        if self.testForFail():
            self.moveDownNext = False
        # if all gaps are small enough
        else:
            if(self.testingMode):
                print("down", inputMove, self.getAllMoves())

            # cut the inputted pole
            self.cutPole(inputMove)

            # add the inputted move and its consequences to the stack
            self.nodeStack.append(self.MoveStatePair(inputMove, self.gaps))

            # see whether this is a complete solution
            if(self.testForSuccess()):
                self.handleSuccess()
                if self.testingMode or self.reccomendedPrints:
                    print("solution found:", self.soln)
            # if it's not, keep going down, starting with the first
            # letter again
            else:
                self.moveDownNext = True
                self.nextDownMove = defaultMove

        self.nodeCost += 1

    # returns true if any gap is too large, false otherwise
    def testForConstraintViolation(self):
        testFailed = False

        lastMove = self.nodeStack[-1].getMove()

        if self.isThisBlockedBySymetry(lastMove) and False:
            testFailed = True

        else:
            # for all gaps in the current state, where j is the assiciated limit
            j = 0
            for i in self.nodeStack[-1].getState():
                # check it's not larger than its maximum size.
                if i >= self.wheels[j]:
                    testFailed = True
                    break
                else:
                    j += 1
        
        self.failureTestCost += len(self.wheels)

        return testFailed

    # removes a day from the stack and tests for a return to the root. if
    # the root has been reached: input is unsolvable. otherwise, will choose
    # the next move.
    def up(self):
        lastMove = 'n'
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

            self.incrementNextDownMove(lastMove)
            self.moveDownNext = True

        return lastMove


    # * * * * * * * * * * * * * * * * *
    # * * * * * OPTIMIZATIONS * * * * *
    # * * * * * * * * * * * * * * * * *

    # functions which support optimizations rather than core functionality

    # increments nextDownMove if a repeated move is attempted
    def preventRepeats(self):
        # repeating the last move is never helpful if there are multiple poles
        # if there are multiple poles, prevent a repeated move.
        lastMove = self.nodeStack[-1].getMove()
        repeatedLastMove = self.nextDownMove == lastMove
        if repeatedLastMove and len(self.wheels) > 1 and lastMove != 'p':
            self.nextDownMove += 1

    # returns true if the chosen element is not blocked by symetry, either
    # because it cannot be, or because the element that could block it has
    # been played
    def isThisBlockedBySymetry(self, inputMove):
        # by default, things are not blocked by symmetry
        blocked = False
        # if the inputMove is in the blocked by symmetry array, the preceeding
        # element should be checked
        if inputMove in self.blockedBySymmetry:
            # if the preceeding element (which must be the blocking
            # element) is untouched, kill this branch

            if self.nodeStack[-1].getState()[inputMove - 1] \
                    == len(self.nodeStack):
                self.moveDownNext = False
                blocked = True

        return blocked

    # constructs an array of all elements which should be blocked by
    # symmetry until the preceeding element is used for the first time
    def constructSymetricBlockArray(self):
        # for all wheels
        i = 0
        while i < len(self.wheels):
            # examin all later wheels
            j = i
            while j < len(self.wheels) - 1:
                # increment first, to prevent any element analysing itself.
                j += 1
                # if a duplicate is found
                if self.wheels[j] == self.wheels[i]:
                    # prevent it from being examined again
                    i = j
                    self.blockedBySymmetry.append(j)
            # go on to the next character
            i += 1

        if len(self.blockedBySymmetry) > 0:
            self.anyBlockedBySymmetry = True

    # sets the minimum solution length - no point in testing any solution
    # shorter than this
    def setNaiveMinSolnLength(self):
        # first check the density, because any pinwheel case that would
        # fail on density grounds has an infinite minimum solution length
        naiveDensity = 0
        for i in self.wheels:
            naiveDensity += 1/i

        if naiveDensity > 1:
            self.solvable = -1
        else:
            # the length of the round robin solution.
            self.minSolnLength = len(self.wheels)
            minNumber = [1]*self.minSolnLength


            actedThisItteration = True
            while actedThisItteration:
                actedThisItteration = False
                for i in range(len(minNumber)):
                    while (minNumber[i]/self.minSolnLength) < \
                            (1/self.wheels[i]):
                        minNumber[i] += 1
                        self.minSolnLength += 1

            #print("minSolnLength = ", self.minSolnLength)


# *****************************************************************
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# * * * * * * * * * * DISCRETE INSTANCE SOLVING * * * * * * * * * *
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# *****************************************************************

class PinwheelSolver(solver_opt):
    "An automated solution finder for one instance of pinwheel scheduling, which \
    uses and produces either a solution or (in testing mode) a witness table \
    (with testing mode disabled, this is generated but not saved)"

    # SDM uses an expanded initiation fn
    def __init__(self, wheels, testingMode, reccomendedPrints):
        self.sharedInitComponents(wheels, testingMode, reccomendedPrints)


    # * * * * * * * * ** * * * * * * * *
    # * * * * * CORE FUNCTIONS * * * * *
    # * * * * * * * * ** * * * * * * * *

    # manages the exploration of the tree with a loop and switch statement.
    # as each move is informed by the last, a recursive process was initially
    # used here - this should work better
    def solve(self):
        self.startTime = time.time()

        # while we don't know whether this pinwheel instance is solvable or not
        while(self.solvable == 0):
            self.preventRepeats()

            # if we want to go down and there's somewhere to go, go down
            if self.moveDownNext and self.nextDownMove < len(self.wheels):
                self.down(self.nextDownMove, 0)
            # otherwise, go up
            else:
                self.up()

        self.solveTimeCost = time.time() - self.startTime
        
        if self.solvable == 1:
            self.sortOutput()

        return self.solvable

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



    # in this method, failure is only possible via constraint vioaltion
    def testForFail(self):
        return self.testForConstraintViolation()

    # in pinwheelSolver nothing needs to be done here - the function is called
    # for use in SchedulableDensityMeasurerer
    def handleSuccess(self):
        pass

    # in sdm a pass is also possible, so a function is used here
    def cutPole(self, poleToCut):
        self.gaps[poleToCut] = 0

    # in sdm a pass is also possible, so a function is used here
    def incrementNextDownMove(self, lastMove):
        self.nextDownMove = lastMove + 1


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


    # * * * * * * * * * * * * * *
    # * * * * * TESTING * * * * *
    # * * * * * * * * * * * * * *

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
        ourPinwheelSolver = PinwheelSolver([6, 4, 2], True, True)
        ourPinwheelSolver.solve()
        print()

        print("should fail:")
        ourPinwheelSolver = PinwheelSolver([4, 3, 2], True, True)
        ourPinwheelSolver.solve()
        # print(ourPinwheelSolver.witnessList)
        print()

        print("should succeed:")
        ourPinwheelSolver = PinwheelSolver([1], True, True)
        ourPinwheelSolver.solve()
        print()

        print("should succeed:")
        ourPinwheelSolver = PinwheelSolver([5, 10, 5, 10, 3], True, True)
        ourPinwheelSolver.solve()
        print()

        print("thing we want:")
        ourPinwheelSolver = PinwheelSolver([3, 4, 4, 14, 14], True, True)
        ourPinwheelSolver.solve()
        print("success: ", ourPinwheelSolver.successTestCost)
        print("failure: ", ourPinwheelSolver.failureTestCost)
        print("nodes:   ", ourPinwheelSolver.nodeCost)

    # this will, by default, be empty - it exists to be populated with 
    # cases that relate to a particular bug/suspected bug
    @staticmethod
    def runBugFixingCases():
        print(" - - - - - bug fixing PinwheelSolver- - - - - ")
        print("problematic input itself:")
        ourPinwheelSolver = PinwheelSolver([14, 14, 11, 11, 11, 10, 10, 9, 5], False, True)
        ourPinwheelSolver.solve()
        print()

        print("potential root 1:")
        ourPinwheelSolver = PinwheelSolver([11, 11, 11, 10, 10, 7, 9, 5], False, True)
        ourPinwheelSolver.solve()
        print()

        print("potential root 2:")
        ourPinwheelSolver = PinwheelSolver([14, 14, 11, 11, 11, 9, 5, 5], False, True)
        ourPinwheelSolver.solve()
        print()


#ourPinwheelSolver = PinwheelSolver([3, 4, 4, 14, 14], True, True)
#ourPinwheelSolver.runBugFixingCases()




# ***********************************************************
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# * * * * * * * * * * DENSITY MEASUREMENT * * * * * * * * * * 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# ***********************************************************

class SchedulableDensityMeasurerer(solver_opt):
    def __init__(self, wheels, testingMode, reccomendedPrints):
        self.sharedInitComponents(wheels, testingMode, reccomendedPrints)

        # the maximum solution length
        self.maxLength = 1
        for i in wheels:
            self.maxLength *= i

        # maximum possible density is 1, if density must be greater than this then the subset is unschedulable
        self.minObservedDensity = 1

        # at the limit, every day could be a cut. - this is in the solution not the sequence.
        self.maxNumberOfCuts = self.maxLength

        # one cut is made in the initiation.
        self.numberOfCuts = 1

        self.naiveDensity = 0
        for i in wheels:
            self.naiveDensity += 1/i

    def solve(self):
        # -1: finished, 0: unscheduled, 1: scheduled 
        # (will repeat until the best solution is found)
        while self.solvable != -1:
            self.preventRepeats()

            moveUpNext = not self.moveDownNext
            # if we want to go down and there's somewhere to go, go down
            if self.moveDownNext:
                if self.nextDownMove == 'p':
                    self.down(self.nextDownMove, 'p')
                elif self.nextDownMove < len(self.wheels):
                    self.down(self.nextDownMove, 'p')

                    # if we've made a cut, increment numberOfCuts
                    self.numberOfCuts += 1
                # if we've tried passing and cutting each pole here to no result, go up
                else:
                    moveUpNext = True
            
            if moveUpNext:
                lastMove = self.up()
                if lastMove != 'p':
                    self.numberOfCuts -= 1

        return self.solvable

    def testForFail(self):
        if self.numberOfCuts > self.maxNumberOfCuts:
            print("IT TOO BIG!!!!!!!!!!!!!!!!!!!")
            exit(1)
        return self.testForConstraintViolation()

    def handleSuccess(self):
        self.solvable = 0
        
        numberOfCuts = len(self.soln)-self.soln.count('p')

        density = numberOfCuts/len(self.soln)
        if density < self.minObservedDensity:
            self.minObservedDensity = density
        print("density =", density)
        print("min observed density =", self.minObservedDensity)
        print("soln: ", self.soln)

        numberOfCutsRepeat = numberOfCuts
        lengthOfRepeatedSoln = len(self.soln)
        while lengthOfRepeatedSoln < self.maxLength:
            numberOfCutsRepeat += numberOfCuts
            lengthOfRepeatedSoln += len(self.soln)

        if numberOfCutsRepeat < self.maxNumberOfCuts:
            self.maxNumberOfCuts = numberOfCutsRepeat

        print("max number of cuts =", self.maxNumberOfCuts)
        print("stack length", len(self.nodeStack))

        print(self.naiveDensity)


    def cutPole(self, poleToCut):
        if poleToCut == 'p':
            pass
        else:
            self.gaps[poleToCut] = 0


    # this is called by up, use it to track number of cuts made
    def incrementNextDownMove(self, lastMove):
        if lastMove == 'p':
            self.nextDownMove = 0
        else:
            self.nextDownMove = lastMove + 1


    def runTestCases():
        
        ourSDM = SchedulableDensityMeasurerer([19,4], False, False)
        ourSDM.solve()
        print("calculated density: ", ourSDM.minObservedDensity)
        
        
        print("naive density: ", ourSDM.naiveDensity)

# ******************************************************
# * * * * * * * * * * * * * ** * * * * * * * * * * * * *
# * * * * * * * * * * APPROXIMATIONS * * * * * * * * * *
# * * * * * * * * * * * * * ** * * * * * * * * * * * * *  
# ******************************************************

# the greedy and cutFastest classes solve the BGT 
# approximately
class Approximation(Solver):
    # init function shared by greedy and cutfastest
    def __init__(self, poleGrowthRates, testingMode):
        # record pole growth rates and sort them in descending order
        self.poleGrowthRates = poleGrowthRates
        self.poleGrowthRates.sort()

        self.numPoles = len(self.poleGrowthRates)

        self.testingMode = testingMode
        if testingMode:
            print("MISSING FEATURE: testing modes for approximations")

        # all poles must have been cut before theres any point in 
        # testing for repetitions. 
        self.minSolnLength = len(self.poleGrowthRates)

        # the heights of all poles, initially 0
        self.poleHeights = []
        for i in range(len(self.poleGrowthRates)):
            self.poleHeights.append(0)

        # each node is a day, with a move state pair recorded in it.
        self.nodeStack = deque()

        # counts comparisons used in success testing, for cost estimation
        self.successTestCost = 0

        # while every possible case can be solved with greedy 
        # (to some quality), this is 0/1 because the same success testing 
        # is used by the DBGT solver, which also needs to be able to say -1
        # when a case is unsolvable
        self.solvable = 0
        self.soln = []

    def solve(self):
        # run until a solution is found
        while self.solvable == 0:
            self.passADay()
            self.testForSuccess()

# approximates by cutting the tallest pole
class GreedySolver(Approximation):
    def passADay(self):
        # let all poles grow.
        for i in range(len(self.poleHeights)):
            self.poleHeights[i] += self.poleGrowthRates[i]

        # cut the tallest pole
        indexOfMax = self.poleHeights.index(max(self.poleHeights))
        self.poleHeights[indexOfMax] = 0

        # record the move state pair
        self.nodeStack.append(self.MoveStatePair(indexOfMax, self.poleHeights.copy()))

    @staticmethod
    def runTestCases():
        print("'Ere we go lads")
        ourGreedySovler = GreedySolver([1, 2, 3], False)
        print(ourGreedySovler.poleGrowthRates)
        print(ourGreedySovler.poleHeights)
        
        ourGreedySovler.solve()
        
        ourGreedySovler.printAllMoveStatePairs()
        print(ourGreedySovler.getSoln())

        print("done son")

# approximates by cutting the fastest growing pole above a certain height
class CutFastestSolver(Approximation):
    def solveWithThreshold(self, threshold):
        dailyGrowth = sum(self.poleGrowthRates)

        self.cutOver = threshold*dailyGrowth
        self.solve()

    def passADay(self):
        # let all poles grow.
        for i in range(len(self.poleHeights)):
            self.poleHeights[i] += self.poleGrowthRates[i]

        # if no pole is over the line, cut nothing
        cutMade = "p"

        # cut the fastest growing pole taller than the threshold*H, if any poles are
        #print(self.poleHeights)
        for i in reversed(range(len(self.poleHeights))):
            if self.poleHeights[i] >= self.cutOver:
                self.poleHeights[i] = 0
                cutMade = i
                break

        # record the move state pair
        self.nodeStack.append(self.MoveStatePair(cutMade, self.poleHeights.copy()))

        # test for infinite repetition
        self.testForFailure()

    def testForFailure(self):
        if max(self.poleHeights) > 100*self.cutOver:
            print("WARNING: threshold unschedulable for instance:", self.poleGrowthRates)
            exit()

    @staticmethod
    def runTestCases():
        ourCutFastestSovler = CutFastestSolver([1,2,3], False)
        ourCutFastestSovler.solveWithThreshold(1)
        print(ourCutFastestSovler.getSoln())
        ourCutFastestSovler.printAllMoveStatePairs()
        print(ourCutFastestSovler.successTestCost)