# *************************************************
# * * * * * * * * * * * * * * * * * * * * * * * * *
# * * * * * * * * * * FORESIGHT * * * * * * * * * *
# * * * * * * * * * * * * * * * * * * * * * * * * *
# *************************************************


from collections import deque
import copy
import numpy as np
import time
from fractions import Fraction

# ***************************************************
# * * * * * * * * * * * * * * * * * * * * * * * * * *
# * * * * * * * * * * DBGT SOLVER * * * * * * * * * *
# * * * * * * * * * * * * * * * * * * * * * * * * * *
# ***************************************************

class solver_foresight:
    "Solver for an instance of the DBGT (decision bamboo garden trimming ",
    "problem) which is equivelant to an instance of pinwheel scheduling. uses ",
    "the foresight method, which explores a tree of urgency states. urgency ",
    "states consider the number of days before a pole will exceed its maximum ",
    "height"

    class MoveStatePair:
        "a node on the exploration tree, consisting of the last move played \
        and the current state of the system"

        # creates a movestate pair, with both relative and absolute recorded
        # for later retrieval
        def __init__(self, rm, am, s):
            if rm.isRelativeMove():
                self.relativeMove = rm
            else:
                print("moveStatePair created wrong: \
                       entered an absolute move as a relative move")
                exit(0)

            # relative and absolute moves both have only the isRelativeMove
            # test: true is for relative, false for absolute and noisy fails
            # for something lacking this method
            if not am.isRelativeMove():
                self.absoluteMove = am
            else:
                print("moveStatePair created wrong: \
                       entered a relative move as an absolute move")
                exit(0)

            self.state = s

        def getRelativeMove(self):
            return self.relativeMove.getRelativeMove()

        def getAbsoluteMove(self):
            return self.absoluteMove.getAbsoluteMove()

        def getState(self):
            return self.state

        # added when the wildcardSovler extension was built: should only be asked by
        # methods from that that class
        def isWildcard(self):
            return False

    class State:
        "an object which stores the urgency state of all the poles in a \
        certain garden, and the forcing number, if this state has any \
        forcing."

        class Pole:
            "a pole in a DBG, with urgency, identity and the ability to cut \
            and grow itself"
            def __init__(self, initialIndex, dBGMaxSeperations):
                # these poles will be repeatedly reordered - this is the
                # initial index, and will not change through this reordering.
                # it identifies this particular pole
                self.initialIndex = initialIndex

                # the urgency when this pole has just been cut, equal to the
                # height of its associated DBGT pole, or the maximum
                # seperation of the associated pinwheel
                self.justCutUrgency = dBGMaxSeperations

                # initial urgency is the just cut value.
                self.urgency = self.justCutUrgency

            # cuts this pole: in an urgency world, this is not as simple
            # as setting the height to zero, you instead have to set the
            # urgency to its maximum value
            def cut(self):
                self.urgency = self.justCutUrgency

            # lets this pole grow for one day, decreasing the number of days
            # before it must be cut
            def grow(self):
                self.urgency -= 1

        # constructs the initial state object - further states will be copies
        # of this state
        def __init__(self, dBGMaxSeperations, testingMode, solver):
            self.testingMode = testingMode
            self.numPoles = len(dBGMaxSeperations)

            if solver.stateCreated:
                print("multiple states created! states should be duplicated \
                      only by copying")
                exit(1)
            else:
                solver.stateCreated = True

            # construct poles array - this is the current system state
            self.poles = []
            for i in range(self.numPoles):
                self.poles.append(self.Pole(i, dBGMaxSeperations[i]))

            # you have to cut one of the poles in the array
            self.forced = len(self.poles)

            # states should only be created by the init function once -
            # on day 1. all later states are copies of their predecessor.
            # thus, this function is called once and the bbs array is
            # reduced in size as time progresses
            self.setBlockedBySymmetryArray()



        # lets all poles grow
        def passADay(self):
            for p in self.poles:
                p.grow()

        # cuts a pole, retaining the ordering of the poles. uses relative pole
        # identity ie position in an array ordered by the dynamic urgency
        # rather than the static just cut urgency
        def cutPole(self, poleToCut):
            # cuts the given pole
            self.poles[poleToCut.getRelativeMove()].cut()

            # pop the cut pole out of the poles array
            justCutPole = self.poles.pop(poleToCut.getRelativeMove())

            # put it back into the poles array, in the correct place
            putItBack = False

            for i in range(self.numPoles - 1):
                if self.poles[i].urgency > justCutPole.urgency:
                    self.poles.insert(i, justCutPole)
                    putItBack = True
                    break

            # if it doesn't go before any element, it goes after all elements
            if not putItBack:
                self.poles.append(justCutPole)

            # remove the just cut pole's neighbour from the list - if it has
            # one
            try:
                self.blockedBySymmetry.remove(justCutPole.initialIndex + 1)
            except ValueError:
                pass

        # prints the urgencies of all poles, in the order of their urgencies
        # (which change)
        def printAllUrgencies(self):
            print("all (pole, urgency) pairs:\t ", end='')

            firstMove = True
            for p in self.poles:
                if firstMove:
                    firstMove = False
                else:
                    print(",", end='')

                print("(", p.initialIndex, ",", p.urgency, ")", end='')
            print()

        # returns all initial index/urgency pairs as a list of tuples
        def getAllUrgencies(self):
            allUrgencies = []
            for p in self.poles:
                allUrgencies.append((p.initialIndex, p.urgency))

            return allUrgencies

        # prints the actual heights of all poles, in order of their index
        def printAbsoluteState(self):
            print("all discretized pole heights:\t (", end='')

            currentIndex = 0
            firstMove = True
            while currentIndex < len(self.poles):
                for p in self.poles:
                    if p.initialIndex == currentIndex:
                        if firstMove:
                            firstMove = False
                        else:
                            print(", ", end='')
                        print(p.justCutUrgency - p.urgency, end='')
                        currentIndex += 1
            print(")")

        # returns true if the test is passed and false if it is failed
        def testForConstraintViolation(self):
            # default value, indicating a pass.
            constraintTestResult = True

            for i in range(self.numPoles):
                # if any pole is too tall, this state fails
                if i > self.poles[i].urgency:
                    constraintTestResult = False
                    break

            # in testing mode, print some information as to what's going on.
            if self.testingMode:
                print()
                if constraintTestResult is True:
                    print("I just passed")
                elif constraintTestResult is False:
                    print("I just failed")
                else:
                    print("there's a problem: constraintTestResult should be ",
                          "true or false, but isn't")
                self.printAllUrgencies()
                print()

            return constraintTestResult

        def testForForcing(self):
            for i in range(self.numPoles):
                # if any pole is forced, record it. note that this will let
                # violating states through but, with failure testing costs
                # vastly dominating success testing costs, this aggressive
                # approach will reduce the cost of failure tests per node
                if i == self.poles[i].urgency:
                    if i < self.forced:
                        self.forced = i
                        break

        # constructs an array of absolute poles which are blocked by symmetry
        def setBlockedBySymmetryArray(self):
            bBS = []

            for i in range(len(self.poles)):
                if i != 0:
                    if self.poles[i-1].justCutUrgency == \
                           self.poles[i].justCutUrgency:
                        bBS.append(i)

            self.blockedBySymmetry = bBS

    class RelativeMove:
        "a move, in the relative world where poles are ordered by their \
        current urgency, this class is used to make absolute and relative \
        moves fail when cross called"

        def __init__(self, move):
            self.move = move

        def getRelativeMove(self):
            return self.move

        def getAbsoluteMove(self):
            print("a relative move was cross called!")
            exit(1)

        @staticmethod
        def isRelativeMove():
            return True

    class AbsoluteMove:
        "a move, in the absolute world where poles are ordered by their \
        origional urgency (aka their growth rate). this class is used to \
        make absolute and relative moves fail when cross called"

        def __init__(self, move):
            self.move = move

        def getRelativeMove(self):
            print("an absolute move was cross called!")
            exit(1)

        def getAbsoluteMove(self):
            return self.move

        @staticmethod
        def isRelativeMove():
            return False

    # initiates the DGBTSolver
    def __init__(self, dBGMaxSeperations, testingMode = False, reccomendedPrints = False, displayMode = False):

        # * * * * * * * * * * INPUT PROCESSING * * * * * * * * * *

        self.reccomendedPrints = reccomendedPrints
        self.testingMode = testingMode
        # display mode is reused for testing.
        self.displayMode = displayMode or testingMode

        # take the inputted max seperations (maximum days of growth per pole)
        # and sort them
        self.dBGMaxSeperations = dBGMaxSeperations
        self.dBGMaxSeperations.sort()

        self.numPoles = len(self.dBGMaxSeperations)

        # * * * * * * * * * * SETUP * * * * * * * * * *

        # states should be duplicated by copying, not by creating new states
        # this variable prevents the state init function from being called
        # twice by the same solver
        self.stateCreated = False

        # the initial state, with every pole at a height of 0 and urgency of
        # u_0
        initialState = self.State(self.dBGMaxSeperations, self.testingMode,
                                  self)

        # setup the node stack
        self.nodeStack = deque()

        # -1 for no, 1 for yes and 0 for we don't know yet
        self.solvable = 0
        self.soln = []

        # naive here means that the density cost of combining poles is
        # not considered
        self.setNaiveMinSolnLength()

        # initialise various cost measures
        self.failureTestCost = 0
        self.successTestCost = 0
        self.nodesVisited = 0
        self.solveTimeCost = 0

        # * * * * * * * * * * DAY 0 AND 1 * * * * * * * * * *

        if self.displayMode:
            print("Discretized poles:\t\t", self.dBGMaxSeperations)

            print()
            print("day 0")
            initialState.printAbsoluteState()
            initialState.printAllUrgencies()
            print()
            print()

        self.firstDay(initialState)

        # needed for exploration
        self.moveDownNext = True
        self.nextDownMove = self.RelativeMove(0)

    # takes and records the first day
    def firstDay(self, initialState):
        firstRelativeMove = self.RelativeMove(0)
        firstAbsoluteMove = self.AbsoluteMove(0)

        # the firt recorded state grows all poles and cuts pole 0
        firstRecordedState = initialState
        firstRecordedState.passADay()
        firstRecordedState.cutPole(firstRelativeMove)

        currentState = firstRecordedState

        if(self.displayMode):
            # the only difference between display and testing modes are
            # whether user input is required for advancement
            if not self.testingMode:
                input()

            print("day 1")
            print("action:\t\t\t\t cut pole 0 in position 0")
            print("reason:\t\t\t\t first cut is a free choice, always choose ", 
                  "to cut 0")
            print("outcome: \t\t\t move successful, continue")
            print()
            initialState.printAbsoluteState()
            firstRecordedState.printAllUrgencies()
            print()
            print("all abs moves:\t\t\t", [0])
            print("all rel moves:\t\t\t", [0])
            print()
            print()

        # record the first day to the first position in the nodeStack
        self.nodeStack.append(self.MoveStatePair(firstRelativeMove,
                                                 firstAbsoluteMove,
                                                 firstRecordedState))
        self.nodesVisited += 1

    # moves down the tree, adding a new node if testing is passed. in this
    # implementation most optimizations happen inside this function.
    def down(self, inputMove):
        # next move won't always be the inputted move.
        nextMove = copy.deepcopy(inputMove)

        # deep copy the state to get a new and completely independent object:
        currentState = copy.deepcopy(self.nodeStack[-1].getState())
        # resets the forced value to prevent memory leaks
        currentState.forced = self.numPoles

        # record the absolute move that was requested.
        nextMoveRelVal = nextMove.getRelativeMove()
        nextMoveAbsVal = currentState.poles[nextMoveRelVal].initialIndex
        nextMoveAbs = self.AbsoluteMove(nextMoveAbsVal)

        # to match the theoretical results, every new node is treated as
        # having the same failure test cost
        self.failureTestCost += self.numPoles
        self.nodesVisited += 1

        # set default action and reason statements, if in display mode
        if self.displayMode:
            action = "cut pole " + str(nextMoveAbsVal)
            action += " in position " + str(nextMoveRelVal)
            if nextMoveRelVal == 0:
                reason = str(nextMoveRelVal)
                reason += " is a valid relative move and the last move was "
                reason += "also valid"
            else:
                reason = str(nextMoveRelVal)
                reason += " is a valid relative move but the last move was "
                reason += "not valid"

        # repetition blocking
        if nextMoveAbs.getAbsoluteMove() == \
                self.nodeStack[-1].getAbsoluteMove():
            # if this letter is a repeat and incrementation is possible,
            # increment the move
            if nextMoveRelVal < self.numPoles - 1:
                newMoveRelVal = nextMoveRelVal + 1
                nextMove = self.RelativeMove(newMoveRelVal)
                nextMoveAbsVal = \
                    currentState.poles[newMoveRelVal].initialIndex
                nextMoveAbs = self.AbsoluteMove(nextMoveAbsVal)

                if self.displayMode:
                    action = "cut pole " + str(nextMoveAbsVal)
                    action += " in position " + str(newMoveRelVal)
                    reason = str(inputMove.getRelativeMove())
                    reason += " is the same as the last move, so the "
                    reason +=  "relative move was incremented"
            # if not, begin backtracking
            else:
                self.moveDownNext = False
                if self.displayMode:
                    action = "none"
                    reason = str(inputMove.getRelativeMove())
                    reason += " is the same as the last move, and there's no "
                    reason += "valid alternative."
                    outcome = "no action taken, begin backtracking"

        # symmetry blocking
        if self.moveDownNext:
            symmetryActedHere = False
            newMoveRelVal = nextMoveRelVal

            while nextMoveAbsVal in \
                    self.nodeStack[-1].getState().blockedBySymmetry:
                newMoveRelVal += 1
                # see if this move is valid
                try:
                    nextMoveAbsVal = \
                        currentState.poles[newMoveRelVal].initialIndex
                    symmetryActedHere = True
                except IndexError:
                    break
                    # no valid solution: see below for consequences

            if symmetryActedHere:
                # if the new move is valid, update both move objects
                # if this move wasn't valid, the above look exited via the
                # exception
                if newMoveRelVal < self.numPoles:
                    nextMove = self.RelativeMove(newMoveRelVal)
                    nextMoveAbs = self.AbsoluteMove(nextMoveAbsVal)

                    if self.displayMode:
                        action = "cut pole " + str(nextMoveAbsVal)
                        action += " in position " + str(newMoveRelVal)
                        reason = str(nextMoveRelVal)
                        reason += "is blocked by symmetry, so another pole "
                        reason += "must be cut before it can be"
                # if no valid new move could be found, symetry rules out
                # continuing at this point
                else:
                    self.moveDownNext = False
                    if self.displayMode:
                        action = "none"
                        reason = str(nextMoveRelVal)
                        reason += "is blocked by symmetry, and no other pole "
                        reason += "can be cut"
                        outcome = "no action taken, begin backtracking"

        # growth constraint testing and testing for forcing
        if self.moveDownNext:
            # let the poles grow
            currentState.passADay()

            # are any poles too tall for all to be allowed as the next cut?
            # compares the urgencies to their indecies
            self.moveDownNext = currentState.testForConstraintViolation()
            if self.displayMode and not self.moveDownNext:
                outcome = "no action taken due to a constraint violation, "
                outcome = "begin backtracking"

            # if the move is valid, test for forcing
            if self.moveDownNext:
                currentState.testForForcing()

            # if a pole is forced, we can only cut that pole, or one more
            # urgent than it is
            if nextMove.getRelativeMove() > currentState.forced:
                self.moveDownNext = False
                if self.displayMode:
                    outcome = "no action taken due to forcing, begin \
                               backtracking"

        # if this move has passed both tests, record it and test for success
        if(self.moveDownNext):
            # cut the inputted pole
            currentState.cutPole(nextMove)

            # add the inputted move and its consequences to the stack
            self.nodeStack.append(self.MoveStatePair(copy.deepcopy(nextMove),
                                                     nextMoveAbs,
                                                     currentState))

            # see whether this is a complete solution
            # if it is, testForSuccess will cause the exiting of the loop in
            # solve() if it's not, nothing happens
            if(self.testForSuccess()):
                if self.displayMode:
                    outcome = "complete solution found!"
            # if it's not, keep going down, starting with the first
            # letter again
            else:
                self.nextDownMove = self.RelativeMove(0)
                if self.displayMode:
                    outcome = "move successful, continue"

        if self.displayMode:
            return[action, reason, outcome]

    # tests whether the solution prefix produced so far is complete.
    # if it is, record it and use self.solvable to communicate completion.
    # if not, nothing happens.
    def testForSuccess(self):
        # there's no reason to test a solution that's too short to possibly
        # be complete - test only sufficiently long solutions.
        if len(self.nodeStack) >= self.minSolnLength:
            # ignore the settling period by starting when we know the settling
            # period is over
            i = self.minSolnLength - 1
            # limit testing to times long enough ago that they might be
            # solutions
            while i < (len(self.nodeStack) - self.minSolnLength):
                # sets were origionally used for this comparison, but the order
                # of identical poles is determined by the previous moves, so
                # this doesn't noticeably improve performance or change the
                # found result
                if self.nodeStack[i].getState().getAllUrgencies() == \
                        self.nodeStack[-1].getState().getAllUrgencies():

                    # if we have a solution, record it and stop looking.
                    self.solvable = 1
                    self.recordSoln(i)
                    break
                # if not, keep looking.
                else:
                    i += 1

                self.successTestCost += self.numPoles

        return self.solvable

    # removes a day from the stack and tests for a return to the root. if
    # the root has been reached: input is unsolvable. otherwise, will choose
    # the next move.
    def up(self):
        # if we've reached the top, this instance is unsolvable
        if len(self.nodeStack) == 1:
            self.solvable = -1

            if self.reccomendedPrints:
                print("input unsolvable due to full tree exploration")
            if self.displayMode:
                action = "stop backtracking"
                reason = "the first move is a free choice, don't have to \
                          delete it"
                outcome = "the discretised garden "
                outcome += str(self.dBGMaxSeperations) + " is unsolvable"

        # if not, pop the last move
        else:
            lastNode = self.nodeStack.pop()

            # prepare for the next move
            lastMove = lastNode.getRelativeMove()
            nextMove = lastMove + 1

            # if this value is legal, it's the next move.
            if nextMove < self.numPoles:
                self.nextDownMove = self.RelativeMove(nextMove)
                self.moveDownNext = True

                if self.displayMode:
                    outcome = "potentially legal moves remain, trying them"

            # otherwise, keep going up.
            else:
                self.moveDownNext = False

                if self.displayMode:
                    outcome = "no legal moves remain, continue backtracking"

            if self.displayMode:
                action = "backtracking, popping " + str(lastMove)
                reason = "see previous day"

        # return if return is needed
        if self.displayMode:
            return [action, reason, outcome]

    # manages the exploration of the tree with a loop and switch statement.
    # as each move is informed by the last, a recursive process was initially
    # used here - this should work better
    def solve(self):
        # start the clock
        startTime = time.time()

        # special case for short inputs
        if len(self.dBGMaxSeperations) == 1:
            if self.dBGMaxSeperations[0] >= 1:
                self.solvable = 1
                return self.solvable

        # while we don't know whether this pinwheel instance is solvable or not
        while(self.solvable == 0):

            # if we want to go down and there's somewhere to go, go down
            if self.moveDownNext:
                if self.displayMode:
                    actionReasonOutcome = self.down(self.nextDownMove)
                else:
                    self.down(self.nextDownMove)

            # otherwise, go up
            else:
                if self.displayMode:
                    actionReasonOutcome = self.up()
                else:
                    self.up()

            if self.displayMode:
                if not self.testingMode:
                    # waits for any input before proceeding
                    input()

                print("day", len(self.nodeStack))
                print()

                nextRelativeMove = self.nextDownMove.getRelativeMove()
                nextAbsoluteMove = self.nodeStack[-1].getAbsoluteMove()
                print("action:\t\t\t\t", actionReasonOutcome[0])
                print("reason:\t\t\t\t", actionReasonOutcome[1])
                print("outcome: \t\t\t", actionReasonOutcome[2])
                print()
                self.nodeStack[-1].getState().printAbsoluteState()
                self.nodeStack[-1].getState().printAllUrgencies()
                print()
                print("all abs moves:\t\t\t", self.getAllAbsoluteMoves())
                print("all rel moves:\t\t\t", self.getAllRelativeMoves())
                print()
                print()

                if self.solvable == -1:
                    print("Garden " + str(self.dBGMaxSeperations) +
                          " is unsolvable!")
                elif self.solvable == 1:
                    print("Garden " + str(self.dBGMaxSeperations) +
                          " is solved!")
                    print("Solution: ", self.soln)

        self.solveTimeCost = time.time() - startTime

        return self.solvable

    # records the absolute moves of the solution to an internal variable
    def recordSoln(self, repetitionStart):
        soln = []
        i = repetitionStart
        while i < len(self.nodeStack) - 1:
            soln.append(self.nodeStack[i].getAbsoluteMove())
            i += 1

        if len(soln) < self.minSolnLength:
            print("WARNING: self reported solution length is wrong")
            exit(0)
        
        self.soln = soln

        if self.reccomendedPrints:
            print("solution found!:\t\t", self.soln)

        self.reorderSolution()

    def reorderSolution(self):
        if self.solvable != 1 or self.soln == []:
            print("warning - attempted to reorder a solution without an existant solution")
            exit(0)

        solnQuality = self.getSolnQuality()
        sortedSolnQuality = solnQuality.copy()
        sortedSolnQuality.sort()
        while solnQuality != sortedSolnQuality:
            for i in range(len(solnQuality)):
                if solnQuality[i+1] < solnQuality[i]:
                    firstSwaplocation = i
                    break

            for i in range(len(self.soln)):
                if self.soln[i] == firstSwaplocation:
                    self.soln[i] = 't'

            for i in range(len(self.soln)):
                if self.soln[i] == firstSwaplocation + 1:
                    self.soln[i] = firstSwaplocation

            for i in range(len(self.soln)):
                if self.soln[i] == 't':
                    self.soln[i] = firstSwaplocation + 1


            solnQuality = self.getSolnQuality()

    # returns all absolute moves taken so far
    def getAllAbsoluteMoves(self):
        allMoves = []

        for n in self.nodeStack:
            absoluteMove = n.getAbsoluteMove()
            allMoves.append(absoluteMove)

        return allMoves

    # returns all relative moves taken so far
    def getAllRelativeMoves(self):
        allMoves = []

        for n in self.nodeStack:
            relativeMove = n.getRelativeMove()
            allMoves.append(relativeMove)

        return allMoves

    # prints all urgency states up to now, for bug testing
    def printAllUrgencyStates(self):
        # allUgencyStates = []

        day = 0
        for n in self.nodeStack:
            urgencyState = n.getState().getAllUrgencies()
            relMove = n.getRelativeMove()
            absMove = n.getAbsoluteMove()
            print("urgency on day", day,
                  "is", urgencyState,
                  "after cutting (rel)", relMove,
                  "aka (abs)", absMove)
            day += 1

    # sets the minimum solution length - no point in testing any solution
    # shorter than this
    def setNaiveMinSolnLength(self):
        # first check the density, because any pinwheel case that would
        # fail on density grounds has an infinite minimum solution length

        self.density = Fraction(0,1)
        for i in self.dBGMaxSeperations:
            self.density += Fraction(1,i)

        # if this garden is unsolvable on density grounds, say so
        if self.density > Fraction(1,1):
            self.solvable = -1
            if self.reccomendedPrints:
                print("input failed due to density")
        # otherwise, set the minimum solution length
        else:
            # start with the length of the round robin solution.
            self.minSolnLength = len(self.dBGMaxSeperations)
            # in the round robin solution, each pole appears once.
            minNumber = [1]*self.minSolnLength

            # when we pass through the loop without acting, we've
            # found a composition that works
            actedThisItteration = True
            while actedThisItteration:
                actedThisItteration = False
                for i in range(len(minNumber)):
                    # while the proportion of this element is too small
                    while (minNumber[i]/self.minSolnLength) < \
                            (1/self.dBGMaxSeperations[i]):
                        # add more of this element
                        minNumber[i] += 1
                        self.minSolnLength += 1
                        actedThisItteration = True

    # find the hardest PWS problem that is solved by this solution.
    def getSolnQuality(self):
        if self.solvable != 1:
            print("unsolved!")
            exit(1)
        else:
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

    # a collection of test cases with known solvability, designed to use
    # different parts of the code and hence to test them
    @staticmethod
    def runTestCases():
        print("should succeed:")
        oursolver_foresight = solver_foresight([6, 4, 2], False, True, False)
        oursolver_foresight.solve()
        print()
        print()

        print("should fail:")
        oursolver_foresight = solver_foresight([4, 3, 2], False, True, False)
        oursolver_foresight.solve()
        print()
        print()

        print("should succeed:")
        oursolver_foresight = solver_foresight([8, 4, 2], False, True, False)
        oursolver_foresight.solve()
        print()
        print()

        print("should succeed:")
        oursolver_foresight = solver_foresight([8, 8, 8, 8, 2], False, True, False)
        oursolver_foresight.solve()
        print()
        print()

        print("should fail:")
        oursolver_foresight = solver_foresight([3, 3, 5, 15], False, True, False)
        oursolver_foresight.solve()
        print()
        print()

    # cases that reproduce a bug - this function will not be cleared out
    # when the bug is fixed, so that they can be used to fix the next bug.
    @staticmethod
    def runBugFixingCases():
        oursolver_foresight = solver_foresight([2], True, True, False)  # 0
        print(oursolver_foresight.solve())

#solver_foresight.runBugFixingCases()