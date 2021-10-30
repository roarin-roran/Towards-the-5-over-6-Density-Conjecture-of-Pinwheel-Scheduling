class pareto_surface_manager:
    class SurfacePoint:
        # adds self.density
        def getDensity(self):
            self.density = 0

            for i in self.recurrenceVector:
                self.density += 1/i

    # a point on the pareto surface, containing a solution and its
    # reoccurence vector.
    class PSPoint(SurfacePoint):
        def __init__(self, soln, recurrenceVector):
            self.soln = soln
            self.recurrenceVector = recurrenceVector

            # ensure that the recurrenceVector is in ascending order
            self.fixPSPointOrder()

            # ensure that the solution is sorted lexographically
            self.chooseSolnStartPoint()

            # adds density based on the rv
            self.getDensity()

        # the recurrence vector of a PSP must ascending, or else the 
        # points have the wrong labels
        def fixPSPointOrder(self):
            while self.recurrenceVector != sorted(self.recurrenceVector):

                # find an element that's out of place
                for i in range(len(self.recurrenceVector)):
                    for j in range(len(self.recurrenceVector) - i - 1):
                        if self.recurrenceVector[i] > self.recurrenceVector[j+i+1]:
                            # swap the characters around

                            firstChar = i
                            secondChar = j+i+1

                            # replace the first character with a placeholder
                            for k in range(len(self.soln)):
                                if self.soln[k] == firstChar:
                                    self.soln[k] = "placeholder"
                            # replace the second character with the first character
                            for k in range(len(self.soln)):
                                if self.soln[k] == secondChar:
                                    self.soln[k] = firstChar
                            # replace the placeholders with secondChar
                            for k in range(len(self.soln)):
                                if self.soln[k] == "placeholder":
                                    self.soln[k] = secondChar

                            # the solution has been updated: update the sqv
                            self.recurrenceVector = pareto_surface_manager.getRecurrenceVectorFor(self.soln)

        # makes all cyclic permutations of a solution, sorts them lexographically and 
        def chooseSolnStartPoint(self):
            allCyclicPerms = []

            for i in range(len(self.soln)):
                allCyclicPerms.append(self.soln[i:] + self.soln[:i])

            allCyclicPerms.sort()

            self.soln = allCyclicPerms[0]

    class FSPoint(SurfacePoint):
        def __init__(self, recurrenceVector):
            self.recurrenceVector = recurrenceVector

            self.getDensity()

    # creates a pareto_surface_manager
    def __init__(self, numberOfPoles):
        self.numberOfPoles = numberOfPoles

        # initialise the pareto surface with the round robin solution
        roundRobin = []
        for i in range(self.numberOfPoles):
            roundRobin.append(i)

        # the roud robin solution has a known recurrence vector
        roundRobinRV = [numberOfPoles]*numberOfPoles


        # initialise the three surfaces using the round robin solution.
        self.paretoSurface = [self.PSPoint(roundRobin, roundRobinRV)]
        self.querySet = [[numberOfPoles - 1]]
        self.failureSurface = []

    def printParetoSurface(self):
        j = 1
        for i in self.paretoSurface:
            print(j, ":\trv", i.recurrenceVector, "\tdensity", round(i.density, 3), "\tsoln", i.soln)
            j += 1

    def printFailureSurface(self):
        j = 1
        for i in self.failureSurface:
            print(j, ":\tdensity", round(i.density, 3), "\t rv", i.recurrenceVector)
            j += 1

    # compares an inputted point with either the pareto surface or failure surface.
    def compareWithSurface(self, inputPoint, pointSatisfiable):
        if pointSatisfiable == False:
            pass
            #print("temporarily adding all points to failure surface")
            #self.failureSurface.append(inputPoint)
            #return

        if pointSatisfiable:
            surface = self.paretoSurface
            compare = pareto_surface_manager.compareTwoPSPoints
        else:
            surface = self.failureSurface
            compare = pareto_surface_manager.compareTwoFSPoints

        inputDominated = False
        
        i = 0
        while i < len(surface):
            comparisonResult = compare(surface[i], inputPoint)

            # old point dominates new point: don't add it and stop comparing it to other things.
            if comparisonResult == 1:
                inputDominated = True
                break
            # new point dominates the old point: remove the old point and keep going
            elif comparisonResult == 2:
                surface.pop(i)
            # both points have their merits: both might be needed
            elif comparisonResult == 3:
                i += 1
            # invalid result: report and exit
            else:
                print("invalid comparison result reached compareWithSurface")
                exit(1)

        # if the input is undominated, add it
        if not inputDominated:
            surface.append(inputPoint)

        # sort the surface - while there's a small cost to doing this every time,
        # it makes bug fixing easier, and isn't going to be significantly expensive
        if pointSatisfiable:
            self.sortParetoSurface()
        else:
            self.sortFailureSurface()

    # compares two pareto surface points - an old one (already on the PS)
    # and a new one (which might be added to the pareto surface)
    # there are four possible outcomes:
    #     1: the old point dominates the new point
    #     2: the new point dominates the old point
    #     3: the two points are different, but neither is dominated
    # PS points must be full length, and therefore must be tight by definition.
    @staticmethod
    def compareTwoPSPoints(oldPSP, newPSP):
        # old/new beats new/old in any respect.
        oldBeatsNew = False
        newBeatsOld = False
        returnValue = -1

        # if they're the same, don't bother comparing the values
        if oldPSP.recurrenceVector == newPSP.recurrenceVector:
            # we want to keep the shorter solution
            if len(oldPSP.soln) <= len(newPSP.soln):
                returnValue = 1
            else:
                returnValue = 2
        else:
            # run through the solutions, comparing them
            for i in range(len(oldPSP.recurrenceVector)):
                if oldPSP.recurrenceVector[i] < newPSP.recurrenceVector[i]:
                    oldBeatsNew = True
                elif oldPSP.recurrenceVector[i] > newPSP.recurrenceVector[i]:
                    newBeatsOld = True

            # if the old point dominates
            if oldBeatsNew and not newBeatsOld:
                returnValue = 1
            # if the new point dominates
            elif newBeatsOld and not oldBeatsNew:
                returnValue = 2
            # if the points both have their merits
            elif newBeatsOld and oldBeatsNew:
                returnValue = 3
            # we should have caught this already
            else:
                print("compareTwoPSPoints had an error - unreachable state reached")

        if returnValue == -1:
            print("comparing two PSPoints had an error - no result found")
            exit(0)
        else:
            return returnValue

    # compares two failure surface points - an old one (already on the FS)
    # and a new one (which might be added to the failure surface)
    # there are three possible outcomes:
    #     1: the old point dominates the new point
    #     2: the new point dominates the old point
    #     3: the two points are different, but neither is dominated
    # the inputs to this function are assumed to be tight - no part of them can be 
    # deleted without giving a solvable query. this is enforced by the query
    # generation method only.
    @staticmethod
    def compareTwoFSPoints(oldFSP, newFSP):
        oldFSRV = oldFSP.recurrenceVector
        newFSRV = newFSP.recurrenceVector
        
        # initialise these: old/new beats new/old in any shared task seperation
        oldBeatsNew = False
        newBeatsOld = False
        returnValue = -1

        # compare points only if they have the same lengths
        if len(newFSRV) != len(oldFSRV):
                #print("tested 6")
                returnValue = 3
        else:
            # compare all RV components
            for i in range(len(newFSRV)):
                # see if any task gaps are dominant
                if oldFSRV[i] > newFSRV[i]:
                    oldBeatsNew = True
                elif oldFSRV[i] < newFSRV[i]:
                    newBeatsOld = True

            if oldBeatsNew and not newBeatsOld:
                returnValue = 1
            elif newBeatsOld and not oldBeatsNew:
                returnValue = 2
            elif newBeatsOld and oldBeatsNew:
                returnValue = 3
            else:
                returnValue = 1
            
        if returnValue > 0:
            return returnValue
        else:
            print("error in compareTwoFSPoints: illegal return value")
            exit(1)

    # accepts a solution and returns the recurence vector for that solution
    @staticmethod
    def getRecurrenceVectorFor(soln):
        numberOfPoles = max(soln) + 1

        recurrenceVector = []

        doubleSoln = soln + soln

        # for all poles
        for i in range(numberOfPoles):
            # measure the max gap between that pole and the next selfSimilar pole
            maxGap = 0
            gap = 0

            # for every part of the solution
            for j in doubleSoln:
                # if we've just found an instance of the pole we're looking for
                if isinstance(j, int):
                    move = j
                else:
                    move = j.move

                if move == i:
                    # see if this is the max
                    if gap > maxGap:
                        maxGap = gap

                    # reset the gap - we just found this pole!
                    gap = 0

                gap += 1

            # test that that task is ever done
            if maxGap == 0:
                print("solution incomplete! task", i, "is never done!")
                exit(1)

            # append the value to the reccurence vector
            recurrenceVector.append(maxGap)
        return recurrenceVector

    # uses a simple insertion sort to put the pareto surface in
    # ascending order of recurrence vector
    def sortParetoSurface(self):
        i = 1
        while i < len(self.paretoSurface):
            j = i
            while j > 0 and self.firstPSPointBigger(j-1,j):
                placeholder = self.paretoSurface[j-1]
                self.paretoSurface[j-1] = self.paretoSurface[j]
                self.paretoSurface[j] = placeholder

                j -= 1

            i += 1

    # returns true if the first of two pareto surface points has a lexographically larger
    # sqv than the second point, else false.
    def firstPSPointBigger(self, firstPSPointPosn, secondPSPointPosn):
        firstPSPoint = self.paretoSurface[firstPSPointPosn]
        secondPSPoint = self.paretoSurface[secondPSPointPosn]

        if firstPSPoint.recurrenceVector > secondPSPoint.recurrenceVector:
            #print("yes")
            return True
        else:
            #print("no")
            return False

    def sortFailureSurface(self):
        i = 1
        while i < len(self.failureSurface):
            j = i
            while j > 0 and self.firstFSPointBigger(j-1,j):
                placeholder = self.failureSurface[j-1]
                self.failureSurface[j-1] = self.failureSurface[j]
                self.failureSurface[j] = placeholder

                j -= 1

            i += 1

    def firstFSPointBigger(self, firstFSPointPosn, secondFSPointPosn):
        firstFSPoint = self.failureSurface[firstFSPointPosn]
        secondFSPoint = self.failureSurface[secondFSPointPosn]

        if firstFSPoint.recurrenceVector > secondFSPoint.recurrenceVector:
            #print("yes")
            return True
        else:
            #print("no")
            return False

    # runs a series of tests that checks that compareTwoFSPoints is working properly
    @staticmethod
    def testCompareTwoFSPoints():
        ourSE = pareto_surface_manager(4)

        allGood = True

        case = 1
        oldFSP = ourSE.FSPoint([3,3,5])
        newFSP = ourSE.FSPoint([2,4,7])
        expectationValue = 3

        if pareto_surface_manager.compareTwoFSPoints(oldFSP, newFSP) != expectationValue:
            print("we have a problem with case", case)
            allGood = False
        else:
            print("no problems with case", case)

        case += 1
        oldFSP = ourSE.FSPoint([2,3,5])
        newFSP = ourSE.FSPoint([2,4,7])
        expectationValue = 1

        if pareto_surface_manager.compareTwoFSPoints(oldFSP, newFSP) != expectationValue:
            print("we have a problem! with case", case)
            allGood = False
        else:
            print("no problems with case", case)

        case += 1
        oldFSP = ourSE.FSPoint([2,3])
        newFSP = ourSE.FSPoint([2,4,7])
        expectationValue = 3

        if pareto_surface_manager.compareTwoFSPoints(oldFSP, newFSP) != expectationValue:
            print("we have a problem! with case", case)
            allGood = False
        else:
            print("no problems with case", case)

        case += 1
        oldFSP = ourSE.FSPoint([2,4,7])
        newFSP = ourSE.FSPoint([2,3,5])
        expectationValue = 2

        if pareto_surface_manager.compareTwoFSPoints(oldFSP, newFSP) != expectationValue:
            print("we have a problem! with case", case)
            allGood = False
        else:
            print("no problems with case", case)

        case += 1
        oldFSP = ourSE.FSPoint([2,4,7])
        newFSP = ourSE.FSPoint([2,3])
        expectationValue = 3

        if pareto_surface_manager.compareTwoFSPoints(oldFSP, newFSP) != expectationValue:
            print("we have a problem! with case", case)
            allGood = False 
        else:
            print("no problems with case", case)

        case += 1
        oldFSP = ourSE.FSPoint([2,3,5])
        newFSP = ourSE.FSPoint([2,3])
        expectationValue = 2

        if pareto_surface_manager.compareTwoFSPoints(oldFSP, newFSP) != expectationValue:
            print("we have a problem! with case", case)
            allGood = False 
        else:
            print("no problems with case", case)

        case += 1
        oldFSP = ourSE.FSPoint([2,3])
        newFSP = ourSE.FSPoint([2,3])
        expectationValue = 1

        if pareto_surface_manager.compareTwoFSPoints(oldFSP, newFSP) != expectationValue:
            print("we have a problem! with case", case)
            allGood = False
        else:
            print("no problems with case", case)

        return allGood

    # runs a series of tests that checks that compareWithSurface is working properly
    # UNFINISHED
    @staticmethod
    def testCompareWithSurface():
        print("I'm not complete, just a solid start - testCompareWithSurface")
        ourSE = pareto_surface_manager(4)
        
        # pareto surface - needs more points and a success condition.

        solns = [[0,1,0,2,0,1,0,3], [0,1,0,2,0,3], [0,1,2,0,1,3]]
        RVs = []

        for i in solns:
            RVs.append(pareto_surface_manager.getRecurrenceVectorFor(i))

        for i in range(len(solns)):
            PSP = ourSE.PSPoint(solns[i], RVs[i])

            ourSE.compareWithSurface(PSP, True)

        print("paretoSurface:")
        ourSE.printParetoSurface()

        # failure surface - also needs more points and a success condition.
        ourSE.compareWithSurface(ourSE.FSPoint([2,4]), False)
        ourSE.compareWithSurface(ourSE.FSPoint([2,3]), False)
        ourSE.compareWithSurface(ourSE.FSPoint([1]), False)

        print("failureSurface: ")
        ourSE.printFailureSurface()

        return ourSE

    # compares a query with the query set, Pareto surface and failure surface
    def addQueryIfUnseen(self, query):
        queryAnswered = False

        # first, test if a query is already in the set. there's no domination here
        # because we don't know if any query in the set is doable or not.
        if query in self.querySet:
            #print("already got", query, "so it wasn't added")
            queryAnswered = True
        else:
            # if the query isn't in the query surface, see if it's satisfied by a
            # point on the pareto surface
            for i in self.paretoSurface:
                #print(i.recurrenceVector, "vs", query)
                
                # default to acceptance, then test if this is true
                pSPointSatisfiesQuery = True
                for j in range(len(query)):
                    if i.recurrenceVector[j] > query[j]:
                        pSPointSatisfiesQuery = False
                        break

                        #print("not this one, because", i.recurrenceVector[j], "is bigger than", query[j])

                if pSPointSatisfiesQuery:
                    #print("not this one, because", query, "is satisfied by the PS")
                    queryAnswered = True
                    break

            # if the query isn't answered by the pareto surface, check the failure surface
            if not queryAnswered:
                #print("still going strong!", query, "is not in the PS")

                # of the query isn't in the pareto surface, see if it's dominated by any
                # point on the failureSurface
                for i in self.failureSurface:
                    if self.compareTwoFSPoints(i, query) == 1:
                        #print("this query,", query, ", is dominated!")
                        queryAnswered = True
                        break

        if not queryAnswered:
            #print("query should be added")
            self.querySet.append(query)
            self.querySet.sort()

    # add all the tests together in one place so they can easily be turned on or off.
    @staticmethod
    def test():
        print("testing compareTwoFSPoints")
        print("not testing everything else - test suites not written.")

        allGood = True

        allGood = allGood and pareto_surface_manager.testCompareTwoFSPoints()

        allGood = allGood and pareto_surface_manager.testCompareWithSurface()

        if allGood:
            print("all's good here chief")
        else:
            print("some tests returned problems")

    # runs whatever we have setup at the moment
    @staticmethod
    def run():
        #add a query trash clearer - if one query produces a result, that result may resolve
        # other queries.

        #ourSE = pareto_surface_manager.testCompareWithSurface()

        #print("querySet:")
        #print(ourSE.querySet)

        #ourSE.addQueryIfUnseen([3, 4, 5])


        ourSE = pareto_surface_manager(5)

        ourPSP = ourSE.PSPoint([4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0], [2, 4, 8, 16, 16])
        ourSE.compareWithSurface(ourPSP, True)
        
        ourPSP = ourSE.PSPoint([4, 0, 1, 2, 0, 1, 3, 0, 1], [3, 3, 9, 9, 9])
        ourSE.compareWithSurface(ourPSP, True)

        ourPSP = ourSE.PSPoint([4, 2, 0, 1, 3, 2, 0, 1, 2, 3, 0, 1], [4, 4, 5, 7, 12])
        ourSE.compareWithSurface(ourPSP, True)        

        ourSE.printParetoSurface()
        
        print("all done pal")


#pareto_surface_manager.test()

#pareto_surface_manager.run()