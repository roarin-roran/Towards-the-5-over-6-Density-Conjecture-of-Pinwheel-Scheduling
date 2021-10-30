import pareto_w_density as pwd
import certificate_checker
from fractions import Fraction

# constants
DENSITY = Fraction(5,6)
# solvers
NAIVE     = 0
OPTIMISED = 1
FORESIGHT = 2
# approximation methods
NO_APPROX   = 0
P_ONE       = 1
P_FIVE      = 2
P_K_MINUS_1 = 3


# * * * * * ** * * * * * *
# * * * * * ** * * * * * *
# * * * * * MAIN * * * * *
# * * * * * ** * * * * * *
# * * * * * ** * * * * * *
if __name__ == '__main__':
	desiredKValues      = [7]
	desiredSolver       = FORESIGHT
	desiredApproxMethod = P_K_MINUS_1
	
	keepCertificates = True
	endlessRunning   = False

	outputFilename = "PWDTimeCosts.csv"

	firstRun = True
	while endlessRunning or firstRun:
		firstRun = False

		for k in desiredKValues:
			ourPWDR = pwd.pareto_w_density(densityLimit = DENSITY,
										   lengthLimit = k, 
				                           solver = desiredSolver, 
				                           approxMethod = desiredApproxMethod, 
				                           keepCerts = keepCertificates,
				                           outputFilename = outputFilename)

			if keepCertificates:
				certificate = ourPWDR.solve()
			else:
				ourPWDR.solve()



	if keepCertificates:
		checker = certificate_checker.certificate_checker(certificate, DENSITY, desiredKValues[-1])
	
		checker.check()



