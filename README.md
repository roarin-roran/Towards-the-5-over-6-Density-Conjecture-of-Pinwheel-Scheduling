this code accompanies the paper "Towards the 5/6-Density Conjecture of Pinwheel Scheduling" published in ALENEX '22

#### BASIC INFO
code was developed and run using python 3.7.9 on Windows, and is stable on this platform.

quoted expected running times are for a Windows computer with a 6 core i5-8400 CPU and 16GB of RAM - performance will vary across platforms, and several operations use a very large amount of RAM.


#### Pareto surfaces (table 1)
to produce any complete Pareto surface (that is, where density = 1), use the run_pareto.py method, by typing:
'python run_pareto.py k'
in your terminal, replacing k with your desired k value. 

expect values <=4 to complete in under a second, k=5 to take several minutes, and k=6 to take over a month (I've yet to get it to complete, but have run it for that long).
The output will be printed to the terminal window in the form of a Pareto surface (explained in the paper) and the analogous *failure surface* - an inclusion minimal set of problems that cannot be solved or dominated. note that low k members of the failure surface like 2,5,7 cannot be the prefix of a solvable instance.


#### Pinwheel oracle
to solve any pinwheel scheduling problem, use the run_pinwheel_solver.py method, by typing:
'python run_pinwheel_solver.py a_0,a_1,a_2...a_k'
in your terminal, replacing a_0,a_1,a_2...a_k with your desired values of a, separated by commas.

expect unsolvable, large k, or large a problems to take substantially longer than solvable, low k, or low a problems - but also to see huge, unpredictable differences between instances.
results will be printed to the terminal window.


#### Pareto surfaces limited by density (figure 3, table 2)
to produce Pareto surfaces limited by density using any pinwheel oracle or approximation method introduced in the paper, open 'run_pwd.py' in a text editor.
fill the list following desiredKValues with all k values you wish to calculate. note that if using P(k-1) as your approximation method, all smaller values will also be calculated.
select your desired oracle by setting desiredSolver = NAIVE, OPTIMISED, or FORESIGHT as desired (note that graph-based solving wasn't implemented here as it's ridiculously slow, and wasn't analyzed this way in the paper).
select your desired approximation method by setting desiredApproxMethod = NO_APPROX, P_ONE, P_FIVE, or P_K_MINUS_ONE as desired.
set 'keepCertificates = True' if you want the certificate for your final run to be tested, and False otherwise (this is quite a slow process)
set 'endlessRunning = True' if you want to calculate the surfaces for the same k values repeatedly for time trials (this will never test a certificate)
use outputFilename to set your desired file for .csv output.

with optimal settings (FORESIGHT, P_K_MINUS_1), expect k=11 to complete in ~10 hours. k=12 can be successfully completed with these settings, but requires over 100GB of RAM at times and between 10 and 24 hours on a high-end server (no time trials could be taken here, as the server we used had a variable task load)
results will be appended to the .csv file that you name above, which will be created if it does not exist. note that these results will not include surfaces - if you want these, you'll have to expand the pareto_w_density.py code to output them to a file yourself as they were huge and we didn't evaluate them manually.

#### Pinwheel solver (figure 2)
to compare the performances of the four pinwheel oracles introduced, use the 'pinwheel_solver_tester.py' method, by typing:
'python pinwheel_solver_tester.py r'
in your terminal, replacing r with which round you wish to reproduce: 1, 2, or 3

each round runs indefinitely, generating random instances and attempting to solve them with the top 5-r methods. 
while the instances are random, the random seed is set, so repeats of the same round of the tournament will use the same instances.
results will be appended to a csv file (one per round), in real-time. 
