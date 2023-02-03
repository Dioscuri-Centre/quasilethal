# quasilethal
The code for the article "Evolutionary adaptation is facilitated by the presence of lethal genotypes" by Viktoryia Blavatska and Bartlomiej Waclaw, arXiv:2301.07599 [q-bio.PE]

# directed_evolution_trajectories_both_models.cpp

This code simulates the stochastic quasispecies model with multiple viable and one non-viable genotypes. It generates a series of random fitness landscapes and simulates the model with and without lethal mutations on each landscape. It outputs the time to best adapted genotype and a sequence of genotypes taken by evolution (evolutionary trajectory).
   
The code has been tested with GNU Compiler (ver. 4.9.2 or above). To compile it, please use
	
g++ directed_evolution_trajectories_both_models.cpp -w -O3 -o evolution.exe
	
then run the executable "evolution.exe" (no command-line parameters required).
	
All parameters are set in the "Parameters" section of the code below. Changing any of the parameters requires the code to be re-compiled.
	
As is, the program will generate data for mu=0.04 and gamma=0.76, for 1000 fitness landscapes, and 20 replicate simulations on each landscape.

# adaptation_time_letal.c"

The code produces the values of adaptation times for quasispecies model as based on averaging on 10000 replicas of randomly constructed fitness landscapes,
with variable parameters L (genotype sequence), mu (mutation rate), gamma (probability to mutate into lethal state)

To compile it, please use 
gcc /adaptation_time_letal.c -o /adaptation_time_lethal.out -lm

By entering the values for L, mu and gamma on flight, one receives the data file "T_L=(L)mu=(mu)gamma=(gamma).dat"
with a value of averaged adaptation time, and the data file "HistogramL=(L)mu=(mu)gamma=(gamma).dat" with a histogram 
(distribution) of adaptation time values.

# fitness_valley_length.c"

The code analyzes all possible pathways available for genotype of length L on a constructed random landscape and choses the ``slowest'' pathways, where the first mutation increases fitness to approx. the same value as the fitness of the final genotype, but subsequent mutations led through a fitness valley. 

To compile it, please use 
gcc /fitness_valley_length.c -o /fitness_valley_length.out -lm

By entering the value for L on flight, one receives the data file "ValleylengthL=(L).dat" with the value corresponding to the minimum length of fitness valley l, obtained by averaging over 10000 random fitness landscapes. 
