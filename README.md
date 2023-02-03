# quasilethal
The code for the article "Evolutionary adaptation is facilitated by the presence of lethal genotypes" by Viktoryia Blavatska and Bartlomiej Waclaw, arXiv:2301.07599 [q-bio.PE]

# directed_evolution_trajectories_both_models.cpp

This code simulates the stochastic quasispecies model with multiple viable and one non-viable genotypes. It generates a series of random fitness landscapes and simulates the model with and without lethal mutations on each landscape. It outputs the time to best adapted genotype and a sequence of genotypes taken by evolution (evolutionary trajectory).
   
The code has been tested with GNU Compiler (ver. 4.9.2 or above). To compile it, please use
	
g++ directed_evolution_trajectories_both_models.cpp -w -O3 -o evolution.exe
	
then run the executable "evolution.exe" (no command-line parameters required).
	
All parameters are set in the "Parameters" section of the code below. Changing any of the parameters requires the code to be re-compiled.
	
As is, the program will generate data for mu=0.04 and gamma=0.76, for 1000 fitness landscapes, and 20 replicate simulations on each landscape.
