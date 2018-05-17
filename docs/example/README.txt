Once SPHERLS is installed and the bin directory is in your path do the following
to create a starting model and run a simulation.

To create starting model and have a look at it's structure.
$ SPHERLSgen
$ plot_profile.py plot_profile.xml --remake-profiles

The first command generates a new starting model from the settings in 
SPHERLSgen.xml, the second command generates a radial profile file 
"T6500_20x1_t00000000_pro.txt" which is an ascii file with column headers etc. 
made from the binary starting model file then it makes a plot of the radial 
profile as specified in the plot_profile.xml settings file.

There are two additional files "SPHERLSgen_first_try.xml" and 
"starting_model_structure_first_try.pdf" which represent my first attempt at 
creating a model. As you can see the gradients at the surface are not well 
resolved and I would consider this a bad starting model. In particular there 
is a jump in temperature around 1e4 K, and also a density inversion around 
there too. This is because of the hydrogen ionization region. One thing to be 
careful of is that you don't want the density to become too low at the surface 
as it becomes more likely to go outside the equation of state tables.

To run the simulation
$ mpirun -n 4 SPHERLS

Note that in SPHERLS.xml I choose the number of surface zones treated 
implicitly based on the starting model structure. Roughly speaking I choose 
this number so that the implicit region ends well below the hydrogen ionization
zone.

$ plot_profile.py plot_profile_output.xml

This creates a number of plots of the structure similar to the one we made for 
the starting model for all the model dumps of the simulation. We only calculate 
10 minutes of simulation time so not much happens over the course of multiple 
plots.
