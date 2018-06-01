You can create a new combined equation of state and opacity file, often referred
to just as an "eos" file using the following command in this directory:

$ eos_interp.py eosNewY25Z01.xml

In the above ".xml" file there are a number of equation of state and opacity 
files listed under the <eos> and <opacity> nodes. This don't often need to be 
changed. More often what is changed are the elements under the 
<interpolatedTables> element near the bottom of the file. Hear you can specify
one ore more <table> elements which define new eos table files. These new
tables contain equaiton-of-state at regular intervals in log_10 density and
log_10 temperature. You can specify the lower values for the temperature and
density as well as the spacing and the number of entries in each varaiable.
You can also specify the composition in hydrogen and metal mass fractions, X and
Z respectively. And also indicate if plots should be made to inspect the table
after it is made. If the <plot> element is set to "True" it will present plots 
of energy, pressure and opacity of the newly created table using matplotlib so 
ensure you have an X11 server running to be able to see the plots. If the <setNans>
element is set to False it will fill in missing data in the tables represented by
"NaNs" with extrapolated data. If you want to ensure that your model only uses
interpolated equation of state and opacity data rather than extrapolated data 
set this to "True", which is the default value, and the "NaNs" will remain in
the final tables and if a NaN is encountered during the calculation it will stop
with a error saying a NaN was encountered.