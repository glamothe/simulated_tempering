#############################################################################
########################## Python Simulated Tempering #######################
##########################            2016            #######################
#############################################################################

Author : Gilles LAMOTHE
References : Nguyen 2013, Park 2007

############################### REQUIREMENTS ################################

python v. 2.7
gromacs v. 5.0.5
python modules : math, numpy, re, shutil, filinput, 
		 os, subprocess, shlex, optparse

################################# ABOUT ST ##################################

Simulated Tempering (ST) consists in doing many successive MD simulations or 
Monte Carlo runs with varying temperatures. This allows to explore more of the
potential energy landscape and better sample the configuration space of a system. 
Temperature weights are important for choosing temperature transitions
throughout the simulation. Previously, weights needed to be determined before
running an ST, which was very costly in time.

This program aims at running an ST with MD runs (using gromacs) according to the
method proposed by the Nguyen 2013 publication. This new ST method is attractive
because it allows for an on-the-fly determination of the temperature weights
which in turn avoids the task of the predetermining the weights. 

################################ HOW TO RUN AN ST ############################

The ST program is run by the STrun.py script. 
Call the script in the shell using python:

	$ python path_to_program/STrun.py [arguments]
-------------------------------------------------------------------------------
A few arguments are required to run this program :

	--gro => .gro file name
	--mdp => .mdp file name
	--top => .top file name
-------------------------------------------------------------------------------
Example :

$ python path_to_program/STrun.py --gro myfile.gro --mdp myfile.mdp --top myfile.top

-------------------------------------------------------------------------------
The following files are required in the working directory: 

--> a gromacs structure file (.gro)
--> a gromacs parameters file (.mdp)
--> a gromacs topology file (.top)
--> the ST parameters file (ST.par)

The ST.par was invented for the purpose of this program. A copy of the file is
provided in the program's directory. Details about the ST parameters are
available in the ST.par file. 

Examples of .gro, .mdp, and .top files are also provided.

########################### MORE ABOUT THE ST.par FILE ########################

The ST parameters are:

--> temp_min : the base temperature of the simulation (K)
--> temp_max : the peak temperature of the simulation (K)
--> temp_interval : difference between two neighbouring temperatures
--> steps_total : number of steps for entire ST (includes all MD runs)
--> steps_md : number of steps for a single MD run
--> steps_trans : number of steps before we attempt a temperature transition
--> steps_save : ST stats are saved every time we run this many steps
--> dt : time interval between two steps (ps)

A separate mdp file will be created for each temperature specified by the ST.par file.

The default values are:

temp_min = 300
temp_max = 500
temp_interval = 4
steps_total = 100000000 ; 200 ns
steps_md = 2000 ; 40 ps
steps_trans = 2000 ; 40 ps
steps_save = 2000 ; 40 ps
dt = 0.002 ; ps

Note1: it is safer if temp_max - temp_min is a multiple of temp_interval
Note2: steps_total > steps_md >= steps_trans and steps_save

####################################  ST OUPUT  #################################

The important output ST files are:

---> *_total.edr   : gromacs energy file (from concatenated MD runs)
---> *_total.xtc   : gromacs trajectory file (from concatenated MD runs)
---> stats.dat     : data about the ST run : steps, temperatures, energy, weights
---> reference.gro : the content of the original .gro file is saved in this file

######################################## TO DO ##################################

--> Smarter error messages
    e.g. Exceptions that tell the user when wrong ST parameter values are inserted in ST.par

--> It would be interesting to implement Monte Carlo simulations as well.

--> Add the option that allows the user to give their own list of temperates instead 
    of the min/max/interval method.

#################################### KNOWN PROBLEM ##############################

When lauching the program on certain machines of the computer cluster, the energy 
was not calculated. It is suspected the issue comes from how those certain machines
handled the piped gromacs commands that are used to fetch the energy values. 



