# -*- coding: utf-8 -*-
"""
Created on Sat Sep 17 22:44:14 2016

@author: gilles
"""

from Temperature_class import *
from MD_class import *
import numpy as np
import re
from shutil import copyfile


class ST:
    """ ST class used to start a Simulated Tempering run
        _____________________________
        Arguments to create instance:
        -----------------------------
        mdp (string) : .mdp file name
        gro (string) : .gro file name
        top (string) : .top file name
        
        Constructor also takes parameters from the ST.param file. 
        
        The parameters taken from ST.param are:
        --> steps_total : number of steps for entire ST (includes all MD runs)
        --> steps_md : number of steps for a signle MD run
        --> steps_trans : number of steps before we attempt a temperature transition
        --> step_save : ST stats are saved everytime we run this many steps
        --> dt : time interval between two steps (ps)
        ___________
        Attributes:
        -----------
        
        steps_total (int)
        steps_md (int)
        steps_trans (int)
        steps_save (int)
        dt (float)
        temp (Temperature object)
        temp_list (list of int) : list of all temperatures used by the ST
    """

    def __init__(self, mdp, gro, top):
        # Set ST parameters        
        self.steps_total,\
        self.steps_md,\
        self.steps_trans,\
        self.steps_save,\
        self.dt = self.init_steps()
        self.temp = Temperature() # create Temperature object
        self.temp_list = range(self.temp.min,
                               self.temp.max + self.temp.interval,
                               self.temp.interval)
        # Create MD object which handles individual MD runs
        self.md = MD(mdp, gro, top, self.steps_md, self.dt, self.temp_list)
        # Copy original structure (gro file) for future RMSD calculations
        copyfile(gro, 'reference.gro')
        # Create and write header in ST output stats file
        self.create_stats_file()

    def init_steps(self):
        """ Fetches ST's step attributes in the ST.par parameters file.
            _____________________________
            Arguments to create instance:
            -----------------------------
            No arguments. Takes parameters from the ST.param file. 
            
            The parameters taken from ST.param are:
            --> steps_total (int): number of steps for entire ST (includes all MD runs)
            --> steps_md (int): number of steps for a signle MD run
            --> steps_trans (int): number of steps before we attempt a temperature transition
            --> step_save (int): ST stats are saved everytime we run this many steps
            --> dt (float): time interval between two steps (ps)
            _______
            Return:
            -------
            steps_total(int), steps_md(int), steps_trans(int), steps_save(int), dt(float)
        """
        # Initialize variables with -1
        steps_total, steps_md, steps_save, steps_trans, dt = [-1]*5
        # Get steps variables from parameters file
        f_param = open('ST.par', 'r')
        for line in f_param:
            if line.startswith('steps_total'):
                steps_total = int(re.findall("\d+", line)[0])
            elif line.startswith('steps_md'):
                steps_md = int(re.findall("\d+", line)[0])
            elif line.startswith('steps_save'):
                steps_save = int(re.findall("\d+", line)[0])
            elif line.startswith('steps_trans'):
                steps_trans = int(re.findall("\d+", line)[0])
            elif line.startswith('dt'):
                dt = float(re.findall("\d+\.\d+", line)[0])
        f_param.close()
        # TODO: Add exception if one of the parameters is not found in file
        # or if there are inconsistancies (ex: steps_trans < steps_md)
        return steps_total, steps_md, steps_trans, steps_save, dt

    def start(self):
        """ Starts a simulated tempereing run
        """
        for step in xrange(self.steps_md, self.steps_total+2, self.steps_md):
            print "##########################################################"
            print "##########################################################"
            # Execute an MD run
            self.md.start(self.temp.current)

            # Increment n for current temperature.
            # (n = nb of times a given temperature was selected for an MD run)
            n = self.temp.get_n(self.temp.current) + 1
            self.temp.set_n(self.temp.current, n)            

            # if the current temperature is the base temperature:
            if self.temp.current == self.temp.min:
                new_acc_Ep = self.md.get_Ep_edr()
            #if the last MD run was the first at this temperature
            elif n == 1:
                # Get average potential energy from last MD run
                new_acc_Ep = self.md.get_Ep_edr()
            else:
                # Get previously accumulated average potential energy
                # associated to current temperature.
                acc_Ep = self.temp.get_Ep(self.temp.current)
                # Get average potential energy from last MD run
                new_Ep = self.md.get_Ep_edr()
                # Get the number of times the current teperature has been
                # selected for an MD run
                n = self.temp.get_n(self.temp.current)
                # Accumulate average potential energy and
                # update accumulated Ep associated to current temperature
                new_acc_Ep = self.accumulate_energy(acc_Ep, new_Ep, n)
            self.temp.set_Ep(self.temp.current, new_acc_Ep)

            # if not the first MD run
            if step > self.steps_md:
                # Concatenate accumulated edr file with edr file from last MD run
                self.md.concat_edr((step-self.steps_md)*self.dt)
                # Concatenate accumulated xtc file with xtc file from last MD run               
                self.md.concat_xtc((step-self.steps_md)*self.dt)

            # Chose next upper or lower temperature.
            # up_down = -1 --> test next lower temp
            # up_down = 1 --> test next upper temp

            # If we're curently at the base temperature:
            if self.temp.current == self.temp.min:
                up_down = 1
            # If we're currently at the peak temperature:
            elif self.temp.current == self.temp.max:
                up_down = -1
            else:
                # randomly choose -1 or 1
                up_down = self.lower_or_upper()

            # Assign the chosen temperature to be tested
            self.temp.update_next(up_down)
            # Recalculate weights for current and next upper temperature            
            self.temp.update_weight(self.temp.current)
            self.temp.update_weight(self.temp.current + self.temp.interval)

            # Save ST stats for this MD run if step is multiple of steps_save
            if step % self.steps_save == 0:
                self.write_stats(step)

            # Perform temperature transition test if
            # step is multiple of steps_trans: 
            if step % self.steps_trans == 0:
                transition = self.temp.test_transition(self.temp.current,
                                                       self.temp.next)
               # if transition test returned True           
                if transition:
                    # Transition to next temperature
                    self.temp.transition()
                    # Convert velocities in gro file produced by last MD run
                    # so they correspond to the new temperature.
                    self.md.conv_velocity(self.temp.conv_coef)

            # Delete gromacs backup files
            self.md.clear_backups()

        # copy last concatenated simulation files
        md.copy_md_output()

    def lower_or_upper(self):
        """Returns -1 or 1 with a 50/50 chance. 
        
        This function is used in order to determin wheather we will attempt
        a temperature transition towards the upper neighbouring temperature (1)
        or towards the lower neighbouring temperature (-1).
        """
        return np.random.choice([-1, 1])

    def accumulate_energy(self, acc_E, new_E, n):
        """ Accumulate average potential energy by merging the newly calculated
            average potential energy (from the last MD run) with the previously
            accumulated average potential energy.    
        ___________
        Attributes:
        -----------
        acc_E (float):  previously accumulated average potential energy
        new_acc (float): newly calculated average potential energy
        n (int): number of average potienal energy values (number of MD runs)
        
        _______
        Return:
        -------
        newly accumulated average potential energy (float)
        """
        return acc_E + (new_E - acc_E)/n

    def create_stats_file(self):
        """ Creates a file that will contain stats of the ST ouput (stats.dat).
        
        This simply created the file and writes the header. 
        Header has: step, time, MD temperature, all temperature weights
        """
        f = open('stats.dat', 'w')
        f.write(('{:>7s}{:>10s}{:>9s}{:>12s}').format('steps', 'time(ns)',
                                                      'temp(K)', 'energy(KJ)'))
        for t in self.temp_list:
            f.write('{:>12s}'.format('weight_' + str(t)))
        f.write('\n')
        f.close()

    def write_stats(self, step):
        """ Writes a line of stats in the ST output stats file (stats.dat)
            for a given step.
        ___________
        Attributes:
        -----------
        step (int): current step
        """
        f = open('stats.dat', 'a')
        f.write(('{:>7d}{:10.3f}{:>9d}{:>12.4g}').format(step,
                                                step*self.dt/1000,
                                                self.temp.current,
                                                self.temp.get_Ep(self.temp.current)))

        for t in self.temp_list:
            f.write(('{:>12.4g}').format(self.temp.get_weight(t)))
        f.write('\n')
        f.close()
