# -*- coding: utf-8 -*-
"""
Created on Sat Sep 17 22:44:15 2016

@author: gilles
"""

import fileinput
import os.path
from subprocess import call, Popen, PIPE
from shutil import copyfile
import shlex


class MD:
    """ MD class used to start a run of molecular dynamics simulation in gromacs (v. 5.0.5)
        _____________________________
        Arguments to create instance:
        -----------------------------
        mdp (string): .mdp file name
        gro (string): .gro file name
        top (string): .top file name
        steps (int): number of steps in the MD run
        dt (float): time between two steps of the MD run
        temp_list (list of int): list of all temperatures used in an ST
        ___________
        Attributes:
        -----------
        mdp (string): .mdp file name
        gro (string): .gro file name
        top (string): .top file name
        
    """
    def __init__(self, mdp, gro, top, steps, dt, temp_list):
        # Filenames:
        self.mdp = mdp
        self.gro = gro
        self.top = top
        # Generate mdp files (one for each temperature)
        self.init_mdp_files(steps, dt, temp_list)

    def init_mdp_files(self, steps, dt, temp_list):
        """Creates different mdp files for different temperatures of an ST. 
            ___________
            Attributes:
            -----------
            steps (int): number of steps in the MD run
            dt (float): time between two steps of the MD run
            temp_list (list of int): list of all temperatures used in an ST        
        """
        self.init_mdp_nsteps(steps)
        self.init_mdp_dt(dt)
        for temp in temp_list:
            self.init_mdp_temp(temp)
            copyfile(self.mdp, str(temp) + '.mdp')

    def init_mdp_nsteps(self, n):
        """ Adds the number of steps for an MD in the main mdp file
            (The main mdp file will then be duplicated with different temperatures.)
            ___________
            Attributes:
            -----------
            n (int): number of steps in the MD run
        """
        start_of_line = 'nsteps'
        newline = 'nsteps = ' + str(n)
        self.line_replace(self.mdp, start_of_line, newline)

    def init_mdp_dt(self, dt):
        """ Adds the time between two MD steps (ps) in the main mdp file
            (The main mdp file will then be duplicated with different temperatures.)
            ___________
            Attributes:
            -----------
            dt (float): time between two MD steps (ps)
        """
        start_of_line = 'dt'
        newline = 'dt = ' + str(dt)
        self.line_replace(self.mdp, start_of_line, newline)

    def init_mdp_temp(self, temp):
        """ Adds the the temperature (K) in the mdp file
            ___________
            Attributes:
            -----------
            temp (int): MD temperature (K)
        """
        start_of_line = 'ref_t'
        newline = 'ref_t = ' + str(temp)
        self.line_replace(self.mdp, start_of_line, newline)

    def line_replace(self, file, start_of_line, newline):
        """ In a given file, replaces a line starting with a given string with a new line. 
            ___________
            Attributes:
            -----------
            file (string): file name
            start_of_line (string)
            newline (string)
        """
        for line in fileinput.input(file, inplace=True): 
            if line.startswith(start_of_line):         
                print newline
            else:
                print line.strip()
        fileinput.close()

    def start(self, temp):
        """ Starts an MD run in gromacs with a given temperature. 
            ___________
            Attributes:
            -----------
            temp (int): MD temperature (K)
        """
        # reassign mdp filename corresponding to right temperature
        self.mdp = str(temp) + '.mdp'
        # if edr file exists, create a dublicate (for concatenation later)
        # because original edr file will be overriden by gromacs MD run
        extention = '.edr'        
        self.check_exists(extention)
        # if xtc file exists, create a dublicate (for concatenation later)
        # because original xtc file will be overriden by gromacs MD run
        extention = '.xtc'
        self.check_exists(extention)        
        # Start MD: gromacs shell commands:
        call('gmx grompp -f ' + self.mdp + ' -c ' \
             + self.gro + ' -p ' + self.top + ' -o ' \
             + self.gro.split('.')[0] + '.tpr -maxwarn 5', shell=True)
        call('export GMX_MAXCONSTRWARN=-1', shell=True)
        call('gmx mdrun -v -deffnm ' + self.gro.split('.')[0], shell=True)

    def check_exists(self, extension):
        """ checks if an md file with the same name as the .gro file but with a
            given extension exists in the current directory. 
            ___________
            Attributes:
            -----------
            extension (string): file extention e.g. '.edr'
        """
        filename = self.gro.split('.')[0] + extension      
        # if file exists        
        if os.path.exists(filename):
            # dublicate file 
            filename_copy = self.gro.split('.')[0] + '_total' + extension
            copyfile(filename, filename_copy)
            
    def copy_md_output():
        """ Copy most recent .edr and .xtc file as _total.edr and _total.xtc
        """
        old = self.gro.split('.')[0] + '.edr'
        new = self.gro.split('.')[0] + '_total.edr'
        copyfile(old, new)
        old = self.gro.split('.')[0] + '.xtc'
        new = self.gro.split('.')[0] + '_total.xtc'
        copyfile(old, new)
            
    def conv_velocity(self, coef):
        """ Convert volecities in .gro file so the correspond to the new temperature. 
            ___________
            Attributes:
            -----------
            coef (float): velocity conversion factor
        """
        new_gro = '' 
        # Parse .gro file:
        f_gro = open(self.gro, 'r')        
        for i, line in enumerate(f_gro):
            if i == 1:
                # number of atoms in gro file
                atoms = int(line)
            elif 1 < i < atoms + 2:
                # fetch velocity coordinates (last 3 items of line)
                items = line.split()
                coords = items[-3:]
                # convert string coordinates to float
                coords = [float(j) for j in coords]
                # apply temperature conversion factor
                coords = [j*coef for j in coords]
                # regenerate line with new velocity
                line = ('{:>8s}'
                        '{:>7.5s}'
                        '{:>5s}'
                        '{:>8.5s}'
                        '{:>8.5s}'
                        '{:>8.5s}'
                        '{:8.4f}'
                        '{:8.4f}'
                        '{:8.4f}\n').format(items[0],
                                            items[1],
                                            items[2],
                                            items[3],
                                            items[4],
                                            items[5],
                                            coords[0],
                                            coords[1],
                                            coords[2])
            new_gro += line
        f_gro.close()
        
        # Write changes in new gro file
        f_new_gro = open(self.gro, 'w')
        f_new_gro.write(new_gro)
        f_new_gro.close()
        
    def get_Ep_edr(self):
        """ Gets average energy potential from the newly generated .edr file.
            Uses the 'gmx energy' command (gromacs v.5.0.5).
        """
        # Get average potential energy from .edr file         
        # bash command with pipe: "echo 'potential' | gmx energy -f md2.edr"
        input_edr = self.gro.split('.')[0] + '.edr'      
        cmd1 = ("gmx energy -f " + input_edr)
        cmd2 = ("echo potential")
        grep = Popen(cmd1.split(), stdin=PIPE, stdout=PIPE)
        ls = Popen(cmd2.split(), stdout=grep.stdin)
        output = grep.communicate()[0]
        ls.wait()
        
        # Iterate through words in output to find potential energy (kJ/mol)
        for i, word in enumerate(output.split()):
            if word == 'Potential':
                # return word following "Potential" and convert it to float
                return float(output.split()[i+1]) * -1 # '* -1' to make positive

# TODO:                
#    def get_RMSD(ref, traj):
#        # RMSD between the trajectory structures and the original structure:
#        # gromacs shell command:
#        call ('gmx rms -f ' + traj + ' -s ' ref + ' -o RMSD_time.xvg')
#        
#        echo protein | gmx cluster -f md2.xtc -s reference.gro -dist RMSD_dist.xvg


    def concat_edr(self, t0):
        """ Concatenates newly genereted .edr file to the previously concatenated 
            .edr file. Uses the 'gmx eneconv' command (gromacs v.5.0.5). 
            ___________
            Attributes:
            -----------
            t0: start time of the newly generated .edr file
        """
        # Adjust start time for the new MD ouput edr file.      
        # And concatenate all .edr files in current directory.
        cmd1 = ("echo '0\n{0}\n'").format(str(t0))        
        input1 = self.gro.split('.')[0] + '_total.edr'
        input2 = self.gro.split('.')[0] + '.edr'
        tmp_output = 'tmp.edr'
        cmd2 = ("gmx eneconv -f {0} {1} -settime -o {2}").format(input1,
                                                                input2,
                                                                tmp_output)
        
        p1 = Popen(shlex.split(cmd1), stdout=PIPE)
        p2 = Popen(shlex.split(cmd2), stdin=p1.stdout)
        p2.wait()
        p1.stdout.close()
        
        copyfile(tmp_output, input2)
        call(['rm', tmp_output])
    
    def concat_xtc(self, t0):
        """ Concatenates newly genereted .xtc file to the previously concatenated 
            .xtc file. Uses the 'gmx trjcat' command (gromacs v.5.0.5). 
            ___________
            Attributes:
            -----------
            t0: start time of the newly generated .xtc file
        """
        # Adjust start time for the new MD ouput xtc file.      
        # And concatenate all .xtc files in current directory.
        cmd1 = ("echo '0\n{0}\n'").format(str(t0))        
        input1 = self.gro.split('.')[0] + '_total.xtc'
        input2 = self.gro.split('.')[0] + '.xtc'
        tmp_output = 'tmp.xtc'
        cmd2 = ("gmx trjcat -f {0} {1} -settime -o {2}").format(input1,
                                                                input2,
                                                                tmp_output)
        
        p1 = Popen(shlex.split(cmd1), stdout=PIPE)
        p2 = Popen(shlex.split(cmd2), stdin=p1.stdout)
        p2.wait()
        p1.stdout.close()
        
        copyfile(tmp_output, input2)
        call(['rm', tmp_output])
        
    def clear_backups(self):
        """ Deletes all gromacs backup files.
        """
        Popen('rm \#*\#', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
