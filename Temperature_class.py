# -*- coding: utf-8 -*-
"""
Created on Sat Sep 17 22:44:15 2016

@author: gilles
"""
import numpy as np
import math
import re

#Boltzmann constant: KJ/(K.mol)
KB = (0.0083144621)


class Temperature:
    """Temperature class used for simulated tempering
        _____________________________
        Arguments to create instance:
        -----------------------------
        No arguments needed. 
        Constructor takes parameters from the ST.param file. 
        
        The parameters taken from ST.param are:
        --> temp_min : the base temperature of the simulation (K)
        --> temp_max : the peak temperature of the simulation (K)
        --> temp_interval : difference between two neighbouring temperatures
        ___________
        Attributes:
        -----------
        min (int)        : base ST temperature
        max (int)        : peak ST temperature
        interval (int)   : difference between two neighbouring temperatures
        current (int)    : current temperature used in the ST
        next (int)       : next temperature to be tested for a temperature transition
        conv_coef (float): velocity conversion factor used to convert a .gro files
        temp_dict        : Global temperature dictionary which contains a 
                           subdictionary associated to each temperature. 
                           Each of these subdictionaries contains:
                             --> beta: 1/KB*temp (KB is Boltzmann constant)
                             --> weight: temperature weight
                             --> Ep: accumulated potential energy
                             --> n: number of simulation runs done at this temperature
    """

    def __init__(self):
        self.min, \
        self.max, \
        self.interval, \
        self.temp_dict = self.init_temperatures_from_file()
        self.current = self.min
        self.next = self.min + self.interval
        self.conv_coef = (self.next/self.current)**0.5

    def init_temperatures_from_file(self):
        """ Initializes a temperature dictionary from file
            __________
            Arguments:
            ----------
            No arguments. Takes parameters from ST.param file.
            
            The parameters taken from ST.param are:
            --> temp_min : the base temperature of the simulation (K)
            --> temp_max : the peak temperature of the simulation (K)
            --> temp_interval : difference between two neighbouring temperatures
            
            _______
            Return:
            -------
            Global temperature dictionary which contains a subdictionary associated
            to each temperature. 
            
            Each of these subdictionaries contains:
                    --> beta: 1/KB*temp (KB is Boltzmann constant)
                    --> weight: temperature weight
                    --> Ep: accumulated average potential energy
                    --> n: number of simulation runs done at this temperature
        """
        # Initialize variables with -1
        min, max, interval = [-1]*3
        # Get temperatures from parameters file
        f_param = open('ST.par', 'r')
        for line in f_param:
            if line.startswith('temp_min'):
                min = int(re.findall("\d+", line)[0])
            elif line.startswith('temp_max'):
                max = int(re.findall("\d+", line)[0])
            elif line.startswith('temp_interval'):
                interval = int(re.findall("\d+", line)[0])
        f_param.close()
        # TODO: throw exception if variables are not found in file 
        temp_dict = self.init_temperatures(min, max, interval)
        return min, max, interval, temp_dict

    def init_temperatures(self, min, max, interval):
        """ Initializes a temperature dictionary from arguments
            __________
            Arguments:
            ----------
            min (int)        : base ST temperature
            max (int)        : peak ST temperature
            interval (int)   : difference between two neighbouring temperatures
            _______
            Return:
            -------
            Global temperature dictionary which contains a subdictionary associated
            to each temperature. 
            
            Each of these subdictionaries contains:
                    --> beta: 1/KB*temp (KB is Boltzmann constant)
                    --> weight: temperature weight
                    --> Ep: accumulated average potential energy
                    --> n: number of simulation runs done at this temperature
        """
        # Temperature dictionary for the starting temperature (min) :
        # Dictionary containing beta, average potential energy, weight,
        # and number of simulations for a given temperature (min in this case).
        # Starting temperature has weight = zero
        t_dict = {'beta': 1/(KB*min),
                  'weight': 0,
                  'Ep': float('NaN'),
                  'n': 0}
        # Adding the first temperature dictionary to the global dictionary t 
        # that will contain the dictionaries of all temperatures 
        # (min to max with a given interval bewteen neighbouring temperatures).
        t = {min: t_dict}

        # Populating the global temperature dictionary:        
        for temp in xrange(min+interval, max+interval, interval):
            # Dictionary containing beta, average potiential energy, weight,
            # and number of simulations for a given temperature, respectively.
            # Starting temperature weight was set to zero.
            # The other temperature weights are set as Not a Number for now.
            t_dict = {'beta': 1/(KB*temp),
                      'weight': 0,
                      'Ep': float('NaN'),
                      'n': 0}
            t[temp] = t_dict

        return t

    def get_beta(self, temp):
        """ Get a given temperature's beta value from the temp_dict attribute
            __________
            Arguments:
            ----------
            temp (int): a temperature
            _______
            Return:
            -------
            beta value associated to given temperature (float)
        """
        return self.temp_dict[temp]['beta']

    def set_beta(self, temp, beta):
        self.temp_dict[temp]['beta'] = beta
        """ Set a given temperature's beta value in the temp_dict attribute
            __________
            Arguments:
            ----------
            temp (int): a temperature
            beta (float): a beta value
        """

    def get_weight(self, temp):
        """ Get a given temperature's weight from the temp_dict attribute
            __________
            Arguments:
            ----------
            temp (int): a temperature
            _______
            Return:
            -------
            weight value associated to given temperature (float)
        """
        return self.temp_dict[temp]['weight']

    def set_weight(self, temp, weight):
        """ Set a given temperature's weight in the temp_dict attribute
            __________
            Arguments:
            ----------
            temp (int): a temperature
            weight (float): temperature's weight
        """
        self.temp_dict[temp]['weight'] = weight

    def get_Ep(self, temp):
        """ Get a given temperature's accumulated average potential energy 
            from the temp_dict attribute
            __________
            Arguments:
            ----------
            temp (int): a temperature
            _______
            Return:
            -------
            Accumulated average potential energy associated to given temperature (float)
        """
        return self.temp_dict[temp]['Ep']

    def set_Ep(self, temp, Ep):
        """ Set a given temperature's accumulated average potential energy
            in the temp_dict attribute
            __________
            Arguments:
            ----------
            temp (int): a temperature
            Ep (float): temperature's accumulated average potential energy
        """
        self.temp_dict[temp]['Ep'] = Ep
        
    def get_n(self, temp):
        """ Get a given temperature's number of simulation runs
            __________
            Arguments:
            ----------
            temp (int): a temperature
            _______
            Return:
            -------
            number of simulation runs associated to given temperature (int)
        """   
        return self.temp_dict[temp]['n']
        
    def set_n(self, temp, n):
        """ Set a given temperature's number of simulation runs 
            in the temp_dict attribute
            __________
            Arguments:
            ----------
            temp (int): a temperature
            n (int): temperature's number of simulation runs
        """
        self.temp_dict[temp]['n'] = n

    def update_next(self, up_down):
        """ Update next temperature to be tested for a temperature transition
            (i.e. the 'next' attribute)
            __________
            Arguments:
            ----------
            up_down (int): can be -1 or 1
            
                if up_down = -1 -->  chose neighbouring lower temperature
                if up_down = 1 --> chose neighbouring upper temperature
        """    
        self.next = self.current + up_down*self.interval

    def update_weight(self, temp2):
        """ Calculate weight for a given temperature
            and set new weight in the temp_dict attribute
            __________
            Arguments:
            ----------
            temp2 (int) : temperature for which we want to calculate the weight
            _______
            Method:
            -------
            Finds weight of temperature n using:
                  ---> beta of temperature n and n-1
                  ---> accumulated average potential energy of temperature n and n-1
                  ---> weight of temperature n-1
            
            According to the formula (Nguyen 2013):
                new_weight2 = weight1 + (b2 - b1)*(Ep2 + Ep1)/2
                
            if temperature n does not yet have a potential energy, use:
                new_weight2 = (b2 - b1)*Ep1/2
        """
        # if temp is base temperature, weight is zero.
        if temp2 == self.min:
            self.set_weight(temp2, 0)
            return
        # if next temperature is withing temperature bounds:
        elif self.min < temp2 <= self.max:
            #continue 
            temp1 = temp2 - self.interval
        #if stepping outside temperature bounds
        else:
            return
        
        weight1 = self.get_weight(temp1)
        b1 = self.get_beta(temp1)
        b2 = self.get_beta(temp2)
        Ep1 = self.get_Ep(temp1)
        Ep2 = self.get_Ep(temp2)
        #if energy for temp2 has not been calculated:
        if math.isnan(Ep2):
            new_weight2 = (b2 - b1)*Ep1/2
        else:
            new_weight2 = weight1 + (b2 - b1)*(Ep2 + Ep1)/2

        self.set_weight(temp2, new_weight2)

    def transition(self):
        """ Applies the temperature transition.
        
            Prepares temp_current attribute for the next md.
            Also recalculates the velocity conversion vactor (conv_coef)
        """
        self.conv_coef = (float(self.next)/self.current)**0.5
        self.current = self.next

    def test_transition(self, temp1, temp2):
        """ Applies the transition test to see if we will change temperatures
            for the next MD run. 
            
            Uses the transition probability.
            _______
            Return:
            -------
            (bool) True if test is successful, False if test failed. 
        """
        p = self.transition_proba(temp1, temp2)
        return np.random.choice([True, False], 1, p=[p, 1-p])[0]

    def transition_proba(self, temp1, temp2):
        """ Caclculates the temperature transition probability from temp1 to temp2.
           
            Uses:
                ---> beta of temperature 1 and 2
                ---> accumulated average potential energy of temperature 1 and 2
                ---> weight of temperature 1
            
            According to the formula (Nguyen 2013):
                proba = exp(-(((beta2 - beta1)*Ep1 - (weight2 - weight1))))
            _______
            Return:
            -------
            (float) temperature transition probability
        """
        b1 = self.get_beta(temp1)
        b2 = self.get_beta(temp2)
        w1 = self.get_weight(temp1)
        w2 = self.get_weight(temp2)
        Ep1 = self.get_Ep(temp1)
        
        try:
            proba = math.exp(-(((b2 - b1)*Ep1 - (w2 - w1))))
            print "######## transition proba: " + str(proba)
        except OverflowError:
            print "######## transition proba: 1"
            if -(((b2 - b1)*Ep1 - (w2 - w1))) > 0:
                proba = 1
            else:
                proba = 0

#        # Write probabilites in file  
#        # TEST: REMOVE LATER
#        f = open('proba.dat', 'a') # TEST TO BE REMOVED
#        f.write(str(proba) + '\n') # TEST TO REMOVE    
#        f.close() # TEST TO REMOVE
        return min(1, proba)