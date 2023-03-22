#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 14:17:29 2023

@author: thibaut
"""

import scipy.interpolate
import numpy as np
from readAndPlot import *

class NumericalStuff:
    def __init__(self):
        pass;
    def dichotomieIntervall(z,u,epsilon):
        # dichotomie with changing interval => we use 1D interpolation to fine the zeros of a function
        # aIni => initial lowest value of inerval
        # bIni => initial highest value of inerval
        #         
        # length => number of sub-interval
        # u      => dichitomie u = 0 : 
        a=0;
        b=0;
        m=0;
        ym=0;
        res = [];
        for i in range(len(u)-1):
            if ( (u[i]*u[i+1] <= 0) and (i > 0 and i+1 < (len(u)) ) ):  
                # Changement de signe  => 0 se situe entre ces deux valeurs => On ne regarde par i=0 et i+1=len(u) => z=0 et z=L car ils sont déjà nulle
                a  = z[i];
                b  = z[i+1];
                m  = (a+b)/2.          # forme une variable milieu => a+b/2
                ym = np.interp(m,z,u)  # Interpole les données : u(m)              
                
                while ( ym > epsilon ):                    
                    # cherche le 0 par interpolation:
                    if(u[i]*ym < 0): # si ui*ym < 0 => le zeros est + proche de u_i
                        b = m;
                    else:            # si ce n'est pas le cas => le zeros est + proche de u_i+1
                        a = m;   
                    # recalcul m et ym                        
                    m = (a+b)/2.
                    ym = np.interp(m,z,u)  # Interpole les données : u(m)  
            
                res.append(m);
          
        
        return np.array(res);
              
        
        
        
        
        
        
        
        
        
                    
                    
                    
            