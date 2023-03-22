#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 12:06:44 2023

@author: thibaut
"""
from readAndPlot import *

class postProcessingFOAM (Reading,Plot):
    
    def __init__ (self):
        Reading.__init__();
        Plot.__init__();
        pass;                
########################### RETURN FIELD (u,Q,vort...)   ###################### 

    def maxVelocity(u):
        uMax      = max(u)     
        return  uMax;
    
    def xyzCoordinate(data,indexX,indexY,indexZ):
        x        = data[:,indexX];       
        y        = data[:,indexY];
        z        = data[:,indexZ];
        return  x,y,z;
    
    def velocityProfil(data,indexX,indexY,indexZ):
        uX        = data[:,indexX];       
        uY        = data[:,indexY];
        uZ        = data[:,indexZ];
        return  uX,uY,uZ;
    
    def vorticityProfil(data,indexX,indexY,indexZ):
        wX        = data[:,indexX];       
        wY        = data[:,indexY];
        wZ        = data[:,indexZ];
        return  wX,wY,wZ;
    
    def criterionQ(data,index):       
        Q         = data[:,index];       
        return Q;
    
    
###############################################################################    


########################## ERROR ##############################################    
    
    def error (path,fileName,
              delimiter,listFolder,D,
              saveName,indexX,indexY,indexZ):     
        
    
            for k in listData: 
                data      = Reading.readingData(path, fileName, delimiter);            
                x         = data[:,indexX];
                z         = data[:,indexZ];
                uz        = data[:,indexUZ];
                error[k]  = np.mean(abs(uz[:,k]-uz[:,k-1]));
                k=k+1;
            return error;

###############################################################################    

######################### BOUNDARY LAYER + RECIRC HEIGHT ######################       
    def boundaryLayerLength(z,L,ur):               
        filmThickness=0.;
        urMax            = max(ur); # valeur max de ur pour z !=0,L
        indexTemp        = np.argwhere(ur == urMax); 
        zTemp2           = z[indexTemp]; # recupere z ou ur==urMax        
        print(zTemp2);
        filmThickness    = abs(L-max(zTemp2))/L; # origine en L 
            
        return filmThickness;
                 
    def recirculationHeight(Q,z,criterion):     
        
        indexQ     = np.argwhere( (Q > criterion) ); # Critère Q > criterion => utiliser un critére Q en fonction de Re
        zTemp      = z[indexQ[:,0]]; # recupere les z correspondant
        if len(zTemp) != 0:
            recircLength = max(zTemp) - min(zTemp);
        else:
            recircLength = 0.;  
               
        return recircLength;
       

###############################################################################    

###################   CHANGING ORIGIN + JET CENTER POSITION ###################
    def changingOrigin(l,r,R):   # changement d'origine de r
        l           = R*(l/2 +1) # espacement entre chaque jet
        r_orig      = r - l;     # calcul le rayon par rapport à la nouvelle origine => en 2*L
        return r_orig;

    def jetCenterPosition(l,R):
        return l*R;

###############################################################################                
        

