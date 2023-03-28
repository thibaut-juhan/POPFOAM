#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 12:06:44 2023

@author: thibaut
"""
from readAndPlot import *
from numericalStuff import *
import scipy

class jetDynamics (Reading,Plot,NumericalStuff):
    
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

    def angleExpansion(uz,z,x,jetSpacing,R):
        
        coordExp          = np.zeros((len(uz[0]),2));
        jetCenter         = jetDynamics.jetCenterPosition(jetSpacing,R);     # center of the jet !  
        
        for l in range (len(uz[0])):             
            uzTemp   = uz[:,l];
            zTemp    = z[:,l];
            xTemp    = x[:,l];
            
            uzTemp   = uzTemp[   xTemp>= jetCenter];
            xTemp    = xTemp[xTemp>= jetCenter];
            
            xCenter  = min(xTemp);
            uzCenter = uzTemp[xTemp == xCenter ]; 
            
            # interpolation => uniquement sur l'extrémité ou il n'y a pas d'interaction avec le 2nd jet
            # interpolation at jetCenterPosition
            uzTemp   = uzTemp - 0.01*uzCenter;
            # recherhe de 0 => f(z)=uz-0.05uz(x=xcenter) => recherche le 0 sur cette fonction
            zerosUz = NumericalStuff.dichotomieIntervall(xTemp, uzTemp, 1.e-9,boundary=False);
            zerosUz = zerosUz[zerosUz > jetCenter]; 
            # recuperer seulement zeros x>xcenter
            if(len(zerosUz)!=0):
                coordExp[l,0]= zerosUz[0];
                coordExp[l,1]= zTemp[0];
            
                
            # recupere uniquement les x > xcenter => on se place sur l'extrémité droit du jet
        
        coordExp = np.delete(coordExp,np.where( coordExp==0 ),axis=0);
        
        xMaxCoordExp = max(coordExp[0:20,0]);
        zMaxCoordExp = max(coordExp[0:20,1]);
        
        xMinCoordExp = min(coordExp[0:20,0]);
        zMinCoordExp = min(coordExp[0:20,1]);
        
        plt.plot(coordExp[0:20,0],coordExp[0:20,1],'o-')
        
        H1              = np.sqrt( (xMaxCoordExp -xMinCoordExp)**2  + (zMaxCoordExp -zMinCoordExp)**2  );
        H2              = abs(zMaxCoordExp -zMinCoordExp);   
        angleExpansion = np.arccos(H2/H1);
        
        return angleExpansion
    
    def jetCoreLength(uz,z,x,jetSpacing,R,lenList,U0):
        uzCenterLine = np.zeros((len(uz[0]),2)); # stock : uz at each slice + corresponding z
        jetCenter    = jetDynamics.jetCenterPosition(jetSpacing,R);     
        for l in range (len(uz[0])):                
            uzTemp   = uz[:,l];
            zTemp    = z[:,l];
            xTemp    = x[:,l];
            f        = scipy.interpolate.interp1d(xTemp,uzTemp); # interpolation
            uzInterp = f(jetCenter);                             # interpolation at jetCenterPosition  
            
            uzCenterLine[l,1]=uzInterp;
            uzCenterLine[l,0]=zTemp[0]; # same z for each zTemp
        
        # realize a dichotomy on uz-0.95uz
        z095 = NumericalStuff.dichotomieIntervall(uzCenterLine[:,0],(uzCenterLine[:,1]-0.95*U0),1e-9,boundary=False);
    
        
        return z095;
###############################################################################                
        

