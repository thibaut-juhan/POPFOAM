#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 14:31:56 2023

@author: thibaut
"""

from postProcessingFOAM import *
from numericalStuff     import *
from scipy import interpolate

path        =  "../dataPostProcessingFOAM/"
fileVel     =  "line_U:Transformed.xy"
fileVelVort =  "line_U:Transformed_vorticity.xy"
fileQ       =  "Q.xy"

R = 5e-4;
L = 3e-3;

listReFolder = os.listdir(path);
listSFolder  = os.listdir( path + listReFolder[0]);

listSFolder  = Reading.sortedFolder(path + listReFolder[0] , listSFolder , "S" , "S" ,1 , len(listSFolder) + 1 );
listReFolder = Reading.sortedFolder(path + listReFolder[0] , listReFolder , "Re_" , "Re_" ,3, len(listSFolder) + 1 );

lenListRe = len(listReFolder); 
lenListS  = len(listSFolder);
count=0;
countRe=0;
epsi=1e-4;

with open("file/inletVelocity","rb") as file:
    inletVelocity = np.loadtxt(file);

############################## POST PROCESS MAIN LOOPS ########################
for i in listReFolder[0:-1]:    
    print(i);
    #### TEST ON OPEN FILE : If it already exists => delete it
    Reading.openingFileTest("file/meanFieldBottomWall_" + str(i) + ".dat" );
    Reading.openingFileTest("file/combinaisonPoint_"    + str(i) + ".dat" );
    Reading.openingFileTest("file/filmThickness_"       + str(i) + ".dat" );
    Reading.openingFileTest("file/mergingPoint_"        + str(i)  + ".dat");
    Reading.openingFileTest("file/fountainFlow_"        + str(i)  + ".dat");
    
    for j in listSFolder:    
        countS=0;
############################### LISTES DE DONNEES ############################# 
        
        pathDataMeanFieldConcentration = path + "/" + str(i) + "/" + str(j) + "/" + "concentrationMeanField/";
        pathDataSliceXZAxial           = path + "/" + str(i) + "/" + str(j) + "/" + "dataAxialSliceXZ/";
        pathDataSliceXZRadial          = path + "/" + str(i) + "/" + str(j) + "/" + "dataRadialSLiceXZ/";
        pathDataCenterJet              = path + "/" + str(i) + "/" + str(j) + "/" + "dataCenterJet/";
        
        listFileMeanConc     = os.listdir(pathDataMeanFieldConcentration); 
        listFileSliceAxial   = os.listdir(pathDataSliceXZAxial); 
        listFileSliceRadial  = os.listdir(pathDataSliceXZRadial);
        listFileCenterJet    = os.listdir(pathDataCenterJet);
       
        listFileSliceAxial  = Reading.sortedFolder(pathDataSliceXZAxial  , listFileSliceAxial  ,"h","singleGraph",10 , 13 );
        listFileSliceRadial = Reading.sortedFolder(pathDataSliceXZRadial , listFileSliceRadial ,"h","singleGraph",10 , 13 );
         
        countS = countS+1;
        lamb = Reading.computeLambda(R, np.float_(j[1:len(j)+1]) )
        
###############################################################################        
        
########################## POST-PROCESS AXIAL DATA#############################
       
        timeFolder = os.listdir(pathDataSliceXZAxial + '/' + str(listFileSliceAxial[0]) + "/" );
        time       = timeFolder[-2];
        data       =  Reading.readingData(pathDataSliceXZAxial + '/' + str(listFileSliceAxial[0]) + "/" +str(time) + "/" , fileVel, ';')
        
        saveX  = np.zeros( (len(data[:,0]) , len(listFileSliceAxial) ) );
        saveZ  = saveX.copy();
        saveUx = saveX.copy();
        saveUz = saveX.copy();
        
        
        for k in listFileSliceAxial: # loop over folder in Re*/S*/data/*
            # loop over runTime data => 
            #print(count);
            timeFolder = os.listdir(pathDataSliceXZAxial + '/' + str(k) + "/" );
            time       = timeFolder[-2];
            data       =  Reading.readingData(pathDataSliceXZAxial + '/' + str(k) + "/" +str(time) + "/" , fileVel, ';')
            
            x,y,z      =  postProcessingFOAM.xyzCoordinate (data,0, 1, 2);
            ux,uy,uz   =  postProcessingFOAM.velocityProfil(data,3 ,4 ,5);
            
            saveX [:,count]  = x;
            saveZ [:,count]  = z;
            saveUx[:,count]  = ux;
            saveUz[:,count]  = uz;
            count = count+1;
        
        coordExp          = np.zeros(( len(listFileSliceAxial),2));
        jetCenterPosition = postProcessingFOAM.jetCenterPosition(np.float_(j[1:len(j)+1]),R);
           
        # calcul la position du centre d'un jet
        
        for l in range (len(saveUz[0])):        
            
            uzTemp  = saveUz[:,l];
            zTemp   = saveZ[:,l];
            xTemp   = saveX[:,l];
            plt.figure();
            plt.plot(xTemp,uzTemp);
            f = scipy.interpolate.interp1d(xTemp,uzTemp);
            uzInterp = f(jetCenterPosition);
            uzTemp  = uzTemp - 0.05*uzInterp; 
            # recherhe de 0 => f(z)=uz-0.95uz(x=xcenter) => recherche le 0 sur cette fonction
            zerosUz = NumericalStuff.dichotomieIntervall(xTemp, uzTemp, 1.e-8);
            zerosUz = zerosUz[zerosUz > jetCenterPosition];
         
            if(len(zerosUz)!=0):
                coordExp[l,0]= zerosUz[0];
                coordExp[l,1]= zTemp[0];
                
            # recupere uniquement les x > xcenter => on se place sur l'extrémité droit du jet
        
        print(coordExp);
        coordExp = np.delete(coordExp,np.where(coordExp==0),axis=0);
        print(coordExp);
    
        
        xMaxCoordExp = max(coordExp[:,0]);
        zMaxCoordExp = max(coordExp[:,1]);
        
        xMinCoordExp = min(coordExp[:,0]);
        zMinCoordExp = min(coordExp[:,1]);
        
        H              = np.sqrt( (xMaxCoordExp -xMinCoordExp)**2  + (zMaxCoordExp -zMinCoordExp)**2  );
        print(H);
        angleExpansion = np.arccos(H/L);
        
        print(angleExpansion);
        
        
        Plot.quiverJetsProfile(
            saveX, saveZ, saveUx, saveUz, 
            inletVelocity[int(countRe)],R, 
            -4, 4, 
            0, 7, 
            r'$\frac{x}{R_{jet}}$', r'$\frac{z}{R_{jet}}$', 
            r'$\frac{U_z}{U_0}$', 6, 
            str(lamb) + r" $Re = $" + str(i[3:len(i)+1])
            ,"plot/" + str(i) + str(lamb) );            
        
        count=0; #reset
        
        # define the jet expansion angle
        
        
        del saveX,saveZ,saveUz,saveUx;
        
###############################################################################       

########################## POST-PROCESS RADIAL DATA############################
        
        filmThickness = np.zeros(len(listFileSliceRadial))
        r             = np.zeros(len(listFileSliceRadial))      
        
        for k in listFileSliceRadial:     
            
            timeFolder = sorted(os.listdir(pathDataSliceXZRadial + '/' + str(k) + "/" ));
            time       = timeFolder[-2];
            
            data                 =  Reading.readingData(pathDataSliceXZRadial + '/' + str(k) + "/" + str(time) + "/" , fileVelVort, ';');
            x,y,z                =  postProcessingFOAM.xyzCoordinate  (data,0, 1, 2);
            ux,uy,uz             =  postProcessingFOAM.velocityProfil (data,3 ,4 ,5);
            filmThickness[count] =  postProcessingFOAM.boundaryLayerLength(z, L, ux);
            r[count]=x[0];
            count = count+1;
            
        count      = 0; #reset
        r          = postProcessingFOAM.changingOrigin(lamb,r,R/2.);
        rChanging  = r[r>0]; # changing origin
        filmThickness = filmThickness[r>0];
        
        dataBL = np.zeros( (len(filmThickness),3) );
        dataBL[:,0] = lamb*np.ones(len(filmThickness));
        dataBL[:,1] = rChanging;
        dataBL[:,2] = filmThickness;
       
###############################################################################   


########################## POST-PROCESS CENTER JETDATA#########################      
 
        timeFolder = os.listdir(pathDataCenterJet + '/' + str( listFileCenterJet[0]) + "/" );
        time       = timeFolder[-2];
        data       =  Reading.readingData(pathDataCenterJet + '/' + str(listFileCenterJet[0]) + "/" + str(time) + "/", fileVel, ';')
    
        x,y,z      =  postProcessingFOAM.xyzCoordinate  (data,0, 1, 2);
        ux,uy,uz   =  postProcessingFOAM.velocityProfil (data,3 ,4 ,5);
    
        maxUz              = postProcessingFOAM.maxVelocity(uz,);
        recombPoint        = z[uz == maxUz ];
    
        if recombPoint==L:
            recombPoint=0;
        
        dataComb        = np.zeros ( (1,2) );
        
        dataComb[0,0]   = lamb;
        dataComb[0,1]   = recombPoint;
        
        
        zerosUz = NumericalStuff.dichotomieIntervall(z,uz,1.e-9);
        
        dataMP = np.zeros((1,2));
        dataFF = np.zeros((1,2));
        
        #plt.plot(z,uz);
        
        if len(zerosUz) != 0: # fountain flow and merging point 
            
            mergingPoint = min(zerosUz);
            fountainFlow = max(zerosUz);
            
            dataMP[0,0] = lamb;
            dataMP[0,1] = mergingPoint;
            
            if mergingPoint == fountainFlow:
                
                dataFF[0,0] = lamb;
                dataFF[0,1] = 0;                
            else:
                dataFF[0,0] = lamb;
                dataFF[0,1] = fountainFlow;              
            
        else:  # pas de zone fountaine / zone de recombinaison 

            dataMP[0,0] = lamb; 
            dataMP[0,1] = 0; # No merging point       

            dataFF[0,0] = lamb;
            dataFF[0,1] = 0; # No fountainFlow  
            
        
###############################################################################  

########################## POST-PROCESS CONCENTRATION FIELD ###################
       
        timeFolder = os.listdir(pathDataMeanFieldConcentration + "/" + 'patchAverage(name=botWall,s)/' )
        data       = Reading.readingDataConcenreation(pathDataMeanFieldConcentration + "/" + 'patchAverage(name=botWall,s)/' + str(timeFolder[-1]) , "/surfaceFieldValue.dat");
    
        t = data[:,0];
        s = data[:,1];
        
        dataWrite = np.zeros( (len(t),3) );
        dataWrite[:,0]  = lamb*np.ones(len(t));
        dataWrite[:,1]  = t;
        dataWrite[:,2]  = s;
        
        hMean = Reading.addHeader("file/meanFieldBottomWall_"+str(i)+".dat",
                          " Re , t , mean concentration on the bot wall");
        
        hComb = Reading.addHeader("file/combinaisonPoint_" + str(i) + ".dat", 
                          "lambda, combinaison point")
        
        hMerge = Reading.addHeader("file/mergingPoint_" + str(i) + ".dat",
                          "merging point");
        
        hFountain = Reading.addHeader("file/fountainFlow_" + str(i) + ".dat"
                          ,"fountainFlow");
        
        hFilm     = Reading.addHeader("file/filmThickness_" + str(i) + ".dat"
                          ,"lamda,boundary layer thickness , radial position");
        
        ########   WRITING DATA  ######## 
        
        with open ( "file/meanFieldBottomWall_"+str(i)+".dat","ab") as file:
            np.savetxt(file,dataWrite,delimiter = " ",header=hMean);
        
        with open ( "file/combinaisonPoint_" + str(i) + ".dat","ab") as file:
            np.savetxt(file,dataComb,delimiter = " ",header=hComb);
            
        with open ( "file/filmThickness_" + str(i) + ".dat","ab") as file:
            np.savetxt(file,dataBL,delimiter = " ",header=hFilm);
                   
        with open ("file/mergingPoint_"      + str(i)  + ".dat","ab") as file:
            np.savetxt(file,dataMP,delimiter = " ",header = hMerge);
        
        with open ("file/fountainFlow_"      + str(i)  + ".dat","ab") as file:
            np.savetxt(file,dataFF,delimiter = " ",header = hFountain);
        
      
###############################################################################  

    countRe=countRe+1;
    

############################ PLOT ############################################

with open("file/ReynoldsNumber","rb") as file:
     Re= np.loadtxt(file);    
 
    
norm = plt.Normalize(np.min(Re), np.max(Re))
cmap = plt.get_cmap('jet')
c = cmap(norm(Re))
 
count=0;
plt.figure(figsize = (12,8) );
for i in listReFolder[0:-1]:
    with open("file/combinaisonPoint_" + str(i) + ".dat","rb") as file:
        dataMean = np.loadtxt(file,delimiter = " ");    
    
    lamb       = dataMean[:,0];
    recomPoint = dataMean[:,1];   
    
    plt.plot(lamb,recomPoint/L,'-o',c = c[count],linewidth = 3.);
    count = count+1;

cb = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap))
cb.ax.tick_params(labelsize=24) 
cb.ax.set_title(r'$Re$',fontsize=30)
plt.xlabel(r"$\lambda$",fontsize = 30.);
plt.ylabel(r"$\frac{z_{comb}}{L}$",fontsize = 30.,rotation=0.,labelpad = 15.);
plt.xticks(fontsize = 24.);
plt.yticks(fontsize = 24.);


marker = [".","s","v","^","<",">","o","P","*","+","h"];
count=0;
plt.figure(figsize = (12,8) );
for i in listReFolder[0:-1]:
    with open("file/meanFieldBottomWall_"+str(i)+".dat","rb") as file:
        dataMean = np.loadtxt(file,delimiter = " ",skiprows = 1);
    
    t = dataMean[:,1];
    s = dataMean[:,2];
    lamb = dataMean[:,0]
    plot = plt.scatter(t[t<2],s[t<2],marker = marker[count],c=lamb[t<2],s=50,cmap ='jet',linewidth = 3.);
    plt.scatter( [],[],marker = marker[count],label = r"$Re = $" + str(i[3:len(i)+1]),c="darkred");
    count = count + 1;
    
plt.legend(bbox_to_anchor=(0.24,0.35),fontsize = 18.,handletextpad=0.005,frameon=False);
cb = plt.colorbar(plot);
cb.ax.tick_params(labelsize=24) 
cb.ax.set_title(r'$\lambda$',fontsize=30)
plt.xlabel(r"$\frac{t}{\tau}$",fontsize = 30.);
plt.ylabel(r"$\frac{\bar{s}}{s_0}$",fontsize = 30.,rotation=0.,labelpad = 15.);
plt.xticks(fontsize = 24.);
plt.yticks(fontsize = 24.);


count=0;
plt.figure(figsize = (12,8));
for i in listReFolder[0:-1]:
    print(i)
    with open("file/mergingPoint_"+str(i)+".dat","rb") as file:
        dataMP = np.loadtxt(file,delimiter =" ",skiprows=1);  # merging point

    plt.plot(dataMP[:,0],dataMP[:,1]/L,'o-',c=c[count]);
    count = count+1

cb = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap))
cb.ax.tick_params(labelsize=24) 
cb.ax.set_title(r'$Re$',fontsize=30)
plt.xlabel(r"$\lambda$",fontsize = 30.);
plt.ylabel(r"$\frac{z_{m}}{L}$",fontsize = 30.,rotation=0.,labelpad = 15.);
plt.xticks(fontsize = 24.);
plt.yticks(fontsize = 24.);


count=0;
plt.figure(figsize = (12,8));
for i in listReFolder[0:-1]:
    print(i)
    with open("file/fountainFlow_"+str(i)+".dat","rb") as file:
        dataFF = np.loadtxt(file,delimiter =" ",skiprows=1);  # merging point

    plt.plot(dataFF[:,0],dataFF[:,1]/L,'o-',c=c[count]);
    count = count+1

cb = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap))
cb.ax.tick_params(labelsize=24) 
cb.ax.set_title(r'$Re$',fontsize=30)
plt.xlabel(r"$\lambda$",fontsize = 30.);
plt.ylabel(r"$\frac{z_{f}}{L}$",fontsize = 30.,rotation=0.,labelpad = 15.);
plt.xticks(fontsize = 24.);
plt.yticks(fontsize = 24.);



###############################################################################  


