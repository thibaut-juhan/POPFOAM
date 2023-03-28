#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 14:31:56 2023

@author: thibaut
"""

from jetDynamics import *
from numericalStuff     import *
import matplotlib.patches as mpatches
import matplotlib.colors  as cl

path        =  "../dataPostProcessingFoam/"
fileVel     =  "line_U:Transformed.xy"
fileVelVort =  "line_U:Transformed_vorticity.xy"
fileQ       =  "Q.xy"

R = 5e-4;
L = 3e-3;

Vcell = np.pi*L*(1.e-2)**2;


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
for i in listReFolder:    
    #### TEST ON OPEN FILE : If it already exists => delete it
    Reading.openingFileTest("file/meanFieldBottomWall_" + str(i) + ".dat" );
    Reading.openingFileTest("file/combinaisonPoint_"    + str(i) + ".dat" );
    Reading.openingFileTest("file/filmThickness_"       + str(i) + ".dat" );
    Reading.openingFileTest("file/mergingPoint_"        + str(i)  + ".dat");
    Reading.openingFileTest("file/fountainFlow_"        + str(i)  + ".dat");
    Reading.openingFileTest("file/lengthCore_"          + str(i) + ".dat" );
    Reading.openingFileTest("file/expansionAngle_"          + str(i) + ".dat" );
    
    
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
        time       = timeFolder[-1];
        
        data       =  Reading.readingData(pathDataSliceXZAxial + '/' + str(listFileSliceAxial[0]) + "/" +str(time) + "/" , fileVel, ';')
        saveX  = np.zeros( (len(data[:,0]) , len(listFileSliceAxial) ) );
        saveZ  = saveX.copy();
        saveUx = saveX.copy();
        saveUz = saveX.copy();
        
        
        for k in listFileSliceAxial: # loop over folder in Re*/S*/data/*
            # loop over runTime data => 
            timeFolder = os.listdir(pathDataSliceXZAxial + '/' + str(k) + "/" );
            time       = timeFolder[-1];
            print(timeFolder);
            
            data       =  Reading.readingData(pathDataSliceXZAxial + '/' + str(k) + "/" +str(time) + "/" , fileVel, ';')
            
            x,y,z      =  jetDynamics.xyzCoordinate (data,0, 1, 2);
            ux,uy,uz   =  jetDynamics.velocityProfil(data,3 ,4 ,5);
            
            saveX [:,count]  = x;
            saveZ [:,count]  = z;
            saveUx[:,count]  = ux;
            saveUz[:,count]  = uz;
            count = count+1;
        
        angleExpansion = jetDynamics.angleExpansion(saveUz,saveZ,saveX,
                                                           np.float_(j[1:len(j)+1]),
                                                           R);
        
        z095 = jetDynamics.jetCoreLength(saveUz, saveZ, saveX,
                                                        np.float_(j[1:len(j)+1]), 
                                                        R,len(listFileSliceAxial),
                                                        inletVelocity[int(countRe)]
                                                        );
        
        print(angleExpansion);
       
        dataLengthCore      = np.zeros((1,2));
        dataLengthCore[0,0] = np.float_(j[1:len(j)+1]);
        dataLengthCore[0,1] = z095;
        
        
        dataExpAngle      = np.zeros((1,2));
        dataExpAngle[0,0] = np.float_(j[1:len(j)+1]);
        dataExpAngle[0,1] = angleExpansion;
        
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
            time       = timeFolder[-1];
            
            data                 =  Reading.readingData(pathDataSliceXZRadial + '/' + str(k) + "/" + str(time) + "/" , fileVelVort, ';');
            x,y,z                =  jetDynamics.xyzCoordinate  (data,0, 1, 2);
            ux,uy,uz             =  jetDynamics.velocityProfil (data,3 ,4 ,5);
            filmThickness[count] =  jetDynamics.boundaryLayerLength(z, L, ux);
            r[count]=x[0];
            count = count+1;
            
        count      = 0; #reset
        r          = jetDynamics.changingOrigin(lamb,r,R/2.);
        rChanging  = r[r>0]; # changing origin
        filmThickness = filmThickness[r>0];
        
        dataBL = np.zeros( (len(filmThickness),3) );
        dataBL[:,0] = lamb*np.ones(len(filmThickness));
        dataBL[:,1] = rChanging;
        dataBL[:,2] = filmThickness;
       
###############################################################################   


########################## POST-PROCESS CENTER JETDATA#########################      
 
        timeFolder = os.listdir(pathDataCenterJet + '/' + str( listFileCenterJet[0]) + "/" );
        time       = timeFolder[-1];
        data       =  Reading.readingData(pathDataCenterJet + '/' + str(listFileCenterJet[0]) + "/" + str(time) + "/", fileVel, ';')
    
        x,y,z      =  jetDynamics.xyzCoordinate  (data,0, 1, 2);
        ux,uy,uz   =  jetDynamics.velocityProfil (data,3 ,4 ,5);
    
        maxUz              = jetDynamics.maxVelocity(uz,);
        recombPoint        = z[uz == maxUz ];
            
        dataComb        = np.zeros ( (1,2) );
        
        dataComb[0,0]   = lamb;
        dataComb[0,1]   = recombPoint;
        
        
        zerosUz = NumericalStuff.dichotomieIntervall(z,uz,1.e-9,boundary=True);
        
        dataMP = np.zeros((1,2));
        dataFF = np.zeros((1,2));
        
        if len(zerosUz) != 0: # fountain flow and merging point exists ! 
            
            mergingPoint = min(zerosUz);
            fountainFlow = max(zerosUz);
            
            
            dataMP[0,0] = lamb;
            dataMP[0,1] = mergingPoint;
            
            if mergingPoint == fountainFlow: 
            # mergingPoint == fountain flow => impact on the top wall / no merging      
                
                dataFF[0,0] = lamb;
                dataFF[0,1] = 0;                # pas d'écoulement fontaine 
           
            else:
                
                dataFF[0,0] = lamb;
                dataFF[0,1] = L-fountainFlow;            
                
        else:  # pas de zone fountaine / zone de recombinaison 
            
            dataMP[0,0] = lamb;
            dataFF[0,0] = lamb;        
            
            if( np.float_(j[1:len(j)+1]) < 3):
                
                dataMP[0,1] = 0; # No merging point
                dataFF[0,1] = 0; # No fountainFlow  
                
            else: # Possibly no interaction => if there is no merging/fountain flow => merging and fountain at L
                dataMP[0,1] = L; 
                dataFF[0,1] = L; 
            
        
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
                          ,"lamdda,boundary layer thickness , radial position");
        hCore     = Reading.addHeader("file/lengthCore_" + str(i) + ".dat"
                          ,"lambda,length core");
        
        hExp     = Reading.addHeader("file/expansionAngle_" + str(i) + ".dat"
                          ,"lambda,expansion angle");
            
            
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
        
        with open ("file/lengthCore_" + str(i) + ".dat","ab" ) as file:
            np.savetxt(file,dataLengthCore,delimiter = " ",header = hCore);
            
        with open ("file/expansionAngle_" + str(i) + ".dat","ab" ) as file:
            np.savetxt(file,dataExpAngle,delimiter = " ",header = hExp);
      
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
for i in listReFolder:
    with open("file/combinaisonPoint_" + str(i) + ".dat","rb") as file:
        dataCB = np.loadtxt(file,delimiter = " ");    
    
    lamb       = dataCB[:,0];
    recomPoint = dataCB[:,1];   
    
    Re=int(i[3:len(i)+1])*np.ones(len(lamb));
    scat=plt.scatter(Re,recomPoint/L,s=100,c=lamb,cmap='jet');

    count = count+1;

cb = plt.colorbar(scat)
cb.ax.tick_params(labelsize=24) 
cb.ax.set_title(r'$\lambda$',fontsize=30)
plt.xlabel(r"$Re$",fontsize = 30.);
plt.ylabel(r"$\frac{z_{cb}}{L}$",fontsize = 30.,rotation=0.,labelpad = 15.);
plt.xticks(fontsize = 24.);
plt.yticks(fontsize = 24.);


marker = [".","s","v","^","<",">","o","P","*","+",'x'];
count=0;
plt.figure(figsize = (12,8) );
for i in listReFolder:
    with open("file/meanFieldBottomWall_"+str(i)+".dat","rb") as file:
        dataMean = np.loadtxt(file,delimiter = " ");
    
    massFlux = 2.*inletVelocity[int(count)]*np.pi*(R**2);
    tau      = Vcell/massFlux;
    t = dataMean[:,1];
    s = dataMean[:,2];
    lamb = dataMean[:,0]
    plot = plt.scatter(t[t<2],s[t<2],marker = marker[count],c=lamb[t<2],s=50,cmap ='jet');
    plt.scatter( [],[],marker = marker[count],label = r"$Re = $" + str(i[3:len(i)+1]),c="darkred");
    count = count + 1;
    
plt.legend(bbox_to_anchor=(0.25,0.35),fontsize = 18.,handletextpad=0.005,frameon=False);
cb = plt.colorbar(plot);
cb.ax.tick_params(labelsize=24) 
cb.ax.set_title(r'$\lambda$',fontsize=30)
plt.xlabel(r"t",fontsize = 30.);
plt.ylabel(r"$\frac{\bar{s}_{Wall} }{s_0}$",fontsize = 30.,rotation=0.,labelpad = 15.);
plt.xticks(fontsize = 24.);
plt.yticks(fontsize = 24.);


count=0;
plt.figure(figsize = (12,8));
for i in listReFolder:
    with open("file/mergingPoint_"+str(i)+".dat","rb") as file:
        dataMP = np.loadtxt(file,delimiter =" ",skiprows=1);  # merging point

    Re=int(i[3:len(i)+1])*np.ones(len(dataMP[:,0]));
    scat=plt.scatter(Re,dataMP[:,1]/L,s=100,c=dataMP[:,0],cmap='jet');
    count = count+1

cb = plt.colorbar(scat)
cb.ax.tick_params(labelsize=24) 
cb.ax.set_title(r'$\lambda$',fontsize=30)
plt.xlabel(r"$Re$",fontsize = 30.);
plt.ylabel(r"$\frac{z_{mp}}{L}$",fontsize = 30.,rotation=0.,labelpad = 15.);
plt.xticks(fontsize = 24.);
plt.yticks(fontsize = 24.);

###################### RECAP FROM THESE OBSERVATIONS ##########################
count=0;
plt.figure(figsize = (12,8));

normalize = cl.Normalize(vmin=0, vmax=1)

for i in listReFolder:
    with open("file/combinaisonPoint_"+str(i)+".dat","rb") as file:
        dataCP = np.loadtxt(file,delimiter =" ");  # merging point
    Re=int(i[3:len(i)+1])*np.ones(len(dataCP[:,0])); 
    scat=plt.scatter(Re,dataCP[:,0],s=100,c=dataCP[:,1]/L,cmap="jet",norm=normalize);
    count = count+1


## add rectangle for each zone :
rect1=mpatches.Rectangle((5,2.1),50,3, 
                        fill = False,
                        color = "purple",
                        linewidth = 3)

rect2=mpatches.Rectangle((55,2.1),50,0.9, 
                        fill = False,
                        color = "purple",
                        linewidth = 3)

rect3=mpatches.Rectangle((55,3.),50,2.1, 
                        fill = False,
                        color = "orange",
                        linewidth = 3)

rect4=mpatches.Rectangle((105,3.),100,2.1, 
                        fill = False,
                        color = "red",
                        linewidth = 3)

rect5=mpatches.Rectangle((105,2.1),100,0.9, 
                        fill = False,
                        color = "orange",
                        linewidth = 3)

plt.gca().add_patch(rect1)
plt.gca().add_patch(rect2)
plt.gca().add_patch(rect3)
plt.gca().add_patch(rect4)
plt.gca().add_patch(rect5)


cb = plt.colorbar(scat)
cb.ax.tick_params(labelsize=24) 
cb.ax.set_title(r'$\frac{z_{CB}}{L}$',fontsize=30)
plt.xlabel(r"$Re$",fontsize = 30.);
plt.ylabel(r"$\lambda$",fontsize = 30.,rotation=0.,labelpad = 15.);
plt.xticks(fontsize = 24.);
plt.yticks(fontsize = 24.);





count=0;
plt.figure(figsize =(12,8) );
for i in listReFolder:
    with open("file/lengthCore_" +str(i)+".dat","rb") as file:
        dataLC = np.loadtxt(file,delimiter =" ",skiprows=1);  # merging point

  
    Re=int(i[3:len(i)+1])*np.ones(len(dataLC[:,0]));
    scat=plt.scatter(Re,dataLC[:,1]/L,s=100,c=dataLC[:,0],cmap='jet');
    count = count+1
cb = plt.colorbar(scat)
cb.ax.tick_params(labelsize=24) 
cb.ax.set_title(r'$\lambda$',fontsize=30)
plt.xlabel(r"$Re$",fontsize = 30.);
plt.ylabel(r"$\frac{l_{c}}{L}$",fontsize = 30.,rotation=0.,labelpad = 15.);
plt.xticks(fontsize = 24.);
plt.yticks(fontsize = 24.);



count=0;
plt.figure(figsize =(12,8) );
for i in listReFolder:
    with open("file/expansionAngle_" +str(i)+".dat","rb") as file:
        dataEA = np.loadtxt(file,delimiter =" ",skiprows=1);  # merging point
  
    Re=int(i[3:len(i)+1])*np.ones(len(dataEA[:,0]));
    scat=plt.scatter(Re,dataEA[:,1]*180/np.pi,s=100,c=dataEA[:,0],cmap='jet');
    count = count+1

cb = plt.colorbar(scat)
cb.ax.tick_params(labelsize=24) 
cb.ax.set_title(r'$\lambda$',fontsize=30)
plt.xlabel(r"$Re$",fontsize = 30.);
plt.ylabel(r"$\Theta$",fontsize = 30.,rotation=0.,labelpad = 15.);
plt.xticks(fontsize = 24.);
plt.yticks(fontsize = 24.);

count=0;
plt.figure(figsize = (12,8));
for i in listReFolder:
    with open("file/fountainFlow_"+str(i)+".dat","rb") as file:
        dataFF = np.loadtxt(file,delimiter =" ",skiprows=1);  # merging point

    Re=int(i[3:len(i)+1])*np.ones(len(dataFF[:,0]));
    scat = plt.scatter(Re,dataFF[:,1]/L,s=100,cmap="jet",c=dataFF[:,0]);
    count = count+1

cb = plt.colorbar(scat)
cb.ax.tick_params(labelsize=24) 
cb.ax.set_title(r'$\lambda$',fontsize=30)
plt.xlabel(r"$Re$",fontsize = 30.);
plt.ylabel(r"$\frac{z_{f}}{L}$",fontsize = 30.,rotation=0.,labelpad = 15.);
plt.xticks(fontsize = 24.);
plt.yticks(fontsize = 24.);

# plot on the same plot : zcp + zmp:
dataPlotMP = np.zeros( (len(dataMP),2,len(listReFolder) -1 ));
dataPlotCB = np.zeros( (len(dataCB),2,len(listReFolder) -1 ));

count=0;
plt.figure(figsize = (12,8)) 
for i in listReFolder:
   with open("file/mergingPoint_"+str(i)+".dat","rb") as file:
        dataMP = np.loadtxt(file,delimiter =" ");  # merging point
   with open("file/combinaisonPoint_" + str(i) + ".dat","rb") as file:
       dataCB = np.loadtxt(file,delimiter = " ");    
       
       
   Re_MP=int(i[3:len(i)+1])*np.ones(len(dataMP[:,0]));
   Re_CB=int(i[3:len(i)+1])*np.ones(len(dataCB[:,0]));
   
   #dataPlotMP[:,:,count] = dataMP[:,:];
   #dataPlotCB[:,:,count] = dataCB[:,:];
    
   scat1=plt.scatter(Re_MP,dataMP[:,1]/L,s=100,marker ='o',c=dataMP[:,0],cmap='jet');   
   scat2=plt.scatter(Re_CB,dataCB[:,1]/L,s=100,marker = 's',c=dataCB[:,0],cmap='jet');
     
   #scat3=plt.scatter(Re_CB,abs(dataCB[:,1]-dataMP[:,1])/L,s=100,marker = 'o',c=dataCB[:,0],cmap='jet');
    
   
   count=count+1;
   
   
plt.scatter([],[],marker = 'o',c = "darkred",label = r"$z_{mp}$"); 
plt.scatter([],[],marker = 's',c = "darkred",label = r"$z_{cb}$"); 
cb = plt.colorbar(scat)
cb.ax.tick_params(labelsize=24) 
cb.ax.set_title(r'$\lambda$',fontsize=30)
plt.xlabel(r"$Re$",fontsize = 30.);
plt.ylabel(r"$\frac{z_{cb}}{L},\frac{z_{mb}}{L}$",fontsize = 30.,rotation=90.,labelpad = 15.);
plt.xticks(fontsize = 24.);
plt.yticks(fontsize = 24.);
plt.legend();









###############################################################################  


