#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 10:48:31 2023

@author: thibaut
"""

######################### CLASS : READING STUFF ###############################


import numpy as np;
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
from matplotlib.colors import Normalize

class Reading:
    
    def _init_(self):  # constructeur 
        pass;
    
    def sortedFolder (path,listFolder,characterDel,characterAdd,startSorted,endSorted): # sortedFolder
        #'                   list of args:                       '
        #'   - path          : path directory of the folder list'        => str
        #'   - listFolder    : list of folder which will be sorted (S*)' => list
        #'   - character     : character added ( if not => ' ') '        => str
        #'   - startdeSorted : started character that we want to keep  ' => str
        #'   - endSorted     : Ended character that we want to keep  '   => str     
        
        chara=[]
        # recupere les 2 derniers caractéres de la liste       
        for i in listFolder:
                    chara.append( i[startSorted:endSorted] );    
                    # get a specific number of character    
        
        for i in range ( len(chara) ):
            if chara[i].isdigit() == False:  
                # si la chaine de caractérer contient une lettre => retourne False
                chara[i] = chara[i].replace(characterDel,"" );
                
        chara = sorted(chara,key=float);
        s     = characterAdd;
        
        list_folder= [ s + x for x in chara];
        
        return list_folder;

        
    def readingData(path,fileName,delimiter): #reading Data
        #'                    list of args:                  '
        #'   - path          : path directory of the file    '        => str
        #'   - fileName      : name of file                  '        => str
        #'   - delimiter     : delimiter between data        '        => str   
        
        with open (path + fileName,"rb") as file:
            data = np.loadtxt(file);
        return data;
    
    def readingDataConcenreation(path,fileName):
        with open (path + fileName,"rb") as file:
            data = np.loadtxt(file,skiprows = 5);
            
        return data;
    
    
    def openingFileTest(nameFile):
        if not os.path.isfile(nameFile):
            print("file doesnt exist")
        else:
            print("file already exists => delete this file ");
            f = open(nameFile,'r+');
            f.truncate(0);
            f.close();
         
    def addHeader(nameFile,header):
        if not os.path.isfile(nameFile):
            header  = header;
        else:
            header = '';
        
        return header;
        
    def computeLambda (D,S):
        #'                    list of args:                  '
        #'   - D             : length                        '        => float/int
        #    - S             :  S list            '        => list 
        l = (S*1e-3/D);
        return l;
       
class Plot:
    def __init__():
        pass;
    def quiverJetsProfile(
            x,z,ux,uz,U0,D,
            xmin,xmax,ymin,ymax,
            xLabel,yLabel,cbarLabel,
            intervall,legend,savePath
            ): ### QUIVER JETS PROFILE : 
        
        #'                   list of args:                       '
        #'   - x          :  x coordinate'        => np.array    '
        #'   - y          :  y coordinate'        => np.array    '
        #'   - z          :  z coordinate'        => np.array    '
        
        
        #'   - ux          :  u in x direction'        => np.array   '
        #'   - uz          :  u in z drection'        => np.array    '
        
        
        #'   - U0          :  mag of initial velocity '        => float   '
        #'   - D           :  length '                         => float     '
        
        #'   - xmin          :  min in x direction  '        => float   '
        #'   - xmax          :  max in x direction  '        => float   '
        #'   - ymin          :  min in y direction  '        => float   '
        #'   - ymax          :  max in y direction '         => float   '
        
        #'   - xLabel         :  x label  '        => str   '
        #'   - yLabel         :  y label  '        => str   '
        #'   - yLimMin        :  min of y limit '        => float   '
        #'   - yLimMax        :  max of y limit'         => float   '
        
        
        
        #'   - legend          :  legend added in plot  '        => str   '
        #'   - savePath        :  savePath  '         => str   '
        
        deltaZ = 0;
        plt.figure(figsize=(12,8));
        for i in range (0,len(x[0]),intervall):
            
            if(i==0):
                deltaZ     = 2*np.abs(z[1,1]-z[2,2]);        
                norm = Normalize();
                norm.autoscale(uz[:,i]/U0)
                colormap = cm.jet
                sm = cm.ScalarMappable(cmap=colormap, norm=norm)
                sm.set_array([])
            
            else:
                deltaZ     = (i+2)*np.abs(z[1,1]-z[2,2]);        
            
            plt.plot(x[:,i]/D,uz[:,i]/U0 + deltaZ/D,c="red",linewidth = 3.);
            plt.quiver(x[:,i]/D,z[:,i]/D,ux[:,i],uz[:,i]/U0,width = 0.005,scale_units = 'xy',scale = 1,color = colormap(norm(uz[:,i]/U0))) #quiver en 0

        cbar = plt.colorbar(sm);  
        cbar.set_label(r'$\frac{U_z}{U_0}$', rotation=0,fontsize = 30,labelpad=15.);
        cbar.ax.tick_params(labelsize = 20)

        plt.vlines(0,ymin=ymin,ymax=ymax,linewidth=3.,linestyle='--',color = "black");
        plt.hlines(6,xmin=xmin,xmax=xmax,linewidth=3.,linestyle='-',color = "black");
        
        plt.xlabel(xLabel,fontsize = 30.)
        plt.ylabel(yLabel,fontsize = 30.,rotation=0,labelpad = 15);
        
        plt.xticks(fontsize = 24.);
        plt.yticks(fontsize = 24.);
        plt.legend(loc ="upper right");
        plt.ylim(ymin,ymax)
        plt.text(-1.08, 6.25, "paroi solide",fontsize = 24)
        plt.title(r"$\lambda = $" + legend,fontsize = 24.)
        plt.savefig(savePath + legend + ".png",bbox_inches = "tight")
        plt.show();
    
    def plotCenterLine(
            z,uz,U0,D
            ,maxZ,maxVel,l,
            saveName,xLabel,yLabel
            ):
        plt.figure(figsize = (12,8));
        plt.plot(z/D,uz/U0,linewidth = 3.,label= l);
        plt.scatter(maxZ/D,maxVel/U0,s=100,c="black",label=r"$z_{cmb}$");
        plt.plot(maxZ/D,maxVel/U0,linewidth=3.,linestyle='--',c="black");
        plt.xlabel(xLabel,fontsize = 30);
        plt.ylabel(yLabel,fontsize = 30,rotation=0.,labelpad=15);
        plt.legend(loc = "upper right",fontsize = 20.)
        plt.xticks(fontsize = 24.);
        plt.yticks(fontsize = 24.);
        plt.axhline(y=0,linestyle = '--',linewidth=3.,c="black");
        plt.axvline(x=0,linestyle = '--',linewidth=3.,c="black");
        plt.savefig(saveName + ".png",bbox_inches = "tight");
        
    def plotFilmThickness(h,r,l,L,
                          loc,savename,name,
                          Re,epsilon,
                          xLabel,yLabel):         
        plt.figure(figsize = (12,8));
        for i in range (len(r[0])):
            indexPos = np.argwhere(r[:,i] >= 0);
            r_pos    = r[indexPos[:,0],i];
            h_pos    = h[indexPos[:,0],i];
            
            indexHcirc = np.argwhere(abs(h_pos[0:-1]-h_pos[1:len(h_pos)])<epsilon);
            h_del = np.delete(h_pos,indexHcirc);
            r_del = np.delete(r_pos,indexHcirc);
            
            plt.scatter(r_del/1.e-2,h_del/L,s=100);
            plt.plot(r_del/1.e-2,h_del/L,linewidth = 3.,label = r"$\lambda = $" + str(l[i]),linestyle='--');
            plt.legend(loc = loc,fontsize = 20.);
            
        plt.xlabel(xLabel,fontsize = 30);
        plt.ylabel(yLabel,fontsize = 30,rotation=0.,labelpad=15);
        plt.xticks(fontsize = 24.);
        plt.yticks(fontsize = 24.);    
        plt.savefig(savename + name + Re,bbox_inches = "tight");
        
        
    def plotRecirculationHeight(h,r,l,
                                L,loc,savename,
                                name,Re,epsilon,
                                xLabel,yLabel):
        
        plt.figure(figsize = (12,8));
        for i in range (len(r[0])):
            indexPos = np.argwhere(r[:,i] >= 0);
            r_pos    = r[indexPos,i];
            h_pos    = h[indexPos,i];
            
            indexHcirc = np.argwhere(abs(h_pos[0:-1]-h_pos[1:len(h_pos)])<epsilon);
            
            h_del = np.delete(h_pos,indexHcirc);
            r_del = np.delete(r_pos,indexHcirc);
        
            plt.scatter(r_del[0:-1]/1.e-2,h_del[0:-1]/L,s=100);
            plt.plot   (r_del[0:-1]/1.e-2,h_del[0:-1]/L,linewidth = 3.,label = r"$\lambda = $" + str(l[i]),linestyle='--');
            plt.legend(loc = loc,fontsize = 20.);
            
        plt.xlabel(xLabel,fontsize = 30);
        plt.ylabel(yLabel,fontsize = 30,rotation=0.,labelpad=15);
        plt.xticks(fontsize = 24.);
        plt.yticks(fontsize = 24.);    
        plt.savefig(savename + name + Re,bbox_inches = "tight");


        


    
    
    
    
    
    
    
    
    
  
        
        
        
