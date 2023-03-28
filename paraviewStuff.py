# trace generated using paraview version 5.10.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

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

path        =  "../dataPostProcessingFoam/"
listReFolder = os.listdir(path);
listSFolder  = os.listdir( path + listReFolder[0]);

listSFolder  = sortedFolder(path + listReFolder[0] , listSFolder , "S" , "S" ,1 , len(listSFolder) + 1 );
listReFolder = sortedFolder(path + listReFolder[0] , listReFolder , "Re_" , "Re_" ,3, len(listSFolder) + 1 );

lenListRe = len(listReFolder); 
lenListS  = len(listSFolder);
count=0;
countRe=0;

for i in listReFolder:
    for j in listSFolder
        # find source
        internalvtu = FindSource( '/Users/thibaut/Desktop/dataPostProcessingFOAM/' + str(i) + '/' + str(j) + 'internal.vtu' )

        # create a new 'Slice'
        slice1 = Slice(registrationName='Slice1', Input=internalvtu)
        slice1.SliceType = 'Plane'
        slice1.HyperTreeGridSlicer = 'Plane'
        slice1.SliceOffsetValues = [0.0]

        # init the 'Plane' selected for 'SliceType'
        slice1.SliceType.Origin = [0.0, 0.0, 0.001500000013038516]

        # init the 'Plane' selected for 'HyperTreeGridSlicer'
        slice1.HyperTreeGridSlicer.Origin = [0.0, 0.0, 0.001500000013038516]

        # Properties modified on slice1.SliceType
        slice1.SliceType.Normal = [0.0, 1.0, 0.0]

        # show data in view
        slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

        # get color transfer function/color map for 'U'
        uLUT = GetColorTransferFunction('U')

        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        slice1Display.ScaleTransferFunction.Points = [-0.013020760379731655, 0.0, 0.5, 0.0, 0.013020760379731655, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        slice1Display.OpacityTransferFunction.Points = [-0.013020760379731655, 0.0, 0.5, 0.0, 0.013020760379731655, 1.0, 0.5, 0.0]

        # hide data in view
        Hide(internalvtu, renderView1)

        # show color bar/color legend
        slice1Display.SetScalarBarVisibility(renderView1, True)

        # update the view to ensure updated data information
        renderView1.Update()

        # get opacity transfer function/opacity map for 'U'
        uPWF = GetOpacityTransferFunction('U')

        # save data
        SaveData( '/Users/thibaut/Desktop/dataPostProcessingFOAM/' + str(i) + "/" + str(j) + 'concentrationPlate/'  + 'concentrationPlate.csv', proxy=slice1, ChooseArraysToWrite=1,
        CellDataArrays=['s'] )

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================
# get layout
layout1 = GetLayout()
#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(2294, 774)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [6.5576844869551e-05, 0.0659658184273601, 0.0014081924011173138]
renderView1.CameraFocalPoint = [6.5576844869551e-05, 0.0, 0.0014081924011173138]
renderView1.CameraViewUp = [1.1102230246251565e-15, 0.0, -1.0]
renderView1.CameraParallelScale = 0.001733367868623468
