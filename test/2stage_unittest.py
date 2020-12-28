import os
from osgeo import gdal

from ARIAtools.unwrapStitching import product_stitch_2stage, product_stitch_overlap, product_stitch_morph
from ARIAtools.shapefile_util import open_shapefile, save_shapefile
from ARIAtools.extractProduct import merged_productbbox
import numpy as np
import pdb







def mkdirs(dirlist):
    '''
        Checks if directory exists, if not creates them
    '''

    if not isinstance(dirlist, list):
        dirlist = [dirlist]
    for dir in dirlist:
        if not os.path.exists(dir):
            os.makedirs(dir)

def unitTest0(files,workdir='unittest0',croptounion='True',outputFormat='Envi',verbose=True):
    
    # selecting the files which will be used in the test
    files = files[0:2]
    
    workdir = os.path.abspath(workdir)
    curdir = os.path.abspath(os.path.curdir)
    
    # stage expected files
    sim_conn_files,sim_unw_files,bounds,prods_TOTbbox,prod_bbox_files = stageExpectedFiles(files,workdir,croptounion)
    
    # overwrite unw-files and connected component file with specific simualtion for unit test 1.
    for file_counter in range(len(sim_unw_files)):
        # extracting phase data (1 band file is the target, so if you are using ENVI make sure to take the right band!)
        
        
        #create simulations
        # doing the phase
        data = gdal.Open(sim_unw_files[file_counter],gdal.GA_ReadOnly)
        data_band = data.GetRasterBand(1)
        phase = data_band.ReadAsArray()
        data = None
        
        # get 0. phase indices
        phase_0s=np.nonzero(phase==0.)
        
        width=int(phase.shape[1]/2)
        length=int(phase.shape[0]/2)
        
        # doing the connected component
        data = gdal.Open(sim_conn_files[file_counter],gdal.GA_ReadOnly)
        data_band = data.GetRasterBand(1)
        connComp = data_band.ReadAsArray()
        data = None
        #get 0/1 conncomp indices
        connComp_0s=np.nonzero(connComp==0)
        connComp_nodatas=np.nonzero(connComp==-1)
        
        # manipulating the connected component
        phase = phase*0.0+np.pi
        phase[0:length,:] = phase[0:length,:]+(2)*np.pi
        phase[length-1:,:] = phase[length-1:,:]+(-4)*np.pi
        phase[phase_0s]=0
        
        connComp = connComp*0 + 1
        connComp[length-1:,:]  = connComp[length-1:,:] +1
        connComp[connComp_0s]=0
        connComp[connComp_nodatas]=-1
        connComp = connComp.astype('int')
        
        
        # writing put new phase
        ds=gdal.Open(sim_unw_files[file_counter],gdal.GA_Update)
        ds.GetRasterBand(1).WriteArray(phase)
        del ds
        cmd = "gdal_translate -of png -scale -ot Byte -q " + sim_unw_files[file_counter] + " " + sim_unw_files[file_counter] + ".png"
        os.system(cmd)
        
        
        # writing out new conncomp
        ds=gdal.Open(sim_conn_files[file_counter],gdal.GA_Update)
        ds.GetRasterBand(1).WriteArray(connComp)
        del ds
        cmd = "gdal_translate -of png -scale -ot Byte -q " + sim_conn_files[file_counter] + " " + sim_conn_files[file_counter] + ".png"
        os.system(cmd)
    
    #run test
    os.chdir(workdir)
    product_stitch_overlap(sim_unw_files,sim_conn_files,prod_bbox_files,bounds,prods_TOTbbox, outFileUnw=os.path.join(workdir,'unwrap'),outFileConnComp= os.path.join(workdir,'conncomp'), mask=mask,outputFormat = outputFormat,verbose=verbose)
    #product_stitch_2stage(sim_unw_files,sim_conn_files,prod_bbox_files,bounds,prods_TOTbbox,outFileUnw=os.path.join(workdir,'unwrap'),outFileConnComp=os.path.join(workdir,'conncomp'), mask=mask, outputFormat=outputFormat,verbose=verbose)
    os.chdir(curdir)

#print('testproduct_stitch_2stage(' + str(sim_unw_files) + ',' + str(sim_conn_files) + ',' + str(prod_bbox_files) + ',' + str(bounds) + ',' + str(prods_TOTbbox) + ', mask=' + str(mask) + ', outputFormat=' + str(outputFormat) + ',verbose=' + str(verbose) +')')


def unitTest1(files,workdir='unittest1',croptounion='True',outputFormat='Envi',verbose=True):
    
    # selecting the files which will be used in the test
    files = files[0:2]

    workdir = os.path.abspath(workdir)
    curdir = os.path.abspath(os.path.curdir)
    
    # stage expected files
    sim_conn_files,sim_unw_files,bounds,prods_TOTbbox,prod_bbox_files = stageExpectedFiles(files,workdir,croptounion)
    
    # overwrite unw-files and connected component file with specific simualtion for unit test 1.
    for file_counter in range(len(sim_unw_files)):
        # extracting phase data (1 band file is the target, so if you are using ENVI make sure to take the right band!)
        
        
        #create simulations
        # doing the phase
        data = gdal.Open(sim_unw_files[file_counter],gdal.GA_ReadOnly)
        data_band = data.GetRasterBand(1)
        phase = data_band.ReadAsArray()
        data = None
        
        # get 0. phase indices
        phase_0s=np.nonzero(phase==0.)
        
        width=int(phase.shape[1]/2)
        length=int(phase.shape[0]/2)
        
        # doing the connected component
        data = gdal.Open(sim_conn_files[file_counter],gdal.GA_ReadOnly)
        data_band = data.GetRasterBand(1)
        connComp = data_band.ReadAsArray()
        data = None
        #get 0/1 conncomp indices
        connComp_0s=np.nonzero(connComp==0)
        connComp_nodatas=np.nonzero(connComp==-1)
        
        # manipulating the connected component
        range_offsets = [0.1 ,-0.1, 0 ,0.2]        # radians
        
        
        
        phase = phase*0.0+range_offsets[file_counter]
        phase[0:length,:] = 2*np.pi +range_offsets[file_counter]
        phase[length-1:,:] = -4*np.pi +range_offsets[file_counter]
        phase[np.int(length/3):np.int(length/3)+600,np.int(width/3):np.int(width/3)+600] = 6*np.pi +range_offsets[file_counter]
        phase[np.int(2*length/3):length,np.int(2*width/3):width] = -12*np.pi+range_offsets[file_counter]
        phase[np.int(length*5/3):2*np.int(length),np.int(width*5/3):2*np.int(width) ] = -2*np.pi+range_offsets[file_counter]

        phase_noise_degrees= 40
        phase = phase + np.random.randn(phase.shape[0],phase.shape[1])*phase_noise_degrees*2*np.pi/360
        phase[phase_0s]=0
        
        connComp = connComp*0 + 1
        connComp[length-1:,:]  = 2
        connComp[np.int(length/3):np.int(length/3)+600,np.int(width/3):np.int(width/3)+600] = 3
        connComp[np.int(2*length/3):length,np.int(2*width/3):width] =  4
        connComp[np.int(length*5/3):2*np.int(length),np.int(width*5/3):2*np.int(width) ] = 5

        connComp[connComp_0s]=0
        connComp[connComp_nodatas]=-1
        connComp = connComp.astype('int')
        
        
        # writing put new phase
        ds=gdal.Open(sim_unw_files[file_counter],gdal.GA_Update)
        ds.GetRasterBand(1).WriteArray(phase)
        del ds
        cmd = "gdal_translate -of png -scale -ot Byte -q " + sim_unw_files[file_counter] + " " + sim_unw_files[file_counter] + ".png"
        os.system(cmd)
        
        
        # writing out new conncomp
        ds=gdal.Open(sim_conn_files[file_counter],gdal.GA_Update)
        ds.GetRasterBand(1).WriteArray(connComp)
        del ds
        cmd = "gdal_translate -of png -scale -ot Byte -q " + sim_conn_files[file_counter] + " " + sim_conn_files[file_counter] + ".png"
        os.system(cmd)
    
    #run test
    os.chdir(workdir)
    product_stitch_morph(sim_unw_files,sim_conn_files,prod_bbox_files,bounds,prods_TOTbbox, outFileUnw=os.path.join(workdir,'unwrap'),outFileConnComp= os.path.join(workdir,'conncomp'), mask=mask,outputFormat = outputFormat,verbose=verbose)
    #product_stitch_2stage(sim_unw_files,sim_conn_files,prod_bbox_files,bounds,prods_TOTbbox,outFileUnw=os.path.join(workdir,'unwrap'),outFileConnComp=os.path.join(workdir,'conncomp'), mask=mask, outputFormat=outputFormat,verbose=verbose)
    os.chdir(curdir)

    #print('testproduct_stitch_2stage(' + str(sim_unw_files) + ',' + str(sim_conn_files) + ',' + str(prod_bbox_files) + ',' + str(bounds) + ',' + str(prods_TOTbbox) + ', mask=' + str(mask) + ', outputFormat=' + str(outputFormat) + ',verbose=' + str(verbose) +')')


def unitTest2(files,workdir='unittest2',croptounion='True',outputFormat='Envi',verbose=True):
    
    # selecting the files which will be used in the test
    files = files[0:2]
    
    workdir = os.path.abspath(workdir)
    curdir = os.path.abspath(os.path.curdir)
    
    # stage expected files
    sim_conn_files,sim_unw_files,bounds,prods_TOTbbox,prod_bbox_files = stageExpectedFiles(files,workdir,croptounion)
    
    # overwrite unw-files and connected component file with specific simualtion for unit test 1.
    for file_counter in range(len(sim_unw_files)):
        # extracting phase data (1 band file is the target, so if you are using ENVI make sure to take the right band!)
        
        
        #create simulations
        # doing the phase
        data = gdal.Open(sim_unw_files[file_counter],gdal.GA_ReadOnly)
        data_band = data.GetRasterBand(1)
        phase = data_band.ReadAsArray()
        data = None
        
        # get 0. phase indices
        phase_0s=np.nonzero(phase==0.)
        
        width=int(phase.shape[1]/2)
        length=int(phase.shape[0]/2)
        
        # doing the connected component
        data = gdal.Open(sim_conn_files[file_counter],gdal.GA_ReadOnly)
        data_band = data.GetRasterBand(1)
        connComp = data_band.ReadAsArray()
        data = None
        #get 0/1 conncomp indices
        connComp_0s=np.nonzero(connComp==0)
        connComp_nodatas=np.nonzero(connComp==-1)
        
        # manipulating the connected component
        range_offsets = [np.pi ,-0.1, 0 ,0.2]        # radians
        
        
        phase = phase*0.0+range_offsets[file_counter]
        phase[0:length,:] = 2*np.pi +range_offsets[file_counter]
        phase[length-1:,:] = -4*np.pi +range_offsets[file_counter]
        phase[np.int(length/3):np.int(length/3)+600,np.int(width/3):np.int(width/3)+600] = 6*np.pi +range_offsets[file_counter]
        phase[np.int(2*length/3):length,np.int(2*width/3):width] = -12*np.pi+range_offsets[file_counter]
        phase[np.int(length*5/3):2*np.int(length),np.int(width*5/3):2*np.int(width) ] = -2*np.pi+range_offsets[file_counter]
        
        phase_noise_degrees= 90
        phase = phase + np.random.randn(phase.shape[0],phase.shape[1])*phase_noise_degrees*2*np.pi/360
        phase[phase_0s]=0
        
        connComp = connComp*0 + 1
        connComp[length-1:,:]  = 2
        connComp[np.int(length/3):np.int(length/3)+600,np.int(width/3):np.int(width/3)+600] = 3
        connComp[np.int(2*length/3):length,np.int(2*width/3):width] =  4
        connComp[np.int(length*5/3):2*np.int(length),np.int(width*5/3):2*np.int(width) ] = 5
        
        connComp[connComp_0s]=0
        connComp[connComp_nodatas]=-1
        connComp = connComp.astype('int')
        
        
        # writing put new phase
        ds=gdal.Open(sim_unw_files[file_counter],gdal.GA_Update)
        ds.GetRasterBand(1).WriteArray(phase)
        del ds
        cmd = "gdal_translate -of png -scale -ot Byte -q " + sim_unw_files[file_counter] + " " + sim_unw_files[file_counter] + ".png"
        os.system(cmd)
        
        
        # writing out new conncomp
        ds=gdal.Open(sim_conn_files[file_counter],gdal.GA_Update)
        ds.GetRasterBand(1).WriteArray(connComp)
        del ds
        cmd = "gdal_translate -of png -scale -ot Byte -q " + sim_conn_files[file_counter] + " " + sim_conn_files[file_counter] + ".png"
        os.system(cmd)
    
    #run test
    os.chdir(workdir)
    product_stitch_morph(sim_unw_files,sim_conn_files,prod_bbox_files,bounds,prods_TOTbbox, outFileUnw=os.path.join(workdir,'unwrap'),outFileConnComp= os.path.join(workdir,'conncomp'), mask=mask,outputFormat = outputFormat,verbose=verbose)

    #product_stitch_2stage(sim_unw_files,sim_conn_files,prod_bbox_files,bounds,prods_TOTbbox,outFileUnw=os.path.join(workdir,'unwrap'),outFileConnComp=os.path.join(workdir,'conncomp'), mask=mask, outputFormat=outputFormat,verbose=verbose)
    os.chdir(curdir)

#print('testproduct_stitch_2stage(' + str(sim_unw_files) + ',' + str(sim_conn_files) + ',' + str(prod_bbox_files) + ',' + str(bounds) + ',' + str(prods_TOTbbox) + ', mask=' + str(mask) + ', outputFormat=' + str(outputFormat) + ',verbose=' + str(verbose) +')')



def stageExpectedFiles(files,workdir='.',croptounion='True'):
    '''
        Normally ARIA tools class has been pre-run, this function will generate the necessery intermediate files
    '''
    
    # generate the working directories
    UnwDir=os.path.join(os.path.abspath(workdir),'unwrappedPhase')
    ConnCompDir=os.path.join(os.path.abspath(workdir),'connectedComponents')
    mkdirs([UnwDir,ConnCompDir])
    
    conn_files = []
    prod_bbox_files = []
    unw_files = []
    file_counter = 0
    for file in files:
        # the datasets which needs to be extracted
        conn_file_in = 'NETCDF:"%s":/science/grids/data/connectedComponents'%(file)
        unw_file_in = 'NETCDF:"%s":/science/grids/data/unwrappedPhase'%(file)
        prod_bbox_file_in ='NETCDF:"%s":productBoundingBox'%(file)
        # staging names of output datasets
        unw_file_out = os.path.join(UnwDir,'part_' + str(file_counter))
        conn_file_out = os.path.join(ConnCompDir,'part_' + str(file_counter))

        # staging unwrapped files
        gdal.BuildVRT(unw_file_out +'.vrt', unw_file_in )
        gdal.Warp(unw_file_out, unw_file_out +'.vrt', options=gdal.WarpOptions(format=outputFormat))
        gdal.BuildVRT(unw_file_out +'.vrt', unw_file_out)
        # staging connected component files
        gdal.BuildVRT(conn_file_out +'.vrt', conn_file_in )
        gdal.Warp(conn_file_out, conn_file_out +'.vrt', options=gdal.WarpOptions(format=outputFormat))
        gdal.BuildVRT(conn_file_out +'.vrt', conn_file_out)
        
        # tracking the staged list
        conn_files.append(conn_file_out)
        unw_files.append(unw_file_out)
        prod_bbox_files.append(prod_bbox_file_in)
    
        # increment the file counter
        file_counter += 1
    
    
    # get bounding boxes
    prods_TOTbbox=os.path.join(workdir, 'productBoundingBox','productBoundingBox.shp')
    mkdirs(os.path.dirname(prods_TOTbbox))

    # Initiate intersection file with first product
    # this is for different scenes
    save_shapefile(prods_TOTbbox, open_shapefile(prod_bbox_files[0], 'productBoundingBox', 1), 'GeoJSON')
    for scene in prod_bbox_files[1:]:
        prods_bbox=open_shapefile(scene, 'productBoundingBox', 1)
        total_bbox=open_shapefile(prods_TOTbbox, 0, 0)
        # Generate footprint for the union of all products
        if croptounion:
            prods_bbox=prods_bbox.union(total_bbox)
        # Generate footprint for the common intersection of all products
        else:
            prods_bbox=prods_bbox.intersection(total_bbox)
        # Check if there is any common overlap
        if prods_bbox.bounds==():
            raise Exception('No common overlap, footprint cannot be generated. Last scene checked: %s'%(scene['productBoundingBox'][0]))
        save_shapefile(prods_TOTbbox, prods_bbox, 'GeoJSON')
    
    bounds=open_shapefile(prods_TOTbbox, 0, 0).bounds

    return conn_files,unw_files,bounds,prods_TOTbbox, prod_bbox_files



def downloadProducts(downloadFiles,downloadDir='.'):
    '''
        Download files if not-existent to the downloadDir location
        Returns the files and their local path.
    '''
    
    # Create downloadDir if needed
    downloadDir = os.path.abspath(downloadDir)
    mkdirs(downloadDir)
    localDir = os.path.abspath(os.path.curdir)
    
    # create mapping of local download file name
    localFiles=[]
    for downloadFile in downloadFiles:
        localFile = os.path.join(downloadDir,os.path.basename(downloadFile))
        localFiles.append(localFile)
    
    for downloadFile in downloadFiles:
        os.chdir(downloadDir)
        if not os.path.exists(localFile):
            os.system('wget ' + downloadFile)
    os.chdir(localDir)

    return localFiles








if __name__ == '__main__':


    #set workdirectory
    workdir=os.path.abspath('unittestsStitching')
    
    # setting variables
    croptounion=True
    mask=None
    outputFormat='ISCE'
    verbose=True

    # pre-staged files
    downloadDate = ['20141116_20141023']
    downloadFiles = ['https://grfn.asf.alaska.edu/door/download/S1-GUNW-D-R-079-tops-20141116_20141023-030926-08430N_06390N-PP-5c6d-v2_0_2.nc',
             'https://grfn.asf.alaska.edu/door/download/S1-GUNW-D-R-079-tops-20141116_20141023-030901-09928N_07891N-PP-7d6c-v2_0_2.nc','https://grfn.asf.alaska.edu/door/download/S1-GUNW-D-R-079-tops-20141116_20141023-030836-11425N_09393N-PP-3df1-v2_0_2.nc',
             'https://grfn.asf.alaska.edu/door/download/S1-GUNW-D-R-079-tops-20141116_20141023-030811-12922N_10895N-PP-a98f-v2_0_2.nc','https://grfn.asf.alaska.edu/door/download/S1-GUNW-D-R-079-tops-20141116_20141023-030746-14420N_12395N-PP-75c0-v2_0_2.nc',
             'https://grfn.asf.alaska.edu/door/download/S1-GUNW-D-R-079-tops-20141116_20141023-030719-16084N_14064N-PP-3545-v2_0_2.nc']


    # download the data into a toplevel folder
    localFiles = downloadProducts(downloadFiles,'products')
    
    # looping over all the different unitTests
    unitTests = ['unitTest1()', 'unitTests2()']
    if not isinstance(unitTests, list):
        unitTests = [unitTests]


    unitTest2(localFiles)
    #unitTest1(localFiles)
    #    unitTest2(localFiles)

    '''
    for unitTest in unitTests:
        # setting up the uwrapped and connected component directory that will be used for stitching
        UnwDir=os.path.join(workdir,unitTest,'unwrappedPhase')
        ConnCompDir=os.path.join(workdir,unitTest,'connectedComponents')
        mkdirs([UnwDir,ConnCompDir])
        
        # setting up the inputs for the unit tests
        input_dict = {}
        input_dict['files']=localFiles
        eval(unitTest,input_dict)
    
        pdb.set_trace()
    '''
    

