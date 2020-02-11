import os
from osgeo import gdal

from ARIAtools.unwrapStitching import testproduct_stitch_2stage
from ARIAtools.shapefile_util import open_shapefile, save_shapefile
from ARIAtools.extractProduct import merged_productbbox
import numpy as np

#set workdirectory and make directories
workdir=os.path.abspath('unittest')
os.mkdir(workdir)
productworkdir=os.path.join(workdir,'products')
os.mkdir(productworkdir)
os.mkdir(os.path.join(workdir,'unwrappedPhase'))
os.mkdir(os.path.join(workdir,'connectedComponents'))
os.mkdir(os.path.join(workdir,'productBoundingBox'))

#set some input variables
croptounion=True
outFileUnw=os.path.join(workdir,'unwrappedPhase/20141116_20141023')
outFileConnComp=os.path.join(workdir,'connectedComponents/20141116_20141023')
mask=None
outputFormat='ISCE'
verbose=True

#get products
os.chdir(productworkdir)
os.system('wget https://aria-products.jpl.nasa.gov/search/dataset/grq_v2.0.2_s1-gunw-released/S1-GUNW-D-R-079-tops-20141116_20141023-030811-12922N_10895N-PP-a98f-v2_0_2/S1-GUNW-D-R-079-tops-20141116_20141023-030811-12922N_10895N-PP-a98f-v2_0_2.nc')
os.system('wget https://aria-products.jpl.nasa.gov/search/dataset/grq_v2.0.2_s1-gunw-released/S1-GUNW-D-R-079-tops-20141116_20141023-030746-14420N_12395N-PP-75c0-v2_0_2/S1-GUNW-D-R-079-tops-20141116_20141023-030746-14420N_12395N-PP-75c0-v2_0_2.nc')
os.system('wget https://aria-products.jpl.nasa.gov/search/dataset/grq_v2.0.2_s1-gunw-released/S1-GUNW-D-R-079-tops-20141210_20141023-030745-14420N_12396N-PP-ec61-v2_0_2/S1-GUNW-D-R-079-tops-20141210_20141023-030745-14420N_12396N-PP-ec61-v2_0_2.nc')
os.system('wget https://aria-products.jpl.nasa.gov/search/dataset/grq_v2.0.2_s1-gunw-released/S1-GUNW-D-R-079-tops-20141210_20141023-030810-12923N_10895N-PP-477c-v2_0_2/S1-GUNW-D-R-079-tops-20141210_20141023-030810-12923N_10895N-PP-477c-v2_0_2.nc')
os.chdir(workdir)


#get more input variables
unw_files=['S1-GUNW-D-R-079-tops-20141116_20141023-030746-14420N_12395N-PP-75c0-v2_0_2.nc','S1-GUNW-D-R-079-tops-20141116_20141023-030811-12922N_10895N-PP-a98f-v2_0_2.nc']
unw_files=[os.path.join(productworkdir,i) for i in unw_files]
conn_files=['NETCDF:"%s":/science/grids/data/connectedComponents'%(i) for i in unw_files]
prod_bbox_files=['NETCDF:"%s":productBoundingBox'%(i) for i in unw_files]
unw_files=['NETCDF:"%s":/science/grids/data/unwrappedPhase'%(i) for i in unw_files]


#get bounding boxes
#bounds=(38.9649142686367, 10.8963088229126, 41.8691085466474, 14.4058259783288)
prods_TOTbbox=os.path.join(workdir, 'productBoundingBox/productBoundingBox.shp')
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

# overwrite unw-files as simulations
sim_unw_files=[]
sim_conn_files=[]
for i in enumerate(unw_files):
# extracting phase data (1 band file is the target, so if you are using ENVI make sure to take the right band!)
    file_unw=outFileUnw+'part%s_phase'%(str(i[0]))
    # building the VRT
    gdal.BuildVRT(file_unw +'.vrt', i[1])
    gdal.Warp(file_unw, file_unw +'.vrt', options=gdal.WarpOptions(format=outputFormat))
    gdal.BuildVRT(file_unw +'.vrt', file_unw)
    sim_unw_files.append(file_unw)
    # extracting connected component data (1 band file is the target!)
    file_conn=outFileConnComp+'part%s_connComp'%(str(i[0]))
    gdal.BuildVRT(file_conn +'.vrt', conn_files[i[0]])
    gdal.Warp(file_conn, file_conn +'.vrt', options=gdal.WarpOptions(format=outputFormat))
    gdal.BuildVRT(file_conn +'.vrt', file_conn)
    sim_conn_files.append(file_conn)

    #create simulations
    # doing the phase     
    data = gdal.Open(file_unw,gdal.GA_ReadOnly)
    data_band = data.GetRasterBand(1)
    phase = data_band.ReadAsArray()
    data = None

    #get 0. phase indices
    phase_0s=np.nonzero(phase==0.)

    width=int(phase.shape[1]/2)
    length=int(phase.shape[0]/2)

    # doing the connected component
    data = gdal.Open(file_conn,gdal.GA_ReadOnly)
    data_band = data.GetRasterBand(1)
    connComp = data_band.ReadAsArray() 
    data = None
    #get 0/1 conncomp indices
    connComp_0s=np.nonzero(connComp==0)
    connComp_nodatas=np.nonzero(connComp==-1)

    # manipulating the connected component
    phase = phase*0.0+np.pi
    phase[0:length,:] = phase[0:length,:]+(2+i[0])*np.pi
    phase[length-1:,:] = phase[length-1:,:]+(-4+i[0])*np.pi
    phase[phase_0s]=0.

    connComp = connComp*0 + 1 
    connComp[length-1:,:]  = connComp[length-1:,:] +1
    connComp[connComp_0s]=0
    connComp[connComp_nodatas]=-1
    connComp = connComp.astype('int')

    # writing put new phase 
    ds=gdal.Open(file_unw,gdal.GA_Update)
    ds.GetRasterBand(1).WriteArray(phase)
    del ds

    # writing out new conncomp 
    ds=gdal.Open(file_conn,gdal.GA_Update)
    ds.GetRasterBand(1).WriteArray(connComp)
    del ds

#run test
testproduct_stitch_2stage(sim_unw_files,sim_conn_files,prod_bbox_files,bounds,prods_TOTbbox, outFileUnw=outFileUnw, outFileConnComp=outFileConnComp, mask=mask, outputFormat=outputFormat,verbose=verbose)
