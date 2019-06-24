#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & David Bekaert
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
import sys
import numpy as np

from osgeo import gdal
gdal.UseExceptions()
#Suppress warnings
gdal.PushErrorHandler('CPLQuietErrorHandler')

# Import functions
from ARIAtools.shapefile_util import open_shapefile
from ARIAtools.unwrapStitching import product_stitch_overlap, product_stitch_2stage

_world_dem = ('https://cloud.sdsc.edu/v1/AUTH_opentopography/Raster/SRTM_GL1_Ellip/SRTM_GL1_Ellip_srtm.vrt')

def createParser():
    '''
        Extract specified product layers. The default is to export all layers.
    '''

    import argparse
    parser = argparse.ArgumentParser(description='Program to extract data and meta-data layers from ARIA standard GUNW products. Program will handle cropping/stiching when needed. By default, the program will crop all IFGs to bounds determined by the common intersection and bbox (if specified)')
    parser.add_argument('-f', '--file', dest='imgfile', type=str,
            required=True, help='ARIA file')
    parser.add_argument('-w', '--workdir', dest='workdir', default='./', help='Specify directory to deposit all outputs. Default is local directory where script is launched.')
    parser.add_argument('-l', '--layers', dest='layers', default=None, help='Specify layers to extract as a comma deliminated list bounded by single quotes. Allowed keys are: "unwrappedPhase", "coherence", "amplitude", "bPerpendicular", "bParallel", "incidenceAngle", "lookAngle","azimuthAngle". If "all" is specified, then all layers are extracted. If blank, will only extract bounding box.')
    parser.add_argument('-d', '--demfile', dest='demfile', type=str,
            default=None, help='DEM file. To download new DEM, specify "Download".')
    parser.add_argument('-p', '--projection', dest='projection', default='WGS84', type=str,
            help='projection for DEM. By default WGS84.')
    parser.add_argument('-b', '--bbox', dest='bbox', type=str, default=None, help="Provide either valid shapefile or Lat/Lon Bounding SNWE. -- Example : '19 20 -99.5 -98.5'")
    parser.add_argument('-m', '--mask', dest='mask', type=str, default=None, help="Provide valid mask file.")
    parser.add_argument('-sm', '--stichMethod', dest='stichMethodType',  type=str, default='overlap', help="Method applied to stich the unwrapped data. Either 'overlap', where product overlap is minimized, or '2stage', where minimization is done on connected components, are allowed methods. Default is 'overlap'.")
    parser.add_argument('-of', '--outputFormat', dest='outputFormat', type=str, default='VRT', help='GDAL compatible output format (e.g., "ENVI", "GTiff"). By default files are generated virtually except for "bPerpendicular", "bParallel", "incidenceAngle", "lookAngle","azimuthAngle", "unwrappedPhase" as these are require either DEM intersection or corrections to be applied')
    parser.add_argument('-croptounion', '--croptounion', action='store_true', dest='croptounion', help="If turned on, IFGs cropped to bounds based off of union and bbox (if specified). Program defaults to crop all IFGs to bounds based off of common intersection and bbox (if specified).")
    parser.add_argument('-verbose', '--verbose', action='store_true', dest='verbose', help="Toggle verbose mode on.")

    return parser

def cmdLineParse(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)

class InterpCube(object):
    '''
        Class to interpolate intersection of cube with DEM
    '''

    def __init__(self, inobj, hgtobj, latobj, lonobj):
        '''
            Init with h5py dataset.
        '''
        self.data = inobj[:]
        self.hgts = hgtobj[:]
        self.offset = None
        self.interp = []
        self.latobj = latobj[:]
        self.lonobj = lonobj[:]

        self.createInterp()

    def createInterp(self):
        '''
            Create interpolators.
        '''
        from scipy.interpolate import RectBivariateSpline
        self.offset = np.mean(self.data)
        for ind, hgt in enumerate(self.hgts):
            self.interp.append( RectBivariateSpline(self.latobj, self.lonobj, self.data[ind]-self.offset))

    def __call__(self, line, pix, h):
        '''
            Interpolate at a single point.
        '''
        from scipy.interpolate import interp1d

        vals  = np.array( [x(line,pix)[0,0] for x in self.interp])
        est = interp1d(self.hgts, vals, kind='cubic')
        return est(h) + self.offset

def prep_dem(demfilename, bbox_file, prods_TOTbbox, proj, arrshape=None, workdir='./', outputFormat='ENVI'):
    '''
        Function which load and export DEM, lat, lon arrays.
        If "Download" flag is specified, DEM will be donwloaded on the fly.
    '''

    # Get bounds of user bbox_file
    bounds=open_shapefile(bbox_file, 0, 0).bounds

    # File must be physically extracted, cannot proceed with VRT format. Defaulting to ENVI format.
    if outputFormat=='VRT':
       outputFormat='ENVI'

    # Download DEM
    if demfilename.lower()=='download':
        demfilename=os.path.join(workdir,'SRTM_3arcsec'+'.dem')
        gdal.Warp(demfilename, '/vsicurl/'+_world_dem, options=gdal.WarpOptions(format=outputFormat, outputBounds=bounds, outputType=gdal.GDT_Int16, width=arrshape[1], height=arrshape[0], dstNodata=-32768.0, srcNodata=-32768.0))
        gdal.Open(demfilename,gdal.GA_Update).SetProjection(proj)
        gdal.Translate(demfilename+'.vrt', demfilename, options=gdal.TranslateOptions(format="VRT")) #Make VRT

    # Load DEM and setup lat and lon arrays
    try:
        demfile = gdal.Warp('', demfilename, options=gdal.WarpOptions(format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=bounds, dstNodata=0))
        demfile.SetProjection(proj)

        ##SS Do we need lon lat if we would be doing gdal reproject using projection and transformation? See our earlier discussions.
        # Define lat/lon arrays for fullres layers
        Latitude=np.linspace(demfile.GetGeoTransform()[3],demfile.GetGeoTransform()[3]+(demfile.GetGeoTransform()[5]*demfile.RasterYSize),demfile.RasterYSize)
        Latitude=np.repeat(Latitude[:, np.newaxis], demfile.RasterXSize, axis=1)
        Longitude=np.linspace(demfile.GetGeoTransform()[0],demfile.GetGeoTransform()[0]+(demfile.GetGeoTransform()[1]*demfile.RasterXSize),demfile.RasterXSize)
        Longitude=np.repeat(Longitude[:, np.newaxis], demfile.RasterYSize, axis=1)
        Longitude=Longitude.transpose()
    except:
        raise Exception('Failed to open user DEM')

    return demfilename, demfile, Latitude, Longitude

def merged_productbbox(product_dict, workdir='./', bbox_file=None, croptounion=False):
    '''
        Extract/merge productBoundingBox layers for each pair and update dict, report common track bbox (default is to take common intersection, but user may specify union), and expected shape for DEM.
    '''

    # Import functions
    from ARIAtools.shapefile_util import save_shapefile

    # If specified workdir doesn't exist, create it
    if not os.path.exists(workdir):
        os.mkdir(workdir)

    # Extract/merge productBoundingBox layers
    for scene in product_dict:
        # Get pair name, expected in dictionary
        pair_name=scene["pair_name"][0]
        outname=os.path.join(workdir, pair_name+'.shp')

        # Create union of productBoundingBox layers
        for frame in scene["productBoundingBox"]:
            prods_bbox=open_shapefile(frame, 'productBoundingBox', 1)
            if os.path.exists(outname):
                union_bbox=open_shapefile(outname, 0, 0)
                prods_bbox=prods_bbox.union(union_bbox)
            save_shapefile(outname, prods_bbox, 'GeoJSON')              ##SS can we track and provide the proj information of the geojson?
        scene["productBoundingBox"]=[outname]

    prods_TOTbbox=os.path.join(workdir, 'productBoundingBox.shp')
    # Initiate intersection file with first product
    # this is for different scenes
    save_shapefile(prods_TOTbbox, open_shapefile(product_dict[0]['productBoundingBox'][0], 0, 0), 'GeoJSON')
    for scene in product_dict[1:]:
        prods_bbox=open_shapefile(scene['productBoundingBox'][0], 0, 0)
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

    # If bbox specified, intersect with common track intersection/union
    if bbox_file is not None:
        user_bbox=open_shapefile(bbox_file, 0, 0)
        total_bbox=open_shapefile(prods_TOTbbox, 0, 0)
        user_bbox=user_bbox.intersection(total_bbox)
        save_shapefile(prods_TOTbbox, user_bbox, 'GeoJSON')

        # Estimate percentage of overlap with bbox
        prods_bbox_area=open_shapefile(prods_TOTbbox, 0, 0).bounds
        prods_bbox_area=(max(prods_bbox_area[0],prods_bbox_area[2]) - min(prods_bbox_area[0],prods_bbox_area[2]))*(max(prods_bbox_area[1],prods_bbox_area[3]) - min(prods_bbox_area[1],prods_bbox_area[3]))
        bbox_area=open_shapefile(bbox_file, 0, 0).bounds
        bbox_area=(max(bbox_area[0],bbox_area[2]) - min(bbox_area[0],bbox_area[2]))*(max(bbox_area[1],bbox_area[3]) - min(bbox_area[1],bbox_area[3]))
        per_overlap=(prods_bbox_area/bbox_area)*100
        if per_overlap<50.:
            print("WARNING: Common track extent only has %d%% overlap with bbox"%per_overlap+'\n')
    else:
        bbox_file=prods_TOTbbox

    # Warp the first scene with the output-bounds defined above
    ds=gdal.Warp('', gdal.BuildVRT('', product_dict[0]['unwrappedPhase'][0]), options=gdal.WarpOptions(format="MEM", outputBounds=open_shapefile(bbox_file, 0, 0).bounds))
    # Get shape of full res layers
    arrshape=[ds.RasterYSize, ds.RasterXSize]
    # Get projection of full res layers
    proj=ds.GetProjection()
    ds = None

    return product_dict, bbox_file, prods_TOTbbox, arrshape, proj

def export_products(full_product_dict, bbox_file, prods_TOTbbox, layers, dem=None, lat=None, lon=None, mask=None, outDir='./',outputFormat='VRT', stichMethodType = 'overlap', verbose=None):
    """
    The function allows for exporting layer and 2D meta-data layers (at the product resolution).
    The function finalize_metadata is called to derive the 2D metadata layer. Dem/lat/lon arrays must be passed for this process.
    The key specifies which layer to extract from the dictionary.
    All products are cropped by the bounds from the input bbox_file, and clipped to the track extent denoted by the input prods_TOTbbox.
    Optionally, a user may pass a mask-file.
    """

    # loading dependencies
    from collections import OrderedDict
    if not layers: return # only bbox

    bounds=open_shapefile(bbox_file, 0, 0).bounds
    # Loop through user expected layers
    for key in layers:
        product_dict=[[j[key] for j in full_product_dict], [j["pair_name"] for j in full_product_dict]]
        workdir=os.path.join(outDir,key)

        # If specified workdir doesn't exist, create it
        if not os.path.exists(workdir):
            os.mkdir(workdir)

        # Iterate through all IFGs
        print("Generating: " + key)
        for i,j in enumerate(product_dict[0]):
            # provide update on which IFG
            print(product_dict[1][i][0])
            outname=os.path.abspath(os.path.join(workdir, product_dict[1][i][0]))

            # Extract/crop metadata layers
            if any(":/science/grids/imagingGeometry" in s for s in j):
                gdal.BuildVRT(outname +'.vrt', j)

                if dem is None:
                    raise Exception('No DEM input specified. Cannot extract 3D imaging geometry layers without DEM to intersect with.')

                # Check if height layers are consistent, and if not exit with error
                if len(set([gdal.Open(i).GetMetadataItem('NETCDF_DIM_heightsMeta_VALUES') for i in j]))==1:
                    gdal.Open(outname+'.vrt').SetMetadataItem('NETCDF_DIM_heightsMeta_VALUES',gdal.Open(j[0]).GetMetadataItem('NETCDF_DIM_heightsMeta_VALUES'))
                else:
                    raise Exception('Inconsistent heights for metadata layer(s) ', j, ' corresponding heights: ', [gdal.Open(i).GetMetadataItem('NETCDF_DIM_heightsMeta_VALUES') for i in j])

                # Pass metadata layer VRT, with DEM filename and output name to interpolate/intersect with DEM before cropping
                finalize_metadata(outname, open_shapefile(bbox_file, 0, 0).bounds, prods_TOTbbox, dem, lat, lon, mask, outputFormat, verbose=verbose)

            # Extract/crop full res layers, except for "unw" and "conn_comp" which requires advanced stiching
            elif key!='unwrappedPhase' and key!='connectedComponents':
                if outputFormat=='VRT' and mask is None:
                    # building the virtual vrt
                    gdal.BuildVRT(outname+ "_uncropped" +'.vrt', j)
                    # building the cropped vrt
                    gdal.Warp(outname+'.vrt', outname+"_uncropped"+'.vrt', options=gdal.WarpOptions(format=outputFormat, cutlineDSName=prods_TOTbbox, outputBounds=bounds))
                else:
                    # building the VRT
                    gdal.BuildVRT(outname +'.vrt', j)

                    # Mask specified, so file must be physically extracted, cannot proceed with VRT format. Defaulting to ENVI format.
                    if outputFormat=='VRT' and mask is not None:
                       outputFormat='ENVI'
                    gdal.Warp(outname, outname+'.vrt', options=gdal.WarpOptions(format=outputFormat, cutlineDSName=prods_TOTbbox, outputBounds=bounds))

                    # Update VRT
                    gdal.Translate(outname+'.vrt', outname, options=gdal.TranslateOptions(format="VRT"))

                    # Apply mask (if specified).
                    if mask is not None:
                        update_file=gdal.Open(outname,gdal.GA_Update)
                        update_file=update_file.GetRasterBand(1).WriteArray(mask*gdal.Open(outname+'.vrt').ReadAsArray())
                        update_file=None

            # Extract/crop "unw" and "conn_comp" layers leveraging the two stage unwrapper
            else:
                # Check if unw phase and connected components are already generated
                if not os.path.exists(os.path.join(outDir,'unwrappedPhase',product_dict[1][i][0])) or not os.path.exists(os.path.join(outDir,'connectedComponents',product_dict[1][i][0])):
                    # extract the inputs required for stitching of unwrapped and connected component files
                    unw_files = full_product_dict[i]['unwrappedPhase']
                    conn_files = full_product_dict[i]['connectedComponents']
                    prod_bbox_files = full_product_dict[i]['productBoundingBoxFrames']
                    # based on the key define the output directories
                    outFileUnw=os.path.join(outDir,'unwrappedPhase',product_dict[1][i][0])
                    outFileConnComp=os.path.join(outDir,'connectedComponents',product_dict[1][i][0])

                    # calling the stiching methods
                    if stichMethodType == 'overlap':
                        product_stitch_overlap(unw_files,conn_files,prod_bbox_files,bounds,prods_TOTbbox, outFileUnw=outFileUnw,outFileConnComp= outFileConnComp,mask=mask,outputFormat = outputFormat,verbose=verbose)
                    elif stichMethodType == '2stage':
                        product_stitch_2stage(unw_files,conn_files,bounds,prods_TOTbbox,outFileUnw=outFileUnw,outFileConnComp= outFileConnComp,mask=mask,outputFormat = outputFormat,verbose=verbose)

    return


def finalize_metadata(outname, bbox_bounds, prods_TOTbbox, dem, lat, lon, mask=None, outputFormat='ENVI', verbose=None):
    '''
        2D metadata layer is derived by interpolating and then intersecting 3D layers with a DEM. Lat/lon arrays must also be passed for this process.
    '''

    # import dependencies
    import scipy

    # Import functions
    from ARIAtools.vrtmanager import renderVRT

    # File must be physically extracted, cannot proceed with VRT format. Defaulting to ENVI format.
    if outputFormat=='VRT':
       outputFormat='ENVI'

    # Check and buffer bounds if <4pix in x/y, which would raise interpolation error
    geotrans=gdal.Open(outname+'.vrt').GetGeoTransform()
    if round((max(bbox_bounds[0::2])-min(bbox_bounds[0::2])),2)<round((abs(geotrans[1])*4),2) or round((max(bbox_bounds[1::2])-min(bbox_bounds[1::2])),2)<round((abs(geotrans[-1])*4),2):
        data_array=gdal.Warp('', outname+'.vrt', options=gdal.WarpOptions(format="MEM"))
    else:
        data_array=gdal.Warp('', outname+'.vrt', options=gdal.WarpOptions(format="MEM", outputBounds=bbox_bounds))

    # Define lat/lon/height arrays for metadata layers
    heightsMeta=np.array(gdal.Open(outname+'.vrt').GetMetadataItem('NETCDF_DIM_heightsMeta_VALUES')[1:-1].split(','), dtype='float32')
    ##SS Do we need lon lat if we would be doing gdal reproject using projection and transformation? See our earlier discussions.
    latitudeMeta=np.linspace(data_array.GetGeoTransform()[3],data_array.GetGeoTransform()[3]+(data_array.GetGeoTransform()[5]*data_array.RasterYSize),data_array.RasterYSize)
    longitudeMeta=np.linspace(data_array.GetGeoTransform()[0],data_array.GetGeoTransform()[0]+(data_array.GetGeoTransform()[1]*data_array.RasterXSize),data_array.RasterXSize)

    # First, using the height/latitude/longitude arrays corresponding to the metadata layer, set-up spatial 2D interpolator. Using this, perform vertical 1D interpolation on cube, and then use result to set-up a regular-grid-interpolator. Using this, pass DEM and full-res lat/lon arrays in order to get intersection with DEM.

    # 2D interpolation
    interp_2d = InterpCube(data_array.ReadAsArray(),heightsMeta,np.flip(latitudeMeta, axis=0),longitudeMeta)
    out_interpolated=np.zeros((heightsMeta.shape[0],latitudeMeta.shape[0],longitudeMeta.shape[0]))

    # 3D interpolation
    for iz, hgt in enumerate(heightsMeta):
        for iline, line in enumerate(latitudeMeta):
            for ipix, pixel in enumerate(longitudeMeta):
                out_interpolated[iz, iline, ipix] = interp_2d(line, pixel, hgt)
    out_interpolated=np.flip(out_interpolated, axis=0)
    # must use 'nearest' interpolation to avoid disprepancies between different crop extents
    interpolator = scipy.interpolate.RegularGridInterpolator((heightsMeta,np.flip(latitudeMeta, axis=0),longitudeMeta), out_interpolated, method='nearest', fill_value=data_array.GetRasterBand(1).GetNoDataValue())
    out_interpolated = interpolator(np.stack((np.flip(dem.ReadAsArray(), axis=0), lat, lon), axis=-1))

    # Save file
    renderVRT(outname, out_interpolated, geotrans=dem.GetGeoTransform(), drivername=outputFormat, gdal_fmt=data_array.ReadAsArray().dtype.name, proj=dem.GetProjection(), nodata=data_array.GetRasterBand(1).GetNoDataValue())

    # Since metadata layer extends at least one grid node outside of the expected track bounds, it must be cut to conform with these bounds.
    # Crop to track extents
    out_interpolated=gdal.Warp('', outname, options=gdal.WarpOptions(format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=bbox_bounds, dstNodata=data_array.GetRasterBand(1).GetNoDataValue())).ReadAsArray()

    # Apply mask (if specified).
    if mask is not None:
        ##SS Just to double check if the no-data is being tracked here. The vrt is setting a no-data value, but here no-data is not used.
        out_interpolated = np.multiply(out_interpolated, mask, out=out_interpolated, where=mask==0)

    update_file=gdal.Open(outname,gdal.GA_Update)
    update_file=update_file.GetRasterBand(1).WriteArray(out_interpolated)

    del out_interpolated, interpolator, interp_2d, data_array, update_file



def main(inps=None):
    '''
        Main workflow for extracting layers from ARIA products
    '''

    from ARIAtools.ARIAProduct import ARIA_standardproduct
    from ARIAtools.shapefile_util import open_shapefile

    print("***Extract Product Function:***")
    # if user bbox was specified, file(s) not meeting imposed spatial criteria are rejected.
    # Outputs = arrays ['standardproduct_info.products'] containing grouped “radarmetadata info” and “data layer keys+paths” dictionaries for each standard product
    # In addition, path to bbox file ['standardproduct_info.bbox_file'] (if bbox specified)
    standardproduct_info = ARIA_standardproduct(inps.imgfile, bbox=inps.bbox, workdir=inps.workdir, verbose=inps.verbose)

    if not inps.layers:
        print ('No layers specified; only creating bounding box shapes')

    elif inps.layers.lower()=='all':
        print('All layers are to be extracted, pass all keys.')
        inps.layers=list(standardproduct_info.products[1][0].keys())
        # Must remove productBoundingBoxes, because it's not a raster layer
        inps.layers.remove('productBoundingBox')
        inps.layers.remove('productBoundingBoxFrames')
        # Must remove pair_name, because it's not a raster layer
        inps.layers.remove('pair_name')

    else:
        inps.layers=list(inps.layers.split(','))
        valid = ["unwrappedPhase", "coherence", "amplitude", "bPerpendicular", "bParallel", "incidenceAngle", "lookAngle","azimuthAngle"]
        for layinp in inps.layers:
            if not layinp in valid:
                raise Exception('Invalid layer: {} entered. Options are any combination of: {} or all'.format(layinp, valid))

    # extract/merge productBoundingBox layers for each pair and update dict,
    # report common track bbox (default is to take common intersection, but user may specify union), and expected shape for DEM.
    standardproduct_info.products[1], standardproduct_info.bbox_file, prods_TOTbbox, arrshape, proj = merged_productbbox(standardproduct_info.products[1], os.path.join(inps.workdir,'productBoundingBox'), standardproduct_info.bbox_file, inps.croptounion)
    # Load mask (if specified).
    if inps.mask is not None:
        inps.mask=gdal.Warp('', inps.mask, options=gdal.WarpOptions(format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=open_shapefile(standardproduct_info.bbox_file, 0, 0).bounds, dstNodata=0))
        inps.mask.SetProjection(proj)
        # If no data value
        if inps.mask.GetRasterBand(1).GetNoDataValue():
            inps.mask=np.ma.masked_where(inps.mask.ReadAsArray() == inps.mask.GetRasterBand(1).GetNoDataValue(), inps.mask.ReadAsArray())
        else:
            inps.mask=inps.mask.ReadAsArray()


    # Download/Load DEM & Lat/Lon arrays, providing bbox, expected DEM shape, and output dir as input.
    if inps.demfile is not None:
        # Pass DEM-filename, loaded DEM array, and lat/lon arrays
        inps.demfile, demfile, Latitude, Longitude = prep_dem(inps.demfile, standardproduct_info.bbox_file, prods_TOTbbox, proj, arrshape=arrshape, workdir=inps.workdir, outputFormat=inps.outputFormat)
    else:
        demfile=None ; Latitude=None ; Longitude=None

    # Extract user expected layers
    export_products(standardproduct_info.products[1], standardproduct_info.bbox_file, prods_TOTbbox, inps.layers, dem=demfile, lat=Latitude, lon=Longitude, mask=inps.mask, outDir=inps.workdir,outputFormat=inps.outputFormat, stichMethodType = inps.stichMethodType, verbose=inps.verbose)
