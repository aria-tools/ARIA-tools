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

def createParser():
    '''
        Extract specified product layers. The default is to export all layers.
    '''

    import argparse
    parser = argparse.ArgumentParser(description='Program to extract data and meta-data layers from ARIA standard GUNW products. Program will handle cropping/stitching when needed. By default, the program will crop all IFGs to bounds determined by the common intersection and bbox (if specified)')
    parser.add_argument('-f', '--file', dest='imgfile', type=str,
            required=True, help='ARIA file')
    parser.add_argument('-w', '--workdir', dest='workdir', default='./', help='Specify directory to deposit all outputs. Default is local directory where script is launched.')
    parser.add_argument('-tp', '--tropo_products', dest='tropo_products', type=str, default=None, help='Path to director(ies) or tar file(s) containing GACOS products.')
    parser.add_argument('-l', '--layers', dest='layers', default=None, help='Specify layers to extract as a comma deliminated list bounded by single quotes. Allowed keys are: "unwrappedPhase", "coherence", "amplitude", "bPerpendicular", "bParallel", "incidenceAngle", "lookAngle","azimuthAngle". If "all" is specified, then all layers are extracted. If blank, will only extract bounding box.')
    parser.add_argument('-d', '--demfile', dest='demfile', type=str,
            default=None, help='DEM file. To download new DEM, specify "Download".')
    parser.add_argument('-p', '--projection', dest='projection', default='WGS84', type=str,
            help='projection for DEM. By default WGS84.')
    parser.add_argument('-b', '--bbox', dest='bbox', type=str, default=None, help="Provide either valid shapefile or Lat/Lon Bounding SNWE. -- Example : '19 20 -99.5 -98.5'")
    parser.add_argument('-m', '--mask', dest='mask', type=str, default=None, help="Path to mask file or 'Download'. File needs to be GDAL compatabile, contain spatial reference information, and have invalid/valid data represented by 0/1, respectively. If 'Download', will use GSHHS water mask")
    parser.add_argument('-at', '--amp_thresh', dest='amp_thresh', default=None, type=str, help='Amplitude threshold below which to mask. Specify "None" to not use amplitude mask. By default "None".')
#    parser.add_argument('-sm', '--stitchMethod', dest='stitchMethodType',  type=str, default='overlap', help="Method applied to stitch the unwrapped data. Either 'overlap', where product overlap is minimized, or '2stage', where minimization is done on connected components, are allowed methods. Default is 'overlap'.")
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
        Function to load and export DEM, lat, lon arrays.
        If "Download" flag is specified, DEM will be donwloaded on the fly.
    '''

    _world_dem = '/vsicurl/https://cloud.sdsc.edu/v1/AUTH_opentopography/Raster/SRTM_GL1_Ellip/SRTM_GL1_Ellip_srtm.vrt'

    # Get bounds of user bbox_file
    bounds=open_shapefile(bbox_file, 0, 0).bounds

    # File must be physically extracted, cannot proceed with VRT format. Defaulting to ENVI format.
    if outputFormat=='VRT':
       outputFormat='ENVI'

    # Download DEM
    if demfilename.lower()=='download':
        demfilename=os.path.join(workdir,'SRTM_3arcsec'+'.dem')
        gdal.Warp(demfilename, _world_dem, options=gdal.WarpOptions(format=outputFormat, outputBounds=bounds, outputType=gdal.GDT_Int16, width=arrshape[1], height=arrshape[0], dstNodata=-32768.0, srcNodata=-32768.0))
        gdal.Open(demfilename,gdal.GA_Update).SetProjection(proj)
        gdal.Translate(demfilename+'.vrt', demfilename, options=gdal.TranslateOptions(format="VRT")) #Make VRT
        print('Downloaded 3 arc-sec SRTM DEM here: '+ demfilename)

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

def prep_mask(product_dict, maskfilename, bbox_file, prods_TOTbbox, proj, amp_thresh=None, arrshape=None, workdir='./', outputFormat='ENVI'):
    '''
        Function to load and export mask file.
        If "Download" flag is specified, GSHHS water mask will be donwloaded on the fly.
    '''

    # Import functions
    from ARIAtools.vrtmanager import renderVRT
    import glob

    _world_watermask = [' /vsizip/vsicurl/http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip/GSHHS_shp/f/GSHHS_f_L1.shp',' /vsizip/vsicurl/http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip/GSHHS_shp/f/GSHHS_f_L2.shp',' /vsizip/vsicurl/http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip/GSHHS_shp/f/GSHHS_f_L3.shp', ' /vsizip/vsicurl/http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip/GSHHS_shp/f/GSHHS_f_L4.shp',' /vsizip/vsicurl/https://osmdata.openstreetmap.de/download/land-polygons-complete-4326.zip/land-polygons-complete-4326/land_polygons.shp']

    # Get bounds of user bbox_file
    bounds=open_shapefile(bbox_file, 0, 0).bounds

    # File must be physically extracted, cannot proceed with VRT format. Defaulting to ENVI format.
    if outputFormat=='VRT':
       outputFormat='ENVI'

    # Download mask
    if maskfilename.lower()=='download':
        maskfilename=os.path.join(workdir,'watermask'+'.msk')

        ###Make coastlines/islands union VRT
        os.system('ogrmerge.py -o ' + os.path.join(workdir,'watermsk_shorelines.vrt') + ''.join(_world_watermask[::2]) + ' -field_strategy Union -f VRT -single')

        ###Make lakes/ponds union VRT
        os.system('ogrmerge.py -o ' + os.path.join(workdir,'watermsk_lakes.vrt') + ''.join(_world_watermask[1::2]) + ' -field_strategy Union -f VRT -single')

        ###Initiate water-mask with coastlines/islands union VRT
        gdal.Rasterize(maskfilename, os.path.join(workdir,'watermsk_shorelines.vrt'), options=gdal.RasterizeOptions(format=outputFormat, outputBounds=bounds, outputType=gdal.GDT_Byte, width=arrshape[1], height=arrshape[0], burnValues=[1], layers='merged'))
        gdal.Open(maskfilename,gdal.GA_Update).SetProjection(proj)
        gdal.Translate(maskfilename+'.vrt', maskfilename, options=gdal.TranslateOptions(format="VRT"))

        ###Must take inverse of lakes/ponds union because of opposite designation (1 for water, 0 for land) as desired (0 for water, 1 for land)
        lake_masks=gdal.Rasterize('', os.path.join(workdir,'watermsk_lakes.vrt'), options=gdal.RasterizeOptions(format='MEM', outputBounds=bounds, outputType=gdal.GDT_Byte, width=arrshape[1], height=arrshape[0], burnValues=[1], layers='merged', inverse=True))
        lake_masks.SetProjection(proj)
        lake_masks=lake_masks.ReadAsArray()

        if amp_thresh:
            ###Make average amplitude mask
            # Iterate through all IFGs
            for i,j in enumerate(product_dict[0]):
                amp_file=gdal.Warp('', j, options=gdal.WarpOptions(format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=bounds))
                amp_file_arr=np.ma.masked_where(amp_file.ReadAsArray() == amp_file.GetRasterBand(1).GetNoDataValue(), amp_file.ReadAsArray())

                # Iteratively update average amplitude file
                # If looping through first amplitude file, nothing to sum so just save to file
                if os.path.exists(os.path.join(workdir,'avgamplitude.msk')):
                    amp_file=gdal.Open(os.path.join(workdir,'avgamplitude.msk'),gdal.GA_Update)
                    amp_file=amp_file.GetRasterBand(1).WriteArray(amp_file_arr+amp_file.ReadAsArray())
                else:
                    renderVRT(os.path.join(workdir,'avgamplitude.msk'), amp_file_arr, geotrans=amp_file.GetGeoTransform(), drivername=outputFormat, gdal_fmt=amp_file_arr.dtype.name, proj=amp_file.GetProjection(), nodata=amp_file.GetRasterBand(1).GetNoDataValue())
                amp_file = coh_val = amp_file_arr = None

            # Take average of amplitude sum
            amp_file=gdal.Open(os.path.join(workdir,'avgamplitude.msk'),gdal.GA_Update)
            arr_mean = amp_file.ReadAsArray()/len(product_dict[0])
            arr_mean = np.where(arr_mean < float(amp_thresh), 0, 1)
            amp_file=amp_file.GetRasterBand(1).WriteArray(arr_mean)
            amp_file = None ; arr_mean = None
            amp_file=gdal.Open(os.path.join(workdir,'avgamplitude.msk')).ReadAsArray()
        else:
            amp_file=np.ones((lake_masks.shape[0],lake_masks.shape[1]))

        ###Update water-mask with lakes/ponds union and average amplitude
        update_file=gdal.Open(maskfilename,gdal.GA_Update)
        update_file=update_file.GetRasterBand(1).WriteArray(update_file.ReadAsArray()*lake_masks*amp_file)
        print('Downloaded water-mask here: '+ maskfilename)
        update_file=None ; lake_masks=None; amp_file = None
        #Delete temp files
        os.remove(os.path.join(workdir,'watermsk_shorelines.vrt')); os.remove(os.path.join(workdir,'watermsk_lakes.vrt'))
        for i in glob.glob(os.path.join(workdir,'avgamplitude.msk*')): os.remove(i)

    # Load mask
    try:
        mask=gdal.Warp('', maskfilename, options=gdal.WarpOptions(format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=bounds, dstNodata=0))
        mask.SetProjection(proj)
        # If no data value
        if mask.GetRasterBand(1).GetNoDataValue():
            mask=np.ma.masked_where(mask.ReadAsArray() == mask.GetRasterBand(1).GetNoDataValue(), mask.ReadAsArray())
        else:
            mask=mask.ReadAsArray()
    except:
        raise Exception('Failed to open user mask')

    return mask

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

def export_products(full_product_dict, bbox_file, prods_TOTbbox, layers, dem=None, lat=None, lon=None, mask=None, outDir='./',outputFormat='VRT', stitchMethodType='overlap', verbose=None):
    """
        Export layer and 2D meta-data layers (at the product resolution).
        The function finalize_metadata is called to derive the 2D metadata layer. Dem/lat/lon arrays must be passed for this process.
        The keys specify which layer to extract from the dictionary.
        All products are cropped by the bounds from the input bbox_file, and clipped to the track extent denoted by the input prods_TOTbbox.
        Optionally, a user may pass a mask-file.
    """

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

            # Extract/crop full res layers, except for "unw" and "conn_comp" which requires advanced stitching
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

                    # calling the stitching methods
                    if stitchMethodType == 'overlap':
                        product_stitch_overlap(unw_files,conn_files,prod_bbox_files,bounds,prods_TOTbbox, outFileUnw=outFileUnw,outFileConnComp= outFileConnComp,mask=mask,outputFormat = outputFormat,verbose=verbose)
                    elif stitchMethodType == '2stage':
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
    # interpolate to interferometric grid
    interpolator = scipy.interpolate.RegularGridInterpolator((heightsMeta,np.flip(latitudeMeta, axis=0),longitudeMeta), out_interpolated, method='linear', fill_value=data_array.GetRasterBand(1).GetNoDataValue())
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


def tropo_correction(full_product_dict, tropo_products, bbox_file, prods_TOTbbox, outDir='./',outputFormat='VRT', verbose=None):
    """
        Perform tropospheric corrections. Must provide valid path to GACOS products.
        All products are cropped by the bounds from the input bbox_file, and clipped to the track extent denoted by the input prods_TOTbbox.
    """

    # Import functions
    from ARIAtools.vrtmanager import renderVRT
    import glob
    import tarfile
    from datetime import datetime
    from ARIAtools.shapefile_util import save_shapefile
    from shapely.geometry import Polygon

    # File must be physically extracted, cannot proceed with VRT format. Defaulting to ENVI format.
    if outputFormat=='VRT':
       outputFormat='ENVI'

    user_bbox=open_shapefile(bbox_file, 0, 0)
    bounds=open_shapefile(bbox_file, 0, 0).bounds
    prods_bbox_area=(max(bounds[0],bounds[2]) - min(bounds[0],bounds[2]))*(max(bounds[1],bounds[3]) - min(bounds[1],bounds[3]))

    product_dict=[[j['unwrappedPhase'] for j in full_product_dict[1]], [j['lookAngle'] for j in full_product_dict[1]], [j["pair_name"] for j in full_product_dict[1]]]
    metadata_dict=[[j['azimuthZeroDopplerStartTime'] for j in full_product_dict[0]], [j['azimuthZeroDopplerEndTime'] for j in full_product_dict[0]], [j['wavelength'] for j in full_product_dict[0]]]
    workdir=os.path.join(outDir,'tropocorrected_products')

    # If specified workdir doesn't exist, create it
    if not os.path.exists(workdir):
        os.mkdir(workdir)

    # Get list of all dates for which standard products exist
    date_list=[]
    for i in product_dict[2]:
        date_list.append(i[0][:8]); date_list.append(i[0][9:])
    date_list=list(set(date_list))

    ### Determine if tropo input is single directory, a list, or wildcard.
    # If list of directories/files
    if len([str(val) for val in tropo_products.split(',')])>1:
        tropo_products=[str(val) for val in tropo_products.split(',')]
    # If single directory or wildcard
    else:
        # If single directory/file
        if os.path.exists(tropo_products):
            tropo_products=[tropo_products]
        # If wildcard
        else:
            tropo_products=glob.glob(os.path.expanduser(os.path.expandvars(tropo_products)))
        # Convert relative paths to absolute paths
        tropo_products=[os.path.normpath(os.path.join(os.getcwd(),i)) if not os.path.isabs(i) else i for i in tropo_products]
    if len(tropo_products)==0:
        raise Exception('No file match found')

    ###Extract tarfiles
    # Setup dictionary to track for products that are to be merged
    tropo_date_dict={}
    for i in date_list: tropo_date_dict[i]=[] ; tropo_date_dict[i+"_UTC"]=[]
    for i,j in enumerate(tropo_products):
        if not os.path.isdir(j):
            untar_dir=os.path.join(os.path.abspath(os.path.join(j, os.pardir)),os.path.basename(j).split('.')[0]+'_extracted')
            if not tarfile.is_tarfile(j):
                raise Exception('Cannot extract %s because it is not a valid tarfile. Resolve this and relaunch'%(j))
            if os.path.exists(untar_dir):
                raise Exception('Cannot extract %s to %s because path already exists. Resolve conflict and relaunch'%(os.path.basename(j),untar_dir))
            print('Extracting GACOS tarfile %s to %s.'%(os.path.basename(j),untar_dir))
            tarfile.open(j).extractall(path=untar_dir)
            tropo_products[i]=untar_dir
        # Loop through each GACOS product file
        for k in glob.glob(os.path.join(tropo_products[i],'*.ztd')):
            tropo_date_dict[os.path.basename(k)[:-4]].append(k)
            tropo_date_dict[os.path.basename(k)[:-4]+"_UTC"].append(os.path.basename(k)[:4]+'-'+os.path.basename(k)[4:6]+'-'+os.path.basename(k)[6:8]+'-'+open(k+'.rsc', 'r').readlines()[-1].split()[1])
            # make corresponding VRT file, if it doesn't exist
            if not os.path.exists(k+'.vrt') and os.path.basename(k)[:-4] in date_list:
                tropo_rsc_dict={}
                for line in open(k+'.rsc', 'r').readlines(): tropo_rsc_dict[line.split()[0]]=line.split()[1]
                gacos_prod=np.fromfile(k, dtype='float32').reshape(int(tropo_rsc_dict['FILE_LENGTH']),int(tropo_rsc_dict['WIDTH']))
                # Save as GDAL file, using proj from first unwrappedPhase file
                renderVRT(k, gacos_prod, geotrans=(float(tropo_rsc_dict['X_FIRST']), float(tropo_rsc_dict['X_STEP']), 0.0, float(tropo_rsc_dict['Y_FIRST']), 0.0, float(tropo_rsc_dict['Y_STEP'])), drivername=outputFormat, gdal_fmt='float32', proj=gdal.Open(os.path.join(outDir,'unwrappedPhase',product_dict[2][0][0])).GetProjection(), nodata=0.)
                gacos_prod = None
                if verbose:
                    print("GACOS product %s successfully converted to GDAL-readable raster"%(k))

    # If multiple GACOS directories, merge products.
    if len(tropo_products)>1:
        print('Stitching/storing GACOS products in %s.'%(tropo_products))
        tropo_products=os.path.join(outDir,'merged_GACOS')
        # If specified merged directory doesn't exist, create it
        if not os.path.exists(os.path.join(outDir,'merged_GACOS')):
            os.mkdir(os.path.join(outDir,'merged_GACOS'))

        for i in tropo_date_dict:
            if 'UTC' not in i:
                outname=os.path.join(outDir,'merged_GACOS',i+'.ztd')
                # building the VRT
                gdal.BuildVRT(outname +'.vrt', tropo_date_dict[i])
                gdal.Warp(outname, outname+'.vrt', options=gdal.WarpOptions(format=outputFormat))
                # Update VRT
                gdal.Translate(outname+'.vrt', outname, options=gdal.TranslateOptions(format="VRT"))
                geotrans=gdal.Open(outname).GetGeoTransform()
                # Create merge rsc file
                with open(outname+'.rsc','w') as merged_rsc:
                    merged_rsc.write('WIDTH %s\n'%(gdal.Open(outname).ReadAsArray().shape[1])) ; merged_rsc.write('FILE_LENGTH %s\n'%(gdal.Open(outname).ReadAsArray().shape[0]))
                    merged_rsc.write('XMIN %s\n'%(0)) ; merged_rsc.write('XMAX %s\n'%(gdal.Open(outname).ReadAsArray().shape[1]))
                    merged_rsc.write('YMIN %s\n'%(0)) ; merged_rsc.write('YMAX %s\n'%(gdal.Open(outname).ReadAsArray().shape[0]))
                    merged_rsc.write('X_FIRST %f\n'%(geotrans[0])) ; merged_rsc.write('Y_FIRST %f\n'%(geotrans[3]))
                    merged_rsc.write('X_STEP %f\n'%(geotrans[1])) ; merged_rsc.write('Y_STEP %f\n'%(geotrans[-1]))
                    merged_rsc.write('X_UNIT %s\n'%('degres')) ; merged_rsc.write('Y_UNIT %s\n'%('degres'))
                    merged_rsc.write('Z_OFFSET %s\n'%(0)) ; merged_rsc.write('Z_SCALE %s\n'%(1))
                    merged_rsc.write('PROJECTION %s\n'%('LATLON')) ; merged_rsc.write('DATUM %s\n'%('WGS84'))
                    merged_rsc.write('TIME_OF_DAY %s\n'%(''.join(tropo_date_dict[i+"_UTC"])))
    else:
        tropo_products=tropo_products[0]

    # Estimate percentage of overlap with tropospheric product
    for i in glob.glob(os.path.join(tropo_products,'*.ztd')):
        # create shapefile
        geotrans=gdal.Open(i).GetGeoTransform()
        bbox=[geotrans[0], geotrans[3]+(gdal.Open(i).ReadAsArray().shape[0]*geotrans[-1]),geotrans[0]+(gdal.Open(i).ReadAsArray().shape[1]*geotrans[1]),geotrans[3]]
        bbox=Polygon(np.column_stack((np.array([bbox[2],bbox[3],bbox[3],bbox[2],bbox[2]]),
                            np.array([bbox[0],bbox[0],bbox[1],bbox[1],bbox[0]]))))
        save_shapefile(i+'.shp', bbox, 'GeoJSON')
        bbox_area=open_shapefile(i+'.shp', 0, 0)
        bbox_area=user_bbox.intersection(bbox_area)
        bbox_area=bbox_area.bounds
        bbox_area=(max(bbox_area[0],bbox_area[2]) - min(bbox_area[0],bbox_area[2]))*(max(bbox_area[1],bbox_area[3]) - min(bbox_area[1],bbox_area[3]))
        per_overlap=(prods_bbox_area/bbox_area)*100
        if per_overlap!=100. and per_overlap!=0.:
            print("WARNING: Common track extent only has %d%% overlap with bbox"%per_overlap+'\n')
        if per_overlap==0.:
            raise Exception('No spatial overlap between tropospheric product %s and defined bounding box. Resolve conflict and relaunch'%(i))

    # Iterate through all IFGs and apply corrections
    for i,j in enumerate(product_dict[0]):
        outname=os.path.join(workdir,product_dict[2][i][0])
        unwname=os.path.join(outDir,'unwrappedPhase',product_dict[2][i][0])
        tropo_reference=os.path.join(tropo_products,product_dict[2][i][0][:8]+'.ztd.vrt')
        tropo_secondary=os.path.join(tropo_products,product_dict[2][i][0][9:]+'.ztd.vrt')
        if os.path.exists(tropo_reference) and os.path.exists(tropo_secondary):
            # Check if tropo products are temporally consistent with IFG
            for k,l in enumerate([tropo_reference, tropo_secondary]):
                # Get ARIA product times
                aria_rsc_dict={}
                aria_rsc_dict['azimuthZeroDopplerStartTime']=[datetime.strptime(os.path.basename(l)[:4]+'-'+os.path.basename(l)[4:6]+'-'+os.path.basename(l)[6:8]+'-'+m[11:], "%Y-%m-%d-%H:%M:%S.%fZ") for m in metadata_dict[0][0]]
                aria_rsc_dict['azimuthZeroDopplerEndTime']=[datetime.strptime(os.path.basename(l)[:4]+'-'+os.path.basename(l)[4:6]+'-'+os.path.basename(l)[6:8]+'-'+m[11:], "%Y-%m-%d-%H:%M:%S.%fZ") for m in metadata_dict[1][0]]
                # Get tropo product UTC times
                tropo_rsc_dict={}
                tropo_rsc_dict['TIME_OF_DAY']=open(l[:-4]+'.rsc', 'r').readlines()[-1].split()[1].split('UTC')[:-1]
                # If stitched tropo product, must account for date change (if applicable)
                if len(tropo_rsc_dict['TIME_OF_DAY'])>1:
                    tropo_rsc_dict['TIME_OF_DAY']=[datetime.strptime(m[:13]+'-'+str(round(float(m[13:])*60)), "%Y-%m-%d-%H-%M") for m in tropo_rsc_dict['TIME_OF_DAY']]
                else:
                    tropo_rsc_dict['TIME_OF_DAY']=[datetime.strptime(os.path.basename(l)[:4]+'-'+os.path.basename(l)[4:6]+'-'+os.path.basename(l)[6:8]+'-'+tropo_rsc_dict['TIME_OF_DAY'][0][:2]+'-'+str(round(float(tropo_rsc_dict['TIME_OF_DAY'][0][2:])*60)), "%Y-%m-%d-%H-%M")]

                # Check and report if tropospheric product falls outside of standard product range
                latest_start = max(aria_rsc_dict['azimuthZeroDopplerStartTime']+[min(tropo_rsc_dict['TIME_OF_DAY'])])
                earliest_end = min(aria_rsc_dict['azimuthZeroDopplerEndTime']+[max(tropo_rsc_dict['TIME_OF_DAY'])])
                delta = (earliest_end - latest_start).total_seconds() + 1
                if delta<0:
                    print("WARNING: tropospheric product was generated %f secs outside of acquisition interval for scene %s in IFG %s"%(abs(delta), os.path.basename(l)[:8], product_dict[2][i][0]))

            # Open unwrappedPhase and mask nodata
            unwphase=gdal.Open(unwname)
            geotrans=unwphase.GetGeoTransform() ; proj=unwphase.GetProjection() ; unwnodata=unwphase.GetRasterBand(1).GetNoDataValue()
            unwphase=unwphase.ReadAsArray()
            unwphase=np.ma.masked_where(unwphase == unwnodata, unwphase)
            # Get reference point and difference phase with it
            ref_point_ind=abs(np.subtract(unwphase,unwphase.mean()))
            ref_point_ind=np.unravel_index(ref_point_ind.argmin(),ref_point_ind.shape)
            unwphase=np.subtract(unwphase,unwphase[ref_point_ind])
            #!#renderVRT(os.path.join(workdir,'ref'+product_dict[2][i][0]), unwphase, geotrans=geotrans, drivername=outputFormat, gdal_fmt='float32', proj=proj, nodata=unwnodata)
            if verbose:
                print("Reference point for IFG %s set to pixel [y=%sy,x=%s], i.e. coord [lat=%s\xb0,lon=%s\xb0]"%(product_dict[2][i][0],str(ref_point_ind[0]),str(ref_point_ind[1]),(geotrans[3]+(ref_point_ind[0]*geotrans[-1])),(geotrans[0]+(ref_point_ind[1]*geotrans[1]))))

            # Open corresponding tropo products and pass the difference
            tropo_product=gdal.Warp('', tropo_reference, options=gdal.WarpOptions(format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=bounds, dstNodata=0.)).ReadAsArray()
            tropo_product=np.ma.masked_where(tropo_product == 0., tropo_product)
            tropo_secondary=gdal.Warp('', tropo_secondary, options=gdal.WarpOptions(format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=bounds, dstNodata=0.)).ReadAsArray()
            tropo_secondary=np.ma.masked_where(tropo_secondary == 0., tropo_secondary)
            tropo_product=np.subtract(tropo_secondary,tropo_product)
            # Difference with reference point
            tropo_product=np.subtract(tropo_product,tropo_product[ref_point_ind])
            if np.ma.is_masked(tropo_product[ref_point_ind]):
                print("WARNING: Must skip IFG %s, because selected reference point [y=%sy,x=%s] is not a valid point in tropospheric product"%(product_dict[2][i][0],str(ref_point_ind[0]),str(ref_point_ind[1])))
                continue

            # Convert troposphere to rad
            tropo_product=np.divide(tropo_product,float(metadata_dict[2][i][0])/(4*np.pi))
            # Account for lookAngle
            lookfile=gdal.Open(os.path.join(outDir,'lookAngle',product_dict[2][i][0])).ReadAsArray()
            lookfile=np.sin(np.deg2rad(np.ma.masked_where(lookfile == 0., lookfile)))
            tropo_product=np.divide(tropo_product,lookfile)

            #Correct phase and save to file
            unwphase=np.subtract(unwphase,tropo_product)
            np.ma.set_fill_value(unwphase, unwnodata); np.ma.set_fill_value(tropo_product, 0.)
            renderVRT(outname+'_tropodiff', tropo_product.filled(), geotrans=geotrans, drivername=outputFormat, gdal_fmt='float32', proj=proj, nodata=0.)
            renderVRT(outname, unwphase.filled(), geotrans=geotrans, drivername=outputFormat, gdal_fmt='float32', proj=proj, nodata=unwnodata)

            del unwphase, tropo_product, tropo_reference, tropo_secondary, lookfile

        else:
            print("WARNING: Must skip IFG %s, because the tropospheric products corresponding to the reference and/or secondary products are not found in the specified folder %s"%(product_dict[2][i][0],tropo_products))


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
        print('No layers specified; only creating bounding box shapes')

    elif inps.layers.lower()=='all':
        print('All layers are to be extracted, pass all keys.')
        inps.layers=list(standardproduct_info.products[1][0].keys())
        # Must remove productBoundingBoxes & pair-names because they are not raster layers
        inps.layers=[i for i in inps.layers if i not in ['productBoundingBox','productBoundingBoxFrames','pair_name']]

    elif inps.tropo_products!='None':
        print('Tropospheric corrections will be applied, making sure at least unwrappedPhase and lookAngle are extracted.')
        if isinstance(inps.layers,str): inps.layers=list(inps.layers.split(','))
        if 'lookAngle' not in inps.layers: inps.layers.append('lookAngle')
        if 'unwrappedPhase' not in inps.layers: inps.layers.append('unwrappedPhase')

    else:
        inps.layers=list(inps.layers.split(','))


    # extract/merge productBoundingBox layers for each pair and update dict,
    # report common track bbox (default is to take common intersection, but user may specify union), and expected shape for DEM.
    standardproduct_info.products[1], standardproduct_info.bbox_file, prods_TOTbbox, arrshape, proj = merged_productbbox(standardproduct_info.products[1], os.path.join(inps.workdir,'productBoundingBox'), standardproduct_info.bbox_file, inps.croptounion)
    # Load or download mask (if specified).
    if inps.mask is not None:
        inps.mask = prep_mask([[j['amplitude'] for j in standardproduct_info.products[1]], [j["pair_name"] for j in standardproduct_info.products[1]]], inps.mask, standardproduct_info.bbox_file, prods_TOTbbox, proj, amp_thresh=inps.amp_thresh, arrshape=arrshape, workdir=inps.workdir, outputFormat=inps.outputFormat)


    # Download/Load DEM & Lat/Lon arrays, providing bbox, expected DEM shape, and output dir as input.
    if inps.demfile is not None:
        # Pass DEM-filename, loaded DEM array, and lat/lon arrays
        inps.demfile, demfile, Latitude, Longitude = prep_dem(inps.demfile, standardproduct_info.bbox_file, prods_TOTbbox, proj, arrshape=arrshape, workdir=inps.workdir, outputFormat=inps.outputFormat)
    else:
        demfile=None ; Latitude=None ; Longitude=None

    # Extract user expected layers
    export_products(standardproduct_info.products[1], standardproduct_info.bbox_file, prods_TOTbbox, inps.layers, dem=demfile, lat=Latitude, lon=Longitude, mask=inps.mask, outDir=inps.workdir, outputFormat=inps.outputFormat, stitchMethodType='overlap', verbose=inps.verbose)

    # Perform GACOS-based tropospheric corrections (if specified).
    if inps.tropo_products:
        tropo_correction(standardproduct_info.products, inps.tropo_products, standardproduct_info.bbox_file, prods_TOTbbox, outDir=inps.workdir, outputFormat=inps.outputFormat, verbose=inps.verbose)
