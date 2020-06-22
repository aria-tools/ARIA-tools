#!/anaconda3/envs/ARIA-tools/bin/python3

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Robert Zinke, Simran Sangha, David Bekaert
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Import modules
import os
import argparse
import pandas as pd
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from osgeo import ogr
from shapely.geometry import Polygon, MultiPolygon



# Parser
def createParser():
    '''
        Use product metadata to assist in area of interest (AOI) creation.
    '''
    parser = argparse.ArgumentParser( description='Preparing preliminary plot of frame extents. First go to the ASF search page, push all SLCs over defined search area to cart, download CSV under the metadata option, and pass the CSV through to this script with the -f flag.')
    parser.add_argument('-f', '--file', dest='imgfile', type=str, required=True,
            help='Full path to CSV file containing SLC frame metadata.')
    parser.add_argument('-w', '--workdir', dest='workdir', default='./',
            help='Specify directory to deposit all outputs. Default is local directory where script is launched.')
    parser.add_argument('-t', '--tracks', dest='tracks', type=str, default='all',
            help='Include only specified track number in results. Can be multiple, separated by spaces. Default : All')
    parser.add_argument('-l', '--lat_bounds', dest='latBounds', type=str, default='-90 90',
            help='Specify a search for only frames that fall within these lat bounds. Default : -90 90')
    parser.add_argument('-s', '--start_date', dest='startDate', type=str, default=None,
            help='Start date. Default : None')
    parser.add_argument('-e', '--end_date', dest='endDate', type=str, default=None,
            help='End date. Default : None')
    parser.add_argument('-x', '--exclude_dates', dest='excludeDates', type=str, default=None,
            help='List of dates to exclude from kml generation. This can be provided as space-separated string in format YYYYMMDD (e.g., \'20180101 20181213 20190428\'), or as a text file with one date to exclude per line. Default : None')
    parser.add_argument('--plot_raw', dest='plotRaw', action='store_true',
            help='Plot raw frames if included in .csv')
    parser.add_argument('--flag_partial_coverage', dest='flagPartialCoverage', action='store_true',
            help='Flag dates that do not cover the full lat/lon extent. This does not remove dates from the lat centers plot, only highlights the dates in red.')
    parser.add_argument('--remove_incomplete_dates', dest='removeIncomplete', action='store_true',
            help='Automatically detect and remove dates that do not entirely fill the given latitude bounds. Note that if lat bounds are left as default, only dates with gaps will be automatically excluded.')
    parser.add_argument('--approximate_AOI', dest='approxAOI', action='store_true',
            help='Create KML of approximate AOI. NOTE: ~20 km or 1 burst must be removed from either end of the AOI--this must be confirmed by the user.')

    return parser

def cmdLineParse(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



# Metadata class
class SentinelMetadata:
    '''
        Class for parsing, filting, and displaying metadata from ASF Vertex csv using Pandas
        functionality. The self.metadata property is a pandas dataframe in which all metadata are
        carried. Columns are originally those from the ASF csv spreadsheet, and some others will
        be added for further sorting. Filtering is done in-place using the dataframe. Additional
        parameters and flags are added without removing metadata entries.
    '''
    # Load data from csv and pre-format
    def __init__(self,imgfile,track,workdir='./',excludeDates=None,plotRaw=False,
        flag_partial_coverage=False,remove_incomplete_dates=False):
        # Record parameters
        self.track=track
        self.minLat=-60
        self.maxLat=60
        self.workdir=workdir
        self.excludeDates=excludeDates
        self.plotRaw=plotRaw
        self.flagPartialCoverage=flag_partial_coverage
        self.removeIncomplete=remove_incomplete_dates

        # Open ASF Vertex csv file and read meta data as Pandas dataframe
        csvfile_name=os.path.abspath(imgfile)
        self.metadata=pd.read_csv(csvfile_name,index_col=False)

        # Pre-formatting
        # Assign datetimes and common dates
        self.__assignFrameID__()
        self.__assignDatetimes__()
        self.__formatExcludeDates__()

        # Pre-filtering
        # Filter by beam mode -- only IW accepted
        self.__filterByTrack__()
        self.__filterByBeamMode__()
        self.__filterByProcessingLevel__()


    # Formatting
    def __assignFrameID__(self):
        '''
            Develop a unique "frame identification code" comprising orbit+path+frame
        '''
        frameIDs=[]
        for frameNdx,frame in self.metadata.iterrows():
            frameProperties=list(frame[['Orbit','Path Number','Frame Number']])
            frameProperties=[str(frameProperty) for frameProperty in frameProperties] # to string
            frameID=''.join(frameProperties) # concatenate to single string
            frameIDs.append(frameID)
        self.metadata['frameID']=frameIDs

    def __assignDatetimes__(self):
        '''
            Add a column with Python datetime objects to represent precisely the date and time at
            which a frame was captured in a format that can be used computationally.
            Additionally, assign a "common" date to which all acquisitions closely correspond. This
            will help with sorting by date, especially in cases where the satellite crosses the
            midnight boundary.
            This function assumes that acqusitions within 24 hours of each other correspond to the
            same date.
        '''
        # Reformat Acquistion Date string into datetime object
        dateFmt='%Y-%m-%dT%H:%M:%S'
        dates=[datetime.strptime(date.split('.')[0],dateFmt) for date in
            self.metadata['Acquisition Date']]
        self.metadata['Datetime']=dates

        # Determine a "Common Date" for the satellite pass
        oneday=timedelta(days=1)
        commonDatetimes=[]
        commonDates=[]
        for date in self.metadata['Datetime']:
            # Find all dates for which time difference is less than 1 day
            timediffs=date-self.metadata['Datetime'] # calculate difference relative to all dates
            dateGroup=(self.metadata['Datetime'][timediffs<oneday]).sort_values() # lt 1 day
            commonDate=dateGroup.to_list()[0] # use first entry of list as common date
            commonDatetimes.append(commonDate) # datetime class
            commonDates.append(commonDate.strftime('%Y%m%d')) # string
        self.metadata['Common Datetime']=commonDatetimes
        self.metadata['Common Date']=commonDates

    def __formatExcludeDates__(self):
        '''
            Determine whether the --exclude_dates input is specified as an input string or a text
            file. If a text file is given, format and transfer the contents to a list attached to the
            input object.
        '''
        # Format user-specified dates to exclude
        # Check whether exclude dates list is given, otherwise, provide an empty list
        if self.excludeDates:
            # Check if the provided input is a file or else a list string
            if os.path.exists(self.excludeDates):
                # Read and format file contents
                with open(self.excludeDates,'r') as excludeFile:
                    lines=excludeFile.readlines()
                    excludeFile.close()
                self.excludeDates=[line.strip('\n') for line in lines] # remove new line formatting
            else:
                self.excludeDates=self.excludeDates.split() # split by spaces
        else:
            self.excludeDates=[]


    # Filtering
    def __filterByTrack__(self):
        '''
            Remove tracks not specified by user. Modify self.tracks parameter to include only those
            tracks.
        '''
        dropIndicesTrack=self.metadata[self.metadata['Path Number']!=self.track].index
        self.metadata.drop(dropIndicesTrack,inplace=True)

        # Record track direction
        self.trackDir=list(self.metadata['Ascending or Descending?'])[0]

        # Short name comprising track number and track direction, e.g., "A41"
        self.trackCode='{}{}'.format(self.trackDir[0],self.track)

    def __filterByBeamMode__(self):
        '''
            Remove frames if beam mode is not IW.
        '''
        dropIndicesIW=self.metadata[self.metadata['Beam Mode']!='IW'].index
        self.metadata.drop(dropIndicesIW,inplace=True)


    def __filterByProcessingLevel__(self):
        '''
            Remove if processing level is not "SLC" or "RAW". Only include "RAW" frames if "SLC" is
            not available for that location and time.
        '''
        # Ensure that only "RAW" and "SLC" frames are included
        dropIndicesProc=self.metadata[(self.metadata['Processing Level']!='SLC') & \
                                    (self.metadata['Processing Level']!='RAW')].index
        self.metadata.drop(dropIndicesProc,inplace=True)

        # SLC frame IDs as integers to compare with RAW frame IDs
        slcIDs=np.array([int(frame['frameID']) for ndx,frame in self.metadata.iterrows() if \
            frame['Processing Level']=='SLC'])

        # Loop through raw frames to check that a similar SLC does not exist, based on frameID
        dropIndicesRAW=[]
        for frameNdx,frame in self.metadata.iterrows():
            if frame['Processing Level']=='RAW':
                # Check if SLC acquisition is within the nearest 4 frames
                ID_differences=abs(int(frame['frameID'])-slcIDs)
                if np.min(ID_differences)<4: dropIndicesRAW.append(frameNdx)
        self.metadata.drop(dropIndicesRAW,inplace=True)

    def filterByDate(self,startDate=None,endDate=None):
        '''
            Clip the data set based on start and end dates if provided. Provide dates in format
            YYYYMMDD.
        '''
        # Convert start and end dates to datetimes
        if startDate is not None: startDate=datetime.strptime(startDate,'%Y%m%d')
        if endDate is not None: endDate=datetime.strptime(endDate,'%Y%m%d')

        # Determine which dates to drop
        dropIndicesDate=[]
        for frameNdx,frame in self.metadata.iterrows():
            # Check start/end date criteria
            if startDate is not None and frame['Datetime']<startDate:
                dropIndicesDate.append(frameNdx)
            if endDate is not None and frame['Datetime']>endDate:
                dropIndicesDate.append(frameNdx)
        self.metadata.drop(dropIndicesDate,inplace=True)

    def filterByLatitude(self,minLat=None,maxLat=None):
        '''
            Remove scenes if they do not meet the specified latitude requirements.
            Loop through all IW SLCs frames to confirm which fall within latitude bounds.
        '''
        # Update latitude bounds
        if minLat is not None: self.minLat=minLat
        if maxLat is not None: self.maxLat=maxLat

        # Filter by latitude
        dropIndicesLats=[]
        for frameNdx,frame in self.metadata.iterrows():
            # All latitude positions
            frameLats=frame[['Near Start Lat','Far Start Lat','Near End Lat','Far End Lat']]

            # Check that frames are within lat bounds
            if frameLats.max() < self.minLat: dropIndicesLats.append(frameNdx)
            if frameLats.min() > self.maxLat: dropIndicesLats.append(frameNdx)
        self.metadata.drop(dropIndicesLats,inplace=True)


    # Spatial checks
    def __addLatExtremes__(self):
        '''
            Add one column each for the minimum and maximum latitude extents. This is designed to ease
            continuity checks.
        '''
        minLats=[]
        maxLats=[]
        for frameNdx,frame in self.metadata.iterrows():
            lats=frame[['Near Start Lat','Far Start Lat','Near End Lat','Far End Lat']].to_numpy()
            minLats.append(lats.min())
            maxLats.append(lats.max())
        self.metadata['minLat']=minLats
        self.metadata['maxLat']=maxLats

    def checkContinuity(self):
        '''
            Check whether an SLC is missing in the track for any given date. Specifically, we want to
            see whether the spacing between frames exceeds a threshold (e.g., one and a half frames).
            If that is the case, we want to assign different colors to the different subsets of frames.
            This information will be carried through the metadata as new columns: "nGaps" (i.e., the
            number of gaps with frames missing), and "subGroup".

            This function also checks whether the frames cover the full latitude extent provided.
            Dates with incomplete coverage will not be automatically excluded, they will simply
            be flagged.

            If the --remove_incomplete_dates option is selected, update the self.excludeDates
            property with the incomplete dates.
        '''
        # Add columns for min/max Lat extremes to each frame
        self.__addLatExtremes__()

        # Add empty columns for gaps, subgroups, extents covered
        self.metadata['gaps']=np.zeros((self.metadata.shape[0],1))
        self.metadata['subgroup']=np.ones((self.metadata.shape[0],1))
        self.metadata['Extent Covered']=np.zeros((self.metadata.shape[0],1),dtype=bool)

        # Use only SLC, not RAW
        SLCindices=self.metadata[self.metadata['Processing Level']=='SLC'].index

        # Dates
        dates=list(set(self.metadata['Common Date']))

        # Loop by date
        for date in dates:
            # Indices of frames matching non-RAW and date
            dateIndices=self.metadata[self.metadata['Common Date']==date].index
            passIndices=set(SLCindices).intersection(dateIndices)

            # Sort tracks south-north and compare latitude extents
            # "satPass" refers to all the acquisitions from a single satellite pass
            satPass=self.metadata.loc[passIndices,:].sort_values(by='Center Lat')

            # Check that maxLat of southern frame > minLat of northern frame
            nAcquisitions=len(passIndices) # number of acquisitions in pass
            currentSubgroup=1
            for i in range(nAcquisitions-1):
                ndxSouth=satPass.index[i]
                ndxNorth=satPass.index[i+1]

                # Action required if N extent of south frame < S extent of north frame
                # i.e., no overlap
                if satPass.loc[ndxSouth,'maxLat']<satPass.loc[ndxNorth,'minLat']:
                    # Update gaps detected for all acquisitions in pass
                    self.metadata.loc[passIndices,'gaps']+=1
                    # Update subgroup
                    self.metadata.loc[satPass.index[i+1],'subgroup']=currentSubgroup+1

            # Check that max/min latitude extents are covered if lat bounds are
            # user-specified
            if (self.maxLat<60) and (self.minLat>-60):
                latExtremesMet=False
                latMax=satPass['maxLat'].to_numpy().max() # max of all frames at this date
                latMin=satPass['minLat'].to_numpy().min() # min of all frames at this date

                if (latMax>self.maxLat) and (latMin<self.minLat):
                    latExtremesMet=True
            else:
                # Ignore lat criteria if bounds are left at default +- 60 degrees
                latExtremesMet=True

            # If there are no gaps and latitude extremes are covered, consider coverage complete
            if (sum(self.metadata.loc[passIndices,'gaps'])==0) and (latExtremesMet==True):
                self.metadata.loc[passIndices,'Extent Covered']=True

        # If remove_incomplete_dates is specified, update self.excludeDates
        if self.removeIncomplete == True:
            for frameNdx,frame in self.metadata.iterrows():
                if frame['Extent Covered']==False:
                    self.excludeDates.append(frame['Common Date'])
            self.excludeDates=list(set(self.excludeDates))
            self.excludeDates.sort()

            # Save updated date list to text file
            exclDatesName='{}_auto-excluded_dates.txt'.format(self.trackCode)
            exclDatesPath=os.path.join(self.workdir,exclDatesName)
            with open(exclDatesPath,'w') as exclDatesFile:
                for date in self.excludeDates:
                    exclDatesFile.write(date+'\n')
                exclDatesFile.close()


    # Plotting
    def plotFrameCenters(self):
        '''
            Plot the center of each frame with attributes.
        '''
        # Spawn figure
        self.Fig=plt.figure(figsize=(80,11))
        self.ax=self.Fig.add_subplot(111)

        # Plot SLC frames only
        SLCindices=self.metadata[self.metadata['Processing Level']=='SLC'].index

        # Loop through each SLC acqusition
        labels=[]
        for frameNdx,frame in self.metadata.loc[SLCindices,:].iterrows():
            # Color based on gaps
            if (self.flagPartialCoverage==True) and (frame['gaps']>0):
                color='r'
                label='Incomplete track'
            else:
                color='k'
                label='Complete track'

            # Color based on spatial extent, if extent is not automatic
            if (self.flagPartialCoverage==True) and (frame['Extent Covered']==False):
                color='r'
                label='Incomplete track'

            # Gray if excluded
            if frame['Common Date'] in self.excludeDates:
                color=(0.6,0.6,0.6)
                label='Excluded date'

            # Plot frame
            self.ax.scatter(frame['Common Datetime'],frame['Center Lat'],s=100,color=color,
                label=label)

        # Plot RAW frames if specified
        if self.plotRaw==True:
            RAWindices=self.metadata[self.metadata['Processing Level']=='RAW'].index

            # Loop through each RAW acquisition
            for frameNdx,frame in self.metadata.loc[RAWindices,:].iterrows():
                self.ax.scatter(frame['Common Datetime'],frame['Center Lat'],s=100,color='m',
                    label='Raw frame')

        # Format x-axis
        dates=list(set(self.metadata['Common Datetime'])); dates.sort()
        datelabels=list(set(self.metadata['Common Date'])); datelabels.sort()
        self.ax.set_xticks(dates)
        self.ax.set_xticklabels(datelabels,rotation=90)
        self.ax.set_xlim([dates[0]-timedelta(days=6),dates[-1]+timedelta(days=6)])
        self.ax.set_xlabel('Date')

        # Highlight dates with only partial coverage
        if self.flagPartialCoverage==True:
            slcIndices=self.metadata[self.metadata['Processing Level']=='SLC'].index
            partialIndices=self.metadata[self.metadata['Extent Covered']==False].index
            partialIndices=set(slcIndices).intersection(partialIndices)
            partialDates=set([date for date in self.metadata.loc[partialIndices,'Common Date']])

            # Change date label to red if only partial coverage
            [self.ax.get_xticklabels()[n].set_color('r') for n,date in enumerate(datelabels) if
                date in partialDates]

        # Format legend
        handles,labels=self.ax.get_legend_handles_labels()
        uniqueLabels = dict(zip(labels,handles))
        self.ax.legend(uniqueLabels.values(),uniqueLabels.keys(),
            bbox_to_anchor=(1.005,1),loc='upper left',borderaxespad=0.)

        # Other formatting
        title='Track {}'.format(self.trackCode)
        self.ax.set_title(title,weight='bold')
        self.ax.set_ylabel('Latitude (degrees)')
        self.ax.margins(x=0)
        self.Fig.tight_layout()

        # Save to file
        figname='{}_lat_extents.eps'.format(self.trackCode)
        self.Fig.savefig(os.path.join(self.workdir,figname))


    # Save spatial extents
    def save2kml(self):
        '''
            Save frame boundaries to kml file - use only SLC frames and non-exluded dates.
        '''
        # Open KML data set
        kmlname='track{}_frames.kml'.format(self.trackCode)
        kmlpath=os.path.join(self.workdir,kmlname)
        DS=ogr.GetDriverByName('LIBKML').CreateDataSource(kmlpath)

        # Filter by SLC only
        slcIndices=self.metadata[self.metadata['Processing Level']=='SLC'].index

        # Loop by date
        dates=list(set(self.metadata.loc[slcIndices,'Common Date']))
        dates=[date for date in dates if date not in self.excludeDates]
        dates.sort() # sort oldest-most recent

        for date in dates:
            dateIndices=self.metadata[self.metadata['Common Date']==date].index
            dateIndices=set(slcIndices).intersection(dateIndices)

            # Create KML layer
            layer=DS.CreateLayer(date,None,ogr.wkbPolygon)
            layer.CreateField(ogr.FieldDefn('id',ogr.OFTInteger)) # add 1 attribute

            # Add frame polygons
            for frameNdx,frame in self.metadata.loc[dateIndices,:].iterrows():
                feat=ogr.Feature(layer.GetLayerDefn())
                feat.SetField('id',frameNdx)
                bbox=Polygon([(float(frame['Near Start Lon']),float(frame['Near Start Lat'])),
                    (float(frame['Far Start Lon']),float(frame['Far Start Lat'])),
                    (float(frame['Far End Lon']),float(frame['Far End Lat'])),
                    (float(frame['Near End Lon']),float(frame['Near End Lat']))])
                geom=ogr.CreateGeometryFromWkb(bbox.wkb)
                feat.SetGeometry(geom)
                feat.SetStyleString("PEN(c:#000000)")
                layer.CreateFeature(feat)
            # Close layer
            layer=feat=geom=None
        # Close KML
        DS=layer=feat=geom=None


    # Suggested AOI based on intersection of all polygons
    def __polygonFromFrame__(self,frame):
        '''
            Create a Shapely polygon from the coordinates of a frame.
        '''
        P=Polygon([(float(frame['Near Start Lon']),float(frame['Near Start Lat'])),
                    (float(frame['Far Start Lon']),float(frame['Far Start Lat'])),
                    (float(frame['Far End Lon']),float(frame['Far End Lat'])),
                    (float(frame['Near End Lon']),float(frame['Near End Lat']))])
        return P

    def __mergeFramesbyDate__(self,date):
        '''
            For a given date, merge all valid (SLC) frames into a polygon.
        '''
        # Use only SLCs
        slcIndices=self.metadata[self.metadata['Processing Level']=='SLC'].index

        # Collect indices of date
        dateIndices=self.metadata[self.metadata['Common Date']==date].index
        dateIndices=set(slcIndices).intersection(dateIndices)

        # Compute polygons
        datePolygons=[]
        for frameNdx,frame in self.metadata.loc[dateIndices,:].iterrows():
            # Convert frame coords to polygon
            datePolygons.append(self.__polygonFromFrame__(frame))

        # Find union of polygons
        dateUnion=datePolygons[0]
        for datePolygon in datePolygons:
            dateUnion=dateUnion.union(datePolygon)

        return dateUnion

    def __simplifyAOI__(self,AOI):
        '''
            Simplify AOI polygon by removing points within several hundred meters of each other.
        '''
        # Assign search distance
        d=0.05 # 0.01 degrees ~ 1 km
        dsq=d**2

        # Convert polygon vertices to "nodes"
        nodes=np.column_stack(AOI.exterior.xy)
        N=nodes.shape[0] # number of nodes

        # Start with first node
        d2all=(nodes[0,0]-nodes[:,0])**2+(nodes[0,1]-nodes[:,1])**2
        nodeCluster=nodes[d2all<dsq,:]
        clusterCenter=np.mean(nodeCluster,axis=0)
        simpleNodes=np.array(clusterCenter).reshape(1,2)

        # Loop through all nodes to find cluster centers
        for n in range(1,N):
            currentNode=nodes[n,:]

            # First, check if a nearby node already exists
            d2simple=(currentNode[0]-simpleNodes[:,0])**2+(currentNode[1]-simpleNodes[:,1])**2
            if np.min(d2simple)>dsq:
                # Compute distance from given node to all other nodes
                d2all=(nodes[:,0]-currentNode[0])**2+(nodes[:,1]-currentNode[1])**2

                # Identify cluster of all nodes within squared search distance
                nodeCluster=nodes[d2all<dsq]
                clusterCenter=np.mean(nodeCluster,axis=0)

                # Add to simple nodes
                simpleNodes=np.vstack([simpleNodes,clusterCenter])

        # Convert simplified nodes back to Shapely Polygon
        simpleNodes=Polygon(simpleNodes)

        return simpleNodes

    def __trimAOIedges__(self,AOI):
        '''
            Remove ~20 km or 1 burst from the N-S edges of the suggested AOI.
            First, remove data centroid. Then, rotate into a semi-vertical orientation.
            Remove ~20 km from top and bottom points.
        '''
        # Assign trim distance
        trimDist=0.2 # degrees

        # Convert polygon vertices to "nodes"
        nodes=np.column_stack(AOI.exterior.xy)

        # Remove centroid
        nodeMeans=np.mean(nodes,axis=0)
        nodes-=nodeMeans # temporarily remove centroid

        # Solve for polygon axes and rotate to horizontal
        U,s,VT=np.linalg.svd(nodes)

        # Rotate into N-S orientation
        R=np.array([[0,-1],[1, 0]]) # 90-degree rotation matrix
        nodes=R.dot(VT).dot(nodes.T).T

        # Remove 20 km from top and bottom
        sortNdx=np.argsort(nodes[:,1])
        nodes[sortNdx[:2],1]+=trimDist
        nodes[sortNdx[-2:],1]-=trimDist

        # Rotate back to original orientation
        nodes=VT.T.dot(R.T).dot(nodes.T).T

        # Add centroid back
        nodes+=nodeMeans

        # Convert back to Shapely Polygon
        AOI=Polygon(nodes)

        return AOI

    def __saveAOI__(self,AOI):
        '''
            Save AOI to kml file.
        '''
        # Open kml data set
        kmlname='track{}_suggestedAOI.kml'.format(self.trackCode)
        kmlpath=os.path.join(self.workdir,kmlname)
        DS=ogr.GetDriverByName('LIBKML').CreateDataSource(kmlpath)

        # Create KML layer
        layer=DS.CreateLayer('suggestedAOI',None,ogr.wkbPolygon)
        layer.CreateField(ogr.FieldDefn('id',ogr.OFTInteger))

        # Create feature
        feat=ogr.Feature(layer.GetLayerDefn())
        feat.SetField('id',0)
        geom=ogr.CreateGeometryFromWkb(AOI.wkb)
        feat.SetGeometry(geom)
        feat.SetStyleString('PEN(c:#000000)')
        layer.CreateFeature(feat)

        # Close layer
        layer=feat=geom=None

        # Close KML
        DS=layer=feat=geom=None

    def intersectionAOI(self):
        '''
            Create a "suggested" AOI based on the intersection of all frame polygons.
            This function does not guarantee a result.
        '''
        # Use only SLCs
        slcIndices=self.metadata[self.metadata['Processing Level']=='SLC'].index

        # Unique dates
        dates=list(set(self.metadata.loc[slcIndices,'Common Date']))
        dates=[date for date in dates if date not in self.excludeDates]

        # Use first date as initial polygon
        dateUnion=self.__mergeFramesbyDate__(dates[0])
        sceneIntersection=dateUnion

        validOverlap=True # assume all frames overlap until detected otherwise
        for date in dates[1:]:
            # Take union of frames in date
            dateUnion=self.__mergeFramesbyDate__(date)

            # Take intersection of dates
            sceneIntersection=sceneIntersection.intersection(dateUnion)

            # Check whether area of intersection remains a single polygon, or if it is broken into
            # multiple pieces
            if isinstance(sceneIntersection,MultiPolygon):
                print('COULD NOT AUTOMATICALLY DETERMINE AOI DUE TO GAPS.')
                print('Area of frame overlap not continuous, beginning on {}'.format(date))
                validOverlap=False
                break

        # Remove ~20km from either end and save suggested AOI
        if validOverlap==True:
            # Simply polygon by removing points within several hundred meters of each other
            sceneIntersection=self.__simplifyAOI__(sceneIntersection)

            # Remove 20 km from edges
            sceneIntersection=self.__trimAOIedges__(sceneIntersection)

            # Save AOI to KML
            self.__saveAOI__(sceneIntersection)


# Main
if __name__ == "__main__":
    '''
        Main workflow for extracting and visualizing meta data.
    '''

    # Inputs
    # Gather inputs from command line
    inps = cmdLineParse(iargs=None)

    # Split latitude values and convert to type float
    inps.latBounds=[float(val) for val in inps.latBounds.split()]


    # Setup
    # Create work directory
    if not os.path.exists(inps.workdir):
        os.mkdir(inps.workdir)


    # Loop by track -- one class instance per track
    # Detect which tracks to analyze
    if inps.tracks=='all':
        # If "all" tracks are specified, need to load file to detect which ones are available
        csv=pd.read_csv(inps.imgfile,index_col=False)
        tracks=list(set(csv['Path Number']))
        csv=None
    else:
        tracks=[int(track) for track in inps.tracks.split()] # split into individual numbers

    # Loop through tracks
    for track in tracks:
        print('Generating outputs for track: {}'.format(track))

        # Instantiate metadata object and load metadata from csv
        track_metadata=SentinelMetadata(imgfile=inps.imgfile,track=track,workdir=inps.workdir,
            excludeDates=inps.excludeDates,
            plotRaw=inps.plotRaw,
            flag_partial_coverage=inps.flagPartialCoverage,
            remove_incomplete_dates=inps.removeIncomplete)

        # Filtering -- remove frames from metadata
        # Clip based on start and end date
        track_metadata.filterByDate(startDate=inps.startDate,endDate=inps.endDate)

        # Filter by latitude bounds
        track_metadata.filterByLatitude(minLat=inps.latBounds[0],maxLat=inps.latBounds[1])

        # Check spatial criteria -- does not remove scenes, only highlights potentially problematic
        # ones
        # Check for gaps
        track_metadata.checkContinuity()

        # Outputs
        # Plot frame centers
        track_metadata.plotFrameCenters()

        # Save to Google Earth KML
        track_metadata.save2kml()

        # Create suggested AOI based on intersection of all polygons
        if inps.approxAOI==True:
            track_metadata.intersectionAOI()


    print('Products generated.')