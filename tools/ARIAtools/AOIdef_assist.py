#!/usr/bin/env python3

import os
import argparse
import datetime
from datetime import datetime, timedelta
import numpy as np
import csv
import itertools
import pandas as pd
from collections import OrderedDict
from osgeo import ogr

import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

from shapely.geometry import Polygon, MultiPolygon

def createParser():
    parser = argparse.ArgumentParser( description='Preparing preliminary plot of frame extents')
    parser.add_argument('-f', '--file', dest='imgfile', type=str, required=True,
            help='Full path to CSV file containing SLC frame metadata.')
    parser.add_argument('-w', '--workdir', dest='workdir', default='./', 
            help='Specify directory to deposit all outputs. Default is local directory where script is launched.')
    parser.add_argument('-m', '--min_frames', dest='min_frames', type=int, default=0,
            help='Minimum percentage of frame overlap required for dates.')
    parser.add_argument('-l', '--lat_bounds', dest='lat_bounds', type=str, default='-60 60',
            help='Specify a search for only frames that fall within these lat bounds.')
    parser.add_argument('-c','--num_connections', dest='num_connections', type=str, default = '2',
            help='Number of interferograms between each date and subsequent dates. -- Default : 2')
    parser.add_argument('-d','--num_days', dest='num_days', type=int, default = None,
            help='If date is within specified days n of a year later, then add pair to list. -- Default : 2')

    return parser

def cmdLineParse(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)

def main(iargs=None):

    inps = cmdLineParse(iargs)

    if not os.path.exists(inps.workdir):
        os.mkdir(inps.workdir)

    inps.lat_bounds=[float(val) for val in inps.lat_bounds.split()]
    csvfile_name=os.path.abspath(inps.imgfile)
    with open(csvfile_name, newline='') as csvfile:
        frame_metadata = csv.reader(csvfile)
        frame_metadata = list(frame_metadata)
        frame_metadata=frame_metadata[1:]

    #Text file delineating simple SNWE box for each track
    f = open(os.path.join(inps.workdir,'track_bboxSNWE.txt'),'w')
    #wkt-format file listing polygon objects for each track
    f.write('ID:         S               N         W              E'+'\n')
    f = open(os.path.join(inps.workdir,'QGIS_latbounds.txt'),'w')
    f.write('Id;POLYGON'+'\n')

    OG_all_scene_lat_dates=[]
    OG_str_all_burst_lats=[]
    OG_str_all_burst_lons=[]
    OG_tot_str_all_burst_lats=[]
    OG_all_scene_frame=[]
    OG_all_pass=[]
    OG_metaindices=[]
    OG_frames=[]
    for i in frame_metadata:
        if i[10]=='SLC':
            if min(float(i[15]),float(i[17]),float(i[19]))>inps.lat_bounds[0] and max(float(i[15]),float(i[17]),float(i[19]))<inps.lat_bounds[1]:
                OG_str_all_burst_lats.append(i[13])
                OG_str_all_burst_lons.append(i[18])
                OG_all_scene_lat_dates.append(i[0][17:25])
                OG_frames.append(i[0])
                OG_metaindices.append(i[5]+i[6]+i[7])
                if i[15]>i[21] and i[15]<i[17]: #descending
                    OG_tot_str_all_burst_lats.append(i[15:23])
                elif i[15]>i[21] and i[15]>i[17]: #ascending glitch
                    OG_tot_str_all_burst_lats.append([i[15],i[16],i[19],i[20],i[17],i[18],i[21],i[22]])
                else: #ascending glitch
                    OG_tot_str_all_burst_lats.append([i[19],i[20],i[21],i[22],i[15],i[16],i[17],i[18]])
                OG_all_scene_frame.append(i[6])
                OG_all_pass.append(i[24][0])

    ###check for raw frames
    for i in frame_metadata:
        if len(set(OG_metaindices).intersection(set([i[5]+i[6]+str(int(i[7])-4),i[5]+i[6]+str(int(i[7])-3),i[5]+i[6]+str(int(i[7])-2),i[5]+i[6]+str(int(i[7])-1),str(int(i[5]+i[6]+i[7])),i[5]+i[6]+str(int(i[7])+1),i[5]+i[6]+str(int(i[7])+2),i[5]+i[6]+str(int(i[7])+3),i[5]+i[6]+str(int(i[7])+4)])))==0 and i[10]=='RAW':
            if min(float(i[15]),float(i[17]),float(i[19]))>inps.lat_bounds[0] and max(float(i[15]),float(i[17]),float(i[19]))<inps.lat_bounds[1]:
                OG_str_all_burst_lats.append(i[13])
                OG_str_all_burst_lons.append(i[18])
                OG_all_scene_lat_dates.append(i[0][17:25])
                OG_frames.append(i[0])
                OG_metaindices.append("RAW")
                if i[15]>i[21] and i[15]<i[17]: #descending
                    OG_tot_str_all_burst_lats.append(i[15:23])
                elif i[15]>i[21] and i[15]>i[17]: #ascending glitch
                    OG_tot_str_all_burst_lats.append([i[15],i[16],i[19],i[20],i[17],i[18],i[21],i[22]])
                else: #ascending glitch
                    OG_tot_str_all_burst_lats.append([i[19],i[20],i[21],i[22],i[15],i[16],i[17],i[18]])
                OG_all_scene_frame.append(i[6])
                OG_all_pass.append(i[24][0])

    OG_all_scene_frame, OG_all_scene_lat_dates, OG_str_all_burst_lats, OG_tot_str_all_burst_lats, OG_all_pass,OG_str_all_burst_lons,OG_metaindices,OG_frames = (list(t) for t in zip(*sorted(zip(OG_all_scene_frame, OG_all_scene_lat_dates, OG_str_all_burst_lats, OG_tot_str_all_burst_lats, OG_all_pass, OG_str_all_burst_lons,OG_metaindices,OG_frames))))

    ###split arrays
    OG_all_scene_lat_dates=[OG_all_scene_lat_dates[OG_all_scene_frame.index(str(k)):len(list(v))+OG_all_scene_frame.index(str(k))] for k,v in itertools.groupby(OG_all_scene_frame, key=lambda x: x[:3])]
    OG_str_all_burst_lons=[OG_str_all_burst_lons[OG_all_scene_frame.index(str(k)):len(list(v))+OG_all_scene_frame.index(str(k))] for k,v in itertools.groupby(OG_all_scene_frame, key=lambda x: x[:3])]
    OG_str_all_burst_lats=[OG_str_all_burst_lats[OG_all_scene_frame.index(str(k)):len(list(v))+OG_all_scene_frame.index(str(k))] for k,v in itertools.groupby(OG_all_scene_frame, key=lambda x: x[:3])]
    OG_tot_str_all_burst_lats=[OG_tot_str_all_burst_lats[OG_all_scene_frame.index(str(k)):len(list(v))+OG_all_scene_frame.index(str(k))] for k,v in itertools.groupby(OG_all_scene_frame, key=lambda x: x[:3])]
    OG_all_pass=[OG_all_pass[OG_all_scene_frame.index(str(k)):len(list(v))+OG_all_scene_frame.index(str(k))] for k,v in itertools.groupby(OG_all_scene_frame, key=lambda x: x[:3])]
    OG_metaindices=[OG_metaindices[OG_all_scene_frame.index(str(k)):len(list(v))+OG_all_scene_frame.index(str(k))] for k,v in itertools.groupby(OG_all_scene_frame, key=lambda x: x[:3])]
    OG_frames=[OG_frames[OG_all_scene_frame.index(str(k)):len(list(v))+OG_all_scene_frame.index(str(k))] for k,v in itertools.groupby(OG_all_scene_frame, key=lambda x: x[:3])]
    OG_all_scene_frame= [list(v) for k,v in itertools.groupby(OG_all_scene_frame, key=lambda x: x[:3])]

    ###Cycle through tracks
    for a in range(0,len(OG_all_scene_frame)):
        OG_str_all_burst_lons[a]=len(OG_str_all_burst_lons[a])*[min(OG_str_all_burst_lons[a])]

    OG_str_all_burst_lons, OG_all_scene_frame, OG_all_scene_lat_dates, OG_str_all_burst_lats, OG_tot_str_all_burst_lats, OG_all_pass, OG_metaindices, OG_frames = (list(t) for t in zip(*sorted(zip(OG_str_all_burst_lons, OG_all_scene_frame, OG_all_scene_lat_dates, OG_str_all_burst_lats, OG_tot_str_all_burst_lats, OG_all_pass, OG_metaindices, OG_frames))))
    alt_colors=[]
    temp_val=0
    temp_val_append=[]
    for i in range(0,len(OG_str_all_burst_lons)):
        temp_val_append.append('T'+str(OG_all_scene_frame[i][0]))
        if len(temp_val_append)>1:
            if temp_val_append[i]!=temp_val_append[i-1]:
               if temp_val==0:
                   temp_val=1
               else:
                   temp_val=0
        if temp_val==0:
            alt_colors.append(len(OG_all_scene_frame[i])*['red'])
        else:
            alt_colors.append(len(OG_all_scene_frame[i])*['blue']) 

    OG_all_scene_frame, OG_all_scene_lat_dates, OG_str_all_burst_lats, OG_tot_str_all_burst_lats, OG_all_pass, alt_colors, OG_str_all_burst_lons, OG_metaindices, OG_frames = (list(t) for t in zip(*sorted(zip(OG_all_scene_frame, OG_all_scene_lat_dates, OG_str_all_burst_lats, OG_tot_str_all_burst_lats, OG_all_pass, alt_colors, OG_str_all_burst_lons, OG_metaindices, OG_frames))))

    
    ###Cycle through tracks
    for a in range(0,len(OG_all_scene_frame)):
        all_scene_lat_dates=OG_all_scene_lat_dates[a]
        all_scene_lon_dates=OG_str_all_burst_lons[a]
        str_all_burst_lats=OG_str_all_burst_lats[a]
        tot_str_all_burst_lats=OG_tot_str_all_burst_lats[a]
        all_scene_frame=OG_all_scene_frame[a]
        all_pass=OG_all_pass[a]
        all_colors=alt_colors[a]
        all_metaindices=OG_metaindices[a]
        all_frames=OG_frames[a]

        f = plt.figure(figsize=(80,11))
        ax=f.add_subplot(111)
        all_scene_lat_dates, str_all_burst_lats, tot_str_all_burst_lats, all_scene_frame, all_metaindices = (list(t) for t in zip(*sorted(zip(all_scene_lat_dates, str_all_burst_lats, tot_str_all_burst_lats, all_scene_frame, all_metaindices))))

        for e, i in reversed(list(enumerate(all_scene_lat_dates))): #check for midnight crossing
            date_dec=datetime.strptime(i, "%Y%m%d")
            date_dec+=timedelta(days=-1)
            date_dec=date_dec.strftime("%Y%m%d")
            if date_dec in all_scene_lat_dates:
                all_scene_lat_dates[e]=date_dec

        array_all_scene_lat_dates=range(len(all_scene_lat_dates))
        unique_all_scene_lat_dates=np.unique(all_scene_lat_dates).tolist()
        dateList=np.unique(all_scene_lat_dates).tolist()
        unique_array_all_scene_lat_dates=range(len(unique_all_scene_lat_dates))
        pd_dates = pd.to_datetime(dateList)
        pd_dates = [datetime(i.year, i.month, i.day) for i in pd_dates]
        ax.set_xticks(pd_dates)

        mostcommonvalue=max(zip((all_scene_lat_dates.count(item) for item in set(all_scene_lat_dates)), set(all_scene_lat_dates)))
        mostcommonvalue_array=[]
        ###get master date overlap
        stack_master_overlap=[[float(tot_str_all_burst_lats[i][0][0]),float(tot_str_all_burst_lats[i][0][-2])] for i, x in enumerate(all_scene_lat_dates) if x == mostcommonvalue[1]]
        stack_master_overlap = list(itertools.chain.from_iterable(stack_master_overlap))
        stack_master_overlap=max(stack_master_overlap)- min(stack_master_overlap)

        #open KML
        ds = ogr.GetDriverByName('LIBKML').CreateDataSource(os.path.join(inps.workdir,'track%s%s_frames.kml'%(all_pass[0],all_scene_frame[0])))

        for i in array_all_scene_lat_dates:
            date_ind=[l for l, e in enumerate(all_scene_lat_dates) if e==all_scene_lat_dates[i]]

            #only initiate KML layer if iterating through first instance of date in list
            if i == date_ind[0]:
                # create KML layer
                layer = ds.CreateLayer(str(all_scene_lat_dates[i]), None, ogr.wkbPolygon)
                layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger)) #Add 1 attribute

            lat_ind=[str_all_burst_lats[i] for i in date_ind]
            lat_ind.sort()
            lat_ind=[str(abs(float(x)-float(lat_ind[i-1]))) for i, x in enumerate(lat_ind)][1:]
            scene_overlap=[float(str_all_burst_lats[i]) for i in date_ind]
            scene_overlap=max(scene_overlap)-min(scene_overlap)
            if all_scene_lat_dates[i]==mostcommonvalue[1] and all_metaindices[i]!='RAW':
                mostcommonvalue_array.append(tot_str_all_burst_lats[i])
            if ((scene_overlap/stack_master_overlap)*100)>inps.min_frames and all(float(i)<1.75 for i in lat_ind)==True and all_metaindices[i]!='RAW':
                ax.plot(pd_dates[unique_all_scene_lat_dates.index(all_scene_lat_dates[i])],float(str_all_burst_lats[i]),'ko',markersize=10, label="valid frame")
                ###only write to KML layer if valid frame
                # Create a new feature (attribute and geometry)
                feat = ogr.Feature(layer.GetLayerDefn())
                feat.SetField('id', i)
                # Make a geometry, from input Shapely object
                bbox = Polygon(np.column_stack((np.array([tot_str_all_burst_lats[i][1],tot_str_all_burst_lats[i][3],tot_str_all_burst_lats[i][7],tot_str_all_burst_lats[i][5],tot_str_all_burst_lats[i][1]]),np.array([tot_str_all_burst_lats[i][0],tot_str_all_burst_lats[i][2],tot_str_all_burst_lats[i][6],tot_str_all_burst_lats[i][4],tot_str_all_burst_lats[i][0]])))) #Pass lons/lats to create polygon
                geom = ogr.CreateGeometryFromWkb(bbox.wkb)
                feat.SetGeometry(geom)
                layer.CreateFeature(feat)

            if ((scene_overlap/stack_master_overlap)*100)<=inps.min_frames and all_metaindices[i]!='RAW':
                ax.plot(pd_dates[unique_all_scene_lat_dates.index(all_scene_lat_dates[i])],float(str_all_burst_lats[i]),'ro',markersize=20, label="invalid frame")
                ax.get_xticklabels()[unique_all_scene_lat_dates.index(all_scene_lat_dates[i])].set_color("red")
                if all_scene_lat_dates[i] in dateList:
                    dateList.remove(all_scene_lat_dates[i])
            if all(float(i)<1.75 for i in lat_ind)==False and all_metaindices[i]!='RAW':
                ax.plot(pd_dates[unique_all_scene_lat_dates.index(all_scene_lat_dates[i])],float(str_all_burst_lats[i]),'ro',markersize=20)
                ax.axvline(x=pd_dates[unique_all_scene_lat_dates.index(all_scene_lat_dates[i])], color='r', linestyle='--')
                ax.get_xticklabels()[unique_all_scene_lat_dates.index(all_scene_lat_dates[i])].set_color("red")
                if all_scene_lat_dates[i] in dateList:
                    dateList.remove(all_scene_lat_dates[i])
            if all_metaindices[i]=='RAW':
                ax.plot(pd_dates[unique_all_scene_lat_dates.index(all_scene_lat_dates[i])],float(str_all_burst_lats[i]),'gx',markersize=25, label="only raw available")
                #track unprocessed raw frames
                fi = open(os.path.join(inps.workdir,'%s%s_rawframes_to_download.txt'%(all_pass[0],all_scene_frame[0])),'a')
                fi.write(all_frames[i])
                fi.write('\n')

            #only close KML layer if iterating through late instance of date in list
            if i == date_ind[-1]:
                # close KML layer
                layer = feat = geom = None
        mostcommonvalue_array.sort()
        mostcommonvalue_array.insert(0,[mostcommonvalue[1]])
        # close KML file
        ds = layer = feat = geom = None

        ###GMT script
        f = open(os.path.join(inps.workdir,'gmtbounds.txt'),'a')
        f.write('#TRACK_'+all_pass[0]+all_scene_frame[0]+'\n')
        f.write('pstext $LLreg $LLproj -G%s -N -O -K << EOF >> $out'%(all_colors[0])+'\n')
        coord_latlon=[str(((float(mostcommonvalue_array[-1][3])+float(mostcommonvalue_array[-1][1]))/2)+0.1),str(((float(mostcommonvalue_array[-1][2])+float(mostcommonvalue_array[-1][0]))/2)+0.1)]
        coord_latlon=[str(round(float(i),1)) for i in coord_latlon]
        if all_pass[0][0]=='D':
            f.write('%s %s 12 -10 0 BC %s'%(coord_latlon[0], coord_latlon[1],all_scene_frame[0])+'\n')
        else:
            f.write('%s %s 12 10 0 BC %s'%(coord_latlon[0], coord_latlon[1],all_scene_frame[0])+'\n')
        f.write('EOF'+'\n'+'\n')

        f.write('psxy  $LLreg $LLproj -L -W2,%s -: -m -O -K << EOF >> $out'%(all_colors[0])+'\n')
        f.write(mostcommonvalue_array[1][6]+" "+mostcommonvalue_array[1][7]+'\n')
        f.write(mostcommonvalue_array[1][4]+" "+mostcommonvalue_array[1][5]+'\n')
        for h in range(1,len(mostcommonvalue_array)):
            f.write(mostcommonvalue_array[h][4]+" "+mostcommonvalue_array[h][5]+'\n')
            f.write(mostcommonvalue_array[h][0]+" "+mostcommonvalue_array[h][1]+'\n')
        f.write(mostcommonvalue_array[-1][0]+" "+mostcommonvalue_array[-1][1]+'\n')
        f.write(mostcommonvalue_array[-1][2]+" "+mostcommonvalue_array[-1][3]+'\n')
        for h in reversed(range(1,len(mostcommonvalue_array))):
            f.write(mostcommonvalue_array[h][2]+" "+mostcommonvalue_array[h][3]+'\n')
            f.write(mostcommonvalue_array[h][6]+" "+mostcommonvalue_array[h][7]+'\n')
        f.write('EOF'+'\n'+'\n'+'\n')

        ###Qgis lat bounds
        f = open(os.path.join(inps.workdir,'QGIS_latbounds.txt'),'a')
        f.write('%s;POLYGON(('%(all_scene_frame[0])+mostcommonvalue_array[1][7]+" "+mostcommonvalue_array[1][6]+" , ")
        for h in range(1,len(mostcommonvalue_array)):
            f.write(mostcommonvalue_array[h][5]+" "+mostcommonvalue_array[h][4]+" , ")
            f.write(mostcommonvalue_array[h][1]+" "+mostcommonvalue_array[h][0]+" , ")
        f.write(mostcommonvalue_array[-1][1]+" "+mostcommonvalue_array[-1][0]+" , ")
        f.write(mostcommonvalue_array[-1][3]+" "+mostcommonvalue_array[-1][2]+" , ")
        for h in reversed(range(1,len(mostcommonvalue_array))):
            f.write(mostcommonvalue_array[h][3]+" "+mostcommonvalue_array[h][2]+" , ")
            if h==1:
                f.write(mostcommonvalue_array[h][7]+" "+mostcommonvalue_array[h][6])
            else:
                f.write(mostcommonvalue_array[h][7]+" "+mostcommonvalue_array[h][6]+" , ")
        f.write(" ))"+'\n')

        ###Bounding box
        f = open(os.path.join(inps.workdir,'track_bboxSNWE.txt'),'a')
        if all_pass[0][0]=='D':
            f.write('%s '%(all_scene_frame[0])+mostcommonvalue_array[1][-2]+" "+mostcommonvalue_array[-1][0]+" "+mostcommonvalue_array[1][-1]+" "+mostcommonvalue_array[-1][1]+"\n")
        else:
            f.write('%s '%(all_scene_frame[0])+mostcommonvalue_array[1][-2]+" "+mostcommonvalue_array[-1][0]+" "+mostcommonvalue_array[-1][1]+" "+mostcommonvalue_array[1][-1]+"\n")

        ###List of pairs
        num_connections = inps.num_connections
        pairs = []
        if num_connections == 'all':
            num_connections = len(dateList) - 1
        else:
            num_connections = int(num_connections)

        pd_dates = pd.to_datetime(dateList)
        num_connections = num_connections
        for i in reversed(range(len(dateList))):
            for j in range(i- num_connections,i):
                if j>=0: ###to prevent from looping to last dates of list again with earlier dates
                    date_ind=[l for l, e in enumerate(all_scene_lat_dates) if e==all_scene_lat_dates[j]]
                    scene_overlap_M=[float(str_all_burst_lats[e]) for e in date_ind]
                    scene_overlap_M=max(scene_overlap_M)-min(scene_overlap_M)
                    date_ind=[l for l, e in enumerate(all_scene_lat_dates) if e==all_scene_lat_dates[i]]
                    scene_overlap_S=[float(str_all_burst_lats[e]) for e in date_ind]
                    scene_overlap_S=max(scene_overlap_S)-min(scene_overlap_S)
                    if scene_overlap_M>0:
                        if j<len(dateList) and ((scene_overlap_S/scene_overlap_M)*100)>inps.min_frames:
                            pairs.append(dateList[j]+'_'+dateList[i])
            #If user specifies addition of annual pairs. See "num_days" argument.
            if inps.num_days:
                try:
                    date_1yr_early = datetime(pd_dates[i].year-1, pd_dates[i].month, pd_dates[i].day)
                except:
                    date_1yr_early = datetime(pd_dates[i].year-1, pd_dates[i].month, pd_dates[i].day-1)
                pot_1yrmaster=pd_dates[pd_dates.get_loc(date_1yr_early,method='nearest')]
                ###If date is within specified days n of a year later, then add pair to list 
                if abs((date_1yr_early-pot_1yrmaster).days)<inps.num_days and str(dateList[pd_dates.get_loc(date_1yr_early,method='nearest')]+'_'+dateList[i]) not in pairs:
                    pairs.append(dateList[pd_dates.get_loc(date_1yr_early,method='nearest')]+'_'+dateList[i])
                try:
                    date_1yr_later = datetime(pd_dates[i].year+1, pd_dates[i].month, pd_dates[i].day)
                except:
                    date_1yr_later = datetime(pd_dates[i].year+1, pd_dates[i].month, pd_dates[i].day-1)
                pot_1yrslave=pd_dates[pd_dates.get_loc(date_1yr_later,method='nearest')]
                ###If date is within specified days n of a year later, then add pair to list 
                if abs((date_1yr_later-pot_1yrslave).days)<inps.num_days and str(dateList[i]+'_'+dateList[pd_dates.get_loc(date_1yr_later,method='nearest')]) not in pairs:
                    pairs.append(dateList[i]+'_'+dateList[pd_dates.get_loc(date_1yr_later,method='nearest')])

        pairs = list(set(pairs)) #remove duplicates
        pairs.sort()
        ###Make list of possible pairs for each track
        f = open(os.path.join(inps.workdir,'%s%s_pairs.txt'%(all_pass[0],all_scene_frame[0])),'w')
        for i in pairs:
            f.write(i)
            f.write('\n')

        # add legend
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0.)

        # make lat extents plot for each track
        ax.set_xlabel('Date',weight='bold')
        ax.set_ylabel('Lat',weight='bold')
        ax.set_title('%s%s lat extents'%(all_pass[0],all_scene_frame[0]),weight='bold')
        ax.margins(x=0)
        myFmt = mdates.DateFormatter('%Y%m%d')
        plt.gca().xaxis.set_major_formatter(myFmt)
        plt.xticks(rotation=90)
        plt.tight_layout()

        # saving the figure
        plt.savefig(os.path.join(inps.workdir,'%s%s_lat_extents.eps'%(all_pass[0],all_scene_frame[0])))
        plt.close()

if __name__ == "__main__":

  # Main engine  
  main()
