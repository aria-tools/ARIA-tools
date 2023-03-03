#!/usr/bin/env bash

export GDAL_HTTP_COOKIEFILE=/tmp//cookies.txt
export GDAL_HTTP_COOKIEJAR=/tmp//cookies.txt
export VSI_CACHE=YES

mkdir test_tmp
cd test_tmp
#product count
ariaDownload.py -t 124 -w products --ifg 20180408_20180502 -v -o count
#product url
ariaDownload.py -t 124 -w products --ifg 20180408_20180502 -v -o url

cd products
wget --no-clobber --load-cookies /tmp/cookies.txt --save-cookies /tmp/cookies.txt --keep-session-cookies https://grfn.asf.alaska.edu/door/download/S1-GUNW-A-R-124-tops-20180502_20180408-043106-21658N_19755N-PP-0dd0-v2_0_1.nc
wget --no-clobber --load-cookies /tmp/cookies.txt --save-cookies /tmp/cookies.txt --keep-session-cookies https://grfn.asf.alaska.edu/door/download/S1-GUNW-A-R-124-tops-20180502_20180408-043040-20161N_18088N-PP-6704-v2_0_1.nc
cd ../

#extraction of 3D metadata layer
ariaExtract.py -f "products/*.nc" -l azimuthAngle -d download
#extraction and mosaicking of coherence data layer
ariaExtract.py -f "products/*.nc" -l coherence
#extraction and mosaicking of unwrappedPhase data layer
ariaExtract.py -f "products/*.nc" -l unwrappedPhase

#time-series prep
ariaTSsetup.py -d download -m download -f "products/*.nc"
#run loading of aria products in mintpy
prep_aria.py -s stack/ -d DEM/SRTM_3arcsec.dem -i incidenceAngle/20180502_20180408.vrt -a azimuthAngle/20180502_20180408.vrt

#download and test NLCD mask
ariaExtract.py -f "products/*.nc" -m 'NLCD'
