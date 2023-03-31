#!/usr/bin/env bash

export GDAL_HTTP_COOKIEFILE=/tmp//cookies.txt
export GDAL_HTTP_COOKIEJAR=/tmp//cookies.txt
export VSI_CACHE=YES

mkdir test_tmp
cd test_tmp
#product count
ariaDownload.py -t 124 -w products --ifg 20180420_20180315 -v -o count
#product url
ariaDownload.py -t 124 -w products --ifg 20180420_20180315 -v -o url

cd products
wget --no-clobber --load-cookies /tmp/cookies.txt --save-cookies /tmp/cookies.txt --keep-session-cookies https://grfn.asf.alaska.edu/door/download/S1-GUNW-A-R-124-tops-20180420_20180315-043105-00157W_00020N-PP-74e7-v2_0_6.nc
wget --no-clobber --load-cookies /tmp/cookies.txt --save-cookies /tmp/cookies.txt --keep-session-cookies https://grfn.asf.alaska.edu/door/download/S1-GUNW-A-R-124-tops-20180420_20180315-043040-00157W_00018N-PP-3f11-v2_0_6.nc
cd ../

#extraction of 3D metadata layer
ariaExtract.py -f "products/*.nc" -l azimuthAngle -d download
#extraction and mosaicking of coherence data layer
ariaExtract.py -f "products/*.nc" -l coherence
#extraction and mosaicking of unwrappedPhase data layer
ariaExtract.py -f "products/*.nc" -l unwrappedPhase

#time-series prep
ariaTSsetup.py -d download -m download -f "products/*.nc" -w 'ts_stack'
#run loading of aria products in mintpy
prep_aria.py -s ts_stack/stack/ -d ts_stack/DEM/glo_30.dem -i ts_stack/incidenceAngle/20180420_20180315.vrt -a ts_stack/azimuthAngle/20180420_20180315.vrt

#download and test NLCD mask
ariaExtract.py -f "products/*.nc" -m 'NLCD'
