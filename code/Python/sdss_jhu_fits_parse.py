#!/usr/bin/env python
""" Parse SDSS sripe82 variable timeseries .fits table files

Call using commandline syntax:
./sdss_jhu_fits_parse.py photo-matches-bin-RA-30to-25N86392-cleaned-MH86338.fit.gz


Printed output is of the form:

-29.9395751953 -0.945328235626 5759 40 1 202 107
-29.9395446777 -0.945354163647 5590 40 1 209 285
-29.9395446777 -0.94540989399 5610 40 1 250 150
-29.9395446777 -0.945373415947 5622 40 1 235 179
-29.9395446777 -0.945381402969 5633 40 1 216 178

where:
RA
Dec
(SDSS object identifier): RUN RERUN CAMCOL FIELD OBJECT

"""


from pyfits import open as pyfits_open
import gzip
import sys

if __name__ == '__main__':

    if len(sys.argv) > 1:
        fpath = sys.argv[1]
    else:
        print #fpath = 'photo-matches-bin-RA-30to-25N86392-cleaned-MH86338.fit.gz'

    hdu = pyfits_open(gzip.open(fpath, 'rb'), mode='readonly')
    #print hdu[1].columns

    for i in xrange(len(hdu[1].data['RA'][0])):
        print hdu[1].data['RA'][0][i], \
              hdu[1].data['DEC'][0][i], \
              hdu[1].data['RUN'][0][i], \
              hdu[1].data['RERUN'][0][i], \
              hdu[1].data['CAMCOL'][0][i], \
              hdu[1].data['FIELD'][0][i], \
              hdu[1].data['OBJECT'][0][i]
        """
              hdu[1].data['OBSMJD'][0][i], \
              hdu[1].data['FLAGS'][0][i], \
              hdu[1].data['PSFMAGS'][0][i], \
              hdu[1].data['PSFERRS'][0][i], \
              hdu[1].data['EXTINCT'][0][i]
        """

        #import pdb; pdb.set_trace()
        #print
