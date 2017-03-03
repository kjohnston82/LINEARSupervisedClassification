#!/usr/bin/env python
"""
Parse relevant columns from SDSS-III photometry fits tables

Shell syntax:
   ./parse_sdss.py /Data/dstarr/src/nersc_paralleldb/data/calibObj-008157-1-star.fits.gz > out


* TODO: objc_type is only useful to include if we will store star and galaxy in the same tables.

* NOTE: I ignore OBJC_FLAG? attribs since their binary masks are harder to query on and the information might be too fine detailed for our use.  With run,rerun,field,camcol,id information we should be able to find this information elsewhere if needed for some sort of forensics.


[(k,sky_hdu[1].header[k]) for k in sky_hdu[1].header.keys()]
...
 ('TTYPE1', 'RUN'),
 ('TTYPE2', 'RERUN'),
 ('TTYPE3', 'CAMCOL'),
 ('TTYPE4', 'FIELD'),
 ('TTYPE21', 'RA'),
 ('TTYPE22', 'DEC'),
...

sky_hdu[1].columns
name = 'RUN'
format = 'I', name = 'RERUN'
format = '3A', name = 'CAMCOL'
format = 'B', name = 'FIELD'
format = 'I', name = 'ID'


sky_hdu[1].data.shape
Out[53]: (73061,)


len(sky_hdu[1].data['RA'])
Out[52]: 73061

"""
import sys
import gzip
from pyfits import open as pyfits_open
from numpy import log, log10, sqrt, array, clip, where

def bin_repr(x, digits=0):
    """ Display binary representation for any oct()'able type.
>>> bin(5) 
'101' 
>>> bin(0x0FFF, 16) 
'0000111111111111' 
>>> bin(5, 8) 
'00000101' 
    """
    oct2bin = ['000','001','010','011','100','101','110','111'] 
    binstring = [oct2bin[int(n)] for n in oct(x)] 
    return ''.join(binstring).lstrip('0').zfill(digits) 



def dered_fluxes(extinction, flux, ivar=None):
    """  Copied from https://code.google.com/p/sdsspy/source/browse/sdsspy/util.py
    supposedly more accurate at the faint end.
    """
    exponent = 0.4*extinction
    flux_correct = (10.0**exponent)*flux

    if ivar is None:
        return flux_correct
    else:
        exponent = -0.8*extinction
        ivar_correct = (10.0**exponent)*ivar

        return flux_correct, ivar_correct


def nmgy2mag(nmgy, ivar=None):
    """ Copied from https://code.google.com/p/sdsspy/source/browse/sdsspy/util.py

    Name:
        nmgy2mag
    Purpose:
        Convert SDSS nanomaggies to a log10 magnitude.  Also convert
        the inverse variance to mag err if sent.  The basic formulat
        is 
            mag = 22.5-2.5*log_{10}(nanomaggies)
    Calling Sequence:
        mag = nmgy2mag(nmgy)
        mag,err = nmgy2mag(nmgy, ivar=ivar)
    Inputs:
        nmgy: SDSS nanomaggies.  The return value will have the same
            shape as this array.
    Keywords:
        ivar: The inverse variance.  Must have the same shape as nmgy.
            If ivar is sent, then a tuple (mag,err) is returned.

    Outputs:
        The magnitudes.  If ivar= is sent, then a tuple (mag,err)
        is returned.

    Notes:
        The nano-maggie values are clipped to be between 
            [0.001,1.e11]
        which corresponds to a mag range of 30 to -5
    """
    nmgy = array(nmgy, ndmin=1, copy=False)

    nmgy_clip = clip(nmgy,0.001,1.e11)

    mag = nmgy_clip.copy()
    mag[:] = 22.5-2.5*log10(nmgy_clip)

    if ivar is not None:

        ivar = array(ivar, ndmin=1, copy=False)
        if ivar.shape != nmgy.shape:
            raise ValueError("ivar must be same shape as input nmgy array")

        err = mag.copy()
        err[:] = 9999.0

        w=where( ivar > 0 )

        if w[0].size > 0:
            err[w] = sqrt(1.0/ivar[w])

            a = 2.5/log(10)
            err[w] *= a/nmgy_clip[w]

        return mag, err
    else:
        return mag


def test_compare_fluxmag_calculation():
    """ This was used to determine which magnitude calculation from flux is correct
    and matches the skyserver.sdss3.org/dr9/   magnitudes for the same source.

    I found that the model magnitude matches (without extinction nor sky flux correction).

    ### NOTE: I do not correct the skyflux since I assume the psfflux has already had sky subtracted, and that skyflux is just for noise statistics purposes.
    #         f = star_hdu[1].data['PSFFLUX'][i][j]# - star_hdu[1].data['SKYFLUX'][i][j]

    ### NOTE: I do not correct using extinction magnitude, since SDSS3 documentation makes it clear the extinction correction is not applied to data products and is left to the user.
    #         f_corr, f_ivar = dered_fluxes(star_hdu[1].data['EXTINCTION'][i][j], star_hdu[1].data['PSFFLUX'][i][j], ivar=star_hdu[1].data['PSFFLUX_IVAR'][i][j])
    #         m, m_err = nmgy2mag(f_corr, ivar=f_ivar)
    """
    import os
    pars = {'photorun_fpath':os.path.abspath('../data/photoRun-007761.fits'),
            'photofield_fpath':os.path.abspath('../data/photoField-007761-1.fits'),
            'calib_star_fpath':os.path.abspath('../data/calibObj-008157-1-star.fits.gz'),
            'calib_gal_fpath':os.path.abspath('../data/calibObj-008157-1-gal.fits.gz'),
            'calib_sky_fpath':os.path.abspath('../data/calibObj-008157-1-sky.fits.gz'),
            }

    #sky_hdu = pyfits_open(gzip.open(pars['calib_sky_fpath'], 'rb'), mode='readonly')
    star_hdu = pyfits_open(gzip.open(pars['calib_star_fpath'], 'rb'), mode='readonly')
    #gal_hdu = pyfits_open(gzip.open(pars['calib_gal_fpath'], 'rb'), mode='readonly')

    #print 'sky_hdu[1].data.shape', sky_hdu[1].data.shape
    print 'star_hdu[1].data.shape', star_hdu[1].data.shape
    #print 'gal_hdu[1].data.shape', gal_hdu[1].data.shape


    single_attribs = ['RUN',
                      'RERUN',
                      'CAMCOL',
                      'FIELD',
                      'ID',
                      'OBJC_TYPE', # all stars are 6, galaxy=3
                      'RA',
                      'DEC',
                      'TMASS_J', # mag
                      'TMASS_J_IVAR', #inverse variance in 1./nmgy^2
                      'TMASS_H', # mag
                      'TMASS_H_IVAR', #inverse variance in 1./nmgy^2
                      'TMASS_K', # mag
                      'TMASS_K_IVAR', #inverse variance in 1./nmgy^2
                      'TMASS_JDATE',
                      'NDETECT', # N detections of obj in SDSS imaging
                      'PM_MATCH', # N obj in USNO-B cat which matched within 1" radius (using pm), if <0: nearest USNO obj matched more than 1 SDSS object
                      'PM_NFIT', # N USNO-B objs +1,
                      'PM_DIST22', # arcsec dist to nn with g<22
                      ]
    #    name = 'EXTINCTION'; format = '5E'
    #    name = 'PSFFLUX'; format = '5E'
    #    name = 'PSFFLUX_IVAR'; format = '5E'

    # TODO: calculate magnitude using nanomaggy equasion and compare with a SIMBAD mag?w

    ### Print the available data for a single object:
    for obj in star_hdu[1].columns:
        k = obj.name
        print k, star_hdu[1].data[k][0]


    #for i in xrange(int(star_hdu[1].data.shape[0])):
    for i in xrange(10):
        row_list = []
        mags = []
        mags_old = []
        mags_model = []
        mags_psf_ext = []
        m_errs = []

        for j in xrange(5):
            f = star_hdu[1].data['PSFFLUX'][i][j]# - star_hdu[1].data['SKYFLUX'][i][j] #I
            m_old = 22.5 - 2.5*log10(f)#/1.e9)
            mags_old.append(m_old)
            m, m_err = nmgy2mag(f, ivar=star_hdu[1].data['PSFFLUX_IVAR'][i][j])
            mags.append(m[0])
            m_errs.append(m_err[0])

            m, m_err = nmgy2mag(star_hdu[1].data['MODELFLUX'][i][j], ivar=star_hdu[1].data['MODELFLUX_IVAR'][i][j])
            mags_model.append(m[0])

            f_corr, f_ivar = dered_fluxes(star_hdu[1].data['EXTINCTION'][i][j], star_hdu[1].data['PSFFLUX'][i][j], ivar=star_hdu[1].data['PSFFLUX_IVAR'][i][j])
            m, m_err = nmgy2mag(f_corr, ivar=f_ivar)
            mags_psf_ext.append(m[0])

        #print '>>', m_errs
            
        #run={run} rerun={rerun} camcol={camcol} field={field} id={id}
        print "i={ind} ra={ra} dec={dec} m_u={m_u:.3f} m_g={m_g:.3f} m_r={m_r:.3f} m_i={m_i:.3f} m_z={m_z:.3f}".format( \
                             ind=i,
                             run=star_hdu[1].data['RUN'][i],
                             rerun=star_hdu[1].data['RERUN'][i],
                             camcol=star_hdu[1].data['CAMCOL'][i],
                             field=star_hdu[1].data['FIELD'][i],
                             id=star_hdu[1].data['ID'][i],
                             ra=star_hdu[1].data['RA'][i],
                             dec=star_hdu[1].data['DEC'][i],
                             m_u=mags[0],
                             m_g=mags[1],
                             m_r=mags[2],
                             m_i=mags[3],
                             m_z=mags[4],
            )
        print '  m_old:', mags_old
        print '  m_ext:', mags_psf_ext
        print '  model:', mags_model
        print '   mags:', mags
        print '  merrs:', m_errs

        for k in single_attribs:
            print "        ", k, star_hdu[1].data[k][i]



def retrieve_star_attrib_csv(star_hdu=None, use_stdout=True):
    """
    """

    single_attribs = ['TMASS_J', # mag
                      'TMASS_H', # mag
                      'TMASS_K', # mag
                      'TMASS_JDATE',
                      'NDETECT', # N detections of obj in SDSS imaging
                      'PM_MATCH', # N obj in USNO-B cat which matched within 1" radius (using pm), if <0: nearest USNO obj matched more than 1 SDSS object
                      'PM_NFIT', # N USNO-B objs +1,
                      'PM_DIST22', # arcsec dist to nn with g<22
                      'OBJC_TYPE', # all stars are 6, galaxy=3
                      ]
    #                  'OBJC_TYPE', # all stars are 6, galaxy=3
    #                  'TMASS_J_IVAR', #inverse variance in 1./nmgy^2
    #                  'TMASS_H_IVAR', #inverse variance in 1./nmgy^2
    #                  'TMASS_K_IVAR', #inverse variance in 1./nmgy^2
    #                  'RUN',
    #                  'RERUN',
    #                  'CAMCOL',
    #                  'FIELD',
    #                  'ID',
    #                  'RA',
    #                  'DEC',
                      

    ### Print the available data for a single object:
    #for obj in star_hdu[1].columns:
    #    k = obj.name
    #    print k, star_hdu[1].data[k][0]
    #import pdb; pdb.set_trace()
    #print


    header_column_order = ['RA','DEC','m_u','m_g','m_r','m_i','m_z','merr_u','merr_g','merr_r','merr_i','merr_z','ext_u','ext_u','ext_r','ext_i','ext_z','RUN','RERUN','CAMCOL','FIELD','ID']
    header_column_order.extend(single_attribs)

    if use_stdout:
        print '#' + ','.join(header_column_order)
    else:
        out_lines = [','.join(header_column_order)]


    #for i in xrange(10):
    for i in xrange(int(star_hdu[1].data.shape[0])):
        mags = []
        m_errs = []

        #for k in single_attribs:
        #    print "        ", k, star_hdu[1].data[k][i]

        for j in xrange(5):
            m, m_err = nmgy2mag(star_hdu[1].data['MODELFLUX'][i][j], ivar=star_hdu[1].data['MODELFLUX_IVAR'][i][j])
            mags.append(m[0])
            m_errs.append(m_err[0])


        #print "i={ind} RA={RA} DEC={DEC} m_u={m_u:.3f} m_g={m_g:.3f} m_r={m_r:.3f} m_i={m_i:.3f} m_z={m_z:.3f} merr_u={merr_u:.3f} merr_g={merr_g:.3f} merr_r={merr_r:.3f} merr_i={merr_i:.3f} merr_z={merr_z:.3f} ext_u={ext_u} ext_g={ext_u} ext_r={ext_r} ext_i={ext_i} ext_z={ext_z} RUN={RUN} RERUN={RERUN} CAMCOL={CAMCOL} FIELD={FIELD} ID={ID}".format( \
        attr_str = "{RA},{DEC},{m_u},{m_g},{m_r},{m_i},{m_z},{merr_u},{merr_g},{merr_r},{merr_i},{merr_z},{ext_u},{ext_u},{ext_r},{ext_i},{ext_z},{RUN},{RERUN},{CAMCOL},{FIELD},{ID}".format( \
                             ind=i,
                             RUN=star_hdu[1].data['RUN'][i],
                             RERUN=star_hdu[1].data['RERUN'][i],
                             CAMCOL=star_hdu[1].data['CAMCOL'][i],
                             FIELD=star_hdu[1].data['FIELD'][i],
                             ID=star_hdu[1].data['ID'][i],
                             RA=star_hdu[1].data['RA'][i],
                             DEC=star_hdu[1].data['DEC'][i],
                             m_u=mags[0],
                             m_g=mags[1],
                             m_r=mags[2],
                             m_i=mags[3],
                             m_z=mags[4],
                             merr_u=m_errs[0],
                             merr_g=m_errs[1],
                             merr_r=m_errs[2],
                             merr_i=m_errs[3],
                             merr_z=m_errs[4],
                             ext_u=star_hdu[1].data['EXTINCTION'][i][0],
                             ext_g=star_hdu[1].data['EXTINCTION'][i][1],
                             ext_r=star_hdu[1].data['EXTINCTION'][i][2],
                             ext_i=star_hdu[1].data['EXTINCTION'][i][3],
                             ext_z=star_hdu[1].data['EXTINCTION'][i][4])
        single_vals_list = [str(star_hdu[1].data[k][i]) for k in single_attribs]
        final_str = attr_str + ',' + ','.join(single_vals_list)
        if use_stdout:
            print final_str
        else:
            out_lines.append(final_str)

    if use_stdout:
        return
    else:
        return out_lines


def retrieve_galaxy_attrib_csv(star_hdu=None, use_stdout=True):
    """
    """

    single_attribs = ['NDETECT', # N detections of obj in SDSS imaging
                      'OBJC_TYPE', # all stars are 6, galaxy=3
                      ]
    #                  'TMASS_J_IVAR', #inverse variance in 1./nmgy^2
    #                  'TMASS_H_IVAR', #inverse variance in 1./nmgy^2
    #                  'TMASS_K_IVAR', #inverse variance in 1./nmgy^2
    #                  'RUN',
    #                  'RERUN',
    #                  'CAMCOL',
    #                  'FIELD',
    #                  'ID',
    #                  'RA',
    #                  'DEC',
                      

    ### Print the available data for a single object:
    #for obj in star_hdu[1].columns:
    #    k = obj.name
    #    print k, star_hdu[1].data[k][0]


    header_column_order = ['RA','DEC','m_u','m_g','m_r','m_i','m_z','merr_u','merr_g','merr_r','merr_i','merr_z','ext_u','ext_u','ext_r','ext_i','ext_z','RUN','RERUN','CAMCOL','FIELD','ID']
    header_column_order.extend(single_attribs)

    if use_stdout:
        print '#' + ','.join(header_column_order)
    else:
        out_lines = [','.join(header_column_order)]


    #for i in xrange(10):
    for i in xrange(int(star_hdu[1].data.shape[0])):
        mags = []
        m_errs = []

        #for k in single_attribs:
        #    print "        ", k, star_hdu[1].data[k][i]

        for j in xrange(5):
            m, m_err = nmgy2mag(star_hdu[1].data['MODELFLUX'][i][j], ivar=star_hdu[1].data['MODELFLUX_IVAR'][i][j])
            mags.append(m[0])
            m_errs.append(m_err[0])


        #print "i={ind} RA={RA} DEC={DEC} m_u={m_u:.3f} m_g={m_g:.3f} m_r={m_r:.3f} m_i={m_i:.3f} m_z={m_z:.3f} merr_u={merr_u:.3f} merr_g={merr_g:.3f} merr_r={merr_r:.3f} merr_i={merr_i:.3f} merr_z={merr_z:.3f} ext_u={ext_u} ext_g={ext_u} ext_r={ext_r} ext_i={ext_i} ext_z={ext_z} RUN={RUN} RERUN={RERUN} CAMCOL={CAMCOL} FIELD={FIELD} ID={ID}".format( \
        attr_str = "{RA},{DEC},{m_u},{m_g},{m_r},{m_i},{m_z},{merr_u},{merr_g},{merr_r},{merr_i},{merr_z},{ext_u},{ext_u},{ext_r},{ext_i},{ext_z},{RUN},{RERUN},{CAMCOL},{FIELD},{ID}".format( \
                             ind=i,
                             RUN=star_hdu[1].data['RUN'][i],
                             RERUN=star_hdu[1].data['RERUN'][i],
                             CAMCOL=star_hdu[1].data['CAMCOL'][i],
                             FIELD=star_hdu[1].data['FIELD'][i],
                             ID=star_hdu[1].data['ID'][i],
                             RA=star_hdu[1].data['RA'][i],
                             DEC=star_hdu[1].data['DEC'][i],
                             m_u=mags[0],
                             m_g=mags[1],
                             m_r=mags[2],
                             m_i=mags[3],
                             m_z=mags[4],
                             merr_u=m_errs[0],
                             merr_g=m_errs[1],
                             merr_r=m_errs[2],
                             merr_i=m_errs[3],
                             merr_z=m_errs[4],
                             ext_u=star_hdu[1].data['EXTINCTION'][i][0],
                             ext_g=star_hdu[1].data['EXTINCTION'][i][1],
                             ext_r=star_hdu[1].data['EXTINCTION'][i][2],
                             ext_i=star_hdu[1].data['EXTINCTION'][i][3],
                             ext_z=star_hdu[1].data['EXTINCTION'][i][4])
        single_vals_list = [str(star_hdu[1].data[k][i]) for k in single_attribs]
        final_str = attr_str + ',' + ','.join(single_vals_list)
        if use_stdout:
            print final_str
        else:
            out_lines.append(final_str)

    if use_stdout:
        return
    else:
        return out_lines


def retrieve_spectra_attrib_csv(star_hdu=None, use_stdout=True):
    """ Adapted from retrieve_star_attrib_csv()

Run using:
python parse_sdss.py spectra /home/dstarr/src/nersc_paralleldb/data/specObj-dr9.fits

    ### specta urls:
    http://data.sdss3.org/datamodel/files/SPECTRO_REDUX/specObj.html
    http://www.sdss3.org/dr9/spectro/catalogs.php

    ### Info on how to match photometry and spectra:
      - http://www.sdss3.org/dr9/algorithms/match.php

      - There are two versions of matching photometry for each
        spectrum: a positional match within 2 arcsec, and a flux-based
        match that tells you which object contributed the most light
        to the spectrum. These two matches differ about 0.5% of the
        time.  bestObjID is the "default" value (best for point sources),
        although fluxObjID may be superior for galaxy redshift use.

# select by specobjid:
http://skyserver.sdss3.org/dr9/en/tools/explore/obj.asp?sid=299489676975171584
# select by objid:
http://skyserver.sdss3.org/dr9/en/tools/explore/obj.asp?id=1237648720142401611


CLASS GALAXY
SUBCLASS 
   Z 0.0212755
   Z_ERR 9.12182e-06
   SPECOBJID     299489676975171584  # 
#BESTOBJID 1237648720142401611     # best positional match of equivalent photometric object
#MJD 51602                         # MJD of observation
#OBJID [756   1   1 206 129]       # photometric objid: run rerun camcol field ID
   PLUG_RA 146.71421                 # J200
   PLUG_DEC -1.0413043

    """
    ### Print the available data for a single object:
    #for obj in star_hdu[1].columns:
    #    k = obj.name
    #    print k, star_hdu[1].data[k][0]

    header_column_order = ['SPECOBJID',
                           'PLUG_RA', 
                           'PLUG_DEC',
                           'Z',
                           'Z_ERR',
                           'CLASS',
                           'SUBCLASS',
                           'MJD',
                           'OBJ_RUN',
                           'OBJ_RERUN',
                           'OBJ_CAMCOL',
                           'OBJ_FIELD',
                           'OBJ_ID']

    if use_stdout:
        print '#' + ','.join(header_column_order)
    else:
        out_lines = [','.join(header_column_order)]

    for i in xrange(int(star_hdu[1].data.shape[0])):
        objid_list = []

        for j in xrange(5):
            obj_rrcfi = star_hdu[1].data['OBJID'][i][j]

        attr_str = '{SPECOBJID},{PLUG_RA},{PLUG_DEC},{Z},{Z_ERR},"{CLASS}","{SUBCLASS}",{MJD},{OBJ_RUN},{OBJ_RERUN},{OBJ_CAMCOL},{OBJ_FIELD},{OBJ_ID}'.format( \
                             SPECOBJID=star_hdu[1].data['SPECOBJID'][i],
                             PLUG_RA=star_hdu[1].data['PLUG_RA'][i],
                             PLUG_DEC=star_hdu[1].data['PLUG_DEC'][i],
                             Z=star_hdu[1].data['Z'][i],
                             Z_ERR=star_hdu[1].data['Z_ERR'][i],
                             CLASS=star_hdu[1].data['CLASS'][i],
                             SUBCLASS=star_hdu[1].data['SUBCLASS'][i],
                             MJD=star_hdu[1].data['MJD'][i],
                             OBJ_RUN=star_hdu[1].data['OBJID'][i][0],
                             OBJ_RERUN=star_hdu[1].data['OBJID'][i][1],
                             OBJ_CAMCOL=star_hdu[1].data['OBJID'][i][2],
                             OBJ_FIELD=star_hdu[1].data['OBJID'][i][3],
                             OBJ_ID=star_hdu[1].data['OBJID'][i][4],
                             )
        print attr_str

    if use_stdout:
        return
    else:
        return out_lines


if __name__ == '__main__':


    if len(sys.argv) > 2:
        fits_type = sys.argv[1]
        fits_fpath = sys.argv[2]
        if fits_type == 'star':
            hdu = pyfits_open(gzip.open(fits_fpath, 'rb'), mode='readonly')
            retrieve_star_attrib_csv(star_hdu=hdu, use_stdout=True)
        elif fits_type == 'spectra':
            hdu = pyfits_open(open(fits_fpath, 'rb'), mode='readonly')
            # /home/dstarr/src/nersc_paralleldb/data/specObj-dr9.fits
            retrieve_spectra_attrib_csv(star_hdu=hdu, use_stdout=True)
        elif fits_type == 'galaxy':
            hdu = pyfits_open(gzip.open(fits_fpath, 'rb'), mode='readonly')
            retrieve_galaxy_attrib_csv(star_hdu=hdu, use_stdout=True)


    if 0:
        import os
        pars = {'photorun_fpath':os.path.abspath('../data/photoRun-007761.fits'),
                'photofield_fpath':os.path.abspath('../data/photoField-007761-1.fits'),
                'calib_star_fpath':os.path.abspath('../data/calibObj-008157-1-star.fits.gz'),
                'calib_gal_fpath':os.path.abspath('../data/calibObj-008157-1-gal.fits.gz'),
                'calib_sky_fpath':os.path.abspath('../data/calibObj-008157-1-sky.fits.gz'),
                }

        #sky_hdu = pyfits_open(gzip.open(pars['calib_sky_fpath'], 'rb'), mode='readonly')
        star_hdu = pyfits_open(gzip.open(pars['calib_star_fpath'], 'rb'), mode='readonly')
        #gal_hdu = pyfits_open(gzip.open(pars['calib_gal_fpath'], 'rb'), mode='readonly')

        #print 'sky_hdu[1].data.shape', sky_hdu[1].data.shape
        #print 'star_hdu[1].data.shape', star_hdu[1].data.shape
        #print 'gal_hdu[1].data.shape', gal_hdu[1].data.shape

        retrieve_star_attrib_csv(star_hdu=star_hdu, use_stdout=True)

        import pdb; pdb.set_trace()
        print
