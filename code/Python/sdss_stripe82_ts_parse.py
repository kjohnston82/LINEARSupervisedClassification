#!/usr/bin/env python
""" Parse SDSS sripe82 full catalog dataset timeseries tsobj*.fits files

NOTE that parse_sdss.py preceeded this script and was intended for
parse calibObj*fits.gz catalog files, not time series data.


CALL on NERSC carver USING:

    module load python/2.7.1 numpy/1.6.1
    python sdss_stripe82_ts_pars.py /global/homes/d/dstarr/paralleldb/SDSS/s82/tsobjSplit_cr_dir

    OR, in PDB mode:
    /usr/common/usg/python/2.7.1-20110310/lib/python2.7/pdb.py sdss_stripe82_ts_parse.py /global/homes/d/dstarr/paralleldb/SDSS/s82/tsobjSplit_cr_dir


PRINTED OUTPUT (to be piped/inserted into a DB):

    COLUMNS:
run camCol rerun field id ra dec mjd cmodel_u cmodel_g cmodel_r cmodel_i cmodel_z cmodel_err_u cmodel_err_g cmodel_err_r cmodel_err_i cmodel_err_z objc_flags objc_flags2

    FOR EXMAPLE:
7110 3 40 97 1 309.848118508 -0.217301662434 54404.0985066 19.8932 18.3616 17.8861 17.7238 17.6106 0.0388644 0.00700434 0.00675674 0.00746204 0.0184178 268566528 36864
7110 3 40 97 4 309.849504998 -0.418446143537 54404.0985172 16.4727 14.8805 14.3652 14.2132 14.1358 0.0061161 0.00307969 0.00350394 0.00381531 0.00408365 -1879048176 16



SOFTWARE PYTHON DEPENDENCIES:

    ### add to ~/.bashrc.ext:
    export PYTHONPATH=/global/homes/d/dstarr/local/lib/python2.7/site-packages:${ggPYTHONPATH}

    module load python/2.7.1 numpy/1.6.1
    source ~/.bashrc.ext
    mkdir -p /global/homes/d/dstarr/local/lib/python2.7/site-packages
    easy_install --prefix=/global/homes/d/dstarr/local pyfits 


    ###This script depends on the sdsspy python package:
    ###   https://code.google.com/p/sdsspy/source/checkout
    hg clone https://code.google.com/p/sdsspy/
    cd sdsspy/
    python setup.py install --prefix=/global/homes/d/dstarr/local
    ### Test it is installed:  (some not loaded warnings are OK)
    python
    import sdsspy



Software COMPILED DEPENDENCY:
    
    wget http://tdc-www.harvard.edu/software/wcstools/wcstools-3.8.7.tar.gz
    tar -xzvf wcstools-3.8.7.tar.gz
    make all
    cp bin/* ~/local/bin/



More information about the timeseries .fits file format:
    http://www.sdss.org/dr7/dm/flatFiles/tsObj.html

TFORM27 = '5E      '          
TTYPE27 = 'psfCounts'          / PSF flux
TUNIT27 = 'asinh mag     '          

TFORM76 = '5E      '          
TTYPE76 = 'counts_deV'         / De Vaucouleurs magnitude fit
TUNIT76 = 'asinh mag '

TFORM84 = '5E      '          
TTYPE84 = 'counts_exp'         / Exponential fit
TUNIT84 = 'asinh mag '

TFORM86 = '5E      '          
TTYPE86 = 'counts_model'       / Better of DeV/Exp magnitude fit
TUNIT86 = 'asinh mag '   

TFORM101= '75E     '          
TTYPE101= 'profMean'           / Mean pixel flux in annulus
TUNIT101= 'maggies/arcsec^2'


FILTERCHARS = ['u','g','r','i','z']
FILTERNUM = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4,
             0:0, 1:1, 2:2, 3:3, 4:4}


"""
import os, sys

old_stderr = sys.stderr
sys.stderr = open('/dev/null','w')

from pyfits import open as pyfits_open
import datetime
#import gzip
from time import strptime
import subprocess
import copy
import numpy
import sdsspy
import flag_mask
import glob

sys.stderr.close()
sys.stderr = old_stderr


def debug_print_fits_rows(hdu):
    """ Just print out fits rows for diagnositcs

   run VALUE=7202     run VALUE=7202    run VALUE=7202    run VALUE=7202    run VALUE=7202    run VALUE=7202
camCol VALUE=6     camCol VALUE=6    camCol VALUE=6    camCol VALUE=6    camCol VALUE=6    camCol VALUE=6   
 rerun VALUE=40     rerun VALUE=40    rerun VALUE=40    rerun VALUE=40    rerun VALUE=40    rerun VALUE=40  
 field VALUE=706    field VALUE=706   field VALUE=706   field VALUE=706   field VALUE=706   field VALUE=706 
parent VALUE=-1    parent VALUE=1    parent VALUE=2    parent VALUE=2    parent VALUE=-1   parent VALUE=5   
    id VALUE=1         id VALUE=2        id VALUE=3        id VALUE=4        id VALUE=5        id VALUE=6   
nchild VALUE=1     nchild VALUE=2    nchild VALUE=0    nchild VALUE=0    nchild VALUE=1    nchild VALUE=0   

    """
    for i_row in xrange(len(hdu[1].data)):
        for i, elem in enumerate(hdu[1].data[i_row]):
            
            try:
                
                extra = ''
                if len(elem) > 5:
                    extra = '!!!!!!!!!'
                print "%20s VALUE=%0.20s LENGTH=%d %s" % (hdu[1].columns[i].name, str(elem[0]), len(elem), extra)


            except:
                print "%20s VALUE=%0.20s" % (hdu[1].columns[i].name, str(elem))
                #print hdu[1].columns[i].name, elem

        import pdb; pdb.set_trace()
        print


def extract_tsObj_data(fpath):
    """ Given an hdu from an opened .fits file, parse the timeseries
    info for seperate sources.
    """
    out_dict = {'run':     [],
                'camCol':  [],
                'rerun':   [],
                'field':   [],
                'id':      [],
                'parent':  [],
                'objc_flags': [], 
                'objc_flags2':  [],
                'ra':      [],
                'dec':     [],
                'counts_model':  [],
                'counts_modelErr':  [],
                }
    """
                'fracPSF': [],
                'psfCounts':  [],
                'psfCountsErr':  [],
                'counts_deV':  [],
                'counts_deVErr':  [],
                'counts_exp':  [],
                'counts_expErr':  [],
    """

    hdu = pyfits_open(fpath, mode='readonly')
    for i_row in xrange(len(hdu[1].data)):
        for k in out_dict.keys():
            out_dict[k].append(hdu[1].data[i_row][k])

    out_dict['counts_model'] = numpy.array(out_dict['counts_model'])
    out_dict['counts_modelErr'] = numpy.array(out_dict['counts_modelErr'])
    #out_dict['psfCounts'] = numpy.array(out_dict['psfCounts'])
    #out_dict['psfCountsErr'] = numpy.array(out_dict['psfCountsErr'])
    #out_dict['counts_deV'] = numpy.array(out_dict['counts_deV'])
    #out_dict['counts_deVErr'] = numpy.array(out_dict['counts_deVErr'])
    #out_dict['counts_exp'] = numpy.array(out_dict['counts_exp'])
    #out_dict['counts_expErr'] = numpy.array(out_dict['counts_expErr'])
    #out_dict['fracPSF'] = numpy.array(out_dict['fracPSF'])
    return out_dict


def extract_fpC_data(fname):
    """ given an un-gzipped fpC...fit file,
    extract data needed for doing SDSS drift time calculation

    This requires xy2sky WCSTOOLS program to be installed.
    
    """
    fimg = pyfits_open(fname, mode='readonly')
    hdr=fimg[0].header
    fimg.close()
    if not hdr.has_key('DATE-OBS') or not hdr.has_key("TAIHMS") or not hdr.has_key("EXPTIME"):
        # damn. not what we were expecting
        raise# continue
    exptime = datetime.timedelta(seconds=float(hdr.get('EXPTIME')))
    tmp = hdr.get('DATE-OBS') + "T" + hdr.get('TAIHMS')

    ## get the approx TAI datetime start

    try:
        im_start_time = datetime.datetime(*strptime(tmp.split(".")[0], "%d/%m/%yT%H:%M:%S")[0:6])
    except:
        im_start_time = datetime.datetime(*strptime(tmp.split(".")[0], "%Y-%m-%dT%H:%M:%S")[0:6])
        im_start_time += datetime.timedelta(microseconds=float('.' + tmp.split(".")[1])*1e6)
        
    #if self.verbose:
    #    print "+ %s image start time is %s " % (filt,str(im_start_time))

    ## now we have the start time for the y=1 pixel
    # get the fiducial RA, DEC
    xfid = hdr.get('naxis1')/2.0
    tmp = "xy2sky -d %s %f %f" % (fname,xfid, 1.0)
    (a,b,c) = os.popen3(tmp)
    a.close()
    c.close()
    ll = b.readlines()
    b.close()
    tmp = ll[0].split()
    (ra1,dec1) = (float(tmp[0]),float(tmp[1]))

    return {'ra1':ra1,
            'dec1':dec1,
            'im_start_time':im_start_time,
            'exptime':exptime,

            }


def parse_times_from_fpC_fits(tsObj_dict, fpc_dirpath=''):
    """
    Adapted from TCP/Software/AstroPhot/sdss_astrophot.py:_apply_astrom()

    INPUT:  tsObj_dict with info:
        field, run
        ra, dec of sources
    OUTPUT:
        list of time_observed for each source

    SOFTWARE DEPENDENCY:
    
    wget http://tdc-www.harvard.edu/software/wcstools/wcstools-3.8.7.tar.gz
    tar -xzvf wcstools-3.8.7.tar.gz
    make all
    cp bin/* ~/local/bin/



### time drift offsets (band, mjd, seconds):
r 54433.4006959 2.92245880266
i 54433.401526 2.98484635386
z 54433.4031859 3.10002294612
u 54433.402357 3.13841571737.
g 54433.4040172 3.27758720706
NOTE: Although 'r' has the bandpass that is at maximum QE and tends to be the brightest wavelengths.
I choose the 'z' fpC...fit file since it is the central filter and thus central time, which is all we care about in the header

    
    """

    run_str = "%06i" % int(tsObj_dict['run'][0])
    field_str = "%04i" % int(tsObj_dict['field'][0])
    filt = 'z' #I choose the 'z' fpC...fit file since it is the central filter and thus central time, which is all we care about in the header
    mjd0 = datetime.datetime(1858,11,17)

    ### TODO: we may want to form the directory path here as well:
    fname = "fpC-" + run_str + "-" + filt + str(tsObj_dict['camCol'][0]) + "-" + field_str + ".fit"    
    download_fpath = fpc_dirpath + "/" + fname

    ### NOTE: downloading 3000 bytes works, 2000 bytes doesnt (incomplete header)  Sticking with 8192
    get_command = "wget -O - http://das.sdss.org/imaging/{run}/{rerun}/corr/{camcol}/{fname}.gz | head -c 3000 | gunzip > {download_fpath}".format(
        fname=fname,
        run=tsObj_dict['run'][0],
        camcol=tsObj_dict['camCol'][0],
        rerun=tsObj_dict['rerun'][0],
        download_fpath=download_fpath,
        )

    ### Unfortunately since we want to use xy2sky, we need to use gunzip in the shell:
    if not os.path.exists(fname):
        p = subprocess.Popen(get_command, shell=True, stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
        sts = os.waitpid(p.pid, 0)

    ### TODO: re-enable try/except to be more robuse to incomplete / errant files:
    try:
        fpC_dict = extract_fpC_data(download_fpath)
    except:
        return [] # calling function will skip this file in this case.
    dist_in_seconds_arr = 3600.0*((numpy.array(tsObj_dict['ra']) - fpC_dict['ra1'])*numpy.cos(fpC_dict['dec1']*numpy.pi/180.)/15.)

    t_mjd_list = []
    for dist_in_seconds in dist_in_seconds_arr:
        t = fpC_dict['im_start_time'] + datetime.timedelta(seconds=dist_in_seconds) + fpC_dict['exptime'] - mjd0
        #tmp = copy.deepcopy(numarray.fromlist([(float(x.days) + float(x.seconds)/(24.0*3600.0) + \
        #    float(x.microseconds)/(24.0*3600*1e6)) for x in t]))
        t_mjd = float(t.days) + float(t.seconds)/(24.0*3600.0) + float(t.microseconds)/(24.0*3600*1e6)
        #print t_mjd, dist_in_seconds
        t_mjd_list.append(t_mjd)

    return t_mjd_list


def calculate_magnitudes(tsObj_dict):
    """ Calculate magnitudes using the fiven ashinh mags (luptitudes).

    nmgy2mag(nmgy, ivar=None)
    lups2nmgy(lups, err=None, band=None)
    make_cmodelmag(objs, doerr=True, dered=False, lups=False)

    objs['devflux']
    objs['expflux']
    objs['fracdev']  or  objs['fracpsf']     fracPSF
    objs['devflux_ivar']
    objs['expflux_ivar']
    """
    out = {}

    nmgy_ivar_arr = sdsspy.lups2nmgy(numpy.array(tsObj_dict['psfCounts']),
                                     err=numpy.array(tsObj_dict['psfCountsErr']))
    mag = sdsspy.nmgy2mag(nmgy_ivar_arr)
    out['mag_psfcounts'] = mag[0]
    out['mag_psfcounts_ivar'] = mag[1]
    

    nmgy_ivar_arr = sdsspy.lups2nmgy(tsObj_dict['counts_model'],
                                     err=tsObj_dict['counts_modelErr'])
    mag = sdsspy.nmgy2mag(nmgy_ivar_arr)
    out['mag_countsmodel'] = mag[0]
    out['mag_countsmodel_ivar'] = mag[1]


    """
    ### Here I'm trying to calculate the deV + Exponential model
    ###    magnitude, but the sdsspy module functions don't seem to
    ###    play nicely with numpy recarray, which they require.
    arr = numpy.array([tsObj_dict['counts_deV'],
                       tsObj_dict['counts_exp'],
                       tsObj_dict['counts_deVErr'],
                       tsObj_dict['counts_expErr'],
                       tsObj_dict['fracPSF']],
                      dtype=[('devflux', tsObj_dict['counts_deV'].dtype),
                             ('expflux', tsObj_dict['counts_exp'].dtype),
                             ('devflux_ivar', tsObj_dict['counts_deVErr'].dtype),
                             ('expflux_ivar', tsObj_dict['counts_expErr'].dtype),
                             ('fracpsf', tsObj_dict['fracPSF'].dtype)])
    rec = arr.view(numpy.recarray)
    model_mags = sdsspy.make_cmodelmag(rec, doerr=True, dered=False, lups=False) # this returns log mags
    """
    return out


def print_sources(tsObj_dict={}, mags_dict={}, mjd_list={}):
    """ combine parents & children for sources

    Then print the results to stdout so it can be inserted into a file or table.

    """
    skip_objc_flags = set(['EDGE', 'BADSKY', 'SATUR']) # CR is ok and shouldnt affect photometry
    skip_objc_flags2 = set(['BAD_MOVING_FIT_CHILD', 'BAD_MOVING_FIT', 'TOO_FEW_GOOD_DETECTIONS', 'LOCAL_EDGE'])
    secondary_objc_flags = set(['CHILD','CR', 'INTERP'])
    secondary_objc_flags2 = set([])

    id_i_dict = {}
    for i, id in enumerate(tsObj_dict['id']):
        id_i_dict[id] = i

    parent_ind_dict = {}
    for i in xrange(len(tsObj_dict['ra'])):
        if tsObj_dict['parent'][i] == -1:
            parent_ind_dict[tsObj_dict['id'][i]] = [(mjd_list[i], i)]

        elif parent_ind_dict.has_key(tsObj_dict['parent'][i]):
            parent_ind_dict[tsObj_dict['parent'][i]].append((mjd_list[i], i))
        else:
            #print "i=%d tsObj_dict['id'][i]=%d tsObj_dict['parent'][i]=%d" % (i, tsObj_dict['id'][i], tsObj_dict['parent'][i])
            ### In this case the child's parent is also a child of another parent...
            id_parent = tsObj_dict['parent'][i]
            i_parent = id_i_dict[id_parent]
            id_parent = tsObj_dict['parent'][i_parent]
            while True:
                if parent_ind_dict.has_key(id_parent):
                    parent_ind_dict[id_parent].append((mjd_list[i], i))
                    break # done, found the parent
                else:
                    i_parent = id_i_dict[id_parent]
                    id_parent = tsObj_dict['parent'][i_parent]

    for parent_id, ind_tup_list in parent_ind_dict.iteritems():
        chosen_indicies = []
        for mjd, i in ind_tup_list:

            objc_flags = flag_mask.get_obj1_flags(tsObj_dict['objc_flags'][i])
            objc_flags2 = flag_mask.get_obj2_flags(tsObj_dict['objc_flags2'][i])

            if len(set(objc_flags) & skip_objc_flags):
                #print '#',
                continue
            elif len(set(objc_flags2) & skip_objc_flags2):
                #print '#',
                continue
            else:
                #print ' ',
                chosen_indicies.append(i)

        n_inds = len(chosen_indicies)
        if n_inds == 0:
            continue # skip this parent since the flags of all of its children failed
        elif n_inds == 1:
            chosen_index = chosen_indicies[0]
        elif n_inds > 1:
            ### Here we tally which observation (index) has the most median magnitude (trying for all 5 bands)
            median_index_tally = {}
            for i in chosen_indicies:
                median_index_tally[i] = 0
                
            for i_band in xrange(5):
                mag_tups = []
                for i in chosen_indicies:
                    mag_tups.append((tsObj_dict['counts_model'][i][i_band], i))
                mag_tups.sort()

                #print '   ', mag_tups, mag_tups[n_inds/2][1], mag_tups[(n_inds/2) - 1][1]
                if (n_inds % 2) == 0:
                    median_index_tally[mag_tups[n_inds/2][1]] += 1
                    median_index_tally[mag_tups[(n_inds/2) - 1][1]] += 1
                else:
                    median_index_tally[mag_tups[n_inds/2][1]] += 1
            tups_list = [(count,index) for index,count in median_index_tally.iteritems()]
            tups_list.sort(reverse=True)

            if tups_list[0][0] != tups_list[1][0]:
                chosen_index = tups_list[0][1] # this is the best median mag observation
            else:
                # There are two observations with compernable number of median occurances.  Now we look at minor SDSS flags to determine
                sort_tups = []
                for i in [tups_list[0][1], tups_list[1][1]]:
                    objc_flags = flag_mask.get_obj1_flags(tsObj_dict['objc_flags'][i])
                    objc_flags2 = flag_mask.get_obj2_flags(tsObj_dict['objc_flags2'][i])
                    sort_tups.append((len((set(objc_flags) & secondary_objc_flags) | (set(objc_flags2) & secondary_objc_flags2)), i))
                sort_tups.sort()  # now [0] is the observation (i) with the fewest number of matching secondary_flags
                ### if there are 2 observations with the same number of middle occurances, this is a pseudo random selection of an index.  Otherwise it is the median observation:
                chosen_index = sort_tups[0][1]
            #print chosen_index
        i = chosen_index
        table = { \
            'run':     tsObj_dict['run'][i],
            'camCol':  tsObj_dict['camCol'][i],
            'rerun':   tsObj_dict['rerun'][i],
            'field':   tsObj_dict['field'][i],
            'parent':  tsObj_dict['parent'][i],
            'id':      tsObj_dict['id'][i],
            'ra':      tsObj_dict['ra'][i],
            'dec':     tsObj_dict['dec'][i],
            'objc_flags':  tsObj_dict['objc_flags'][i],
            'objc_flags2': tsObj_dict['objc_flags2'][i],
            'mjd': mjd,
            'cmodel_u':     tsObj_dict['counts_model'][i][0],
            'cmodel_g':     tsObj_dict['counts_model'][i][1],
            'cmodel_r':     tsObj_dict['counts_model'][i][2],
            'cmodel_i':     tsObj_dict['counts_model'][i][3],
            'cmodel_z':     tsObj_dict['counts_model'][i][4],
            'cmodel_err_u':     tsObj_dict['counts_modelErr'][i][0],
            'cmodel_err_g':     tsObj_dict['counts_modelErr'][i][1],
            'cmodel_err_r':     tsObj_dict['counts_modelErr'][i][2],
            'cmodel_err_i':     tsObj_dict['counts_modelErr'][i][3],
            'cmodel_err_z':     tsObj_dict['counts_modelErr'][i][4],
            }
        #THIS WORKS ON anathem:# print "{run} {camCol} {rerun} {field} {id} {ra} {dec} {mjd} {cmodel_u:05.5f} {cmodel_g:05.5f} {cmodel_r:05.5f} {cmodel_i:05.5f} {cmodel_z:05.5f} {cmodel_err_u:05.5f} {cmodel_err_g:05.5f} {cmodel_err_r:05.5f} {cmodel_err_i:05.5f} {cmodel_err_z:05.5f} {objc_flags} {objc_flags2}".format(**table)
        print "{run} {camCol} {rerun} {field} {id} {ra} {dec} {mjd} {cmodel_u} {cmodel_g} {cmodel_r} {cmodel_i} {cmodel_z} {cmodel_err_u} {cmodel_err_g} {cmodel_err_r} {cmodel_err_i} {cmodel_err_z} {objc_flags} {objc_flags2}".format(**table)

    """
                'objc_flags':  tsObj_dict['objc_flags'][i],
                'objc_flags2': tsObj_dict['objc_flags2'][i],
            'objc_flags':  objc_flags,
            'objc_flags2': objc_flags2,


                'fracPSF':     tsObj_dict['fracPSF'][i],
                'counts_model':     tsObj_dict['counts_model'][i][0],
                'psfCounts':        tsObj_dict['psfCounts'][i][0],
                'mag_countsmodel':     mags_dict['mag_countsmodel'][i][0],
                'mag_countsmodel_ivar':mags_dict['mag_countsmodel_ivar'][i][0],
                'mag_psfcounts':       mags_dict['mag_psfcounts'][i][0],
                'mag_psfcounts_ivar':  mags_dict['mag_psfcounts_ivar'][i][0],
    """



if __name__ == '__main__':
    """ 
    wget http://das.sdss.org/imaging/7202/40/calibChunks/6/tsField-007202-6-40-0706.fit
    wget http://das.sdss.org/imaging/7202/40/objcs/6/fpObjc-007202-6-0706.fit
    wget http://das.sdss.org/imaging/7202/40/corr/6/fpC-007202-u6-0706.fit.gz
    """
    if len(sys.argv) > 1:
        dirpath = sys.argv[1]
    else:
        dirpath = '/global/homes/d/dstarr/paralleldb/SDSS/s82/tsobjSplit_cr_dir'

    fpaths = glob.glob("%s/tsObj*.fit" % (dirpath))

    fpc_dirpath = dirpath.replace('tsobj','fpc')
    if not os.path.exists(fpc_dirpath):
        os.mkdir(fpc_dirpath)
    
    for fpath in fpaths:
        #fpath = '/home/dstarr/src/starvars/tsObj-007202-6-40-0706.fit'
        tsObj_dict = extract_tsObj_data(fpath)
        #mags_dict = calculate_magnitudes(tsObj_dict)

        old_stderr = sys.stderr
        sys.stderr = open('/dev/null','w')
        mjd_list = parse_times_from_fpC_fits(tsObj_dict, fpc_dirpath=fpc_dirpath)
        sys.stderr.close()
        sys.stderr = old_stderr

        if len(mjd_list) > 0:
            print_sources(tsObj_dict=tsObj_dict, mjd_list=mjd_list) # , mags_dict=mags_dict

    # TODO: Check that the generated m(t) timeseries matches exsisting lightcurve catalogs?

    """
    with '/home/dstarr/src/starvars/tsField-007202-6-40-0706.fit' #
    'MJD(TAI)' look at sdss_astrophot.py:L454 to see how to use
    wcstools:xy2sky to convert xy displacement to ra,dec disp and then
    calculate actual time info.
    #tsField_data = extract_tsField_data(fpath)
    """
