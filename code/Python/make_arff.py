#!/usr/bin/env python
import os
import glob
import numpy as np
import cStringIO

def make_combined_200k_arff(combo_arff_fpath = "/home/dstarr/src/hakeem_20130309/starvars/200k_combo.arff"):
    """
    """
    files = glob.glob("/home/dstarr/src/hakeem_20130309/200k_arffs/*/*.arff")
    
    combo_arff_lines = []
    srcid_list = []
    for i_fpath, fpath in enumerate(files):
        lines = open(fpath).readlines()
        in_header = True
        i_header = 0
        header_list = []
        for raw_line in lines:
            line = raw_line.strip()
            if len(line) == 0:
                continue
            if line == '@data':
                in_header = False
                continue
            elif ((line[0] == '%')):
                ### comment
                continue #print line
            elif line[:10] == '@ATTRIBUTE':
                ### in header
                print 
                header_list.append(line)
                i_header += 1
                #print line
                continue
            if not in_header:
                ### data instance lines
                combo_arff_lines.append(line)
                srcid_list.append(int(line[:line.find(',')]))

        if i_fpath == 0:
            ### first arff file
            first_header_list = header_list
        else:
            if first_header_list != header_list:
                print "mismatch"
                import pdb; pdb.set_trace()
                print
            #import pdb; pdb.set_trace()
            #print
    fp = open(combo_arff_fpath, 'w')
    for e in first_header_list:
        fp.write("%s\n" % (e))
    fp.write("@data\n")
    for e in combo_arff_lines:
        fp.write("%s\n" % (e))
    fp.close()

    ### Sanity check:
    n_set = len(set(srcid_list))
    n_list =  len(srcid_list)
    print "n_set=%d n_list=%d" % (n_set, n_list)


def parse_supplemental_data():
    """
#sdss_ra sdss_dec ra dec objectID objtype mag_median rExt noSATUR uMod gMod rMod iMod zMod uErr gErr rErr iErr zErr J H K JErr HErr KErr stdev rms chi2pdf nObs skew kurt p1 phi1 p2 phi2 p3 phi3 spec RRc Rchi2pdf
    
    """
    train_summary_fpath = "/home/dstarr/src/hakeem_20130309/starvars/masterMain.dat.txt"
    data_tups = np.genfromtxt(train_summary_fpath,
                              delimiter=" ",
                              filling_values=np.nan,
                              dtype=None,
                              names=True)

    return data_tups


def parse_200k_arff(combo_arff_fpath=""):
    """

    """
    lines = open(combo_arff_fpath).readlines()
    in_header = True
    i_header = 0
    header_list = []
    feature_name_list = []
    i_data = 0
    for i, raw_line in enumerate(lines):
        line = raw_line.strip()
        if len(line) == 0:
            continue
        if line == '@data':
            in_header = False
            continue
        elif ((line[0] == '%')):
            ### comment
            continue #print line
        elif line[:10] == '@ATTRIBUTE':
            ### in header
            print 
            header_list.append(line)
            elems = line.split(' ')
            feature_name_list.append(elems[1])
            i_header += 1

            #print line
            continue
        if not in_header:
            ### data instance lines
            i_data = i
            break # done

    data_str = '#' + ','.join(feature_name_list) + '\n' + ''.join(lines[i_data:])
    fp = cStringIO.StringIO(data_str)
    data_tups = np.genfromtxt(fp,
                              delimiter=",",
                              filling_values=np.nan,
                              dtype=None,
                              names=True,
                              missing_values="?")
    return data_tups


def combine_suplemental_with_200k(data_200k=None, data_sup=None, combined_csv_fpath=''):
    """
    """        
    id_index_sup = {}
    for i, srcid_sup in enumerate(data_sup['objectID']):
        id_index_sup[srcid_sup] = i

    jh_list = []
    hk_list = []

    for srcid in data_200k['source_id']:
        i_sup = id_index_sup[srcid]
        jh_list.append(data_sup['J'][i_sup] - data_sup['H'][i_sup])
        hk_list.append(data_sup['H'][i_sup] - data_sup['K'][i_sup])

    header_name_sublist = []
    for attrib_name in data_200k.dtype.names:
        if attrib_name in ['color_diff_bj', 'color_diff_rj', 'color_diff_vj','color_bv_extinction']:
            continue
        header_name_sublist.append(attrib_name)

    header = '#' + ",".join(header_name_sublist)
    csv_row_list = [header]
    for i in xrange(len(data_200k['source_id'])):
        row_list = []
        for attrib_name in header_name_sublist:
            if attrib_name == 'color_diff_jh':
                row_list.append(str(jh_list[i]))
            elif attrib_name == 'color_diff_hk':
                row_list.append(str(hk_list[i]))
            elif data_200k[attrib_name][i] == True:
                row_list.append('?')
            else:
                row_list.append(str(data_200k[attrib_name][i]))
        line = ','.join(row_list)
        csv_row_list.append(line)


    fp = open(combined_csv_fpath, 'w')
    fp.write('\n'.join(csv_row_list))
    fp.close()
    import pdb; pdb.set_trace()
    print



if __name__ == '__main__':


    combo_arff_fpath = "/home/dstarr/src/hakeem_20130309/starvars/200k_combo.arff" 
    final_combined_csv_fpath = "/home/dstarr/src/hakeem_20130309/starvars/200k_final_combo.csv" 

    ### This takes smaller .arff files and combines into a single large .arff:
    #make_combined_200k_arff(combo_arff_fpath=combo_arff_fpath)

    #### This parses the larger arff file into a data structure:
    data_200k = parse_200k_arff(combo_arff_fpath=combo_arff_fpath)


    data_sup = parse_supplemental_data()

    combine_suplemental_with_200k(data_200k=data_200k, data_sup=data_sup,
                                  combined_csv_fpath=final_combined_csv_fpath)

