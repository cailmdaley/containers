'''
Given a list of arguments and descriptors creates the python and idl functions for store and load the hdf5 files.
Right now it assumes a flat hdf5 file structure
'''



import numpy as np
import string

'''
Notes on format:

R = required


arg_name: R the name of the argument to use
is_arr:  R  true or false, if it is an array set to true if scalar set to false


'''

rw_python_file_name = '../../../mapspectra/python/bi_inout.py'
USE_H5PY = False


bispec_in_args = [ {'arg_name':'map',                'is_arr':True},
                   {'arg_name':'reso_arcmin',        'is_arr':False},
                   {'arg_name':'apod_mask',          'is_arr':True},
                   {'arg_name':'winf_grid',          'is_arr':True},
                   #{'arg_name':'winf_thresh',        'is_arr':False},
                   {'arg_name':'ps_grid',            'is_arr':True},

                   {'arg_name':'cmb_ells',           'is_arr':True},
                   {'arg_name':'cmb_cls',            'is_arr':True},
                   #{'arg_name':'ellvec',             'is_arr':True},
                   {'arg_name':'field_name',         'is_arr':False},
                   {'arg_name':'bands',               'is_arr':True},
                   {'arg_name':'pt_clust_mask',      'is_arr':True},

                   #{'arg_name':'test_bispec_out',    'is_arr':True},
                   #{'arg_name':'test_weight_out',    'is_arr':True},

                 ]


bispec_out_args = [ {'arg_name':'bispec',           'is_arr':True},
                    {'arg_name':'bispec_weight',    'is_arr':True},
                    {'arg_name':'field_name',       'is_arr':False},
                    {'arg_name':'bands',            'is_arr':True},
                    {'arg_name':'grid_bands',       'is_arr':True},
                    {'arg_name':'bispec_tag',        'is_arr':True},
                    {'arg_name':'ell_vec',        'is_arr':True},
                    {'arg_name':'delta_ell',        'is_arr':False}
                 ]



bispec_field_out_args = [ {'arg_name':'bispec',           'is_arr':True},
                          {'arg_name':'bispec_variance',    'is_arr':True},
                          {'arg_name':'bispec_weight_used', 'is_arr':True},
                          {'arg_name':'field_name',       'is_arr':False},
                          {'arg_name':'bands',            'is_arr':True},
                          {'arg_name':'grid_bands',       'is_arr':True},
                          {'arg_name':'bispec_tag',        'is_arr':True},
                          {'arg_name':'ell_vec',        'is_arr':True},
                          {'arg_name':'delta_ell',        'is_arr':False}
                          ]


debug_info_args = [ {'arg_name':'itfgrid',          'is_arr':True},
                    {'arg_name':'fftm_real',          'is_arr':True},
                    {'arg_name':'fftm_imag',          'is_arr':True},
                    {'arg_name':'ellgrid',          'is_arr':True},
                 ]





def get_idl_file_name(desc_name, is_store = True, idl_folder = 'io_funcs/'):
    type_desc = 'store' if is_store else 'load'
    return '%s/%s_%s_hdf5.pro' % (idl_folder, type_desc, desc_name )



def get_idl_store(desc_name, f_args):
    out_str = 'pro store_%s_hdf5, out_file, ' % desc_name
    for arg in f_args:
        out_str += arg['arg_name'] + ', '
    out_str = out_str[:-2] + '\n  fid = H5F_CREATE(out_file)\n'
    for arg in f_args:
        if arg['is_arr']:
            t_str = 'array'
        else:
            t_str = 'scalar'
        out_str += "  simple_store_hdf5_%s, fid, %s, '%s'\n" %(t_str, arg['arg_name'], arg['arg_name']  )
    out_str += '  H5F_CLOSE,fid\nend\n\n'
    return out_str

def get_idl_load(desc_name, f_args):
    out_str = 'pro load_%s_hdf5, in_file, ' %desc_name
    for arg in f_args:
        out_str += arg['arg_name'] + ', '
    out_str = out_str[:-2] + '\n  fid = H5F_OPEN(in_file)\n'
    for arg in f_args:
        out_str += "  %s = simple_load_hdf5( fid,  '%s')\n" %( arg['arg_name'], arg['arg_name']  )
        if  not arg['is_arr']:
            out_str += '%s = %s[0]\n' %(arg['arg_name'], arg['arg_name'])
    out_str += '  H5F_CLOSE,fid\nend\n\n'
    return out_str






def get_python_load_h5py(desc_name, f_args):
    out_str = "def load_%s_hdf5(in_file):\n" % (desc_name)
    out_str += "    fid = h5py.File(in_file, 'r')\n"
    out_str += "    out_dic = {}\n"
    for arg in f_args:
        out_str += "    out_dic['%s'] = fid['/%s'].value\n"%(arg['arg_name'], arg['arg_name'])
    out_str += "    fid.close()\n"
    out_str += "    return out_dic\n\n"        
    return out_str

def get_python_store_h5py(desc_name, f_args):
    out_str = "def store_%s_hdf5(out_file, " %desc_name
    for arg in f_args:
        out_str += '%s, '%arg['arg_name']
    out_str = out_str[:-2]+ "):\n"
    out_str += "    fid = h5py.File(out_file, 'w')\n"
    for arg in f_args:
        out_str += "    fid['/%s'] = %s\n"%(arg['arg_name'], arg['arg_name'])
    out_str += "    fid.close()\n\n"
    return out_str




def get_python_load_tables(desc_name, f_args):
    out_str = "def load_%s_hdf5(in_file):\n" % (desc_name)
    out_str += "    fid = tables.openFile(in_file, 'r')\n"
    out_str += "    out_dic = {}\n"
    for arg in f_args:
        out_str += "    out_dic['%s'] = fid.getNode('/%s').read()\n"%(arg['arg_name'], arg['arg_name'])
    out_str += "    fid.close()\n"
    out_str += "    return out_dic\n\n"        
    return out_str

def get_python_store_tables(desc_name, f_args):
    out_str = "def store_%s_hdf5(out_file, " %desc_name
    for arg in f_args:
        out_str += '%s, '%arg['arg_name']
    out_str = out_str[:-2]+ "):\n"
    out_str += "    fid = tables.openFile(out_file, 'w')\n"
    for arg in f_args:
        out_str += "    fid.createArray('/','%s', %s)\n"%(arg['arg_name'], arg['arg_name'])
    out_str += "    fid.close()\n\n"
    return out_str


if USE_H5PY:
    python_rw_header = '''import h5py
import numpy as np
'''
    get_python_store = get_python_store_h5py
    get_python_load = get_python_load_h5py

else:
    python_rw_header = '''import tables
import numpy as np
'''
    get_python_store = get_python_store_tables
    get_python_load = get_python_load_tables
    

'''    
fid = h5py.File('test_store.h5', 'r')
val = fid['/map'].value


fid['/test2'] = np.array([3,2], dtype = 'float64')
'''

if __name__ == '__main__':

    f = open(get_idl_file_name('debug_info', is_store = True), 'w')
    f.write(get_idl_store('debug_info', debug_info_args))
    f.close()


    f = open(get_idl_file_name('bispec_input', is_store = True), 'w')
    f.write(get_idl_store('bispec_input', bispec_in_args))
    f.close()

    f = open(get_idl_file_name('bispec_input', is_store = False), 'w')
    f.write(get_idl_load('bispec_input', bispec_in_args))
    f.close()


    f = open(get_idl_file_name('bispec_output', is_store = True), 'w')
    f.write(get_idl_store('bispec_output', bispec_out_args))
    f.close()

    f = open(get_idl_file_name('bispec_output', is_store = False), 'w')
    f.write(get_idl_load('bispec_output', bispec_out_args))
    f.close()



    f = open(get_idl_file_name('bispec_field_output', is_store = True), 'w')
    f.write(get_idl_store('bispec_field_output', bispec_field_out_args))
    f.close()

    f = open(get_idl_file_name('bispec_field_output', is_store = False), 'w')
    f.write(get_idl_load('bispec_field_output', bispec_field_out_args))
    f.close()




    f = open(rw_python_file_name, 'w')
    f.write(python_rw_header)
    f.write( get_python_store('bispec_input', bispec_in_args))
    f.write(get_python_load('bispec_input', bispec_in_args))

    #f.write( get_python_store('fake_bispec_input', fake_bispec_in_args))
    #f.write(get_python_load('fake_bispec_input', fake_bispec_in_args))

    f.write( get_python_store('bispec_output', bispec_out_args))
    f.write(get_python_load('bispec_output', bispec_out_args))

    f.write( get_python_store('bispec_field_output', bispec_field_out_args))
    f.write(get_python_load('bispec_field_output', bispec_field_out_args))


    f.write( get_python_store('debug_info', debug_info_args))
    f.write(get_python_load('debug_info', debug_info_args))


    f.write('''
if __name__ == '__main__':
    print(load_bispec_input_hdf5('../test_store.h5'))
''')


    f.close()
