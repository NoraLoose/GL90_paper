import yaml
import numpy as np
from copy import deepcopy
import xarray as xr

class YAMLParser():
    """A simple class to read a YAML into a dict of dicts.
    This doesn't need to be a class, but maybe we'll want more customizations later.

    Note:
        To read in numbers written in scientific notation, they must have a "." to
        be recognized as a float rather than a string. For example:

            - tikhonov_parameter : 4.e-5

        will be read in as expected. However

            - tikhonov_parameter : 4e-5

        will be read in as a string.

    Note:
        Any blank input will be received as None type. Any string "None"
        will be converted to None type.

    Note:
        Any tab character will cause errors.
    """

    def __init__(self):
        pass


    def read(self, filename):
        """Read a yaml file and return a dict of dicts

        Args:
            filename (str): e.g. '/path/to/config.yaml'

        Returns:
            params (dict): with structure like the yaml file
        """
        params = {}
        with open(filename,'r') as file:

            doc = yaml.safe_load(file)

            for key,val in doc.items():
                if isinstance(val[0], dict):
                    params[key] = {}
                    for item in val:
                        params[key].update(item)

                    # --- Now go back and unpack any nested dicts
                    for k2,v2 in params[key].items():
                        this_param = params[key][k2]
                        if isinstance(this_param, list):
                            if np.all([isinstance(x,dict) for x in this_param]):
                                tmp = deepcopy(this_param)
                                params[key][k2] = {}
                                for item in tmp:
                                    params[key][k2].update(item)
                else:
                    params[key] = val


        # --- Convert strings with value "None" to be None type
        #     Do this for each dict in this big dict
        params = _convert_none(params)
        return params


    def write(self, params, filename='config.yaml'):
        """Dump a copy of the config file so you remember what you did!

        Args:
            params (dict): a massive dictionary with all options
            filename (str): path to dump config in
        """
        # first need to undo our nested dict operations
        out = {}
        for key,val in params.items():
            if isinstance(val,dict):
                out[key] = []
                for k2,v2 in val.items():
                    out[key].append({k2:v2})
            else:
                out[key] = val

        with open(filename,'w') as file:
            yaml.dump(out, stream=file)


def _convert_none(this_dict):
    for key,val in this_dict.items():
        if isinstance(val,dict):
            this_dict[key] = _convert_none(val)
        elif val == "None":
            this_dict[key] = None
    return this_dict


from xgcm import Grid 
def make_grid(st, ds):
    coords = {'X': {'center': 'xh', 'outer': 'xq'},
                'Y': {'center': 'yh', 'outer': 'yq'},
                #'Z': {'center': 'zl', 'outer': 'zi'} 
             }
    metrics = {('X',):['dxCu','dxCv','dxT','dxBu'],
               ('Y',):['dyCu','dyCv','dyT','dyBu']
              }

    if 'zi' in ds:
        st['zi'] = ds['zi']
        st['zl'] = ds['zl']
        coords = {'X': {'center': 'xh', 'outer': 'xq'},
                 'Y': {'center': 'yh', 'outer': 'yq'},
                 'Z': {'center': 'zl', 'outer': 'zi'} 
             }

    grid = Grid(st, coords=coords, periodic=['X'])

    st['dxT'] = grid.interp(st.dxCu,'X')
    st['dyT'] = grid.interp(st.dyCv,'Y')
    st['dxBu'] = grid.interp(st.dxCv,'X')
    st['dyBu'] = grid.interp(st.dyCu,'Y',boundary='fill')

    grid = Grid(st, coords=coords, periodic=['X'], metrics=metrics)
    return grid


###### read parameterized runs
import glob, os
from os.path import exists

def read_parameterized_runs(exps, read_from_time_list=True, read_snapshots=False):
    chunks = {'time': 1}
    for exp in exps:
        print(exp)
        os.chdir('%s' %exps[exp]['dir'])
        if read_from_time_list:
            file_list = []
            for time in exps[exp]['times']:
                filename = '%s/averages_%08d.nc' % (exps[exp]['dir'], time)
                file_list.append(filename)
            ds = xr.open_mfdataset(file_list , decode_times=False, chunks=chunks, combine='by_coords', parallel=True)
        else:
            time_list = []
            for file in glob.glob('averages*'):
                time_list.append(int(file[9:17]))
            found = False
            while found == False:
                max_time = max(time_list)
                filename = '%s/averages_%08d.nc' % (exps[exp]['dir'], max_time)
                #print(filename)
                ds = xr.open_dataset(filename, decode_times=False, chunks=chunks)
                if len(ds.time) == 0:
                    time_list.remove(max_time)
                else:
                    found = True

        exps[exp]['ds'] = ds
        
        ## snapshots
        if read_snapshots:
            if read_from_time_list:
                filename = '%s/snapshots_%08d.nc' % (exps[exp]['dir'], exps[exp]['time'])
                sn = xr.open_dataset(filename , decode_times=False, chunks=chunks)
            else:
                time_list = []
                for file in glob.glob('snapshots*'):
                    time_list.append(int(file[10:18]))
                if len(time_list):
                    found = False
                    while found == False:
                        max_time = max(time_list)
                        filename = '%s/snapshots_%08d.nc' % (exps[exp]['dir'], max_time)
                        sn = xr.open_dataset(filename, decode_times=False, chunks=chunks)
                        if len(sn.time) == 0:
                            time_list.remove(max_time)
                        else:
                            found = True
            exps[exp]['sn'] = sn

        stats = xr.open_dataset('%s/ocean.stats.nc' % (exps[exp]['dir']), decode_times=False)
        for start, end in zip(exps[exp]['cut_start'], exps[exp]['cut_end']):
            stats0 = stats.isel(Time=slice(None, start))
            stats1 = stats.isel(Time=slice(end, None))
            stats = xr.concat([stats0, stats1], dim='Time')
        exps[exp]['os'] = stats
            
            
        filename = '%s//static.nc' % (exps[exp]['dir'])
        if exists(filename):
                st = xr.open_dataset(filename, decode_times=False)
                grid = make_grid(st, ds)
                exps[exp]['st'] = st
                exps[exp]['grid'] = grid

    return exps
        
###### read unparameterized runs

def read_stats(run):
    chunks = {'time': 1}
    filename = '%s/ocean_stats.nc' %(run)
    stats = xr.open_dataset(filename, decode_times=False)
    return stats

def read_unparameterized_runs(exps, read_from_time_list=True, read_snapshots=False):
    chunks = {'time': 1}
    for exp in exps:
        if read_from_time_list:
            file_list = []
            for time in exps[exp]['times']:
                filename = '%s/averages_%08d.nc' % (exps[exp]['dir'], time)
                file_list.append(filename)
            ds = xr.open_mfdataset(file_list , decode_times=False, chunks=chunks, combine='by_coords', parallel=True)
        else:
            ds = xr.open_mfdataset('%s/averages_*.nc' %exps[exp]['dir'], decode_times=False, chunks=chunks, combine='by_coords')
            ds = ds.isel(time=slice(-100, None))
        exps[exp]['ds'] = ds
        
        # time averaged diagnostics    
        filename = '%s/time_averaged_diags_500days' %exps[exp]['dir_work']
        file_exists = exists(filename)
        if file_exists:
                dst = xr.open_zarr(filename, decode_times=False, chunks=chunks)
                exps[exp]['dst'] = dst

        ## snapshots
        if read_snapshots:
            if exps[exp]['time'] is not None:
                filename = '%s/snapshots_%08d.nc' %(exps[exp]['dir'], exps[exp]['time']) 
                ds = xr.open_dataset(filename, decode_times=False, chunks=chunks)
            else: # for 1/32 degree case
                ds = xr.open_mfdataset('%s/snapshots_*.nc' %exps[exp]['dir'], decode_times=False, chunks=chunks, combine='by_coords')
                ds = ds.isel(time=slice(-100, None))
            exps[exp]['sn'] = ds

        stats = read_stats(exps[exp]['dir_orig'])
        exps[exp]['os'] = stats
            
            
        filename = '%s//static.nc' % (exps[exp]['dir'])
        st = xr.open_dataset(filename, decode_times=False)
        grid = make_grid(st, ds)
        exps[exp]['st'] = st
        exps[exp]['grid'] = grid
        
        degree = exps[exp]['degree']

    return exps

def read_offline_runs(exps):
    chunks = {'time': 1}
    time_list = [2402, 2502, 2602, 2702, 2802]
    for exp in exps:
        #ds
        for i, time in zip(range(len(time_list)), time_list):
            filename = '%s/averages_%08d_filtered_fac%i' %(exps[exp]['dir_scratch'], time, exps[exp]['filter_fac']) 
            ds_tmp = xr.open_zarr(filename, decode_times=False, chunks=chunks)
    
            if i == 0:
                ds = ds_tmp
            else:
                ds = xr.combine_nested([ds, ds_tmp], concat_dim='time', combine_attrs='drop_conflicts') 
        exps[exp]['ds'] = ds

        # time averaged diagnostics    
        filename = '%s/time_averaged_diags_fac%i_500days' %(exps[exp]['dir_work'], exps[exp]['filter_fac'])
        file_exists = exists(filename)
        if file_exists:
                dst = xr.open_zarr(filename, decode_times=False, chunks=chunks)
                exps[exp]['dst'] = dst

        #energy diagnostics
        for i, time in zip(range(len(time_list)), time_list):
            filename = '%s/bleck_cycle_%08d_fac%i' %(exps[exp]['dir_scratch'], time, exps[exp]['filter_fac']) 
            ds_tmp = xr.open_zarr(filename, decode_times=False, chunks=chunks)
    
            if i == 0:
                ds = ds_tmp
            else:
                ds = xr.combine_nested([ds, ds_tmp], concat_dim='time', combine_attrs='drop_conflicts') 
        exps[exp]['ds_energy'] = ds
 
        filename = '%s/bleck_cycle_fac%i_500days.nc' % (exps[exp]['dir_work'], exps[exp]['filter_fac'])
        ds = xr.open_dataset(filename, decode_times=False)

        filename = '%s/MKE2MPE_TWA_inaccurate_fac%i_500days.nc' % (exps[exp]['dir_work'], exps[exp]['filter_fac'])
        file_exists = exists(filename)
        if file_exists:
            ds2 = xr.open_dataset(filename, decode_times=False)
            ds['MKE_to_MPE_TWA_inaccurate'] = ds2['MKE_to_MPE_TWA_inaccurate']

        exps[exp]['dst_energy'] = ds
        #st and grid
        filename = '/glade/campaign/univ/unyu0004/NeverWorld2/nw2_0.03125deg_N15_baseline_hmix20/static.nc'
        st = xr.open_dataset(filename, decode_times=False)
        grid = make_grid(st, ds)
        exps[exp]['st'] = st
        exps[exp]['grid'] = grid

    return exps

