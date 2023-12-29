#!/usr/bin/env python3

import argparse
import os
import sys
import numpy as np
from pathlib import Path
from datetime import datetime as dt
from mintpy.objects import ifgramStack
from mintpy.cli import view

def create_parser():
    parser = argparse.ArgumentParser(description=
        'Plot ifgs per temporal baseline')
    parser.add_argument('-i', '--ifgs_stack', type=str,
        help='Input Mintpy ifgramStack.h5', dest='ifg_file')
    parser.add_argument('-g', '--geom_stack', default=None, type=str,
        help='Input Mintpy geometryGeo.h5', dest='geom_file')
    parser.add_argument('-o', '--output', default='./ifgs_plots',
                        help='Output dir', type=str, dest='out_dir')
    parser.add_argument('-w', '--wrap', action='store_true',
        help='Plot wrapped', dest='wrap')
    parser.add_argument('-wr', '--wrap_range', nargs=2, 
        metavar=('VMIN', 'VMAX'), type=float, default=[-3.14, 3.14],
        help='Display limits for wrapped plotting', dest='wrap_range')
    parser.add_argument('-s', '--start', dest='start', type=int, 
        help='Start temp baseline for plotting in days', default=-1)
    parser.add_argument('-e', '--end', dest='end', type=int, 
        help='End temp baseline for plotting in days', default=-1)
    parser.add_argument('--all', action='store_true',
        help='Plot unwrap and wrapped phase w 3.14, 10, 20', dest='all')
    parser.add_argument('-v', '--verbose', action='store_true',
        help='Display print messages', dest='verbose')
    return parser

def fix_datetime(x):
    return dt.strftime(dt.strptime(str(np.squeeze(x)), '%Y%m%d'),
                '%Y-%m-%d')

def plot_stack(ifgram_file, outdir, 
               start=-1, end=-1, 
               geom_file=None, 
               wrap=False, wrap_range=[-3.14, 3.14],
               verbose=True):
    try:
        stack_obj = ifgramStack(Path(ifgram_file).resolve())
        stack_obj.open()
    except:
        raise ValueError(f'Cannot open {ifgram_file}')
    
    # Initialize outdir
    out = Path(outdir).resolve()

    # Get date list
    date12List = stack_obj.get_date12_list(dropIfgram=False)
    dates12array = np.vstack([d.split('_') for d in date12List])
    date1, date2 = dates12array[:,0], dates12array[:,1]
    stack_obj.close()

    # Get temporal baselines
    date1d = np.vstack(np.apply_along_axis(fix_datetime, 
                0, np.atleast_2d(date1))).astype('datetime64[D]')
    date2d = np.vstack(np.apply_along_axis(fix_datetime,
                0, np.atleast_2d(date2))).astype('datetime64[D]')

    t_baseline = (date2d-date1d).astype(int).flatten()
    unique_tbase, tcount = np.unique(t_baseline, return_counts=True)
    
    # Plot only specific temporal baseline range
    if start == -1: start = unique_tbase[0]
    if end == -1: end = unique_tbase[-1] 
    
    # Prepare Mintpy view.py command parameters
    cmd =  [f'{ifgram_file}']
    if geom_file is not None:
        cmd += ['-d',f'{geom_file}']
    if wrap is True:
        cmd += ['--wrap', '-c' ,'cmy',  
                '--wrap-range', str(wrap_range[0]), str(wrap_range[1])]
        flag = f'wrap_{np.int16(np.abs(wrap_range[0]))}'
        out = out / flag
    else:
        flag = 'unwrap'
        out = out / flag

    cmd += ['--zm', '--nrows', '2', '--ncols', '5',
               '--save', '--nodisplay', '-n']

    # Create output
    out.mkdir(parents=True, exist_ok=True)

    # Plot
    for tbase in unique_tbase:
        if (int(tbase) >= int(start)) & (int(tbase) <= int(end)):
            pairs = (np.where(t_baseline == tbase)[0])
            plist = [str(p) for p in pairs.tolist()]
            if verbose: 
                print(f'Plotting ifgs w {tbase} temp baseline')
            view.main(cmd + plist + ['-o', str(out / f'ifg_{flag}_{tbase}d.png')])

def main(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    
    if inps.all is True:
        wrap_range = [0, 10, 20]
        for wr in wrap_range:
            wrap_flag = False if wr == 0 else True
            # Plot
            plot_stack(inps.ifg_file, inps.out_dir, 
               inps.start, inps.end, 
               inps.geom_file, 
               wrap_flag, [-wr, wr],
               inps.verbose) 
    else:           
        # plot
        plot_stack(inps.ifg_file, inps.out_dir, 
                inps.start, inps.end, 
                inps.geom_file, 
                inps.wrap, inps.wrap_range,
                inps.verbose)

if __name__ == '__main__':
    main(sys.argv[1:])

