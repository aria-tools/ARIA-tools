#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime 
from mintpy.objects import ifgramStack

def create_parser():
    parser = argparse.ArgumentParser(description=
        'Plot ifgs per temporal baseline')
    parser.add_argument('-i', '--ifgs_stack', type=str,
        help='Input Mintpy ifgramStack.h5', dest='ifgram_file')
    parser.add_argument('-n', default=[],  nargs='+',
                        help='List of ifg to drop', type=str, dest='dropList')
    parser.add_argument('-o', '--output', default='./date12_list.txt',
                        help='Output dir', type=str, dest='output')
    parser.add_argument('-s', '--start', dest='start', type=int, 
        help='Start temporal baseline in days', default=None)
    parser.add_argument('-e', '--end', dest='end', type=int, 
        help='End temporal baseline  in days', default=None)
    parser.add_argument('--invert', action='store_true',
        help='Invert selection', dest='invert_flag')
    parser.add_argument('-v', '--verbose', action='store_true',
        help='Display print messages', dest='verbose')
    return parser

def _drop_ix(x, df):
    # Option 1: '20150921_20160319'
    if len(x.split('_')) > 1:
        ix = df[df.date12.isin([x])].index
        if len(ix) ==0:
            print(f'    Skip: {x}')
            ix = None
        else:
            ix = ix[0]
    # Option 2: string index, e.g. 0
    else:
        ix = int(x)
    return ix

def get_stack_ix(ifgram_file, dropList=[], 
                 start=None, end=None,
                 invert_flag=False, verbose=False):
    try:
        stack_obj = ifgramStack(Path(ifgram_file).resolve())
        stack_obj.open()
    except:
        raise ValueError(f'Cannot open {ifgram_file}')

    # Get dates
    date12List = stack_obj.get_date12_list(dropIfgram=False)
    dates12array = np.vstack([d.split('_') for d in date12List])
    
    # Create Dataframe
    cols = ['date12', 'date1', 'date2']
    df = pd.DataFrame(np.c_[date12List, dates12array], columns=cols)

    # Add temporal baseline
    date_str2dt = lambda x: datetime.strptime(x, '%Y%m%d')
    d1 = df.date1.apply(date_str2dt)
    d2 = df.date2.apply(date_str2dt)
    df['btemporal'] = (d2 - d1).dt.days

    # Get drop_index could be both: index and value
    #  e.g., 0, '20150921_20160319'
    dropList = [_drop_ix(date, df) for date in dropList]
    dropList = [date for date in dropList if date is not None]

    # Plot only specific temporal baseline range
    if start is None: start = df.btemporal.min()
    if end is None: end = df.btemporal.max()

    # Filter based on start and end
    df = df[(df.btemporal >= start) & (df.btemporal <= end)]

    # Get the drop flag
    flag = ~df.index.isin(dropList)

    if invert_flag:
        flag = ~flag

    if verbose:
        print(f'Number of kept ifgs: {np.sum(flag)}')
        print(f'Number of dropped ifgs: {np.sum(~flag)}')
        print(df[~flag])
        print(f'Drop indexes: {df[~flag].index.values}')

    return df[flag], df[~flag]

def main(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    
    print(inps)
    df, drop_df = get_stack_ix(inps.ifgram_file, 
                      inps.dropList, 
                      inps.start,
                      inps.end,
                      inps.invert_flag, 
                      inps.verbose)

    # Initialize outdir
    out = Path(inps.output).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    # Overwrite
    if out.exists(): out.unlink()

    # Get index values for output
    print(f'Writing referenceFile to {out}')
    out_df = df.copy()
    out_df['index'] = df.index.values
    f = open(str(out.parent / 'drop_ix.txt'), 'a')
    f.write(f'# Drop ix: {drop_df.index.values}\n')
    f.close()
    out_df.to_csv(str(out), columns=['date12', 'index'],
                sep=' ', header=False, index=False)
    

if __name__ == '__main__':
    main(sys.argv[1:])