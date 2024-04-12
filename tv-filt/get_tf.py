# Reading four args from command line:
# arg1: dir, target directory; arg2: nsamp_col, NSAMP; arg3: target_col, measured time col; ar4:ns_per_samp
# Reading all csvs in dir, calculating t=met_time - nsamp * ns_per_samp, write results to a csv file.

import os
import numpy as np
import pandas as pd


def read_csvs(dir_path, scol, mcol, nsps):
    dfs = []
    files = os.listdir(dir_path)
    for i in range(0, len(files)):
        df = pd.read_csv(dir_path + '/' + files[i], header=None)
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)
    arr = df[mcol].to_numpy() - df[scol].to_numpy() * nsps
    np.savetxt('./tf.csv', arr, fmt='%d' , delimiter='\n')


if __name__ == '__main__':
    import sys
    csv_dir = sys.argv[1]
    nsamp_col = int(sys.argv[2])
    data_col = int(sys.argv[3])
    ns_per_samp = float(sys.argv[4])
    read_csvs(csv_dir, nsamp_col, data_col, ns_per_samp)
    print("Done.")