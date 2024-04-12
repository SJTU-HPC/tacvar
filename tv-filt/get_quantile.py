# Reading three args from command line:
# arg1: dir, target directory; arg2: col, target column; arg3: quant, quantile.
# If a fourth argument is given, then it is the run time per vkern, print the number of vkerns.
# Reading all csvs in dir, calculate the median of the target column, and print the result.

import os
import numpy as np
import pandas as pd


def read_csvs(dir_path, col, q=0.5):
    dfs = []
    files = os.listdir(dir_path)
    for i in range(0, len(files)):
        df = pd.read_csv(dir_path + '/' + files[i], header=None)
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)
    arr = df[col].to_numpy()
    return np.quantile(arr, q)


if __name__ == '__main__':
    import sys
    csv_dir = sys.argv[1]
    data_col = int(sys.argv[2])
    quant = float(sys.argv[3])
    if len(sys.argv) == 4:
        print(read_csvs(csv_dir, data_col, quant))
    elif len(sys.argv) == 5:
        nspk = float(sys.argv[4])
        ns = read_csvs(csv_dir, data_col, quant)
        print(int(ns / nspk))