# Reading three args from command line:
# arg1: dir, target directory; arg2: col, target column
# Reading all csvs in dir, aggregate target column and dump to a csv file.

import os
import numpy as np
import pandas as pd


def read_csvs(dir_path, col):
    dfs = []
    files = os.listdir(dir_path)
    for i in range(0, len(files)):
        df = pd.read_csv(dir_path + '/' + files[i], header=None)
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)
    arr = df[col].to_numpy()
    np.savetxt('./met.csv', arr, fmt='%d' , delimiter='\n')


if __name__ == '__main__':
    import sys
    csv_dir = sys.argv[1]
    data_col = int(sys.argv[2])
    read_csvs(csv_dir, data_col)
    print("Done.")