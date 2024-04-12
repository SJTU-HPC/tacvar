# binw = max (10, (mid - min) / 50)

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
    gap = np.quantile(arr, 0.5) - np.quantile(arr, 0)
    binw = max(gap / 50, 10)
    binw = int(binw / 10) * 10
    return binw


if __name__ == '__main__':
    import sys
    csv_dir = sys.argv[1]
    data_col = int(sys.argv[2])
    print(read_csvs(csv_dir, data_col))
