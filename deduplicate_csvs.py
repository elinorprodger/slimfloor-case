"""
deduplicate_csvs.py

Removes duplicate rows from all results_*.csv files caused by running
the optimisation scripts multiple times without deleting the CSV first.

The save function in plot_datasets.py uses mode='a' (append), so each
rerun adds another full copy of the results.  This script deduplicates
in-place, keeping only the first occurrence of each unique row.

Run once:
    python deduplicate_csvs.py
"""

import glob
import os
import pandas as pd

csv_files = sorted(glob.glob('results_*.csv'))

if not csv_files:
    print('No results_*.csv files found.')
else:
    for path in csv_files:
        df = pd.read_csv(path)
        before = len(df)
        df_clean = df.drop_duplicates()
        after = len(df_clean)
        removed = before - after
        if removed > 0:
            df_clean.to_csv(path, index=False)
            print(f'{os.path.basename(path)}: {before} rows -> {after} rows '
                  f'({removed} duplicates removed)')
        else:
            print(f'{os.path.basename(path)}: {before} rows — already clean')

print('\nDone.')
