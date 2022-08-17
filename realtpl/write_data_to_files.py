from dataclasses import dataclass
import pandas as pd
import os


def write_csv(df: pd.DataFrame, fp: dataclass, output_dir: str):
    path = (os.path.join(output_dir, fp.name, 'data'))
    os.makedirs(path, exist_ok=True)

    for kind, dff in df.groupby('kind'):
        dff.to_csv(os.path.join(path, str(kind) + '.csv'),
                   sep='\t', index=False)
