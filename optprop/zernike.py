"""
This module contains functions required to compute the Zernike polynomials
"""

import pandas as pd
import pathlib
import os

package_dir = pathlib.Path(__file__).parent.resolve()
zernike_path = os.path.join(package_dir,"data","zernike_index.csv")

data = pd.read_csv(zernike_path, index_col=False, header=None)
data = data.to_numpy(dtype=int)

n = data[:,0]
m = data[:,1]
i = data[:,2]