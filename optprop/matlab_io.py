"""
This module contains functions that interface with MATLAB
"""

from scipy.io import loadmat
import os.path as osp
from pathlib import Path
import numpy as np
import pandas as pd

def mat_to_ndarray(mat_path: str | Path) -> np.ndarray:
    """
    mat_to_ndarray takes MATLAB workspace variable (.mat file) and converts to ndarray assuming it is a matrix.
    
    Args:
        path: Path to .mat file
    Returns:
        data_array: Outputs a ndarray of the .mat file
    """
    fn,ext = osp.splitext(osp.basename(osp.normpath(mat_path)))
    data_mat = loadmat(mat_path)
    data_array = data_mat[fn]
    return np.asarray(data_array)

def mat_to_csv(mat_path: str | Path, csv_path: str | Path):
    """
    mat_to_csv takes MATLAB workspace variable (.mat file) and converts to a csv file assuming it is a matrix.
    
    Args:
        path: Path to .mat file
    """
    fn,ext = osp.splitext(osp.basename(osp.normpath(mat_path)))
    data_mat = loadmat(mat_path)
    data_array = data_mat[fn]
    df = pd.DataFrame(data_array)
    df.to_csv(csv_path,index=False,header=False)