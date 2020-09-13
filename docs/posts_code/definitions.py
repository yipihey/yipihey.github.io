#%%
from scipy import stats
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.style.use('./Modified5308.mplstyle') 
figtransparent = 0

import scipy.spatial as sptl
import scipy.integrate

# %%
def peaked(y):
    """For all values of y<0.5 return y. For y>0.5 return 1-y. 
       This is useful for peaked emperical CDFs. """
    yp = y
    yp[y>0.5] = 1-yp[y>0.5]
    return yp

# %%
