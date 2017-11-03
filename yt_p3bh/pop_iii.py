"""
Pop III star functions.



"""

#-----------------------------------------------------------------------------
# Copyright (c) Britton Smith <brittonsmith@gmail.com>.  All rights reserved.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.units.yt_array import \
    YTArray, \
    YTQuantity

# Woosley & Weaver (1995)
_mprog = YTArray([12, 13, 15, 18, 20, 22,
                  25, 30, 35, 40, 140], "Msun")
_mrem = YTArray([1.35, 1.28, 1.53, 3.40, 4.12, 1.49,
                 6.36, 8.17, 12.8, 16.6, 65], "Msun")
def get_remnant_mass(mass):
    """
    Return black hole mass for a given Pop III stellar mass.
    """
    global _mprog, _mrem
    cc = mass < _mprog[-1]
    mrem = np.zeros_like(mass)
    if cc.any():
        i = np.digitize(mass[cc], _mprog) - 1
        i = np.clip(i, 0, _mprog.size-2)
        slope = (_mrem[i+1] - _mrem[i]) / \
          (_mprog[i+1] - _mprog[i])
        mrem[cc] = slope * (mass[cc] - _mprog[i]) + _mrem[i]
    pisn = mass >= _mprog[-1]
    if pisn.any():
        mrem[pisn] = (13. / 24) * \
          (mass[pisn] - YTQuantity(20, "Msun"))
    return mrem

# Schaerer (2002)
_ms_m = YTArray([5, 9, 15, 25, 40, 60, 80,
                 120, 200, 300, 400, 500], "Msun")
_ms_t = YTArray([6.190e7, 2.022e7, 1.040e7, 6.459e6, 3.864e6,
                 3.464e6, 3.012e6, 2.521e6, 2.204e6, 2.047e6,
                 1.974e6, 1.899e6], "yr")
def get_MS_lifetime(mass):
    """
    Return Pop III main-sequence lifetime for a given Pop III
    stellar mass.
    """
    global _ms_m, _ms_t
    i = np.digitize(mass, _ms_m) - 1
    i = np.clip(i, 0, _ms_m.size-2)
    slope = (_ms_t[i+1] - _ms_t[i]) / \
      (_ms_m[i+1] - _ms_m[i])
    t = slope * (mass - _ms_m[i]) + _ms_t[i]
    return t
