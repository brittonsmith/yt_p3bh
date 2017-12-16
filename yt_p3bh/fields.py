"""
field stuff



"""

#-----------------------------------------------------------------------------
# Copyright (c) Britton Smith <brittonsmith@gmail.com>.  All rights reserved.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.utilities.physical_constants import G

def _mdot_bh_norm(field, data):
    val = data["gas", "density"].astype(np.float64) / \
          data["gas", "sound_speed"].astype(np.float64)**3
    return val * np.pi * G**2

def _bh_to_edd_norm(field, data):
    D = data.ds.quan(2.2e-8, "1/yr")
    return data["gas", "bondi_hoyle_accretion_rate"] / D

def add_fields(ds):
    ds.add_field(("gas", "bondi_hoyle_accretion_rate"),
                 function=_mdot_bh_norm,
                 units="1/(Msun*yr)",
                 sampling_type="cell")

    ds.add_field(("gas", "BH_to_Eddington"),
                 function=_bh_to_edd_norm,
                 units="1/Msun",
                 sampling_type="cell")
