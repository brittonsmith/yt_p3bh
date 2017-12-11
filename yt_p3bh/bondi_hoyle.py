"""
Bondi-hoyle accretion functions



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

def bondi_hoyle_accretion_rate(data, current_mass):
    """
    Calculate Bondi-hoyle accretion rate for a data container
    and a black hole mass.
    """
    ds = data.ds
    v_rel = ds.arr([(data["particle_velocity_%s" % ax] -
                     data["velocity_%s" % ax])
                     for ax in "xyz"])
    v_rel_mag = np.sqrt((v_rel**2).sum(axis=0))
    v_BH = ds.arr([data["sound_speed"], v_rel_mag]).max(axis=0)

    mdot = ((np.pi * data["density"] * G**2 * current_mass**2) /
            v_BH**3).to("Msun/yr")
    return mdot
