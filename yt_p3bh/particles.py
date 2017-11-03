"""
particle stuff



"""

#-----------------------------------------------------------------------------
# Copyright (c) Britton Smith <brittonsmith@gmail.com>.  All rights reserved.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt_p3bh.pop_iii import \
    get_MW_lifetime

def _pop3_black_hole(pfilter, data):
    ct = data["creation_time"]
    stars = (ct > 0)
    if not stars.any():
        return stars

    Zcr = data.ds.parameters['PopIIIMetalCriticalFraction']
    Zfr = data["metallicity_fraction"]
    stars &= (Zfr < Zcr)
    if not stars.any():
        return stars

    pt = data["particle_type"]
    stars &= ((pt == 1) | (pt == 5))
    if not stars.any():
        return stars

    # mass is multiplied by 1e-20 when main-sequence lifetime is over
    mass = data["particle_mass"].to("Msun")
    mass[mass < 1e-9] *= 1e20
    stars &= (((mass >= 25) & (mass <= 140)) | (mass >= 260))
    if not stars.any():
        return stars

    t_ms = get_MS_lifetime(mass[stars])
    age = (data.ds.current_time - ct[stars]).to(t_ms.units)
    stars[stars] &= (age > t_ms)
    return stars
