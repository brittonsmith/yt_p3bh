"""
Pop III black hole analysis functions.



"""

#-----------------------------------------------------------------------------
# Copyright (c) Britton Smith <brittonsmith@gmail.com>.  All rights reserved.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt import \
    add_field, \
    add_particle_filter

from yt_p3bh.pop_iii import \
    get_remnant_mass, \
    get_MS_lifetime

from yt_p3bh.fields import \
    _mdot_bh_norm, \
    _bh_to_edd_norm

add_field("bondi_hoyle_accretion_rate",
          function=_mdot_bh_norm,
          units="1/(Msun*yr)",
          sampling_type="cell")
add_field("BH_to_Eddington",
          function=_bh_to_edd_norm,
          units="1/Msun",
          sampling_type="cell")

from yt_p3bh.particles import \
     _pop3_black_hole

add_particle_filter(
    "black_hole", function=_pop3_black_hole, filtered_type="all",
    requires=["creation_time", "metallicity_fraction",
              "particle_type", "particle_mass"])
