"""
Create halo catalog with various black hole stats.

Usage: python black_hole_halo_catalog.py <simulation dataset>

Must run black_hole_growth_stats.py first.
"""

import numpy as np
import os
import sys
import yt
yt.enable_parallelism()

from yt.extensions.astro_analysis.halo_analysis.api import *

from yt.extensions.p3bh import *

def count_black_holes(halo):
    dds = halo.halo_catalog.data_ds
    hds = halo.halo_catalog.halos_ds

    sp = halo.data_object
    return sp["black_hole", "particle_index"].size
add_quantity("n_black_holes", count_black_holes)

def black_hole_stats(halo, bhds):
    dds = halo.halo_catalog.data_ds
    hds = halo.halo_catalog.halos_ds

    fields = ["relative_growth", "absolute_growth", "max_growth_rate"]
    units = ["", "Msun", "Msun/yr"]
    calcs = ["max", "mean"]

    sp = halo.data_object
    h_bhids = sp["black_hole", "particle_index"].astype(np.int64)

    bhids = bhds.data["particle_index"].astype(np.int64)
    common_ids = np.intersect1d(h_bhids, bhids)
    dsids = np.array([np.where(bhids == bhid)[0][0] for bhid in common_ids])
    for field, unit in zip(fields, units):
        for calc in calcs:
            quant = "%s_%s" % (field, calc)
            if quant not in halo.halo_catalog.quantities:
                halo.halo_catalog.quantities.append(quant)
            if common_ids.size == 0:
                halo.quantities[quant] = dds.quan(0, unit)
            else:
                halo.quantities[quant] = getattr(bhds.data[field][dsids], calc)()
add_callback("black_hole_stats", black_hole_stats)

if __name__ == "__main__":
    dds = yt.load(sys.argv[1])
    add_particle_filters(dds)
    hds = yt.load("rockstar_halos/halos_%s.0.bin" % dds.directory)
    bhds = yt.load("black_hole_growth_stats.h5")

    ad = hds.all_data()
    # cr = ad.cut_region(["obj['particle_mass'].to('Msun') > 1e8"])

    data_dir = "halo_catalogs"
    hc = HaloCatalog(data_ds=dds, halos_ds=hds, # data_source=cr,
                     output_dir=os.path.join(data_dir, dds.directory))
    hc.add_callback("sphere")
    hc.add_quantity("n_black_holes")
    hc.add_callback("black_hole_stats", bhds)
    hc.create()
