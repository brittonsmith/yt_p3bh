"""
Do clump finding on the BH_to_Eddington field for each halo in the merger tree.

Usage: python tree_accretion_clumps.py <simulation filename> <merger tree filename> <BH_to_Eddington threshold>
Example: python tree_accretion_clumps.py rs_normal_bg1.h5 halo_2170858/tree_2170858/tree_2170858.h5 1e-4
"""

import gc
import numpy as np
import os
import sys
import yt
yt.enable_parallelism()
import ytree

from yt.analysis_modules.level_sets.api import *
from yt.funcs import \
    ensure_dir
from yt.utilities.physical_constants import G

from yt.utilities.logger import \
    ytLogger
ytLogger.setLevel(20)

from yt.extensions.p3bh import *
from yt.extensions.p3bh.merger_tree_analysis import \
    get_existing_datasets

def _bulk_velocity(clump, **kwargs):
    bv = clump.quantities.bulk_velocity(**kwargs)
    return "Bulk velocity: %s.", bv
add_clump_info("bulk_velocity", _bulk_velocity)

def _total_volume(clump):
    v = clump.data["index", "cell_volume"].sum()
    return "Total volume: %s.", v
add_clump_info("total_volume", _total_volume)

if __name__ == "__main__":
    es = yt.load(sys.argv[1])
    exi, fns = get_existing_datasets(es, "%s")

    pfields = [("black_hole", field) for field in
               ["particle_mass", "particle_type", "particle_index",
                "creation_time", "metallicity_fraction",
                "particle_velocity_x", "particle_velocity_y",
                "particle_velocity_z"]]
    cfield = "BH_to_Eddington"

    a = ytree.load(sys.argv[2])
    cfield_min = float(sys.argv[3])

    data_dir = os.path.join(a.filename.split("/")[-3], "clumps_%.0e" % cfield_min)
    last_snap_idx = int(a[0]["Snap_idx"])
    idx_offset = fns.size - last_snap_idx - 1

    for i, fn in yt.parallel_objects(enumerate(fns)):
        if not os.path.exists(fn):
            continue
        ds = yt.load(fn)

        c_min = ds.quan(cfield_min, "1/Msun")

        output_dir = os.path.join(data_dir, ds.directory)
        ensure_dir(output_dir)
        setup_ds = True


        search_idx = i - idx_offset
        sel = 'tree["tree", "Snap_idx"] == %d' % search_idx
        my_halos = a.select_halos(sel, fields=["Snap_idx"])

        yt.mylog.info("Finding clumps for %d halos." % my_halos.size)
        for halo in my_halos:
            if halo["phantom"] > 0:
                continue
            if halo["n_black_holes"] <= 0:
                continue
            output_filename = os.path.join(output_dir, "halo_%d_clumps.h5" % halo.uid)
            if os.path.exists(output_filename):
                continue

            c_max = halo["BH_to_Eddington"]
            if c_max < c_min:
                continue

            if setup_ds:
                add_fields(ds)
                add_particle_filters(ds)
                setup_ds = False

            center = ds.arr(a.arr([halo["position_%s" % ax] for ax in "xyz"]).to("unitary"))
            radius = ds.quan(halo["virial_radius"].to("kpccm"))
            sphere = ds.sphere(center, radius)

            yt.mylog.info("Extrema: %s %s." % (c_min, c_max))
            step = 1.1 * c_max / c_min

            master_clump = Clump(sphere, cfield)

#             import pdb ; pdb.set_trace()
#             cr = sphere.cut_region(["obj['BH_to_Eddington'] > %f" % c_min])
#             master_clump = Clump(cr, cfield)

            master_clump.clear_clump_info()
            master_clump.add_info_item("total_cells")
            master_clump.add_info_item("total_volume")
            master_clump.add_info_item("cell_mass")
            master_clump.add_info_item("center_of_mass")
            master_clump.add_info_item("bulk_velocity")

            find_clumps(master_clump, c_min, c_max, step)
            master_clump.save_as_dataset(
                filename=output_filename,
                fields=["density", "BH_to_Eddington", "bondi_hoyle_accretion_rate",
                        "velocity_x", "velocity_y", "velocity_z"] + pfields)

            master_clump.data.clear_data()
            del master_clump
            sphere.clear_data()
            del sphere
            gc.collect()

        del ds
        del my_halos
        gc.collect()
