"""
Create a ytree merger tree for a specific tree and add some analysis fields
to be used later.  Specifically, we add fields for the number of black holes
and the max value of the BH_to_Eddington field in each halo.  These have been
calculated by the black_hole_halo_catalog.py script, but can also be gotten
from the simulation datasets by setting from_halo_catalogs to False.

Usage: python get_analysis_merger_tree.py <simulation filename>
"""

import numpy as np
import os
import sys
import yt
import ytree

from yt.funcs import \
    ensure_dir

from yt.utilities.logger import \
    ytLogger
ytLogger.setLevel(20)

from yt.extensions.p3bh import *

from yt.extensions.p3bh.merger_tree_analysis import \
    get_existing_datasets, \
    interpolate_phantoms

def get_my_tree():
    """
    Halo with the most black holes from the normal region.
    """

    ### Get the tree with the most black holes.
    ds = yt.load("halo_catalogs/RD0077/RD0077.0.h5")
    a = ytree.load("rockstar_halos/trees/tree_0_0_0.dat")

    i = ds.r["n_black_holes"].argmax()
    hid = ds.r["particle_identifier"][i].d

    t = a[a["Orig_halo_ID"] == hid][0]
    uid = t["uid"]

    data_dir = "halo_%d" % uid
    ensure_dir(data_dir)

    fn = t.save_tree(filename=os.path.join(data_dir, "tree_%d" % uid))
    return fn

def get_my_tree_rare():
    """
    Halo with the most black holes from the normal region.
    """

    ### Get the tree with the most black holes.
    ds = yt.load("halo_catalogs_nosub/RD0041/RD0041.0.h5")
    a = ytree.load("rockstar_halos/trees/tree_0_0_0.dat")

    # halo with highest relative growth
    gsort = ds.r["relative_growth_max"].argsort()[::-1]
    hid = ds.r["particle_identifier"][gsort][0].d

    # # halo with second highest relative growth
    # gsort = ds.r["relative_growth_max"].argsort()[::-1]
    # hid = ds.r["particle_identifier"][gsort][1].d

    t = a[a["Orig_halo_ID"] == hid][0]
    uid = t["uid"]

    data_dir = "halo_%d" % uid
    ensure_dir(data_dir)

    fn = t.save_tree(filename=os.path.join(data_dir, "tree_%d" % uid))
    return fn

if __name__ == "__main__":
    ### Get the target tree here.
    # a_fn = "halo_2170858/tree_2170858/tree_2170858.h5"
    # if not os.path.exists(a_fn):
    #     get_my_tree()
    a_fn = get_my_tree_rare()

    a = ytree.load(a_fn)

    uid = a["uid"]
    data_dir = "halo_%d" % uid
    ensure_dir(data_dir)
    
    es = yt.load(sys.argv[1])
    exi, fns = get_existing_datasets(es, "%s")

    last_snap_idx = int(a[0]["Snap_idx"])
    idx_offset = fns.size - last_snap_idx - 1
    fields = ["n_black_holes", "BH_to_Eddington"]
    field_units = ["", "1/Msun"]
    for f, u in zip(fields, field_units):
        if f not in a.field_list:
            yt.mylog.info("Initializing analysis field: %s." % f)
            a.add_analysis_field(f, u)
            for t in a:
                for halo in t["tree"]:
                    halo[f] = -1

    from_halo_catalogs = True

    for i, fn in enumerate(fns):
        if not os.path.exists(fn):
            continue

        if from_halo_catalogs:
            hfn = os.path.dirname(fn)
            ds = yt.load("halo_catalogs/%s/%s.0.h5" % (hfn, hfn))

        else:
            ds = yt.load(fn)
            setup_ds = True

        search_idx = i - idx_offset
        sel = 'tree["tree", "Snap_idx"] == %d' % search_idx
        my_halos = a.select_halos(sel, fields=["Snap_idx"])

        pbar = yt.get_pbar("Assigning fields for halos", my_halos.size)
        for halo in my_halos:
            pbar.update(1)
            if halo["phantom"] > 0:
                continue
            if halo["BH_to_Eddington"] >= 0:
                continue
            if halo["n_black_holes"] == 0:
                continue

            if from_halo_catalogs:
                my_halo = np.where(halo["Orig_halo_ID"] == ds.r["particle_identifier"])[0][0]
                for field, units in zip(fields, field_units):
                    halo[field] = ds.r[field][my_halo].to(units)

            else:
                if setup_ds:
                    add_fields(ds)
                    add_particle_filters(ds)
                    setup_ds = False

                center = ds.arr(a.arr([halo["position_%s" % ax] for ax in "xyz"]).to("unitary"))
                radius = ds.quan(halo["virial_radius"].to("kpccm"))
                sphere = ds.sphere(center, radius)

                n_bh = sphere["black_hole", "particle_mass"].size
                halo["n_black_holes"] = n_bh
                if n_bh == 0:
                    sphere.clear_data()
                    del sphere
                    continue

                fex = sphere.quantities.extrema("BH_to_Eddington")
                halo["BH_to_Eddington"] = fex[1]

                sphere.clear_data()
                del sphere

        pbar.finish()

        if not from_halo_catalogs and not setup_ds:
            fn = a.save_arbor(filename=a.filename)
            del a
            a = ytree.load(fn)

        del ds
        del my_halos

    # do some interpolation for phantom halos
    interpolate_phantoms(a, fields)
    a.save_arbor(filename=a.filename)
