"""
Gather Pop III star particles and fields necessary for computing Bondi-Hoyle accretion rates.

Usage: python gather_black_holes.py <simulation filename>
Example: python gather_black_holes rs_normal_bg1.h5
"""
from collections import defaultdict
import numpy as np
import os
import sys
import yt
yt.enable_parallelism()
from yt.funcs import \
    ensure_dir

def find_stars(ds, filename, min_level=4):
    fields=["particle_mass", "particle_index", "particle_type",
            "particle_position_x", "particle_position_y", "particle_position_z",
            "particle_velocity_x", "particle_velocity_y", "particle_velocity_z",
            "creation_time", "metallicity_fraction"]
    bfields = ["density", "temperature", "sound_speed",
               "velocity_x", "velocity_y", "velocity_z"]
    data = defaultdict(list)

    Zcr = ds.parameters['PopIIIMetalCriticalFraction']

    ns = 0
    pbar = yt.get_pbar("Reading grids", ds.index.grids.size)
    for i, grid in enumerate(ds.index.grids):
        if ds.index.grid_levels[i][0] >= min_level:
            ct = grid["creation_time"]
            stars = (ct > 0)
            if not stars.any():
                grid.clear_data()
                continue

            Zfr = grid["metallicity_fraction"]
            stars &= (Zfr < Zcr)
            if not stars.any():
                grid.clear_data()
                continue

            pt = grid["particle_type"]
            stars &= ((pt == 1) | (pt == 5))
            if not stars.any():
                grid.clear_data()
                continue

            # mass is multiplied by 1e-20 when main-sequence lifetime is over
            mass = grid["particle_mass"].to("Msun")
            mass[mass < 1e-9] *= 1e20
            stars &= (((mass >= 25) & (mass <= 140)) | (mass >= 260))
            if stars.any():
                ns += stars.sum()
                for field in fields:
                    data[field].append(grid[field][stars])
            grid.clear_data()
        pbar.update(i)
    pbar.finish()

    ndata = {}
    if len(data["particle_mass"]) > 0:
        for field in fields:
            a = ds.arr(np.empty(ns), data[field][0].units)
            ip = 0
            for chunk in data[field]:
                a[ip:ip+chunk.size] = chunk
                ip += chunk.size
            ndata[field] = a

        yt.mylog.info("Getting %d point field values for %s." % (ndata["particle_mass"].size, str(ds)))

        p = ds.arr([ndata["particle_position_%s" % ax] for ax in "xyz"])
        p = np.rollaxis(p, 1)
        bdatal = ds.find_field_values_at_points(bfields, p)
        bdata = dict((field, bdatal[i]) for i, field in enumerate(bfields))
        ndata.update(bdata)

    con_args = ["center", "left_edge", "right_edge"]
    extra_attrs = dict((field, getattr(ds, "domain_%s" % field))
                       for field in con_args)
    extra_attrs["con_args"] = con_args
    extra_attrs["data_type"] = "yt_data_container"
    extra_attrs["container_type"] = "region"
    extra_attrs["dimensionality"] = 3
    ftypes = dict((field, "star") for field in fields + bfields)
    yt.save_as_dataset(ds, filename, ndata,
                       field_types=ftypes, extra_attrs=extra_attrs)

if __name__ == "__main__":
    data_dir = "black_holes_bh"
    ensure_dir(data_dir)

    es = yt.load(sys.argv[1])
    fns = es.data["filename"]

    for fn in yt.parallel_objects(fns):
        if not os.path.exists(fn):
            continue
        ds = yt.load(fn)

        star_file = os.path.join(data_dir, "%s.h5" % ds)
        if not os.path.exists(star_file):
            find_stars(ds, star_file)
