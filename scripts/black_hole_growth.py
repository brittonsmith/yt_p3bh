"""
Compute mass and mass growth histories for black hole particles.

Usage: python black_hole_growth.py <simulation dataset filename>

Need to have directory called "black_holes_bh" containing black hole
particles for each dataset.
"""
from collections import defaultdict
import numpy as np
import os
import sys
import yt

from yt.utilities.logger import \
    ytLogger
ytLogger.setLevel(40)

from yt.extensions.p3bh import *

if __name__ == "__main__":
    black_holes = {}

    sim_fn = os.path.abspath(sys.argv[1])
    es = yt.load(sim_fn)
    data_dir = os.path.dirname(sim_fn)
    for i, fn in enumerate(es.data["filename"]):
        if isinstance(fn, bytes):
            fn = fn.decode("utf")
        fn = os.path.basename(fn)
        dsfn = os.path.join(data_dir, "black_holes_bh", "%s.h5" % fn)
        if not os.path.exists(dsfn):
            continue

        ds = yt.load(dsfn)
        if not ds.field_list:
            continue
        pids = ds.r["particle_index"].d.astype(int)
        mbh = ds.arr(np.zeros(pids.size), "Msun")
        mstar = ds.r["particle_mass"].to("Msun")
        mstar[mstar < 1e-9] *= 1e20
        mrem = get_remnant_mass(mstar)
        tms = get_MS_lifetime(mstar)
        print ("%s (z = %f, t = %s): calculating %d accretion rates." %
               (fn, ds.current_redshift, ds.current_time.to("Myr"), mbh.size))
        for i, pid in enumerate(pids):
            if pid not in black_holes:
                age = ds.current_time - ds.r["creation_time"][i]
                if age < tms[i]:
                    continue
                black_holes[pid] = defaultdict(list)
                black_holes[pid]["mass"].append(mrem[i])
            else:
                t_el = ds.current_time - black_holes[pid]["time"][-1]
                m_new = black_holes[pid]["mass"][-1] + \
                  black_holes[pid]["mdot"][-1] * t_el
                black_holes[pid]["mass"].append(m_new)
    
            black_holes[pid]["time"].append(ds.current_time)
            mbh[i] = black_holes[pid]["mass"][-1]

        mdot = bondi_hoyle_accretion_rate(ds.r, mbh)

        for i, pid in enumerate(pids):
            if pid not in black_holes:
                continue
            black_holes[pid]["mdot"].append(mdot[i])

    for pid in black_holes:
        if black_holes[pid]["time"][-1] >= ds.current_time:
            continue
        t_el = ds.current_time - black_holes[pid]["time"][-1]
        m_new = black_holes[pid]["mass"][-1] + \
          black_holes[pid]["mdot"][-1] * t_el
        black_holes[pid]["mass"].append(m_new)
        black_holes[pid]["time"].append(ds.current_time)

    data = {}
    for bh in black_holes:
        for field in black_holes[bh]:
            key = ("p_%d_%s" % (bh, field))
            data[key] = es.arr(black_holes[bh][field])

    ftypes=dict((field, "star") for field in data)
    yt.save_as_dataset(es, os.path.join(data_dir, "black_holes_bh.h5"),
                       data, field_types=ftypes)
