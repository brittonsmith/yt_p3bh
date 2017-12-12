"""
Calculate some growth statistics for all black hole particles.
This makes a file that is faster to load than the one with the
full growth histories.

Usage: python black_hole_growth_stats.py
"""
from collections import defaultdict
import yt

if __name__ == "__main__":
    ds = yt.load("black_holes_bh.h5")

    data = defaultdict(list)
    
    for field in ds.field_list:
        if "mass" not in field[1]:
            continue

        bhid = int(field[1].split("_")[1])
        data["particle_index"].append(bhid)
        mfield = ds.data[("star", "p_%d_mass" % bhid)]
        gfield = ds.data[("star", "p_%d_mdot" % bhid)]
        data["relative_growth"].append(mfield[-1] / mfield[0])
        data["absolute_growth"].append(mfield[-1] - mfield[0])
        data["max_growth_rate"].append(gfield.max())

    for field in data:
        data[field] = ds.arr(data[field])

    yt.save_as_dataset(ds, "black_hole_growth_stats.h5", data)
