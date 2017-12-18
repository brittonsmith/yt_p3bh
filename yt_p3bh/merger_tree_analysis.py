import numpy as np
import os
import yt

class NodeArray(np.ndarray):
    def __new__(cls, input_array):
        obj = np.asarray(input_array).view(cls)
        return obj

    def __del_(self):
        if hasattr(self, "cache"):
            self.cache = []

    def __getitem__(self, key):
        if not isinstance(key, str):
            return super(NodeArray, self).__getitem__(key)
        if self.size == 0:
            return None
        if not hasattr(self, "cache"):
            self.cache = {}
        if key not in self.cache:
            d = np.empty(self.size)
            q = self[0][key]
            if hasattr(q, "units"):
                d = yt.YTArray(d, q.units, registry=q.units.registry)
            for i, n in enumerate(self):
                d[i] = n[key]
            self.cache[key] = d
        return self.cache[key]

def get_descendent_line(halo):
    dline = [halo]
    while hasattr(dline[-1], "descendent"):
        dline.append(dline[-1].descendent)
    return NodeArray(dline)

def get_datasets(es, dt_min):
    t = es.data["time"].to("Myr")
    print ("Minimum timestep : %s." % dt_min)
    d = [0]
    for i in range(1, t.size):
        if (t[i] - t[d[-1]]) > dt_min or i == t.size - 1:
            d.append(i)
    d = np.array(d)
    print ("Total datasets: %d." % d.size)
    files = es.data["filename"][d].astype(str)
    return d, np.array(files)

def get_existing_datasets(es, key):
    fns = es.data["filename"].astype(str)
    existing = np.array(
        [i for i in range(fns.size)
         if os.path.exists(key % fns[i])])
    print ("Total datasets: %d." % existing.size)
    return existing, fns[existing]

def interpolate_phantoms(a, fields):
    """
    Do linear interpolation of fields for phantom halos.
    """

    phantoms = []
    for t in a:
        for halo in t["tree"]:
            if halo["phantom"] > 0:
                phantoms.append(halo)
                continue

    pbar = yt.get_pbar("Fixing phantoms", len(phantoms))
    for halo in phantoms:
        pbar.update(1)
        real_anc = halo["prog", "phantom"][1:] == 0
        if not (real_anc).any():
            anc = None
        else:
            anc = halo["prog"][1:][real_anc][0]
        dline = get_descendent_line(halo)
        real_desc = dline["phantom"][1:] == 0
        if not (real_desc).any():
            desc = None
        else:
            desc = dline[1:][real_desc][0]
        del dline

        if anc is None and desc is None:
            continue

        for f in fields:
            if anc is None:
                halo[f] = desc[f]
            elif desc is None:
                halo[f] = anc[f]
            else:
                slope = (desc[f] - anc[f]) / \
                  (desc["scale_factor"] - anc["scale_factor"])
                halo[f] = slope * (halo["scale_factor"] -
                                   anc["scale_factor"]) + \
                  anc[f]
    pbar.finish()
