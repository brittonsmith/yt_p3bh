from collections import defaultdict
import glob
import numpy as np
import os
import yt

from yt.utilities.logger import ytLogger as \
    mylog
llevel = mylog.level

def get_distribution(vals, n=101):
    id_sort = np.argsort(vals)
    dist = vals[id_sort[np.clip(np.round(np.linspace(0, 1, n) *
                                         id_sort.size).astype(int),
                                         0, id_sort.size-1)]]
    return dist

if __name__ == "__main__":
    es = yt.load("rs_normal_bg1.h5")
    fns = es.data["filename"].astype(str)

    value = 1e-3
    halo_dir = "halo_2170858"
    # halo_dir = "halo_2171203"
    data_dir = os.path.join(halo_dir, "clumps_%.0e_inner" % value)
    ofn = os.path.join(halo_dir, "bh_clump_distance_edge_%.0e_inner.h5" % value)

    contained_info = {}

    my_storage = {}
    for sto, (i, fn) in yt.parallel_objects(enumerate(fns), storage=my_storage):
        my_dmin = []
        contained = 0
        cfns = glob.glob(os.path.join(data_dir, os.path.dirname(fn), "*.h5"))
        pbar = yt.get_pbar("%s (z = %f) - Calculating distances" %
                           (os.path.dirname(fn), es.data["redshift"][i]), len(cfns))
        for cfn in cfns:
            clump_dmin = []
            mylog.setLevel(40)
            ds = yt.load(cfn)
            mylog.setLevel(llevel)

            if not hasattr(ds, "tree") or \
               ds.tree.children is None or \
               len(ds.tree.children) == 0:
                del ds
                pbar.update(1)
                continue

            bhps = ds.tree["black_hole", "particle_position"]
            bhids = ds.tree["black_hole", "particle_index"].d.astype(np.int64)
            if bhps.size == 0:
                del ds
                pbar.update(1)
                continue

            contained += sum([c["black_hole", "particle_mass"].size
                              for c in ds.leaves])

            for ib, bhp in enumerate(bhps):
                min_sep = None
                for my_clump in ds.leaves:
                    cbhids = my_clump["black_hole", "particle_index"].d.astype(np.int64)
                    ### check manually instead
                    # if bhids[ib] in cbhids:
                    #     min_sep = 0.0 * cpos.uq
                    #     break                    
                    # else:
                    if True:
                        cpos = np.rollaxis(
                            ds.arr([my_clump["grid", ax] for ax in "xyz"]), 1)
                        dx = np.rollaxis(
                            ds.arr([my_clump["grid", "d%s" % ax] for ax in "xyz"]), 1)
                        diff = np.abs(bhp - cpos)
                        within = diff - 0.5 * dx
                        np.clip(within, 0, np.inf, out=within)
                        d = within.sum(axis=1)
                        if min_sep is None:
                            min_sep = d.min()
                        else:
                            min_sep = min(min_sep, d.min())
                        if (min_sep <= 0):
                            if bhids[ib] not in contained_info:
                                contained_info[bhids[ib]] = defaultdict(list)
                            bhinfo = contained_info[bhids[ib]]
                            bhinfo["mass"].append(
                                ds.tree["black_hole", "particle_mass"][ib])
                            bhinfo["accretion_rate"].append(
                                my_clump["grid", "bondi_hoyle_accretion_rate"][d.argmin()])
                            bhinfo["BH_to_Eddington"].append(
                                my_clump["grid", "BH_to_Eddington"][d.argmin()])
                            bhinfo["redshift"].append(ds.current_redshift)
                            bhinfo["dt"].append(es.data["time"][min(i+1, fns.size-1)] -
                                                es.data["time"][i])
                            print ("In clump %s (cell %d): %d." %
                                   (str(my_clump), d.argmin(), bhids[ib]))
                            break

                clump_dmin.append(min_sep.to("pc"))
            # save dataset with bh-clump distances
            if cfn == os.path.join(data_dir, "%s/RD0077/halo_2170858_clumps.h5"):
                distance_fn = "%s_distances.h5" % cfn[:-3]
                distance_data = \
                  {"particle_index": bhids,
                   "clump_distances": ds.arr(clump_dmin)}
                yt.save_as_dataset(ds, distance_fn, distance_data)
            my_dmin.extend(clump_dmin)
            del ds
            pbar.update(1)
        pbar.finish()

        sto.result_id = i
        if my_dmin:
            my_dmin = yt.YTArray(my_dmin)
            print (my_dmin.min(), my_dmin.max())

            sto.result = [contained, get_distribution(my_dmin)]

    if yt.is_root():
        ids = []
        all_dmin = []
        all_contained = []
        for my_id, my_data in sorted(my_storage.items()):
            if my_data is None:
                continue
            ids.append(my_id)
            all_contained.append(my_data[0])
            all_dmin.append(my_data[1])
        ids = np.array(ids)
        all_contained = np.array(all_contained)
        all_dmin = yt.YTArray(all_dmin)
        data = dict((field[1], es.data[field][ids])
                    for field in es.field_list)
        data["distance_distribution"] = all_dmin
        data["black_holes_in_clumps"] = all_contained
        yt.save_as_dataset(es, ofn, data)

        cont_data = {}
        ftypes = {}
        for bhid in contained_info:
            for field in contained_info[bhid]:
                cont_data[("p_%d" % bhid, field)] = es.arr(contained_info[bhid][field])
                ftypes[("p_%d" % bhid, field)] = "p_%d" % bhid

        yt.save_as_dataset(es, "black_holes_contained.h5", cont_data,
                           field_types=ftypes)
