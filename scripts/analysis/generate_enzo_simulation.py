"""
Generate an hdf5 file describing available datasets for this simulation.

Usage: python generate_enzo_simulation.py <parameter file> <output file>
Example: python generate_enzo_simulation.py RD0003/RedshiftOutput0003 rs_normal_bg1.h5
Example: python generate_enzo_simulation.py simulation.enzo rs_normal_bg1.h5
"""
import numpy as np
import sys
import yt

def save_simulation(es, filename=None):
    def to_arr(my_list):
        if hasattr(my_list[0], "units"):
            f = es.arr
        else:
            f = np.array
        return f(my_list)
    fields = ["filename", "time"]
    ex_keys = ["box_size", "initial_time", "final_time"]
    if es.cosmological_simulation:
        fields.append("redshift")
        ex_keys.extend(["initial_redshift", "final_redshift"])
    data = dict((field, to_arr([d[field] for d in es.all_outputs]))
                for field in fields)
    for i in range(data["filename"].size):
        if data["filename"][i].startswith("./"):
            data["filename"][i] = data["filename"][i][2:]
    if filename is None:
        filename = str(es)
        filename = filename[:filename.rfind(".")] + ".h5"
    extra_attrs = dict((field, getattr(es, field))
                       for field in ex_keys)
    yt.save_as_dataset(es, filename, data, extra_attrs=extra_attrs)

if __name__ == "__main__":
    es = yt.simulation(sys.argv[1], "Enzo", find_outputs=True)
    save_simulation(es, filename=sys.argv[1])
