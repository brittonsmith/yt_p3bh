"""
Keep all fields of input halo catalog, but filter out subhalos.

Usage: python halo_catalog_nosub.py <halo catalog>
"""

import numpy as np
import os
import sys
import yt
yt.enable_parallelism()

from yt.extensions.astro_analysis.halo_analysis.api import *

def replicate_catalog(hds, hc, field_type="halos"):
    hfields = getattr(hds.fields, field_type)
    for field in hfields:
        if field.name in hds.field_list and \
           field.name[1] not in hc.quantities:
            hc.add_quantity(field.name[1], field_type=field_type, prepend=True)

if __name__ == "__main__":
    hds = yt.load(sys.argv[1])

    data_dir = "halo_catalogs_nosub"
    hc = HaloCatalog(
        data_ds=None, halos_ds=hds,
        output_dir=os.path.join(data_dir, os.path.basename(hds.directory)))
    replicate_catalog(hds, hc)
    hc.add_filter("not_subhalo")
    hc.create()
