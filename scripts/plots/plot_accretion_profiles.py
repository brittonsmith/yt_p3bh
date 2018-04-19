import copy
import glob
import h5py
from matplotlib import pyplot, cm
from matplotlib.ticker import \
     FuncFormatter, FixedLocator, StrMethodFormatter, NullFormatter
import matplotlib as mpl
import numpy as np
import os
import yt
from yt.visualization.color_maps import *
from yt.units.yt_array import YTQuantity, YTArray
from yt.utilities.cosmology import Cosmology

from yt.utilities.logger import \
    ytLogger
ytLogger.setLevel(40)

from grid_figure import GridFigure

mpl.rcParams['axes.unicode_minus'] = False

def _z_from_t(t, pos):
    co = Cosmology(omega_matter=0.266, 
                   omega_lambda=0.734, 
                   hubble_constant=0.71)
    return "%d" % np.round(co.z_from_t(co.quan(t, "Myr")))

def plot_profile(my_fig, my_axes, data_dir, title):
    my_files = glob.glob(os.path.join(data_dir, "profiles/*.h5"))
    my_files = my_files[5:]
    ds = yt.load(my_files[0])
    z_i = ds.current_redshift
    ds = yt.load(my_files[-1])
    z_f = ds.current_redshift
    del ds

    x_min = None
    x_max = None

    pbar = yt.get_pbar("Plotting profiles", len(my_files))
    for my_fn in my_files:
        ds = yt.load(my_fn)

        my_x = ds.profile.x
        v_tot = ds.profile["cell_volume"].sum()
        my_y = ds.profile["cell_volume"]
        my_y = my_y[::-1].cumsum()[::-1]
        my_y /= v_tot
        if x_min is None:
            x_min = ds.profile.x_bins.min()
        x_min = min(x_min, ds.profile.x_bins.min())
        if x_max is None:
            x_max = ds.profile.x_bins.max()
        x_max = max(x_max, ds.profile.x_bins.max())

        my_z = ds.current_redshift
        c_i = (z_i - my_z) / (z_i - z_f)
        my_axes.loglog(my_x, my_y, color=cmap(c_i), alpha=0.8)
        pbar.update(1)
    pbar.finish()

    my_pos = np.array(my_axes.get_position())
    my_cax = my_fig.add_cax(my_axes, "right", buffer=0.01)
    norm = mpl.colors.Normalize(vmin=z_i, vmax=z_f)
    cb1 = mpl.colorbar.ColorbarBase(my_cax, cmap=cmap,
                                    norm=norm,
                                    orientation='vertical')
    cb1.set_label("z", fontsize=fontsize)
    cb1.solids.set_edgecolor("face")
    if z_i - z_f > 8:
        step = -2
    else:
        step = -1
    cb1.set_ticks(np.arange(int(z_i), int(z_f), step))

    my_axes.set_xlim(5e-7, 10)
    tx = 10**(np.log10(5e-7 * 10) / 2)
    ty = 1e-2
    my_axes.text(tx, ty, title, horizontalalignment="center")
    for xmaj in np.logspace(-6, 1, 8):
        my_axes.axvline(x=xmaj, color="black", alpha=0.2, linestyle=":",
                        zorder=1)
    my_axes.xaxis.set_ticks(np.logspace(-6, 0, 4))
    my_axes.xaxis.set_ticks(np.logspace(-6, 0, 7), minor=True)
    my_axes.xaxis.set_minor_formatter(NullFormatter())

    my_axes.set_ylim(1e-16, 1)
    my_axes.yaxis.set_ticks(np.logspace(-16, 0, 5))
    my_axes.yaxis.set_ticks(np.logspace(-16, -1, 16), minor=True)
    my_axes.yaxis.set_minor_formatter(NullFormatter())
    for ymaj in np.logspace(-16, -1, 16):
        my_axes.axhline(y=ymaj, color="black", alpha=0.2, linestyle=":",
                        zorder=1)

    return my_cax

if __name__ == "__main__":
    fontsize = 16
    cmap = yt_colormaps["algae"]
    my_fig = GridFigure(3, 1, figsize=(6, 9),
                        left_buffer=0.18, right_buffer=0.16,
                        bottom_buffer=0.08, top_buffer=0.02,
                        vertical_buffer=0)

    data_dirs = ["Rarepeak_LWB", "normal_BG1", "void_BG1"]
    titles = ["rare peak", "normal", "void"]

for i, my_axes in enumerate(my_fig):
    my_cax = plot_profile(my_fig, my_axes, data_dirs[i], titles[i])

    if i == len(my_fig) - 1:
        my_axes.xaxis.set_label_text(
            "Bondi-Hoyle to Eddington ratio [1 / M$_{\\odot}$]", fontsize=fontsize)
    else:
        my_axes.xaxis.set_ticklabels([])
        tl = my_axes.yaxis.get_ticklabels()
        tl[0].set_visible(False)

    my_axes.yaxis.set_label_text("f (>V)", fontsize=fontsize)

    tick_labels = my_axes.xaxis.get_ticklabels() + \
      my_axes.yaxis.get_ticklabels() + my_cax.yaxis.get_ticklabels()
    for tl in tick_labels:
        tl.set_size(fontsize)

pyplot.savefig("figures/accretion_profiles.pdf")
