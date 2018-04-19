import copy
import glob
import h5py
from matplotlib import pyplot, cm
from matplotlib.ticker import FuncFormatter
import matplotlib as mpl
import numpy as np
import os
import yt
from yt.visualization.color_maps import yt_colormaps
from yt.utilities.cosmology import Cosmology

mpl.rcParams['axes.unicode_minus'] = False

from grid_figure import GridFigure

def _z_from_t(t, pos):
    co = Cosmology(omega_matter=0.266, 
                   omega_lambda=0.734, 
                   hubble_constant=0.71)
    return "%d" % np.round(co.z_from_t(co.quan(t, "Myr")))

def plot_distribution(my_axes, x, y, step=1, show_dist=True,
                      **kwargs):
    ihalf = y.shape[1] // 2
    alpha = kwargs.get("alpha", 0.9)
    color = kwargs.get("color", "red")
    kwargs.update({"alpha": alpha, "color": color})
    my_axes.plot(x, y[:, ihalf], **kwargs)
    if not show_dist:
        return
    my_alpha = step / ihalf
    for i in range(step, ihalf+1, step):
        my_axes.fill_between(x, y1=y[:, ihalf-i], y2=y[:, ihalf+i],
                             color=color, alpha=my_alpha, linewidth=0)
def plot_legend(my_axes, items, step=1):
    my_pos = np.array(my_axes.get_position())
    lfontsize = 8
    lax = pyplot.axes((my_pos[1][0] - 0.11, my_pos[0][1] + 0.03,
                       my_fig.panel_width * 0.1,
                       my_fig.panel_height * 0.2))

    ihalf = 50
    my_alpha = step / ihalf
    ones = np.ones(2)
    for i, item in enumerate(items):
        color, label = item
        lax.plot([i-0.4, i+0.4], [50, 50], color=color, alpha=0.9)
        for offset in range(step, ihalf+1, step):
            lax.fill_between([i-0.4, i+0.4],
                             y1=ihalf-offset,
                             y2=ihalf+offset,
                             alpha=my_alpha, color=color,
                             linewidth=0)
    lax.set_ylim(0, 100)
    lax.yaxis.tick_left()
    lax.yaxis.set_ticks(np.arange(0, 101, 25))
    lax.yaxis.set_label_text("%", fontsize=lfontsize)
    lax.yaxis.set_label_position("left")
    lax.yaxis.labelpad = 0
    lax.xaxis.tick_bottom()
    lax.xaxis.set_ticks(np.arange(1)-0.2)
    #lax.set_xticklabels(lax.xaxis.get_majorticklabels(), rotation=315)
    lax.xaxis.set_ticklabels([r""])
    for tick in lax.xaxis.get_ticklines():
            tick.set_visible(False)
    for ticklabel in lax.xaxis.get_majorticklabels():
        ticklabel.set_visible(True)
        ticklabel.set_horizontalalignment("left")
        ticklabel.set_verticalalignment("top")
        ticklabel.set_fontsize(10)
    for ticklabel in lax.yaxis.get_majorticklabels():
        ticklabel.set_fontsize(lfontsize)

if __name__ == "__main__":
    halo_dir = "halo_2170858"
    # halo_dir = "halo_2171203"
    ds = yt.load(os.path.join(halo_dir, "bh_clump_distance_edge_1e-03.h5"))
    ds_i = yt.load(os.path.join(halo_dir, "bh_clump_distance_edge_1e-03_inner_0.25.h5"))

    fontsize = 12
    my_fig = GridFigure(1, 1, top_buffer=0.12, bottom_buffer=0.13,
                        left_buffer=0.14, right_buffer=0.04,
                        horizontal_buffer=0.05, vertical_buffer=0.1,
                        figsize=(6, 4.5))
    my_axes = my_fig[0]
    my_axes.set_yscale("log")
    plot_distribution(my_axes, ds.data["time"].to("Myr"),
                      ds.data["distance_distribution"].to("pc"), step=10)
    plot_distribution(my_axes, ds_i.data["time"].to("Myr"),
                      ds_i.data["distance_distribution"].to("pc"), step=3,
                      show_dist=False, color="black", linestyle="--",
                      alpha=0.5)
    plot_legend(my_axes, [("red", "")], step=10)

    xlim = (200, 400)
    my_axes.xaxis.set_ticks(np.arange(200, 401, 50))
    my_axes.xaxis.set_ticks(np.arange(180, 401, 10), minor=True)
    my_axes.xaxis.set_label_text("t [Myr]")
    my_axes.set_xlim(*xlim)

    co = Cosmology(omega_matter=0.266, 
                   omega_lambda=0.734, 
                   hubble_constant=0.71)

    tx = my_axes.twiny()
    tx.set_xlim(*xlim)
    tx.xaxis.tick_top()
    z_ticks_in_t = co.t_from_z(np.array([19, 18, 17,
                                         16, 15,
                                         14, 13, 12])).in_units("Myr")
    tx.xaxis.set_ticks(z_ticks_in_t.d)
    tx.xaxis.set_major_formatter(FuncFormatter(_z_from_t))
    tx.set_xlim(xlim)
    tx.xaxis.set_label_text("z", fontsize=fontsize)

    for y in np.logspace(0, 3, 4):
        my_axes.axhline(y=y, color="black", linestyle=":", alpha=0.1)

    my_axes.yaxis.set_label_text("distance to clump [pc]")
    my_axes.set_ylim(.1, 1e4)
    ty = my_axes.twinx()
    ty.set_yscale("log")
    ty.set_ylim(my_axes.get_ylim())
    ty.set_yticklabels([])
    ty.tick_params(direction="in", axis="y", which="both")

    ofn = os.path.join(halo_dir, "distance_distribution_edge.pdf")
    print ("Saving to %s." % ofn)
    pyplot.savefig(ofn)
