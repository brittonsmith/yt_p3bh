from matplotlib.ticker import FuncFormatter
from matplotlib import pyplot, ticker
import numpy as np
# import seaborn as sns
import yt

from grid_figure import GridFigure
from yt.utilities.cosmology import Cosmology
from yt.visualization.color_maps import \
    yt_colormaps

def _z_from_t(t, pos):
    co = Cosmology(omega_matter=0.266, 
                   omega_lambda=0.734, 
                   hubble_constant=0.71)
    return "%d" % np.round(co.z_from_t(co.quan(t, "Myr")))

def get_distribution(vals, n=101):
    if vals.size == 0:
        return np.zeros(n)
    id_sort = np.argsort(vals)
    dist = vals[id_sort[np.clip(np.round(np.linspace(0, 1, n) *
                                         id_sort.size).astype(int),
                                         0, id_sort.size-1)]]
    return dist

if __name__ == "__main__":
    my_fig = GridFigure(1, 1, figsize=(6, 4.5),
                        left_buffer=0.17, right_buffer=0.04,
                        bottom_buffer=0.13, top_buffer=0.11,
                        vertical_buffer=0)

    # palette = sns.color_palette(palette="colorblind")
    palette = \
      [(0.0, 0.4470588235294118, 0.6980392156862745),
       (0.0, 0.6196078431372549, 0.45098039215686275),
       (0.8352941176470589, 0.3686274509803922, 0.0),
       (0.8, 0.4745098039215686, 0.6549019607843137),
       (0.9411764705882353, 0.8941176470588236, 0.25882352941176473),
       (0.33725490196078434, 0.7058823529411765, 0.9137254901960784)]

    cmap = yt_colormaps["arbre"]

    fn = "black_holes_p2p.h5"
    ds = yt.load(fn)
    co = Cosmology(omega_matter=ds.omega_matter,
                   omega_lambda=ds.omega_lambda,
                   hubble_constant=ds.hubble_constant)

    t_nosf = ds.quan(387.11424659890685, "Myr")

    my_axes = my_fig[0]
    my_axes.set_yscale("log")

    tall = []
    grall = []
    nbh = len(ds.particle_types) - 1
    i = 0

    ptypes = [ptype for ptype in ds.particle_types
              if ptype != "all"]
    
    for i, p in enumerate(ptypes):
        t = ds.data[p, "time"].to("Myr")
        m = ds.data[p, "mass"].to("Msun")
        gr = ds.data[p, "mdot"].to("Msun/yr")
        tall.extend(t[:-1])
        grall.extend(gr)
        my_axes.plot(t[:-1], gr, color=cmap(i/(nbh-1)), alpha=0.7)
        isf = np.where(t > t_nosf)[0][0]
        print (m[0], m[isf], m[-1], (m[isf] / m[0]) - 1, (m[-1] / m[0]) - 1,
               t[-1] - t[0])

    my_axes.xaxis.set_label_text("t [Myr]")
    my_axes.yaxis.set_label_text("$\dot{\\rm M}$ [M$_{\\odot}$ / yr]")
    my_axes.axvline(x=t_nosf, color="red", linestyle="--")
    edd = ds.quan(16.6 * 2.2e-8, "Msun/yr")
    my_axes.axhline(y=edd, color="orange", linewidth=5, alpha=0.5)

    tall = ds.arr(tall)
    grall = ds.arr(grall)
    data = []
    tbins = np.linspace(100, 500, 81)
    bi = np.digitize(tall, tbins) - 1
    for ti in range(tbins.size):
        data.append(get_distribution(grall[bi == ti], 3))
    data = np.rollaxis(np.array(data), 1)
    pdata = data[1]
    imin = np.where(pdata > 0)[0].min()
    imax = np.where(pdata > 0)[0].max() + 1
    my_axes.plot(tbins[imin:imax], pdata[imin:imax],
                 color="black", linewidth=3, alpha=0.75)

    xlim = (125, 500)
    my_axes.set_xlim(*xlim)
    my_axes.xaxis.set_ticks(np.linspace(125, 500, 16), minor=True)
    my_axes.xaxis.set_minor_formatter(ticker.NullFormatter())
    tx = my_axes.twiny()
    tx.set_xlim(*xlim)
    tx.xaxis.tick_top()
    z_ticks_in_t = co.t_from_z(np.array([24, 20, 17,
                                         15, 13, 12,
                                         11, 10])).in_units("Myr")
    tx.xaxis.set_ticks(z_ticks_in_t.d)
    tx.xaxis.set_major_formatter(FuncFormatter(_z_from_t))
    tx.set_xlim(xlim)
    tx.xaxis.set_label_text("z")

    my_axes.yaxis.set_ticks(np.logspace(-20, -6, 8))
    my_axes.yaxis.set_ticks(np.logspace(-20, -7, 14), minor=True)
    my_axes.yaxis.set_minor_formatter(ticker.NullFormatter())
    my_axes.set_ylim(1e-20, 1e-6)
    for my_y in np.logspace(-20, -8, 7):
        my_axes.axhline(y=my_y, color="black", linestyle=":",
                        linewidth=1.0, alpha=0.2, zorder=-100)

    pyplot.savefig("figures/p2p_growth.pdf")
