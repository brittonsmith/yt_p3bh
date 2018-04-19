from matplotlib import pyplot, ticker
import numpy as np
# import seaborn as sns
import yt

from grid_figure import GridFigure

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
                        left_buffer=0.15, right_buffer=0.04,
                        bottom_buffer=0.14, top_buffer=0.04,
                        vertical_buffer=0)

    # palette = sns.color_palette(palette="colorblind")
    palette = \
      [(0.0, 0.4470588235294118, 0.6980392156862745),
       (0.0, 0.6196078431372549, 0.45098039215686275),
       (0.8352941176470589, 0.3686274509803922, 0.0),
       (0.8, 0.4745098039215686, 0.6549019607843137),
       (0.9411764705882353, 0.8941176470588236, 0.25882352941176473),
       (0.33725490196078434, 0.7058823529411765, 0.9137254901960784)]

    labels = ["rare peak", "normal", "void", "p2p"]
    ds_list = [
        yt.load("Rarepeak_LWB/bh_growths.h5"),
        yt.load("normal_BG1/bh_growths.h5"),
        yt.load("void_BG1/bh_growths.h5"),
        yt.load("pop2prime/bh_growths.h5"),
    ]
    data = []

    my_axes = my_fig[0]

    my_axes.set_xscale("log")
    my_axes.set_yscale("log")

    for ip, ds in enumerate(ds_list):
        fdata = ds.data["mdot"] / ds.data["mass"] / \
          ds.quan(2.2e-8, "1/yr")
        fdata.convert_to_units("")
        data.append(fdata)
    dmin = np.floor(np.log10(min([d.min() for d in data]))) - 1
    dmax = np.ceil(np.log10(max([d.max() for d in data]))) + 1
    bpd = 4
    bins = np.logspace(dmin, dmax, (dmax - dmin) * bpd + 1)

    for ip, ds in enumerate(ds_list):
        hist, hbins = np.histogram(data[ip], bins=bins)
        hist = hist / hist.sum()
        hist = hist[::-1].cumsum()[::-1]
        #hist[hist == 0] = 1e-10
        chist = hist[::-1].cumsum()[::-1]
        my_axes.step(bins[:-1], hist, color=palette[ip],
                     where="post", alpha=0.9, label=labels[ip])
        # my_axes.step(bins[:-1], chist, color=palette[ip],
        #              where="post", alpha=0.9, label=labels[ip])

    my_axes.legend(loc=(0.3, 0.0), frameon=False,
                   markerfirst=False,
                   markerscale=5)
    my_axes.set_xlim(1e-19, 1e1)
    my_axes.xaxis.set_ticks(np.logspace(-16, 0, 5))
    my_axes.xaxis.set_ticks(np.logspace(-19, 0, 20), minor=True)
    my_axes.xaxis.set_minor_formatter(ticker.NullFormatter())
    for my_x in np.logspace(-16, 0, 5):
        my_axes.axvline(x=my_x, color="black", linestyle=":",
                        linewidth=1.0, alpha=0.2, zorder=-100)

    # my_axes.yaxis.set_ticks(np.logspace(-24, -8, 5))
    # my_axes.yaxis.set_ticks(np.logspace(-24, -8, 17), minor=True)
    # my_axes.yaxis.set_minor_formatter(ticker.NullFormatter())
    my_axes.set_ylim(5e-6, 1.0)
    for my_y in np.logspace(-5, -2, 4):
        my_axes.axhline(y=my_y, color="black", linestyle=":",
                        linewidth=1.0, alpha=0.2, zorder=-100)

    my_axes.xaxis.set_label_text("$\dot{\\rm M}_{\\rm B-H}\ /\ \dot{\\rm M}_{\\rm Edd}$")
    my_axes.yaxis.set_label_text("f")

    pyplot.savefig("figures/growth_cdf.pdf")
