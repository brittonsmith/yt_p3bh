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
                        left_buffer=0.17, right_buffer=0.02,
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

    labels = ["rare peak", "normal", "void"]
    ds_list = [
        yt.load("Rarepeak_LWB/black_hole_halos.h5"),
        yt.load("normal_BG1/black_hole_halos.h5"),
        yt.load("void_BG1/black_hole_halos.h5")
    ]

    m_halo = [ds.data["particle_mass"].to("Msun") for ds in ds_list]
    r_growth = [(ds.data["relative_growth"] - 1) /
                ds.data["growth_time"].to("yr") for ds in ds_list]

    my_axes = my_fig[0]

    my_axes.set_xscale("log")
    my_axes.set_yscale("log")
    for ip in range(len(ds_list)):
        my_filter = r_growth[ip] > 0
        mass_data = m_halo[ip][my_filter]
        growth_data = r_growth[ip][my_filter]
        bins = np.logspace(6, 10, 17)
        bi = np.digitize(mass_data, bins) - 1
        data = []
        for b in range(bins.size - 1):
            data.append(get_distribution(growth_data[bi == b], n=3))
        pbins = (bins[:-1] + bins[1:]) / 2
        data = np.rollaxis(np.array(data), 1)
        imin = np.where(data[1] > 0)[0].min()
        imax = np.where(data[1] > 0)[0].max()
        my_axes.scatter(m_halo[ip], r_growth[ip], color=palette[ip],
                        alpha=0.15, s=5, marker=".")
        my_axes.plot(pbins[imin:imax+1], data[1][imin:imax+1],
                     color=palette[ip], 
                     alpha=0.7, label=labels[ip])

    my_axes.legend(loc="lower right", frameon=False,
                   markerfirst=False,
                   markerscale=5, handletextpad=0.5)
    my_axes.set_xlim(1e6, 4e9)

    my_axes.yaxis.set_ticks(np.logspace(-24, -8, 5))
    my_axes.yaxis.set_ticks(np.logspace(-24, -8, 17), minor=True)
    my_axes.yaxis.set_minor_formatter(ticker.NullFormatter())
    my_axes.set_ylim(1e-24, 1e-8)

    my_axes.xaxis.set_label_text("M$_{\\rm halo}$ [M$_{\\odot}$]")
    # my_axes.yaxis.set_label_text(
    #     "$\dot{\\rm M}_{\\rm ave}$ / M$_{\\rm i}$ [1 / yr]")
    # my_axes.yaxis.set_label_text(
    #     "(M$_{\\rm f}$ - M$_{\\rm i}$) / M$_{\\rm i}$ / age [1 / yr]")
    my_axes.yaxis.set_label_text("specific growth rate [1/yr]")

    for my_x in np.logspace(7, 9, 3):
        my_axes.axvline(x=my_x, color="black", linestyle=":",
                        linewidth=0.75, alpha=0.2, zorder=-100)

    for my_y in np.logspace(-24, -8, 5):
        my_axes.axhline(y=my_y, color="black", linestyle=":",
                        linewidth=0.75, alpha=0.2, zorder=-100)

    pyplot.savefig("figures/m_halo-r_growth_median.pdf")
