from matplotlib import pyplot
import numpy as np
# import seaborn as sns
import yt

from grid_figure import GridFigure

def get_distribution(vals, n=101):
    id_sort = np.argsort(vals)
    dist = vals[id_sort[np.clip(np.round(np.linspace(0, 1, n) *
                                         id_sort.size).astype(int),
                                         0, id_sort.size-1)]]
    return dist

if __name__ == "__main__":
    my_fig = GridFigure(3, 1, figsize=(5, 7),
                        left_buffer=0.14, right_buffer=0.18,
                        bottom_buffer=0.08, top_buffer=0.02,
                        vertical_buffer=0)

    # palette = sns.color_palette(palette="colorblind")
    palette = \
      [(0.0, 0.4470588235294118, 0.6980392156862745),
       (0.0, 0.6196078431372549, 0.45098039215686275),
       (0.8352941176470589, 0.3686274509803922, 0.0),
       (0.8, 0.4745098039215686, 0.6549019607843137),
       (0.9411764705882353, 0.8941176470588236, 0.25882352941176473),
       (0.33725490196078434, 0.7058823529411765, 0.9137254901960784)]

    labels = ["rare peak\nz = 15",
              "normal\nz = 11.6",
              "void\nz = 9.9"]
    ds_list = [
        yt.load("Rarepeak_LWB/halo_catalogs_nosub/RD0041/RD0041.0.h5"),
        yt.load("normal_BG1/halo_catalogs_nosub/RD0077/RD0077.0.h5"),
        yt.load("void_BG1/halo_catalogs_nosub/RD0094/RD0094.0.h5")
    ]

    m_halo = [ds.r["particle_mass"].to("Msun") for ds in ds_list]
    n_bh = [ds.r["n_black_holes"] for ds in ds_list]

    for i, my_axes in enumerate(my_fig):
        my_axes.set_xscale("log")
        tx = my_axes.twinx()
        for j in range(i+1, i+len(ds_list)+1):
            ip = j % len(ds_list)
            if ip == i:
                s = 0.5
                alpha = 0.9
                label = labels[ip]
            else:
                continue
                s = 5
                alpha = 0.2
                label = ""
            my_axes.scatter(m_halo[ip], n_bh[ip], color=palette[ip],
                            s=s, alpha=alpha, label=label,
                            rasterized=True)

            # plot median
            bins = np.logspace(6, 10, 9)
            bi = np.digitize(m_halo[ip], bins) - 1
            data = []
            mean = []
            f = []
            myds = ds_list[ip]
            for b in range(bins.size - 1):
                sel = bi == b
                if sel.any():
                    data.append(get_distribution(n_bh[ip][bi == b], n=3))
                    f.append((n_bh[ip][bi == b] > 0).sum() /
                             n_bh[ip][bi == b].size)
                    if ip == 0:
                        mean.append(myds.r["n_pop_3"][bi == b].mean())
                else:
                    data.append(np.zeros(3))
                    f.append(0.)
                    if ip == 0:
                        mean.append(0)
            data = np.rollaxis(np.array(data), 1)
            pbins = (bins[:-1] + bins[1:]) / 2
            median = data[1]
            imax = np.where(median > 0)[0].max() + 1
            my_axes.plot(pbins[:imax], median[:imax],
                         color="black", alpha=0.7, linestyle="-")
            # if ip == 0:
            #     my_axes.plot(pbins[:imax], mean[:imax],
            #                  color="red", alpha=0.7, linestyle="-")

            tx.loglog(pbins[:imax], f[:imax], color="red", alpha=0.75)

        tx.set_ylim(8.912509381337459e-05, 1.0)
        tx.yaxis.tick_right()
        tx.yaxis.set_label_text("f (N > 0)")

        my_axes.legend(loc=(-0.05, 0.7), frameon=False,
                       markerscale=5, handletextpad=-0.25,
                       scatteryoffsets = [1.35])
        my_axes.set_xlim(1e6, 4e9)
        my_axes.set_ylim(-1, 80)

        my_axes.yaxis.set_ticks(np.linspace(0, 80, 5))
        my_axes.yaxis.set_ticks(np.linspace(0, 75, 16), minor=True)

        if i < len(ds_list) - 1:
            my_axes.xaxis.set_ticklabels([])
            tl = my_axes.yaxis.get_ticklabels()
            tl[0].set_visible(False)
            tl = tx.yaxis.get_majorticklabels()
            tl[2].set_visible(False)

        if i == len(ds_list) - 1:
            my_axes.xaxis.set_label_text("M$_{\\rm halo}$ [M$_{\\odot}$]")
        my_axes.yaxis.set_label_text("N$_{\\rm black\ holes}$")

        for my_x in np.logspace(7, 9, 3):
            my_axes.axvline(x=my_x, color="black", linestyle=":",
                            linewidth=0.75, alpha=0.2, zorder=-100)

        for my_y in np.linspace(20, 60, 3):
            my_axes.axhline(y=my_y, color="black", linestyle=":",
                            linewidth=0.75, alpha=0.2, zorder=-100)

    pyplot.savefig("figures/m_halo-n_bh.pdf")
