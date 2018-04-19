from matplotlib import pyplot, ticker
import numpy as np
# import seaborn as sns
from scipy import stats
import yt
from yt.utilities.cosmology import Cosmology

from grid_figure import GridFigure

from yt.extensions.p3bh import *

if __name__ == "__main__":
    my_fig = GridFigure(2, 2, figsize=(7, 5.5),
                        left_buffer=0.12, right_buffer=0.04,
                        bottom_buffer=0.12, top_buffer=0.02,
                        vertical_buffer=0.12, horizontal_buffer=0.16)

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
        yt.load("Rarepeak_LWB/black_hole_growth_stats.h5"),
        yt.load("normal_BG1/black_hole_growth_stats.h5"),
        yt.load("void_BG1/black_hole_growth_stats.h5")
    ]
    age_ds_list = [
        yt.load("Rarepeak_LWB/black_holes_bh/RedshiftOutput0041.h5"),
        yt.load("normal_BG1/black_holes_bh/RedshiftOutput0077.h5"),
        yt.load("void_BG1/black_holes_bh/RedshiftOutput0094.h5")
    ]
    z_f = [15, 11.6, 9.9]
    co = Cosmology(hubble_constant=ds_list[0].hubble_constant,
                   omega_matter=ds_list[0].omega_matter,
                   omega_lambda=ds_list[0].omega_lambda)

    for i, my_axes in enumerate(my_fig):

        if i == 0:
            my_axes.set_yscale("log")
            for j, ds in enumerate(age_ds_list):
                mass = ds.data["particle_mass"].to("Msun")
                mass[mass < 1e-3] *= 1e20
                t_ms = get_MS_lifetime(mass)
                data = co.z_from_t(ds.data["creation_time"] + t_ms)
                hist, bins = np.histogram(
                    data, bins=np.arange(9, 32, 0.5))
                my_axes.step(bins[:-1], hist, color=palette[j],
                             where="post", alpha=0.9, label=labels[j])

                my_axes.set_xlim(31, 9)

                my_axes.xaxis.set_ticks(np.arange(30, 9, -5))
                my_axes.xaxis.set_ticks(np.arange(31, 8, -1), minor=True)

                my_axes.xaxis.set_label_text("z$_{\\rm form}$")

        elif i == 1:
            my_axes.set_xscale("log")
            my_axes.set_yscale("log")
            for j, ds in enumerate(ds_list):
                data = ds.data["relative_growth"] - 1
                data[data == 0] = 1e-19
                hist, bins = np.histogram(
                    data, bins=np.logspace(-20, 0, 41))
                my_axes.step(bins[:-1], hist, color=palette[j],
                             where="post", alpha=0.9, label=labels[j])

            my_axes.set_xlim(1e-20, 1)

            my_axes.xaxis.set_ticks(
                np.concatenate([np.array([1.78e-19]), np.logspace(-15, 0, 4)]))
            my_axes.xaxis.set_ticks(np.logspace(-16, 0, 17), minor=True)
            my_axes.xaxis.set_minor_formatter(ticker.NullFormatter())
            my_fig.figure.canvas.draw()
            labels = [item.get_text() for item in my_axes.get_xticklabels()]
            labels[0] = "0"
            my_axes.set_xticklabels(labels)
            my_axes.xaxis.set_label_text("M$_{\\rm f}$ / M$_{\\rm i}$ - 1")

        # elif i == 2:
        #     my_axes.set_xscale("log")
        #     my_axes.set_yscale("log")
        #     for j, ds in enumerate(ds_list):
        #         data = ds.data["final_mass"]
        #         hist, bins = np.histogram(
        #             data, bins=np.logspace(0, 3, 61))
        #         my_axes.step(bins[:-1], hist, color=palette[j], label=labels[j],
        #                      where="post")

        #         my_axes.set_xlim(5, 300)
        #         my_axes.xaxis.set_label_text("M$_{\\rm bh}$ [M$_{\\odot}$]")
        #         my_axes.legend(loc="lower left", fontsize=10, frameon=False,
        #                        numpoints=3, handletextpad=0.5)

        elif i == 2:
            my_axes.set_xscale("log")
            my_axes.set_yscale("log")
            for j, ds in enumerate(ds_list):
                data = (ds.data["absolute_growth"] /
                        ds.data["growth_time"]).to("Msun/yr")
                data[data != data] = 1e-25
                data[data == 0] = 1e-25
                hist, bins = np.histogram(
                    data, bins=np.logspace(-27, -5, 45))
                my_axes.step(bins[:-1], hist, color=palette[j],
                             where="post", alpha=0.9)

            my_axes.set_xlim(1e-26, 1e-5)
            my_axes.xaxis.set_ticks(
                np.concatenate([np.array([1.78e-25]), np.logspace(-20, -5, 4)]))
            my_axes.xaxis.set_ticks(np.logspace(-23, -6, 18), minor=True)
            my_axes.xaxis.set_minor_formatter(ticker.NullFormatter())
            my_fig.figure.canvas.draw()
            labels = [item.get_text() for item in my_axes.get_xticklabels()]
            labels[0] = "0"
            my_axes.set_xticklabels(labels)
            my_axes.xaxis.set_label_text("ave. growth rate [M$_{\\odot}$/yr]")

        elif i == 3:
            my_axes.set_xscale("log")
            my_axes.set_yscale("log")
            for j, ds in enumerate(ds_list):
                data = ds.data["max_growth_rate_eddington"] / \
                  ds.quan(2.2e-8, "1/yr")
                hist, bins = np.histogram(
                    data, bins=np.logspace(-18, 1, 41))
                my_axes.step(bins[:-1], hist, color=palette[j],
                             where="post", alpha=0.9)

                my_axes.set_xlim(1e-18, 10)
                my_axes.xaxis.set_ticks(np.logspace(-15, 0, 4))
                my_axes.xaxis.set_ticks(np.logspace(-18, 1, 20), minor=True)
                my_axes.xaxis.set_minor_formatter(ticker.NullFormatter())

                my_axes.xaxis.set_label_text("max growth rate / Eddington")

        my_axes.yaxis.set_label_text("N$_{\\rm black\ holes}$")
        my_axes.yaxis.set_ticks(np.kron(np.logspace(0, 4, 5), np.arange(1, 10)), minor=True)
        my_axes.yaxis.set_minor_formatter(ticker.NullFormatter())
        my_axes.set_ylim(0.5, 2000)
        for my_y in np.logspace(0, 3, 4):
            my_axes.axhline(y=my_y, color="black", linestyle=":",
                            linewidth=0.75, alpha=0.2, zorder=-100)

    pyplot.savefig("figures/black_hole_growth.pdf")
