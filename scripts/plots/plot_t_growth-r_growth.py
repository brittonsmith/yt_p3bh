from matplotlib import pyplot, ticker
import numpy as np
# import seaborn as sns
from scipy import stats
import yt

from grid_figure import GridFigure

if __name__ == "__main__":
    my_fig = GridFigure(3, 1, figsize=(4.5, 7),
                        left_buffer=0.22, right_buffer=0.02,
                        bottom_buffer=0.09, top_buffer=0.02,
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
        yt.load("Rarepeak_LWB/black_hole_growth_stats.h5"),
        yt.load("normal_BG1/black_hole_growth_stats.h5"),
        yt.load("void_BG1/black_hole_growth_stats.h5")
    ]

    t_growth = [ds.data["growth_time"].to("Myr") for ds in ds_list]
    r_growth = [ds.data["relative_growth"] - 1 for ds in ds_list]

    for i, my_axes in enumerate(my_fig):
        my_axes.set_xscale("linear")
        my_axes.set_yscale("log")
        ip = i
        if ip == i:
            s = 0.5
            alpha = 0.7
            label = labels[ip]
        else:
            s = 3
            alpha = 0.2
            label = ""
        my_filter = r_growth[ip] > 0
        my_axes.scatter(t_growth[ip][my_filter], r_growth[ip][my_filter],
                        color=palette[ip],
                        s=s, alpha=alpha, label=label,
                        rasterized=True)

        my_axes.legend(loc="lower right", frameon=False,
                       numpoints=3, markerfirst=False,
                       markerscale=5, handletextpad=0)

        my_axes.set_xlim(0, 350)
        my_axes.set_ylim(1e-16, 0.2)

        my_axes.xaxis.set_ticks(np.linspace(0, 300, 4))
        my_axes.xaxis.set_ticks(np.linspace(0, 350, 15), minor=True)
        my_axes.xaxis.set_minor_formatter(ticker.NullFormatter())

        my_axes.yaxis.set_ticks(np.logspace(-15, 0, 4))
        my_axes.yaxis.set_ticks(np.logspace(-16, 0, 17), minor=True)
        my_axes.yaxis.set_minor_formatter(ticker.NullFormatter())

        if i < len(ds_list) - 1:
            my_axes.xaxis.set_ticklabels([])
        else:
            my_axes.xaxis.set_label_text("black hole age [Myr]")

        if i > 0:
            tl = my_axes.yaxis.get_ticklabels()
            tl[-1].set_visible(False)


        my_axes.yaxis.set_label_text("M$_{\\rm f}$ / M$_{\\rm i}$ - 1")

        for my_x in np.linspace(0, 300, 4):
            my_axes.axvline(x=my_x, color="black", linestyle=":",
                            linewidth=0.75, alpha=0.2, zorder=-100)

        for my_y in np.logspace(-15, 0, 4):
            my_axes.axhline(y=my_y, color="black", linestyle=":",
                            linewidth=0.75, alpha=0.2, zorder=-100)

    pyplot.savefig("figures/t_growth-r_growth.pdf")
