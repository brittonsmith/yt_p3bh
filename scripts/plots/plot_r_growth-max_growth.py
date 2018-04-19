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

    r_growth = [ds.data["relative_growth"] - 1 for ds in ds_list]
    max_growth = [ds.data["max_growth_rate"] for ds in ds_list]

    for i, my_axes in enumerate(my_fig):
        my_axes.set_xscale("log")
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
        my_axes.scatter(r_growth[ip][my_filter], max_growth[ip][my_filter],
                        color=palette[ip],
                        s=s, alpha=alpha, label=label,
                        rasterized=True)

        x_r = np.log10(r_growth[ip][my_filter])
        y_r = np.log10(max_growth[ip][my_filter])
        slope, intercept, r_value, p_value, std_err = \
          stats.linregress(x_r, y_r)
        p_x = np.array([1e-16, 0.2])
        my_axes.plot(p_x, 10**intercept*p_x**slope, color=palette[3],
                     alpha=0.9, label="$\\alpha$ = %.3f" % slope)

        handles, my_labels = my_axes.get_legend_handles_labels()
        my_axes.legend(handles[::-1], my_labels[::-1],
                       loc="upper left", frameon=False,
                       numpoints=3,
                       markerscale=5, handletextpad=0)

        my_axes.set_xlim(1e-16, 0.2)
        my_axes.set_ylim(1e-23, 5e-6)

        my_axes.xaxis.set_ticks(np.logspace(-14, -2, 4))
        my_axes.xaxis.set_ticks(np.logspace(-16, -1, 16), minor=True)
        my_axes.xaxis.set_minor_formatter(ticker.NullFormatter())

        my_axes.yaxis.set_ticks(np.logspace(-22, -6, 5))
        my_axes.yaxis.set_ticks(np.logspace(-22, -6, 17), minor=True)
        my_axes.yaxis.set_minor_formatter(ticker.NullFormatter())

        if i < len(ds_list) - 1:
            my_axes.xaxis.set_ticklabels([])

        if i == len(ds_list) - 1:
            my_axes.xaxis.set_label_text(
                "(M$_{\\rm f}$ / M$_{\\rm i}$ - 1)$_{\\rm max}$")
        my_axes.yaxis.set_label_text(
            "$\\dot{\\rm M}_{\\rm max}$ [M$_{\\odot}$ / yr]")

        for my_x in np.logspace(-16, -1, 4):
            my_axes.axvline(x=my_x, color="black", linestyle=":",
                            linewidth=0.75, alpha=0.2, zorder=-100)

        for my_y in np.logspace(-22, -6, 5):
            my_axes.axhline(y=my_y, color="black", linestyle=":",
                            linewidth=0.75, alpha=0.2, zorder=-100)

    pyplot.savefig("r_growth-max_growth.pdf")
