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

    labels = ["rare peak", "normal", "void"]
    ds_list = [
        yt.load("Rarepeak_LWB/top_growers.h5"),
        yt.load("normal_BG1/top_growers.h5"),
        yt.load("void_BG1/top_growers.h5")
    ]

    my_axes = my_fig[0]

    my_axes.set_xscale("log")
    my_axes.set_yscale("log")
    for ip in range(len(ds_list)):
        my_axes.plot([1e-4], [2], color=palette[ip], label=labels[ip])
        myds = ds_list[ip]
        pids = []
        for field in myds.field_list:
            pid = int(field[1].split("_")[1])
            if pid not in pids:
                pids.append(pid)

        for pid in pids:
            t_field = myds.data["p_%d_time" % pid].to("Myr")
            m_field = myds.data["p_%d_mass" % pid]
            g_field = myds.data["p_%d_mdot" % pid].to("Msun/yr")
            g_edd = myds.quan(2.2e-8, "1/yr") * m_field
            f_edd = g_field / g_edd
            my_sedd = f_edd >= 0.25
            x_field = t_field - t_field[0]
            y_field = m_field / m_field[0] - 1
            print (m_field[0])
            if (my_sedd).any():
                print (f_edd[my_sedd], np.where(my_sedd)[0])
                if (my_sedd).sum() > 1:
                    raise NotImplementedError                
                sedd = np.where(my_sedd)[0][0]
                if sedd > 0:
                    my_axes.plot(x_field[0:sedd+1], y_field[0:sedd+1],
                                 color=palette[ip], alpha=0.6)
                if sedd < t_field.size - 1:
                    my_axes.plot(x_field[sedd:sedd+2], y_field[sedd:sedd+2],
                                 color=palette[ip], alpha=0.6,
                                 linestyle="--")
                if sedd < t_field.size - 2:
                    my_axes.plot(x_field[sedd+1:], y_field[sedd+1:],
                                 color=palette[ip], alpha=0.6)
            else:
                my_axes.plot(x_field, y_field,
                             color=palette[ip], alpha=0.6)
            pgr = (m_field - m_field[0]) / (m_field[-1] - m_field[0])
            i_gr = np.where((pgr > 0.9))[0].min()
            my_axes.scatter([(t_field - t_field[0])[i_gr]],
                            (m_field / m_field[0] - 1)[i_gr],
                            s=100, marker="*", alpha=0.9,
                            color=palette[ip])
            my_axes.scatter([(t_field - t_field[0])[-1]],
                            (m_field / m_field[0] - 1)[-1],
                            s=100, marker=".", alpha=0.9,
                            color=palette[ip])
    my_axes.legend(loc="lower right", frameon=False,
                   markerfirst=False,
                   markerscale=5)
    my_axes.set_xlim(1, 1.5e2)
    for my_x in np.logspace(1, 1, 1):
        my_axes.axvline(x=my_x, color="black", linestyle=":",
                        linewidth=1.0, alpha=0.2, zorder=-100)

    # my_axes.yaxis.set_ticks(np.logspace(-24, -8, 5))
    # my_axes.yaxis.set_ticks(np.logspace(-24, -8, 17), minor=True)
    # my_axes.yaxis.set_minor_formatter(ticker.NullFormatter())
    my_axes.set_ylim(1e-3, 0.2)
    for my_y in np.logspace(-2, -1, 2):
        my_axes.axhline(y=my_y, color="black", linestyle=":",
                        linewidth=1.0, alpha=0.2, zorder=-100)

    my_axes.xaxis.set_label_text("black hole age [Myr]")
    my_axes.yaxis.set_label_text(
        "M(t) / M$_{\\rm i}$ - 1")

    pyplot.savefig("figures/top_growers.pdf")
