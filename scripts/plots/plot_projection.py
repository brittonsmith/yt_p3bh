from matplotlib import pyplot, colors
import numpy as np
import os
import yt

from grid_figure import GridFigure
from yt.visualization.color_maps import \
    yt_colormaps

from yt.extensions.p3bh import *

def plot_length_bar(my_axes, ds, length, pwidth, color="white"):
    width = (ds.parameters["right_edge"] -
             ds.parameters["left_edge"])[0]
    lsize = (pwidth * length / width).to("")
    x0 = pwidth / 16
    y0 = pwidth - x0
    my_axes.plot([x0, x0 + lsize], [y0, y0],
                 color=color, linewidth=2, zorder=100)

    lstring = "%d %s" % (length.d, length.units)
    my_axes.text(x0 + lsize / 2, y0 - .5 * x0, "%s" % lstring,
                 horizontalalignment="center",
                 color=color)

def plot_particles(my_axes, pds, bhds, proj_axis,
                   color="black", s=9, cmap=None):
    if not isinstance(color, str):
        mynorm = colors.LogNorm(color.min(), color.max())
        mynorm = colors.Normalize(color.min(), color.max())
    else:
        mynorm = None

    axi = "xyz".find(proj_axis)
    image_x = "xyz"[(axi + 1) % 3]
    image_y = "xyz"[(axi + 2) % 3]
    bwidth = pds.parameters["right_edge"] - pds.parameters["left_edge"]
    px = pshape[0] * (bhds.r["particle_position_%s" % image_x] -
                      pds.parameters["left_edge"][0]) / bwidth[0]
    py = pshape[0] * (bhds.r["particle_position_%s" % image_y] -
                      pds.parameters["left_edge"][1]) / bwidth[1]
    image = my_axes.scatter(px, py, c=color, s=s, 
                            norm=mynorm, cmap=cmap, edgecolor="black", linewidth=0.5)
    return image

if __name__ == "__main__":
    my_fig = GridFigure(2, 1, figsize=(5, 8), square=True,
                        left_buffer=0.02, right_buffer=0.2,
                        vertical_buffer=0.02)
    cbar_buffer = 0.02
    cbar_length = 0.95
    cbar_width = 0.05
    data_dir = "normal_BG1/halo_2170858"
    # data_dir = "normal_BG1/halo_2171203"
    dsfn = "RedshiftOutput0077"
    # data_dir = "Rarepeak_LWB/halo_1367222/"
    # data_dir = "Rarepeak_LWB/halo_1372918"
    # dsfn = "RedshiftOutput0041"

    proj_axis = "x"
    for i, my_axes in enumerate(my_fig):
        if i == 0:
            pds = yt.load(os.path.join(
                data_dir, "projections_%s" % proj_axis, "%s_proj_frb.h5" % dsfn))
            pdata = pds.data["density"].d
            pshape = pdata.shape
            norm = colors.LogNorm(pdata.min(), pdata.max())
            image = my_axes.imshow(pdata, cmap=yt_colormaps["arbre"],
                                   norm=norm, aspect="equal",
                                   interpolation="nearest")
            my_cax = my_fig.add_cax(my_axes, "right", buffer=cbar_buffer,
                                    length=cbar_length, width=cbar_width)
            cbar = pyplot.colorbar(image, orientation="vertical",
                                   cax=my_cax)
            my_cax.yaxis.set_label_text("$\\rho$ [g / cm$^{\\rm 3}$]")
            plot_length_bar(my_axes, pds, pds.quan(2, "kpc"), pdata.shape[0])

            sqw = pshape[0] / 4
            sqcenter = pshape[0] / 2
            x0 = sqcenter - sqw / 2
            my_axes.plot([x0, x0 + sqw], [x0, x0], color="black")
            my_axes.plot([x0, x0 + sqw], [x0 + sqw, x0 + sqw], color="black")
            my_axes.plot([x0 + sqw, x0 + sqw], [x0, x0 + sqw], color="black")
            my_axes.plot([x0, x0], [x0, x0 + sqw], color="black")
            # my_axes.plot([0, x0], [pshape[1], x0 + sqw], color="black", linestyle="--")
            # my_axes.plot([pshape[0], x0+sqw], [pshape[1], x0 + sqw], color="black", linestyle="--")
            # my_axes.set_xlim(0, pshape[0])
            # my_axes.set_ylim(pshape[1], 0)

        if i == 1:
            pds = yt.load(os.path.join(
                data_dir, "projections_0.25_%s" % proj_axis, "%s_proj_frb.h5" % dsfn))
            pdata = pds.data["BH_to_Eddington"].d
            pshape = pdata.shape
            norm = colors.LogNorm(1e-6, 1e-1)
            image = my_axes.imshow(pdata, 
                                   cmap="Greens",
                                   norm=norm, aspect="equal",
                                   interpolation="nearest")
            my_cax = my_fig.add_cax(my_axes, "right", buffer=cbar_buffer,
                                    length=cbar_length, width=cbar_width)
            cbar = pyplot.colorbar(image, orientation="vertical",
                                   cax=my_cax)
            my_cax.yaxis.set_label_text(
                "($\\dot{\\rm m}_{\\rm B-H}$ / $\\dot{\\rm m}_{\\rm Edd}$) / m$_{\\rm bh}$ [1 / M$_{\\odot}$]")
            plot_length_bar(my_axes, pds, pds.quan(500, "pc"), pdata.shape[0],
                            color="black")

            bhds = yt.load(os.path.join(
                data_dir, "projections_0.25_%s" % proj_axis, "%s_sphere.h5" % dsfn))

            ###############
            ### super hack!
            # get the clump/bh distances and use them for color
            disds = yt.load("normal_BG1/halo_2170858/clumps_1e-03/RD0077/halo_2170858_clumps_distances.h5")
            disdict = dict([(my_id, my_dis) for (my_id, my_dis) in
                            zip(disds.data["particle_index"].astype(np.int64),
                                disds.data["clump_distances"])])
            my_distances = bhds.arr([disdict[myid] for myid in
                                     bhds.data["particle_index"].d.astype(np.int64)])
            pimage = plot_particles(
                my_axes, pds, bhds, proj_axis,
                color=my_distances, s=12, cmap=pyplot.cm.inferno_r)
            my_pcax = my_fig.add_cax(my_axes, "bottom", buffer=cbar_buffer,
                                     length=cbar_length, width=cbar_width)
            pcbar = pyplot.colorbar(pimage, orientation="horizontal",
                                   cax=my_pcax, ticks=[100, 400, 700])
            minorticks = pcbar.norm(np.linspace(100, 800, 8))
            my_pcax.xaxis.set_ticks(minorticks, minor=True)
            
            my_pcax.xaxis.set_label_text("bh-clump distance [pc]")
            ###############
            #                           cmap=pyplot.cm.winter)
            #plot_particles(my_axes, pds, bhds, proj_axis, color="black", s=9)
    
        my_axes.xaxis.set_visible(False)
        my_axes.yaxis.set_visible(False)

    ofn = os.path.join(data_dir, "%s_proj_%s.pdf" % (dsfn, proj_axis))
    print ("Saving: %s." % ofn)
    pyplot.savefig(ofn, bbox_inches="tight")
                           
