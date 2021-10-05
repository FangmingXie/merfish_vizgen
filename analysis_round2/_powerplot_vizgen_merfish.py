
import numpy as np
import datashader as ds
import colorcet

import sys
sys.path.insert(0, '/cndd2/fangming/projects/SingleCellRoutines')
from __init__plots import *
import utils
import powerplot


def fig_plot_gene_insitu_routine_unified_colorbar(
    thedatagmat, samples, x, y, hue, 
    samples_annot=dict(),
    scale_paras=dict(pxl_scale=20),
    cmap=colorcet.cm.blues,
    vmin=0, 
    vmax=0,
    nx=3,
    ny=3,
    figsize=(9*3,6*3),
    output='',
    close=False,
    ):
    """
    """

    fig, axs = plt.subplots(ny, nx, figsize=figsize)
    if nx == 1 and ny == 1:
        flataxs = [axs]
    else:
        flataxs = axs.flat

    for i, (ax, sample) in enumerate(zip(flataxs, samples)):
        data = thedatagmat[thedatagmat['sample']==sample]
        if len(samples_annot) > 0:
            title = samples_annot[sample]
        else:
            title = sample

        if i == 0:
            configs = dict(
                arrows=True,
                scalebar=True,
                vmin=vmin,
                vmax=vmax,
                )
        else:
            configs = dict(
                arrows=False,
                scalebar=True,
                vmin=vmin,
                vmax=vmax,
                )
        powerplot.plot_gene_insitu_routine(ax, data, x, y, hue, scale_paras, cmap, title, **configs)

    # colorbar
    cax = fig.add_axes([0.25, 0.1, 0.1, 0.01])
    powerplot.add_colorbar_unified_colorbar(fig, cax, vmin=vmin, vmax=vmax, cmap=cmap, orientation='horizontal')

    fig.subplots_adjust(wspace=-0.2)
    fig.suptitle(hue, y=0.93)
    
    if output:
        utils.savefig(fig, output)
    if close:
        plt.close()
    return fig

def fig_plot_gene_insitu_routine(
    thedatagmat, samples, x, y, hue, 
    samples_annot=dict(),
    scale_paras=dict(pxl_scale=20),
    cmap=colorcet.cm.blues,
    vmaxp=99,
    nx=3,
    ny=3,
    figsize=(9*3,6*3),
    output='',
    close=False,
    ):
    """
    """
    fig, axs = plt.subplots(ny, nx, figsize=figsize)
    if nx == 1 and ny == 1:
        flataxs = [axs]
    else:
        flataxs = axs.flat

    for i, (ax, sample) in enumerate(zip(flataxs, samples)):
        data = thedatagmat[thedatagmat['sample']==sample]
        if len(samples_annot) > 0:
            title = samples_annot[sample]
        else:
            title = sample

        if i == 0:
            configs = dict(
                arrows=True,
                scalebar=True,
                vmaxp=vmaxp,
                )
        else:
            configs = dict(
                arrows=False,
                scalebar=True,
                vmaxp=vmaxp,
                )
        powerplot.plot_gene_insitu_routine(ax, data, x, y, hue, scale_paras, cmap, title, **configs)

    # colorbar
    cax = fig.add_axes([0.25, 0.1, 0.1, 0.01])
    powerplot.add_colorbar(fig, cax, cmap=cmap, orientation='horizontal')

    fig.subplots_adjust(wspace=-0.2)
    fig.suptitle(hue, y=0.93)
    
    if output:
        utils.savefig(fig, output)
    if close:
        plt.close()
    return fig

def fig_plot_gene_umap_routine(
    thedatagmat, x, y, hue, 
    scale_paras=dict(npxlx=300),
    cmap=colorcet.cm.blues,
    vmaxp=99,
    figsize=(6,6),
    output='',
    close=False,
    ):
    """
    """
    title = hue
    configs = dict(
        arrows=True,
        )

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    powerplot.plot_gene_umap_routine(
        ax, thedatagmat, x, y, hue, scale_paras, cmap, title, vmaxp=vmaxp, **configs)
    # colorbar
    cax = fig.add_axes([0.5, 0, 0.2, 0.02])
    powerplot.add_colorbar(fig, cax, cmap=cmap, orientation='horizontal')
    
    if output:
        utils.savefig(fig, output)
    if close:
        plt.close()
    else:
        plt.show()

def fig_plot_cluster_insitu_routine(
    thedatagmat, samples, x, y, hue,
    clstcolors_obj,
    samples_annot=dict(),
    scale_paras=dict(npxlx=300),
    nx=3,
    ny=3,
    figsize=(9*3,6*3),
    cbar_fontsize=15,
    showticks=False,
    roi=[],
    suptitle='colored by cluster',
    close=False,
    output='', 
    ): 
    """
    roi - a list of 4 numbers, [xmin, xmax, ymin, ymax], empty [] means plotting everything
    
    """
    # agg data for each sample
    fig, axs = plt.subplots(ny, nx, figsize=figsize)
    for i, (ax, sample) in enumerate(zip(axs.flat, samples)):
        if len(samples_annot) > 0:
            title = samples_annot[sample]
        else:
            title = sample

        if i == 0:
            arrows=True
            scalebar=True
        else:
            arrows=False
            scalebar=True
        
        # agg data for each sample
        data = thedatagmat[thedatagmat['sample']==sample]
        aggdata, ps, cluster_labels = powerplot.agg_count_cat(data, x, y, hue, scale_paras, clip_max=0, reduce=True)
        assert np.all(clstcolors_obj.labels == cluster_labels)
        
        powerplot.imshow_routine(
            ax, aggdata+0.5, cmap=clstcolors_obj.cmap, norm=clstcolors_obj.norm)
        if arrows:
            powerplot.add_arrows(ax, 'in situ')
        if scalebar:
            bar_length = 1000 # (micron)
            powerplot.add_scalebar(ax, ps.npxlx-ps.len2pixel(bar_length), ps.npxlx, '1 mm')
        if showticks:
            ax.axis('on')
        if len(roi) > 0 and len(roi) == 4:
            ax.set_xlim(roi[:2])
            ax.set_ylim(roi[2:4])
            
        ax.set_title(title)
        
    clstcolors_obj.add_colorbar(fig, fontsize=cbar_fontsize)
    fig.suptitle(suptitle, y=0.93)

    fig.subplots_adjust(wspace=-0.2)
    if output:
        utils.savefig(fig, output)
    if close:
        plt.close()
    else:
        plt.show(
        )
    return ps
        
def fig_plot_cluster_umap_routine(
    thedatagmat, x, y, hue,
    clstcolors_obj,
    scale_paras=dict(npxlx=300),
    title='colored by cluster',
    figsize=(6,6),
    close=False,
    output="",
    ):
    # plot all clusters UMAP
    aggthedatagmat_umap, _, _ = powerplot.agg_count_cat(thedatagmat, x, y, hue, scale_paras, clip_max=0, reduce=True)

    # plot 
    fig, ax = plt.subplots(figsize=figsize)
    powerplot.imshow_routine(
            ax, aggthedatagmat_umap+0.5, cmap=clstcolors_obj.cmap, norm=clstcolors_obj.norm)
    powerplot.add_arrows(ax, 'UMAP')
    clstcolors_obj.add_colorbar(fig)
    ax.set_title(title)
    if output:
        utils.savefig(fig, output)
    if close:
        plt.close()
    else:
        plt.show()
    plt.show()
