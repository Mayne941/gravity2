from matplotlib.colors import LinearSegmentedColormap as LSC
import matplotlib.pyplot as plt
import numpy as np

def get_blue_cmap():
    return LSC('MyBlues',
               {
                    'red':  ((0., 0., 0.), (1., 1., 1.)),
                    'green':  ((0., 0., 0.), (1., 1., 1.)),
                    'blue':  ((0., 1., 1.), (1., 1., 1.))
                },
                1024
    )

def get_red_cmap():
    return LSC('MyReds',
               {
                    'red':  ((0., 1., 1.), (1., 1., 1.)),
                    'green':  ((0., 0., 0.), (1., 1., 1.)),
                    'blue':  ((0., 0., 0.), (1., 1., 1.))
                },
                1024
    )

def get_purple_cmap():
    return LSC('MyPurples',
               {
                    'red':  ((0., 0.5, 0.5), (1., 1., 1.)),
                    'green':  ((0., 0., 0.), (1., 1., 1.)),
                    'blue':  ((0., 1., 1.), (1., 1., 1.))
                },
                1024
    )

def get_hmap_params(n_viruses, n_pphmms=99, is_square=True,virus_dendrogram=None,do_logscale=False):
    hmap_params = {}
    '''General'''
    hmap_params['Outer_margin'] = 0.5
    hmap_params['FontSize'] = 6
    hmap_params['dpi']=100 # TODO How high can we get this without massive delays?
    if not is_square:
        if n_pphmms > 3000:
            hmap_params['Heatmap_width'] = float(60)
            hmap_params['Heatmap_height'] = float(15)
        elif n_pphmms > 500 and n_pphmms <= 3000:
            hmap_params['Heatmap_width'] = float(30)
            hmap_params['Heatmap_height'] = float(15)
        elif n_pphmms > 100 and n_pphmms <= 500:
            hmap_params['Heatmap_width'] = float(15)
            hmap_params['Heatmap_height'] = float(10)
        else:
            hmap_params['Heatmap_width'] = float(10)
            hmap_params['Heatmap_height'] = float(8)
    else:
        hmap_params['Heatmap_width'] = float(12)
        hmap_params['Heatmap_height'] = hmap_params['Heatmap_width']
    hmap_params['TaxoLable_space'] = 1.00

    hmap_params['CBar_Heatmap_gap'] = 0.05
    hmap_params['CBar_width'] = hmap_params['Heatmap_width']
    hmap_params['CBar_height'] = 0.50
    hmap_params['CBarLable_space'] = 0.25

    hmap_params['Dendrogram_width'] = hmap_params['Heatmap_width']/3
    hmap_params['Dendrogram_height'] = hmap_params['Heatmap_height']
    hmap_params['Dendrogram_Heatmap_gap'] = 0.1

    hmap_params['ScaleBar_Dendrogram_gap'] = hmap_params['CBar_Heatmap_gap']
    hmap_params['ScaleBar_width'] = hmap_params['Dendrogram_width']
    hmap_params['ScaleBar_height'] = hmap_params['CBar_height']
    hmap_params['ScaleBarLable_space'] = hmap_params['CBarLable_space']

    '''Dynamic sizing'''
    hmap_params['linewidth_major'] = 0.4
    hmap_params['linewidth_minor'] = 0.2

    if n_viruses < 100:
        hmap_params['FontSize'] = 14

    elif n_viruses >= 200 and n_viruses <400:
        '''Reduce draw parameter size if large n'''
        hmap_params['FontSize'] = hmap_params['FontSize']/2
        hmap_params['linewidth_major'] = hmap_params['linewidth_major']/2
        hmap_params['linewidth_minor'] = hmap_params['linewidth_minor']/2
        # hmap_params['dpi']=hmap_params['dpi']*1.5

    elif n_viruses >= 400 and n_viruses <600:
        hmap_params['FontSize'] = hmap_params['FontSize']/8
        hmap_params['linewidth_major'] = hmap_params['linewidth_major']/4
        hmap_params['linewidth_minor'] = hmap_params['linewidth_minor']/4
        # hmap_params['dpi']=hmap_params['dpi']*1.5 ## RM TODO << FLAVI SET IS LESS THAN 600!!

    elif n_viruses >= 600:
        hmap_params['FontSize'] = 0.5
        hmap_params['linewidth_major'] = hmap_params['linewidth_major']/10
        hmap_params['linewidth_minor'] = hmap_params['linewidth_minor']/10
        # hmap_params['dpi']=hmap_params['dpi']*1.5
        hmap_params['Heatmap_width'] = float(24)

    hmap_params['Fig_width'] = hmap_params['Outer_margin'] + hmap_params['Dendrogram_width'] + hmap_params['Dendrogram_Heatmap_gap'] + \
        hmap_params['Heatmap_width'] + hmap_params['TaxoLable_space'] + hmap_params['Outer_margin']
    hmap_params['Fig_height'] = hmap_params['Outer_margin'] + hmap_params['CBarLable_space'] + hmap_params['CBar_height'] + \
        hmap_params['CBar_Heatmap_gap'] + hmap_params['Heatmap_height'] + hmap_params['TaxoLable_space'] + hmap_params['Outer_margin']

    hmap_params['ax_Dendrogram_L'] = hmap_params['Outer_margin']/hmap_params['Fig_width']
    hmap_params['ax_Dendrogram_B'] = (hmap_params['Outer_margin'] + hmap_params['ScaleBarLable_space'] +
                        hmap_params['ScaleBar_height'] + hmap_params['ScaleBar_Dendrogram_gap'])/hmap_params['Fig_height']
    hmap_params['ax_Dendrogram_W'] = hmap_params['Dendrogram_width']/hmap_params['Fig_width']
    hmap_params['ax_Dendrogram_H'] = hmap_params['Dendrogram_height']/hmap_params['Fig_height']

    hmap_params['ax_ScaleBar_L'] = hmap_params['Outer_margin']/hmap_params['Fig_width']
    hmap_params['ax_ScaleBar_B'] = (hmap_params['Outer_margin'] + hmap_params['ScaleBarLable_space'])/hmap_params['Fig_height']
    hmap_params['ax_ScaleBar_W'] = hmap_params['ScaleBar_width']/hmap_params['Fig_width']
    hmap_params['ax_ScaleBar_H'] = hmap_params['ScaleBar_height']/hmap_params['Fig_height']

    hmap_params['ax_Heatmap_L'] = (hmap_params['Outer_margin'] + hmap_params['Dendrogram_width'] +
                    hmap_params['Dendrogram_Heatmap_gap'])/hmap_params['Fig_width']
    hmap_params['ax_Heatmap_B'] = (hmap_params['Outer_margin'] + hmap_params['CBarLable_space'] +
                    hmap_params['CBar_height'] + hmap_params['CBar_Heatmap_gap'])/hmap_params['Fig_height']
    hmap_params['ax_Heatmap_W'] = hmap_params['Heatmap_width']/hmap_params['Fig_width']
    hmap_params['ax_Heatmap_H'] = hmap_params['Heatmap_height']/hmap_params['Fig_height']

    hmap_params['ax_CBar_L'] = (hmap_params['Outer_margin'] + hmap_params['Dendrogram_width'] +
                    hmap_params['Dendrogram_Heatmap_gap'])/hmap_params['Fig_width']
    hmap_params['ax_CBar_B'] = (hmap_params['Outer_margin'] + hmap_params['CBarLable_space'])/hmap_params['Fig_height']
    hmap_params['ax_CBar_W'] = hmap_params['CBar_width']/hmap_params['Fig_width']
    hmap_params['ax_CBar_H'] = hmap_params['CBar_height']/hmap_params['Fig_height']

    '''Plot heat map'''
    fig = plt.figure(figsize=(hmap_params['Fig_width'], hmap_params['Fig_height']), dpi=hmap_params['dpi'])

    '''Draw Dendrogram'''
    ax_dendrogram = fig.add_axes(
        [hmap_params['ax_Dendrogram_L'], hmap_params['ax_Dendrogram_B'], hmap_params['ax_Dendrogram_W'], hmap_params['ax_Dendrogram_H']], frame_on=False, facecolor="white")

    '''Dendrogram scale bar'''
    ax_ScaleBar = fig.add_axes(
        [hmap_params['ax_ScaleBar_L'], hmap_params['ax_ScaleBar_B'], hmap_params['ax_ScaleBar_W'], hmap_params['ax_ScaleBar_H']], frame_on=False, facecolor="white")


    if do_logscale:
        ScaleBarTicks = [10**0, 10**0.25, 10**0.50, 10**0.70, 10**0.9, 10**1]
        ScaleBarTickLabels = [0, 0.25, 0.50, 0.70, 0.9, 1.0]

        ax_ScaleBar		.plot([10**0, 10**1], [0, 0], 'k-')
        for Tick in ScaleBarTicks:
            ax_ScaleBar.plot([Tick, Tick], [-0.05, 0.05], 'k-')

        ax_ScaleBar		.set_xlim([10**1, 10**0])
        ax_ScaleBar		.set_xlabel('10^Distance', rotation=0, size=hmap_params['FontSize']+2)

    else:
        ScaleBarTickLabels = ScaleBarTicks = [0, 0.25, 0.5, 0.75, 1]
        ax_ScaleBar		.plot([0,1], [0, 0], 'k-')
        for Tick in ScaleBarTicks:
            ax_ScaleBar.plot([Tick, Tick], [-0.05, 0.05], 'k-')

        ax_ScaleBar		.set_xlim([1, 0])
        ax_ScaleBar		.set_xlabel('Distance', rotation=0, size=hmap_params['FontSize']+2)

    ax_ScaleBar		.set_xticks(ScaleBarTicks)
    ax_ScaleBar		.set_xticklabels(
        list(map(str, ScaleBarTickLabels)), rotation=0, size=hmap_params['FontSize'])

    ax_ScaleBar		.xaxis.set_label_position('bottom')
    ax_ScaleBar		.tick_params(top=False,
                                bottom=False,
                                left=False,
                                right=False,
                                labeltop=False,
                                labelbottom=True,
                                labelleft=False,
                                labelright=False,
                                direction='out')

    '''Heatmap'''
    ax_Heatmap = fig.add_axes(
        [hmap_params['ax_Heatmap_L'], hmap_params['ax_Heatmap_B'], hmap_params['ax_Heatmap_W'], hmap_params['ax_Heatmap_H']], frame_on=True, facecolor="white")

    return hmap_params, fig, ax_dendrogram, ax_Heatmap

def construct_hmap_lines(ax_Heatmap, len_x, LineList_major, LineList_minor, hmap_params, ClassLabelList_x, ClassLabelList_y, TickLocList):
    '''Draw grouping major & minor lines on heatmap'''
    if len_x < 200: ## TODO TEST
        for l in LineList_major:
            ax_Heatmap.axvline(l, color='k', lw=hmap_params['linewidth_major'])
            ax_Heatmap.axhline(l, color='k', lw=hmap_params['linewidth_major'])

    if len_x < 200:
        for l in LineList_minor:
            ax_Heatmap.axhline(l, color='gray', lw=hmap_params['linewidth_minor'])
            ax_Heatmap.axvline(l, color='gray', lw=hmap_params['linewidth_minor'])

    ax_Heatmap			.set_xticks(TickLocList)
    ax_Heatmap			.set_xticklabels(
        ClassLabelList_x, rotation=90, fontsize=hmap_params['FontSize'])

    ax_Heatmap			.set_yticks(TickLocList)
    ax_Heatmap			.set_yticklabels(
        ClassLabelList_y, rotation=0, fontsize=hmap_params['FontSize'])

    ax_Heatmap			.tick_params(top=True,
                            bottom=False,
                            left=False,
                            right=True,
                            labeltop=True,
                            labelbottom=False,
                            labelleft=False,
                            labelright=True,
                            direction='out')
    return ax_Heatmap

def construct_wide_hmap_lines(ax_Heatmap, len_x, LineList_major, LineList_minor, hmap_params, ClassLabelList_x, ClassLabelList_y, TickLocList_y, data):
    for l in LineList_major:
        ax_Heatmap.axhline(l, color='k', lw=hmap_params['linewidth_major']*2)
    for l in LineList_minor:
        ax_Heatmap.axhline(l, color='gray', lw=hmap_params['linewidth_minor'])

    x_lines = np.arange(-1,len_x) + 0.5
    TickLocList_x = np.arange(len_x)
    if len_x < 200:
        '''No lines or labels for massive graphs'''
        for l in x_lines:
            ax_Heatmap.axvline(l, color='gray', lw=hmap_params['linewidth_minor'])
        ax_Heatmap			.set_xticks(TickLocList_x)
        ax_Heatmap			.set_xticklabels(
            ClassLabelList_x, rotation=90, size=hmap_params['FontSize'])
    else:
        # profile_xlables = np.round(np.nanmean(data, axis=0), 2)
        # profile_xlables.sort()
        # ax_Heatmap			.set_xticks(TickLocList_x)
        # ax_Heatmap			.set_xticklabels(
        #     profile_xlables, rotation=90, size=hmap_params['FontSize'])
        ...

    ax_Heatmap			.set_yticks(TickLocList_y)
    ax_Heatmap			.set_yticklabels(
        ClassLabelList_y, rotation=0, size=hmap_params['FontSize'])

    ax_Heatmap			.tick_params(top=False,
                            bottom=False,
                            left=False,
                            right=True,
                            labeltop=True,
                            labelbottom=False,
                            labelleft=False,
                            labelright=True,
                            direction='out')
    return ax_Heatmap
