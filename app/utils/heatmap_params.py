from matplotlib.colors import LinearSegmentedColormap as LSC
import matplotlib.pyplot as plt

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

def get_hmap_params(n_viruses):
    hmap_params = {}
    '''General'''
    hmap_params['Outer_margin'] = 0.5
    hmap_params['FontSize'] = 6
    hmap_params['dpi']=600
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

    if n_viruses >= 200:
        '''Reduce draw parameter size if large n'''
        hmap_params['FontSize'] = hmap_params['FontSize']/2
        hmap_params['linewidth_major'] = hmap_params['linewidth_major']/2
        hmap_params['linewidth_minor'] = hmap_params['linewidth_minor']/2
        hmap_params['dpi']=hmap_params['dpi']*1.5


    elif n_viruses >= 400:
        hmap_params['FontSize'] = hmap_params['FontSize']/4
        hmap_params['linewidth_major'] = hmap_params['linewidth_major']/4
        hmap_params['linewidth_minor'] = hmap_params['linewidth_minor']/4
        hmap_params['dpi']=hmap_params['dpi']*3

    elif n_viruses >= 600:
        hmap_params['FontSize'] = 1
        hmap_params['linewidth_major'] = hmap_params['linewidth_major']/8
        hmap_params['linewidth_minor'] = hmap_params['linewidth_minor']/8
        hmap_params['dpi']=hmap_params['dpi']*3
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
    ax_ScaleBar		.plot([0, 1], [0, 0], 'k-')
    ScaleBarTicks = [0, 0.25, 0.5, 0.75, 1]
    for Tick in ScaleBarTicks:
        ax_ScaleBar.plot([Tick, Tick], [-0.05, 0.05], 'k-')

    ax_ScaleBar		.set_xlim([1, 0])
    ax_ScaleBar		.set_xticks(ScaleBarTicks)
    ax_ScaleBar		.set_xticklabels(
        list(map(str, ScaleBarTicks)), rotation=0, size=hmap_params['FontSize'])
    ax_ScaleBar		.set_xlabel('Distance', rotation=0, size=hmap_params['FontSize']+2)
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

def construct_hmap_lines(ax_Heatmap, LineList_major, LineList_minor, hmap_params, ClassLabelList_x, ClassLabelList_y, TickLocList):
    '''Draw grouping major & minor lines on heatmap'''
    for l in LineList_major:
        ax_Heatmap.axvline(l, color='k', lw=hmap_params['linewidth_major'])
        ax_Heatmap.axhline(l, color='k', lw=hmap_params['linewidth_major'])

    for l in LineList_minor:
        ax_Heatmap.axvline(l, color='gray', lw=hmap_params['linewidth_minor'])
        ax_Heatmap.axhline(l, color='gray', lw=hmap_params['linewidth_minor'])
    ax_Heatmap			.set_xticks(TickLocList)
    ax_Heatmap			.set_xticklabels(
        ClassLabelList_x, rotation=90, size=hmap_params['FontSize'])

    ax_Heatmap			.set_yticks(TickLocList)
    ax_Heatmap			.set_yticklabels(
        ClassLabelList_y, rotation=0, size=hmap_params['FontSize'])

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
