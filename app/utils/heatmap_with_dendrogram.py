
from Bio import Phylo
from matplotlib.colors import LinearSegmentedColormap as LSC
import matplotlib.pyplot as plt

### RM < NEEDS LOTS MORE WORK
def heatmap_with_dendrogram(Bootstrap,
                            BootstrappedVirusDendrogramFile,
                            VirusDendrogramFile,
                            TaxoLabelList_AllVirus,
                            ):
    print("\tConstruct GRAViTy heat map with dendrogram")
    #-------------------------------------------------------------------------------
    #Load the tree
    #-------------------------------------------------------------------------------
    if Bootstrap == True:
        VirusDendrogram		= Phylo.read(BootstrappedVirusDendrogramFile, "newick")
    else:
        VirusDendrogram		= Phylo.read(VirusDendrogramFile, "newick")
    
    #Determine virus order
    #-------------------------------------------------------------------------------
    _ 			= VirusDendrogram.ladderize(reverse = True)
    OrderedTaxoLabelList	= [Clade.name for Clade in VirusDendrogram.get_terminals()]
    VirusOrder		= [TaxoLabelList_AllVirus.index(TaxoLabel) for TaxoLabel in OrderedTaxoLabelList]
    
    #Re-order the distance matrix
    #-------------------------------------------------------------------------------
    OrderedDistMat	= DistMat[VirusOrder][:,VirusOrder]
    
    #Remove clade support values that are < Heatmap_DendrogramSupport_Cutoff
    #-------------------------------------------------------------------------------
    N_InternalNodes = len(VirusDendrogram.get_nonterminals())
    for InternalNode_i in range(N_InternalNodes):
        if VirusDendrogram.get_nonterminals()[InternalNode_i].confidence >= Heatmap_DendrogramSupport_Cutoff:
            VirusDendrogram.get_nonterminals()[InternalNode_i].confidence = round(VirusDendrogram.get_nonterminals()[InternalNode_i].confidence, 2)
        else:
            VirusDendrogram.get_nonterminals()[InternalNode_i].confidence=""
    
    #Colour terminal branches: reference virus's branch is blue, unclassified virus's branch is red
    #-------------------------------------------------------------------------------
    N_Viruses = N_RefViruses + N_UcfViruses
    for Virus_i in range(N_Viruses):
        Taxolabel = VirusDendrogram.get_terminals()[Virus_i].name
        if Taxolabel in TaxoLabelList_RefVirus:
            VirusDendrogram.get_terminals()[Virus_i].color="blue"
        elif Taxolabel in TaxoLabelList_UcfVirus:
            VirusDendrogram.get_terminals()[Virus_i].color="red"
    
    #Labels, label positions, and ticks
    #-------------------------------------------------------------------------------
    TaxoGroupingList_AllVirus	= TaxoGroupingList_RefVirus.tolist() + TaxoAssignmentList
    Taxo2ClassDict		= {TaxoLabel: TaxoGrouping for TaxoLabel, TaxoGrouping in zip(TaxoLabelList_AllVirus, TaxoGroupingList_AllVirus)}
    ClassDendrogram		= copy(VirusDendrogram)
    for Clade in ClassDendrogram.find_clades(terminal=True):
        Clade.name = Taxo2ClassDict[Clade.name]
    
    ClassLabelList = []
    LineList = [-1]
    TerminalNodeList = [TerminalNode for TerminalNode in ClassDendrogram.get_terminals()]
    while len(TerminalNodeList)!=0:
        FarLeftNode = TerminalNodeList[0]
        for Clade in ([ClassDendrogram]+ClassDendrogram.get_path(FarLeftNode)):
            DescendantNodeList = Clade.get_terminals()
            DescendantClassLabelList = list(set([c.name for c in DescendantNodeList]))
            if len(DescendantClassLabelList)==1:
                ClassLabelList.append(DescendantClassLabelList[0])
                LineList.append(LineList[-1]+len(DescendantNodeList))
                TerminalNodeList = TerminalNodeList[len(DescendantNodeList):]
                break
    
    ClassLabelList	= np.array(ClassLabelList)
    LineList	= np.array(LineList) + 0.5
    TickLocList	= np.array(list(map(np.mean, list(zip(LineList[0:-1],LineList[1:])))))
    
    #Heat map colour indicators
    #-------------------------------------------------------------------------------
    IndicatorMat_RefVirus	= np.tile([True]*N_RefViruses + [False]*N_UcfViruses, N_Viruses).reshape(N_Viruses,-1)
    IndicatorMat_RefVirus	= IndicatorMat_RefVirus * IndicatorMat_RefVirus.T
    IndicatorMat_RefVirus	= IndicatorMat_RefVirus[VirusOrder][:,VirusOrder]
    
    IndicatorMat_UcfVirus	= np.tile([False]*N_RefViruses + [True]*N_UcfViruses, N_Viruses).reshape(N_Viruses,-1)
    IndicatorMat_UcfVirus	= IndicatorMat_UcfVirus * IndicatorMat_UcfVirus.T
    IndicatorMat_UcfVirus	= IndicatorMat_UcfVirus[VirusOrder][:,VirusOrder]
    
    IndicatorMat_CrossGroup	= ~IndicatorMat_RefVirus * ~IndicatorMat_UcfVirus
    
    #Masked OrderedDistMat
    #-------------------------------------------------------------------------------
    OrderedDistMat_RefVirus = np.ma.masked_where(~IndicatorMat_RefVirus, OrderedDistMat)
    OrderedDistMat_UcfVirus = np.ma.masked_where(~IndicatorMat_UcfVirus, OrderedDistMat)
    OrderedDistMat_CrossGroup = np.ma.masked_where(~IndicatorMat_CrossGroup, OrderedDistMat)
    
    #Colour map construction
    #-------------------------------------------------------------------------------
    BlueColour_Dict = {
    'red'  :  ((0., 0., 0.), (1., 1., 1.)),
    'green':  ((0., 0., 0.), (1., 1., 1.)),
    'blue' :  ((0., 1., 1.), (1., 1., 1.))
    }
    
    RedColour_Dict = {
    'red'  :  ((0., 1., 1.), (1., 1., 1.)),
    'green':  ((0., 0., 0.), (1., 1., 1.)),
    'blue' :  ((0., 0., 0.), (1., 1., 1.))
    }
    
    PurpleColour_Dict = {
    'red'  :  ((0., 0.5, 0.5), (1., 1., 1.)),
    'green':  ((0., 0., 0.), (1., 1., 1.)),
    'blue' :  ((0., 1., 1.), (1., 1., 1.))
    }
    
    MyBlues	= LSC('MyBlues', BlueColour_Dict, 1024)
    MyReds = LSC('MyReds', RedColour_Dict, 1024)
    MyPurples = LSC('MyPurples', PurpleColour_Dict, 1024)
    
    #Plot configuration
    #-------------------------------------------------------------------------------
    Heatmap_width		= float(12)
    Heatmap_height		= Heatmap_width
    TaxoLable_space		= 1.00
    
    CBar_Heatmap_gap	= 0.05
    CBar_width		= Heatmap_width
    CBar_height		= 0.50
    CBarLable_space		= 0.25
    
    Dendrogram_width	= Heatmap_width/3
    Dendrogram_height	= Heatmap_height
    Dendrogram_Heatmap_gap	= 0.1
    
    ScaleBar_Dendrogram_gap	= CBar_Heatmap_gap
    ScaleBar_width		= Dendrogram_width
    ScaleBar_height		= CBar_height
    ScaleBarLable_space	= CBarLable_space
    
    Outer_margin		= 0.5
    FontSize		= 6
    
    Fig_width		= Outer_margin + Dendrogram_width + Dendrogram_Heatmap_gap + Heatmap_width + TaxoLable_space + Outer_margin
    Fig_height		= Outer_margin + CBarLable_space + CBar_height + CBar_Heatmap_gap + Heatmap_height + TaxoLable_space + Outer_margin
    
    ax_Dendrogram_L		= Outer_margin/Fig_width
    ax_Dendrogram_B		= (Outer_margin + ScaleBarLable_space + ScaleBar_height + ScaleBar_Dendrogram_gap)/Fig_height
    ax_Dendrogram_W		= Dendrogram_width/Fig_width
    ax_Dendrogram_H		= Dendrogram_height/Fig_height
    
    ax_ScaleBar_L		= Outer_margin/Fig_width
    ax_ScaleBar_B		= (Outer_margin + ScaleBarLable_space)/Fig_height
    ax_ScaleBar_W		= ScaleBar_width/Fig_width
    ax_ScaleBar_H		= ScaleBar_height/Fig_height
    
    ax_Heatmap_L		= (Outer_margin + Dendrogram_width + Dendrogram_Heatmap_gap)/Fig_width
    ax_Heatmap_B		= (Outer_margin + CBarLable_space + CBar_height + CBar_Heatmap_gap)/Fig_height
    ax_Heatmap_W		= Heatmap_width/Fig_width
    ax_Heatmap_H		= Heatmap_height/Fig_height
    
    ax_CBar_L		= (Outer_margin + Dendrogram_width + Dendrogram_Heatmap_gap)/Fig_width
    ax_CBar_B		= (Outer_margin + CBarLable_space)/Fig_height
    ax_CBar_W		= CBar_width/Fig_width
    ax_CBar_H		= CBar_height/Fig_height
    
    #Plot the heat map
    #-------------------------------------------------------------------------------
    fig			= plt.figure(figsize = (Fig_width, Fig_height), dpi = 300)
    
    ax_Dendrogram		= fig.add_axes([ax_Dendrogram_L, ax_Dendrogram_B, ax_Dendrogram_W, ax_Dendrogram_H], frame_on = False, facecolor = "white")
    Phylo			.draw(VirusDendrogram, label_func = lambda x: "", do_show = False,  axes = ax_Dendrogram)
    VirusDendrogramDepth	= max([v for k,v in VirusDendrogram.depths().items()])
    ax_Dendrogram		.set_xlim([(VirusDendrogramDepth - 1), VirusDendrogramDepth])
    ax_Dendrogram		.set_ylim([N_Viruses+0.5, 0.5])
    ax_Dendrogram		.set_axis_off()
    
    ax_ScaleBar		= fig.add_axes([ax_ScaleBar_L, ax_ScaleBar_B, ax_ScaleBar_W, ax_ScaleBar_H], frame_on = False, facecolor = "white")
    ax_ScaleBar		.plot([0,1],[0,0],'k-')
    ScaleBarTicks		= [0, 0.25, 0.5, 0.75, 1]
    for Tick in ScaleBarTicks:
        ax_ScaleBar.plot([Tick, Tick],[-0.05, 0.05],'k-')
    
    ax_ScaleBar		.set_xlim([1, 0])
    ax_ScaleBar		.set_xticks(ScaleBarTicks)
    ax_ScaleBar		.set_xticklabels(list(map(str, ScaleBarTicks)), rotation = 0, size = FontSize)
    ax_ScaleBar		.set_xlabel('Distance', rotation = 0, size = FontSize+2)
    ax_ScaleBar		.xaxis.set_label_position('bottom')
    ax_ScaleBar		.tick_params(	top = 'off',
                        bottom = 'off',
                        left = 'off',
                        right = 'off',
                        labeltop = 'off',
                        labelbottom = 'on',
                        labelleft = 'off',
                        labelright = 'off',
                        direction = 'out')
    
    ax_Heatmap		= fig.add_axes([ax_Heatmap_L, ax_Heatmap_B, ax_Heatmap_W, ax_Heatmap_H], frame_on = True, facecolor = "white")
    ax_Heatmap		.imshow(OrderedDistMat_RefVirus, cmap = MyBlues, aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
    ax_Heatmap		.imshow(OrderedDistMat_UcfVirus, cmap = MyReds, aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
    ax_Heatmap		.imshow(OrderedDistMat_CrossGroup, cmap = MyPurples, aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
    for l in LineList:
        ax_Heatmap.axvline(l, color = 'k', lw = 0.2)
        ax_Heatmap.axhline(l, color = 'k', lw = 0.2)
    
    ax_Heatmap		.set_xticks(TickLocList)
    ax_Heatmap		.set_xticklabels(ClassLabelList, rotation = 90, size = FontSize)
    ax_Heatmap		.set_yticks(TickLocList)
    ax_Heatmap		.set_yticklabels(ClassLabelList, rotation = 0, size = FontSize)
    ax_Heatmap		.tick_params(	top = 'on',
                        bottom = 'off',
                        left = 'off',
                        right = 'on',
                        labeltop = 'on',
                        labelbottom = 'off',
                        labelleft = 'off',
                        labelright = 'on',
                        direction = 'out')
    
    ax_CBar_RefVirus	= fig.add_axes([ax_CBar_L, ax_CBar_B + 2*ax_CBar_H/3, ax_CBar_W, ax_CBar_H/3], frame_on = True, facecolor = "white")
    ax_CBar_RefVirus	.imshow(np.linspace(0, 1, 1025).reshape(1,-1), cmap = MyBlues, aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
    ax_CBar_RefVirus	.set_yticks([0.0])
    ax_CBar_RefVirus	.set_yticklabels(["Ref viruses"], rotation = 0, size = FontSize + 2)
    ax_CBar_RefVirus	.tick_params(	top = 'off',
                        bottom = 'off',
                        left = 'off',
                        right = 'on',
                        labeltop = 'off',
                        labelbottom = 'off',
                        labelleft = 'off',
                        labelright = 'on',
                        direction = 'out')
    ax_CBar_UcfVirus	= fig.add_axes([ax_CBar_L, ax_CBar_B + 1*ax_CBar_H/3, ax_CBar_W, ax_CBar_H/3], frame_on = True, facecolor = "white")
    ax_CBar_UcfVirus	.imshow(np.linspace(0, 1, 1025).reshape(1,-1), cmap = MyReds, aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
    ax_CBar_UcfVirus	.set_yticks([0.0])
    ax_CBar_UcfVirus	.set_yticklabels(["Ucf viruses"], rotation = 0, size = FontSize + 2)
    ax_CBar_UcfVirus	.tick_params(	top = 'off',
                        bottom = 'off',
                        left = 'off',
                        right = 'on',
                        labeltop = 'off',
                        labelbottom = 'off',
                        labelleft = 'off',
                        labelright = 'on',
                        direction = 'out')
    ax_CBar_CrossGroup	= fig.add_axes([ax_CBar_L, ax_CBar_B + 0*ax_CBar_H/3, ax_CBar_W, ax_CBar_H/3], frame_on = True, facecolor = "white")
    ax_CBar_CrossGroup	.imshow(np.linspace(0, 1, 1025).reshape(1,-1), cmap = MyPurples, aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
    ax_CBar_CrossGroup	.set_yticks([0.0])
    ax_CBar_CrossGroup	.set_yticklabels(["Ref VS Ucf viruses"], rotation = 0, size = FontSize + 2)
    ax_CBar_CrossGroup	.set_xticks(np.array([0, 0.25, 0.50, 0.75, 1])*(1025)-0.5)
    ax_CBar_CrossGroup	.set_xticklabels(['0', '0.25', '0.50', '0.75', '1'], rotation = 0, size = FontSize)
    ax_CBar_CrossGroup	.set_xlabel("Distance", rotation = 0, size = FontSize + 2)
    ax_CBar_CrossGroup	.tick_params(	top = 'off',
                        bottom = 'on',
                        left = 'off',
                        right = 'on',
                        labeltop = 'off',
                        labelbottom = 'on',
                        labelleft = 'off',
                        labelright = 'on',
                        direction = 'out')
    
    #Save the plot to file
    HeatmapWithDendrogramFile = self.VariableShelveDir_UcfVirus+"/HeatmapWithDendrogram.RefVirusGroup=%s.IncompleteUcfRefGenomes=%s.Scheme=%s.Method=%s.p=%s.pdf"%(RefVirusGroup, str(int(IncludeIncompleteGenomes_UcfVirus))+str(int(IncludeIncompleteGenomes_RefVirus)), SimilarityMeasurementScheme, Dendrogram_LinkageMethod, p)
    plt.savefig(HeatmapWithDendrogramFile, format = "pdf")