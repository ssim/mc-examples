# load defaults for different journals
def loadDefaultPlotParams( journal=None ):
    """
    load default parameters for plots:
     - if journal=None, load color-blind friendly colors and dashes
     - if journal='mnras' or 'aa', set corresponding font and figure sizes
    """
    import matplotlib.pyplot as plt
    from itertools import cycle
    # colors suitable for color-blind users
    plt.rcParams['axes.color_cycle'] = [
      (0,0,0), # black
      ( 86./255.,180./255.,233./255.), #sky blue
      (213./255., 94./255.,  0./255.), #vermilion
      (  0./255.,158./255.,115./255.), #bluish green
      (204./255.,121./255.,167./255.), #reddish purple
      (240./255.,228./255., 66./255.), #yellow
      (  0./255.,114./255.,178./255.), #blue
      (230./255.,159./255.,  0./255.), #orange
      ]
#      (0.8,0.4,0),    # vermillion (cold)
#      (0.35,0.7,0.9), # sky blue
#      (0.9,0.6,0),    # orange (cold)
#      (0,0.6,0.5),    # bluish green
#      (0.8,0.6,0.7),  # reddish purple
#      (0,0.45,0.7),   # blue
#      (0.95,0.9,0.25),# yellow (cold)
    dashes = [
      [1,0],
      [7,3],
      [2,2],
      [7,3,2,3],
      [7,3,2,3,2,3],
      [7,3,7,3,2,3],
      [7,3,7,3,2,3,2,3],
      [7,3,2,3,2,3,2,3],
      ]
    dashcycler = cycle(dashes)
    # use in plots: plot( ..., dashes=next(dashcycler))

    if ( journal != None ):
        # set fonts
        plt.rcParams['font.family'] = 'serif'
        #plt.rcParams['font.serif'] = 'txfont'
        #plt.rcParams['font.serif'] = 'Times'
        plt.rcParams['font.serif'] = ''
        plt.rcParams['font.size'] = 8
        plt.rcParams['axes.labelsize'] = 8
        plt.rcParams['legend.fontsize'] = 8
        plt.rcParams['legend.handlelength'] = 2.5
        plt.rcParams['xtick.labelsize'] = 8
        plt.rcParams['ytick.labelsize'] = 8
        plt.rcParams['text.usetex'] = True
        plt.rcParams['text.latex.preamble'] = r'\usepackage{times},\usepackage{siunitx}'
        plt.rcParams['legend.frameon'] = False
        plt.rcParams['axes.formatter.use_mathtext'] = True
        plt.rcParams['lines.markersize'] = 3
        plt.rcParams['lines.linewidth'] = 0.5
        plt.rcParams['lines.markeredgewidth'] = 0.5
        plt.rcParams['patch.linewidth'] = 0.5
        plt.rcParams['axes.linewidth'] = 0.5
        plt.rcParams['figure.subplot.bottom'] = 0.1
        plt.rcParams['figure.subplot.wspace'] = 0.1
        plt.rcParams['figure.subplot.hspace'] = 0.1

        # set page width according to latex templates
        inches_per_pt = 1.0/72.27
        if ( journal == 'mnras' ):
            # MNRAS
            pagewidth = 504.0*inches_per_pt
            columnwidth = 240.0*inches_per_pt
        elif ( journal == 'aa' ):
            # A&A
            pagewidth = 523.5307*inches_per_pt
            columnwidth = 255.76535*inches_per_pt
        elif ( journal == 'black' ):
            # A&A
            pagewidth = 523.5307*inches_per_pt
            columnwidth = 255.76535*inches_per_pt
            # black background, for presentations
            #fgcolor = 'w'
            #fgcolor = (1.0, 0.596, 0)
            fgcolor = '#FF9800'
            plt.rcParams['text.color'] = fgcolor
            plt.rcParams['xtick.color'] = fgcolor
            plt.rcParams['ytick.color'] = fgcolor
            plt.rcParams['axes.edgecolor'] = fgcolor
            plt.rcParams['axes.labelcolor'] = fgcolor
            plt.rcParams['figure.edgecolor'] = 'k'
            plt.rcParams['figure.facecolor'] = 'k'
            plt.rcParams['savefig.edgecolor'] = 'k'
            plt.rcParams['savefig.facecolor'] = 'k'
            plt.rcParams['axes.facecolor'] = 'k'

            plt.rcParams['axes.color_cycle'] = [
              (1,1,1), # white instead of black!
              (0.8,0.6,0.7),  # reddish purple
              (0.95,0.9,0.25),# yellow (cold)
              (0.35,0.7,0.9), # sky blue
              (0.8,0.4,0),    # vermillion (cold)
              (0.9,0.6,0),    # orange (cold)
              (0,0.6,0.5),    # bluish green
              (0,0.45,0.7),   # blue
              ]

            # larger fonts for posters
            #plt.rcParams['font.family'] = 'serif'
            ##plt.rcParams['font.serif'] = 'txfonts'
            #plt.rcParams['font.size'] = 15
            #plt.rcParams['axes.labelsize'] = 15
            #plt.rcParams['legend.fontsize'] = 12
            #plt.rcParams['xtick.labelsize'] = 12
            #plt.rcParams['ytick.labelsize'] = 12
            #plt.rcParams['text.usetex'] = False

        return dict(pagewidth=pagewidth, columnwidth=columnwidth, dashcycler=dashcycler, dashes=dashes, colors=plt.rcParams['axes.color_cycle'])
    else:
        return dict(dashcycler=dashcycler, dashes=dashes, colors=plt.rcParams['axes.color_cycle'])
