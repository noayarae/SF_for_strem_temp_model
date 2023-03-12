### This is a sub-code to plot areas of SF contribution

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def test_plot1():
    from matplotlib.patches import Patch
    
    #red_patch = mpatches.Patch(facecolor="green", edgecolor="red",label='The red data')
    h1 = mpatches.Rectangle((0.25, 0.25), 0.5, 0.5,edgecolor="red", facecolor=[0.5, 1.0, 1.0, 1],
                     hatch=' ... ', label=r'$\pi_e\ free$')
    h2 = mpatches.Rectangle((0.25, 0.25), 0.5, 0.5, edgecolor="red", facecolor=[0.5, 1.0, 1.0, 1],
                     hatch='xxxx', label=r'$\pi_e\ free$')
    h3 = mpatches.Rectangle((0.25, 0.25), 0.5, 0.5,edgecolor="red", facecolor=[0.5, 1.0, 1.0, 1],
                     hatch=' //OO ', label=r'$\pi_e\ free$', zorder = 0)

    # draw hatch
    #ax1.bar(range(1, 5), range(1, 5), color='none', edgecolor='red', hatch="/", lw=1., zorder = 0)
    # draw edge
    h3b = mpatches.Rectangle((0.25, 0.25), 0.5, 0.5, color='b', edgecolor='k',
                             label=r'$\pi_e\ free$', zorder=1, lw=2.)


    
    leg = plt.legend(handles=[h1,h2,h3,h3b],
                     labelspacing = 0.2 # space between labels (0.8 by default)
                     )
    print (leg.get_patches())

    leg.get_patches()[0].set_height(20)
    leg.get_patches()[1].set_height(15)
    leg.get_patches()[2].set_height(10)

    '''
    for patch in leg.get_patches():
        print patch
        patch.set_height(20) # Set high of hatch-rectangle
        patch.set_y(-6)     # Set gap or distance of hatch-rectangle from the top
        #patch.set_x(-20)     # Set gap or distance of hatch-rectangle from the right
        #'''

    plt.show()#'''
    '''
    import matplotlib.lines as mlines
    blue_line = mlines.Line2D([],[], color='blue', marker='*',
                          markersize=15, label='Blue stars')
    plt.legend(handles=[blue_line])
    plt.show()#'''

def test_plot2():
    import matplotlib.pyplot as plt
    from matplotlib.path import Path
    import matplotlib.patches as patches

    verts1 = [
        (0., 0.), # left, bottom
        (0., 1.), # left, top
        (1., 1.), # right, top
        (1., 0.), # right, bottom
        (0., 0.), # ignored
        ]

    codes1 = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY, ]

    path = Path(verts1, codes1)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    p1 = patches.PathPatch(path, color='green', fill=False, hatch='xx', lw=2)
    p2 = patches.PathPatch(path, color='k', fill=False, lw=2)
    ax.add_patch(p1)
    ax.add_patch(p2)
    ax.set_xlim(-2,2)
    ax.set_ylim(-2,2)
    plt.show()



#def plot_areas (sf_yr_top, sf_yr_R_tree, sf_yr_L_tree, sf_yr_all):
def plot_areas (sb, sf_top, sf_R, sf_L, sf_all, strm_id, sf_id):
    sf_cases = ["(No riparian)", "(Current riparian)", "(Max. riparian)",
                "(Effic. riparian)","(Effic2. riparian)"]
    
    ### ---- Plot Shade factor over the year ---- ###
    ### Data to plot
    x_pts = np.array(list(range(1,366)))
    top_pts = np.array(sf_top)
    Rtree_pts = np.array(sf_R)
    Ltree_pts = np.array(sf_L)
    all_pts = np.array(sf_all)

    ar1 = top_pts
    ar2 = np.add(ar1, Rtree_pts)
    ar3 = np.add(ar2, Ltree_pts)

    ### Basic setings of the plot
    fig = plt.figure(figsize=(8, 4), facecolor='white') # Figure size, background color
    ### Setings of the 4 lines: Top, Rtree, Ltree, All
    #p1, = plt.plot(x_pts, top_pts, color='b', label='Top',dashes=[4,2], linestyle = ':', linewidth = '2')
    #p2, = plt.plot(x_pts, Rtree_pts, color='g', label='R-tree',dashes=[10,3], linestyle = '--', linewidth = '2')
    #p3, = plt.plot(x_pts, Ltree_pts, color='r', label='L-tree',dashes=[12,2,3,2,3,2], linestyle = '-.', linewidth = '2')
    #p1, = plt.plot(x_pts, ar1, color='b', label='Top',dashes=[4,2], linestyle = ':', linewidth = '2')
    #p2, = plt.plot(x_pts, ar2, color='g', label='R-tree',dashes=[10,3], linestyle = '--', linewidth = '2')
    #p3, = plt.plot(x_pts, ar3, color='r', label='L-tree',dashes=[12,2,3,2,3,2], linestyle = '-.', linewidth = '2')
    p1, = plt.plot(x_pts, ar1, color='b', label='Top',linestyle = '-', linewidth = '1')
    p2, = plt.plot(x_pts, ar2, color='#808080', label='R-tree', linestyle = '-', linewidth = '1')
    p3, = plt.plot(x_pts, ar3, color='r', label='L-tree', linestyle = '-', linewidth = '1')
    #p4, = plt.plot(x_pts, all_pts, color='k', label='All', linestyle = '-', linewidth = '2')
    p4, = plt.plot(x_pts, all_pts, color='k', label='All', linestyle = '-', linewidth = '2.5')
    h1 = mpatches.Rectangle((0.25, 0.25), 0.5, 0.5,edgecolor="black", facecolor="#ffffff",
                     hatch=' ... ', label=r'$\pi_e\ free$') # Right
    h2 = mpatches.Rectangle((0.25, 0.25), 0.5, 0.5, edgecolor="red", facecolor="#ffffff",
                     hatch='xxxx', label=r'$\pi_e\ free$') # Left
    h3 = mpatches.Rectangle((0.25, 0.25), 0.5, 0.5,edgecolor="blue",
                            facecolor="#ffffff", # background color
                            hatch=' //OO ', label=r'$\pi_e\ free$') # Topography

    #plt.fill_between(x_pts, 0, ar1, color='#539ecd') # Blue - topography
    plt.fill_between(x_pts, 0, ar1, color="none", hatch=" //OO ", edgecolor="b", linewidth=0.0)
    #plt.fill_between(x_pts, ar1,ar2, color='#ff9980') # Red - Left #ff9980
    plt.fill_between(x_pts, ar1,ar2, color="none", hatch=" .... ", edgecolor="k", linewidth=0.0)
    #plt.fill_between(x_pts, ar2,ar3, hatch="X", color='#00b300') # Green - Right #00b300
    plt.fill_between(x_pts, ar2,ar3, color="none", hatch=" xxx ", edgecolor="r", linewidth=0.0)

    plt.xticks(np.arange(0, 370, 20))
    plt.xlim([0, 370])
    plt.ylim([0, 1.48])
    plt.grid()
    title = 'Contribution of topography and riparian vegetation on the shade factor \n Sub-basin '
    title = title + str(sb+1)+'; Azimuth='+ str(strm_id)+'; '+sf_cases[sf_id] 
    plt.subplots_adjust(top=0.88) # Set space for title
    fig.suptitle(title, fontsize=12)
    plt.xlabel('Day in the year', fontsize=10)
    plt.ylabel('Shade Factor (SF)', fontsize=10)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=10)
    
    ### Text for legend
    txt1 = "Shade factor by topographic"
    txt2 = "Shade-factor by Right bank riparian vegetation"
    txt3 = "Shade-factor by Left bank riparian vegetation"
    txt4 = "Shade-factor by all obstacles"
    th1 = "Shade-factor by Right bank riparian vegetation"
    th2 = "Shade-factor by Left bank riparian vegetation"
    th3 = "Shade-factor by topography"

    #leg = plt.legend([p1, p2, p3, p4,h1,h2,h3], [txt1, txt2,txt3,txt4, th1, th2, th3], loc='upper left',  ## Old one
    leg = plt.legend([p4,h1,h2,h3], [txt4, th1, th2, th3], loc='upper left',
                     scatterpoints = 1,
                     prop={'size': 9}, # font size in 
                     bbox_to_anchor=(0.01,0.999), # coordinates of legend box (0,0 is left bottom)
                     handlelength = 4.5, # length of lines in legend
                     labelspacing = 0.75,   # space between labels
                     handletextpad = 0.5  # gap between marks and labels  0.8 by default
                     )
    #leg.get_lines()[0].set_linewidth(2.5) # Set line-width in legend for top
    #leg.get_lines()[1].set_linewidth(2.5) # Set line-width in legend for R
    #leg.get_lines()[2].set_linewidth(2.5) # Set line-width in legend for L
    leg.get_lines()[0].set_linewidth(3.5) # Set line-width in legend for all

    leg.get_patches()[0].set_linewidth(1.0) # Set line-width in legend for all
    leg.get_patches()[0].set_height(12)
    #leg.get_patches()[0].set_y(-5) 

    leg.get_patches()[1].set_linewidth(0.5) # Set line-width in legend for all
    leg.get_patches()[1].set_height(12)
    #leg.get_patches()[1].set_y(-10)

    leg.get_patches()[2].set_linewidth(0.5) # Set line-width in legend for all
    leg.get_patches()[2].set_height(12)
    #leg.get_patches()[2].set_y(-15)

    plt.show()
    #'''

if __name__ == "__main__":
    # Input data
    plot = test_plot1()
    #plot = test_plot2()

    print ("donee...")











