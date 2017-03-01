### Works with Python 2.7
# The gridspec module in matplotlib is necessary and may have to be installed separately

### Import other packages
from matplotlib import *
import matplotlib.pylab as pylab
from pylab import *
from numpy import *
from scipy import *
import numpy as np
import matplotlib.gridspec as gs
from matplotlib.ticker import MaxNLocator

def plotresY(res, priors=[0,0,0,0,0], startn=0, color1='blue', color2='gray', color3='red', ls1=':', lsty='-', lw1=1.5, ms1=0.15, mk='ko'):
    """
    This is a function that will take a single population output results file and plot the sampling history and resulting PDFs
    of the following, assuming that the results file is formatted as a matric with columns of age, metallicity, distance, 
    extinction, and helium.

    Call signature: plotresY(res)
    
    where res is the results file; this is typically loaded in a separate line as: res=loadtxt('/path/cluster_results.res',skiprows=1).

    Optional keyword arguments:
    
    =========   =======================================================
    Keyword     Description
    =========   =======================================================
    priors :    default=[0,0,0,0,0] in the order of age, metallicity, distance, extinction, and helium. When the default is assumed,
                no lines for the priors are plotted. Whenever a value other than zero is specified, the prior for that variable is plot
                as a line.
    startn :    default=0; this is the starting iteration value for plotting the variables.
    color1 :    default='blue'; this is the color of the median and 90% bayesian interval lines in the PDF plots
    color2 :    default='gray'; this is the color of the histograms in the PDF plots
    color3 :    default='red'; this is the color of the prior lines across the entire plot
    ls1 :       default=':', a dotted line; this is the linestyle of the prior lines and bayesian interval lines
    lsty :      default='-', a solid line; this is the linestyle of the median bayesian interval line
    lw1 :       default=1.5; this is the linewidth of the all lines in the plot
    ms1 :       default=0.15; this is the marker size of the sampling history points
    mk :        default='ko', black points; this is the marker format of the sampling history points

    Additional options for plotting colors, sizes, formats, etc. can be found in pyplot.plot() and related descriptions. Colors may take
    any standard HTML string descriptor (re: http://www.w3schools.com/html/html_colornames.asp).

    =========
    Example
    =========
    ## loading the data beforehand:
    results=loadtxt('NGC288.single.res',skiprows=1)
    ## plotting the data
    plotresY(results, startn=1001, priors=[10.06,-1.07,15,0.08,0.],color2='greenyellow')
    ## saving and showing the plot
    savefig('NGC288.png',dpi=200)
    plt.show()

    """
    ############  set up figure  ############
    fig = matplotlib.pylab.figure(figsize=(13, 12.5))#**figprops)                                                  # New figure
    fig.subplots_adjust(wspace=0.,hspace=0.)                                              # Tunes the subplot layout
    
    ############  set up data  ############
    age=[];feh=[];dist=[];Av=[];post=[];Y=[]
    for i in arange(startn,len(res)):
        age.append(res[i][0])
        feh.append(res[i][2])
        dist.append(res[i][3])
        Av.append(res[i][4])
        Y.append(res[i][1])
    
     ############  set up axes column 1  ############
    ax = plt.subplot2grid((5, 6), (0, 0),colspan=5)
    bx = plt.subplot2grid((5, 6), (1, 0),colspan=5)
    cx = plt.subplot2grid((5, 6), (2, 0),colspan=5)
    dx = plt.subplot2grid((5, 6), (3, 0),colspan=5)
    ex = plt.subplot2grid((5, 6), (4, 0),colspan=5)

    ############  plot sampling history of variables in first column  ############
    variables=[age,feh,dist,Av,Y]
    varnames=['log(age)','[Fe/H]','Dist Mod','A$_{V}$','Y','Iteration']
    xx=arange(0,len(age))
    for i,axes in enumerate([ax,bx,cx,dx,ex]):
        axes.plot(xx,variables[i],mk,ms=ms1)
        axes.set_ylabel(varnames[i])
    ex.set_xlabel('Iteration')

    ############  find and print to terminal the 90% CI for parameters  ############
    age_p=round(np.percentile(age,5),4),round(np.percentile(age,95),4)
    feh_p=round(np.percentile(feh,5),4),round(np.percentile(feh,95),4)
    dist_p=round(np.percentile(dist,5),4),round(np.percentile(dist,95),4)
    Av_p=round(np.percentile(Av,5),4),round(np.percentile(Av,95),4)
    Y_p=round(np.percentile(Y,5),4),round(np.percentile(Y,95),4)
    
    varnames2=['log(age)  ','[Fe/H]    ','Dist Mod  ','A$_{V}$   ','Y         ']
    print '              median      90% Confidence Interval'
    for i,varnames2 in enumerate(varnames2):
        print varnames2, '  ', round(median(variables[i]),4), '   ', round(np.percentile(variables[i],5),4),'  ',round(np.percentile(variables[i],95),4)
    
    ############  plot the PDFs to the side  ############
    ax2 = plt.subplot2grid((5, 6), (0, 5))
    bx2 = plt.subplot2grid((5, 6), (1, 5))
    cx2 = plt.subplot2grid((5, 6), (2, 5))
    dx2 = plt.subplot2grid((5, 6), (3, 5))
    ex2 = plt.subplot2grid((5, 6), (4, 5))
    
    CIs=[age_p,feh_p,dist_p,Av_p,Y_p]
    for i,axes in enumerate([ax2,bx2,cx2,dx2,ex2]):
        axes.hist(variables[i],color=color2,bins=15, orientation='horizontal', normed=True)
        axes.axhline(median(variables[i]),color=color1,linestyle=lsty,linewidth=lw1)
        axes.axhline(CIs[i][0],color=color1,linestyle=ls1,linewidth=lw1)
        axes.axhline(CIs[i][1],color=color1,linestyle=ls1,linewidth=lw1)

    ############  format axes labels and tickmarks  ############
    ex2.set_xlabel('PDF')
    for axis in (ax2,bx2,cx2,dx2,ex2):
        axis.yaxis.tick_right()
    for axis in (ax,bx,cx,dx,ax2,bx2,cx2,dx2,ex2):
        setp(axis.get_xticklabels(),visible=False)
    for axis in (ax,bx,cx,dx,ex,ax2,bx2,cx2,dx2,ex2):
        axis.get_yaxis().get_major_formatter().set_useOffset(False)
        axis.yaxis.set_major_locator(MaxNLocator(nbins=7))#prune='both'))
    for axis in [ax,bx,cx,dx,ex]:
        axis.set_xlim(xmax=len(age)+1)


    ############  plot prior lines, if defined as non-zero  ############
    axess=[ax,ax2,bx,bx2,cx,cx2,dx,dx2,ex,ex2];j=0
    for i in arange(0,len(axess),2):
        if priors[j]!=0:
            axess[i].axhline(y=priors[j],xmin=0,xmax=1,color=color3,linestyle=ls1,lw=lw1)
            axess[i+1].axhline(y=priors[j],xmin=0,xmax=1,color=color3,linestyle=ls1,lw=lw1)
        j=j+1
    
    return
