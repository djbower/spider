#!/usr/bin/env python

import spider_utils as su
import matplotlib.pyplot as plt
import numpy as np

#====================================================================
def main():

    figw = 4.7747/2.0
    figh = 4.7747/2.0

    fig_o = su.FigureData( 1, 1, figw, figh, 'radius' )

    ax0 = fig_o.ax

    data_a = np.loadtxt( 'radius.dat', unpack=True )

    arcsinhscale = 1.0E3
    time_fmt = su.MyFuncFormatter( arcsinhscale )

    time = data_a[0] 
    time /= 1.0E6 # to Myr
    time = time_fmt.ascale( time )
    radius = data_a[1]
    radius_10mb = data_a[2]
    radius_1mb = data_a[3]

    ref = radius[0]
    radius = (radius-ref)/ref*100.0
    radius_10mb = (radius_10mb-ref)/ref*100.0
    radius_1mb = (radius_1mb-ref)/ref*100.0

    h1, = ax0.plot( time, radius, label='Surface' )
    h2, = ax0.plot( time, radius_10mb, label='10 mb' )
    h3, = ax0.plot( time, radius_1mb, label='1 mb' )
    handles_l = [h1,h2,h3]

    xticks = [1.0E-6,1.0E-2,1.0E-1,1.0E0,1.0E1,1.0E2]
    yticks = [-5,-4,-3,-2,-1,0,1,2,3,4]

    title = 'Evolution of radius'
    fig_o.set_myaxes(  ax0, title=title, xlabel='Time (Myr)',
        ylabel='Radius (\% change)', xticks=xticks, xfmt=time_fmt, yticks=yticks, yrotation=90 )
    fig_o.set_mylegend( ax0, handles_l, loc='upper right', TITLE='Radius' )
    #ax0.yaxis.set_label_coords(-0.14,0.5)

    fig_o.savefig(1)

#====================================================================

if __name__ == "__main__":

    main()
    plt.show()
