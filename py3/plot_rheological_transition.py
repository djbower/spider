#!/usr/bin/env python

import spider_utils as su
import logging
import matplotlib.pyplot as plt
import numpy as np
import argparse

logger = logging.getLogger(__name__)

#====================================================================
def figure7( times ):

    '''Extract magma ocean depth for partitioning of moderately
       siderophile elements'''

    # TODO: this is an intermediary plot that has been superseded
    # by fig. 8.  It can probably be removed soon, once I have the
    # desired plotting for fig. 8.

    width = 4.7747
    height = 4.7747
    fig_o = su.FigureData( 2, 2, width, height, times=times )

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]

    handle_l = []

    time = fig_o.time[0]
    myjson_o = su.MyJSON( 'output/{}.json'.format(time) )
    xx_pres_b = myjson_o.get_scaled_field_values_internal('pressure_b')
    xx_pres_b *= 1.0E-9

    for nn, time in enumerate( fig_o.time ):
        # read json
        myjson_o = su.MyJSON( 'output/{}.json'.format(time) )

        color = fig_o.get_color( nn )
        # use melt fraction to determine mixed region
        MIX = myjson_o.get_mixed_phase_boolean_array( 'basic_internal' )

        label = fig_o.get_legend_label( time )

        # regime
        yy = myjson_o.get_scaled_field_values_internal( 'regime_b' )
        ax0.plot( xx_pres_b, yy, '--', color=color )
        handle, = ax0.plot( xx_pres_b*MIX, yy*MIX, '-', color=color, label=label )
        handle_l.append( handle )

        # regime based on melt fraction
        yy = myjson_o.get_scaled_field_values_internal( 'phi_b' )
        regime = np.copy(yy)
        regime[:] = 0
        regime[yy>=0.4] = 1.0
        regime[yy<0.4] = 2.0

        ax1.plot( xx_pres_b, regime, '--', color=color )
        handle, = ax1.plot( xx_pres_b*MIX, regime*MIX, '-', color=color, label=label )
        handle_l.append( handle )

        # temperature
        #xx, yy = fig_o.get_xy_data( 'temp', time )
        #ax1.plot( xx, yy, '--', color=color )
        #ax1.plot( xx*MIX, yy*MIX, '-', color=color )
        # melt fraction
        #xx, yy = fig_o.get_xy_data( 'phi', time )
        #ax2.plot( xx, yy, '-', color=color )


    xticks = [0,50,100,135]

    # titles and axes labels, legends, etc
    title = '(a) Regime from dynamics'
    yticks = [1.0,2.0]
    yvals = ['liquid','solid']
    fig_o.set_myaxes( ax0, title=title, xlabel='$P$ (GPa)', ylabel='regime', xticks=xticks, yticks=yticks )
    ax0.set_yticklabels(yvals)
    ax0.set_ylim([0.9,2.1])
    #ax0.yaxis.set_label_coords(-0.075,0.59)

    title = '(b) Regime from melt fraction'
    yticks = [1.0,2.0]
    yvals = ['liquid','solid']
    fig_o.set_myaxes( ax1, title=title, xlabel='$P$ (GPa)', ylabel='regime', xticks=xticks, yticks=yticks )
    ax1.set_yticklabels(yvals)
    ax1.set_ylim([0.9,2.1])
    #ax1.yaxis.set_label_coords(-0.075,0.59)

    fig_o.savefig(7)

#====================================================================
def figure8():

    width = 4.7747/2*3
    height = 4.7747/2
    fig_o = su.FigureData( 1, 3, width, height, 'rheological_transition' )

    ax0 = fig_o.ax[0]
    ax1 = fig_o.ax[1]
    ax2 = fig_o.ax[2]

    handle_l = []

    fig_o.time = su.get_all_output_times()

    data_l = []

    time = fig_o.time[0]
    myjson_o = su.MyJSON( 'output/{}.json'.format(time) )
    xx_pres_b = myjson_o.get_scaled_field_values('pressure_b')
    xx_pres_b *= 1.0E-9

    for nn, time in enumerate( fig_o.time ):

        print('time=',time)
        myjson_o = su.MyJSON( 'output/{}.json'.format(time) )

        yy1 = myjson_o.get_scaled_field_values( 'regime_b' )
        yy2 = 1.5

        # always need temperature
        y_temp = myjson_o.get_scaled_field_values( 'temp_b' )

        # find the pressure(s) that correspond(s) to the cross-over point
        pres = su.find_xx_for_yy( xx_pres_b, yy1, yy2 )
        ind = su.get_first_non_zero_index( pres )

        if ind is not None:
            data_l.append( (xx_pres_b[ind], y_temp[ind], time) )

        # break loop once the rheological transition has reached
        # the surface
        if ind == 1:
            break

    # reshape data
    # column order: pressure, temperature, time
    data_a = np.reshape( data_l, (-1,3) )

    # set time offset to zero, such that t0 defines the time at which the
    # rheological transition is at the CMB
    data_a[:,2] -= data_a[0,2]

    polydegree = 3 # hard-coded

    # just fit the part where the rheological transition is advancing
    # through the magma ocean
    pres_coeff = np.polyfit( data_a[:,2], data_a[:,0], polydegree ) # pressure
    pres_poly = np.poly1d( pres_coeff )
    pres_fit = pres_poly( data_a[:,2] )
    temp_coeff = np.polyfit( data_a[:,2], data_a[:,1], polydegree ) # temperature
    temp_poly = np.poly1d( temp_coeff )
    temp_fit = temp_poly( data_a[:,2] )

    tmax = np.max( data_a[:,2] )
    tmin = np.min( data_a[:,2] )

    vmax = np.max( data_a[:,2] )
    vmin = np.min( data_a[:,2] )

    print( '---------------------------------------------------' )
    print( 'rheological transition advancing through the mantle' )
    print( '---------------------------------------------------' )
    print( 'fitting coefficients for %(polydegree)sth order polynomial' % vars() )
    print( 'between t_min= %(tmin)s yrs and t_max= %(tmax)s yrs' % vars() )
    print( 'pres_coeff=', pres_coeff )
    print( 'temp_coeff=', temp_coeff )

    fig_o.ax[0].scatter( data_a[:,2], data_a[:,0], c=data_a[:,2], cmap='inferno', vmin=vmin, vmax=vmax, s=3.0  )
    h1, = fig_o.ax[0].plot( data_a[:,2], pres_fit, linestyle='-', color='green', label=r'Fit' )

    fig_o.ax[1].scatter( data_a[:,2], data_a[:,1], c=data_a[:,2], cmap='inferno', vmin=vmin, vmax=vmax, s=3.0 )
    h2, = fig_o.ax[1].plot( data_a[:,2], temp_fit, linestyle='-', color='green', label=r'Fit' )

    handle_scatter = fig_o.ax[2].scatter( data_a[:,0], data_a[:,1], c=data_a[:,2], cmap='inferno', vmin=vmin, vmax=vmax, s=3.0 )
    handle_fit, = fig_o.ax[2].plot( pres_fit, temp_fit, color='green', linestyle='-', label=r'Fit' )

    hf_l = [handle_fit]

    # titles and axes labels, legends, etc
    title = r'$P(t)$, GPa'
    fig_o.set_myaxes( ax0, title=title, ylabel='$P$', xlabel='Time (yrs)' )
    ax0.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    title = r'$T(t)$, K'
    fig_o.set_myaxes( ax1, title=title, ylabel='$T$', xlabel='Time (yrs)' )
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    title = r'$P-T-t$'
    xticks = [0,50,100,135]
    yticks = [1000,2000,3000,4000,5000]
    fig_o.set_myaxes( ax2, title=title, ylabel='$T$', yticks=yticks, xlabel='$P$ (GPa)', xticks=xticks )
    #ax2.yaxis.set_label_coords(-0.15,0.575)
    fig_o.set_mylegend( ax0, hf_l, loc='upper right', ncol=1, TITLE="" )
    fig_o.set_mylegend( ax1, hf_l, loc='upper right', ncol=1, TITLE="" )
    fig_o.set_mylegend( ax2, hf_l, loc='upper right', ncol=1, TITLE="" )
    cbar = fig_o.fig.colorbar( handle_scatter, format='%.0e' )#, ticks=[0,200,400,600,800,1000,2000] )
    cbar.set_label('Time (yrs)')

    fig_o.savefig(8)

#====================================================================
def main():

    # arguments (run with -h to summarize)
    #parser = argparse.ArgumentParser(description='SPIDER plotting script')
    #parser.add_argument('-t', '--times', type=str, help='Comma-separated (no spaces) list of times');
    #args = parser.parse_args()

    #figure7( args.times )
    figure8()
    plt.show()

#====================================================================

if __name__ == "__main__":

    main()
