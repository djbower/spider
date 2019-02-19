#!/usr/bin/env python

import spider_utils as su
import matplotlib.pyplot as plt
import numpy as np
import argparse

#====================================================================
def plot_radius_evolution():

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
def plot_interior_atmosphere( times ):

    figw = 4.7747 / 2.0
    figh = 4.7747

    fig_o = su.FigureData( 2, 1, figw, figh, 'interior_atmosphere', times, units='Myr' )
    plt.subplots_adjust(hspace=0.0)

    ax0 = fig_o.ax[0]
    ax1 = fig_o.ax[1]
    ax0.axes.xaxis.set_visible(False)

    time = fig_o.time[0] # first timestep since liquidus and solidus
                         # are time-independent

    myjson_o = su.MyJSON( 'output/{}.json'.format(time) )

    TIMEYRS = myjson_o.data_d['nstep']

    xx_pres = myjson_o.get_dict_values_internal(['data','pressure_b'])
    xx_pres *= 1.0E-9

    xx_pres_s = myjson_o.get_dict_values(['data','pressure_s'])
    xx_pres_s *= 1.0E-9

    # shade grey between liquidus and solidus
    yy_liqt = myjson_o.get_dict_values_internal(['data','liquidus_temp_b'])
    yy_solt = myjson_o.get_dict_values_internal(['data','solidus_temp_b'])
    ax1.fill_betweenx( xx_pres, yy_liqt, yy_solt, facecolor='grey', alpha=0.35, linewidth=0 )

    # dotted lines of constant melt fraction
    #for xx in range( 0, 11, 2 ):
    #    yy_b = xx/10.0 * (yy_liq - yy_sol) + yy_sol
    #    if xx == 0:
    #        # solidus
    #        ax0.plot( yy_b, xx_pres, '-', linewidth=0.5, color='black' )
    #    elif xx == 3:
    #        # typically, the approximate location of the rheological transition
    #        ax0.plot( yy_b, xx_pres, '-', linewidth=1.0, color='white')
    #    elif xx == 10:
    #        # liquidus
    #        ax0.plot( yy_b, xx_pres, '-', linewidth=0.5, color='black' )
    #    else:
    #        # dashed constant melt fraction lines
    #        ax0.plot( yy_b, xx_pres, '--', linewidth=1.0, color='white' )

    handle_l = [] # handles for legend

    for nn, time in enumerate( fig_o.time ):
        # read json
        myjson_o = su.MyJSON( 'output/{}.json'.format(time) )

        color = fig_o.get_color( nn )

        # atmosphere structure
        hatm_interp1d = myjson_o.get_atm_struct_depth_interp1d()
        height_point = hatm_interp1d( 10.0E-3 ) / 1.0E3 # 10 mb
        tatm_interp1d = myjson_o.get_atm_struct_temp_interp1d()
        temp_point = tatm_interp1d( 10.0E-3 ) # 10 mb
        ax0.plot( temp_point, height_point, color=color, marker='o' )

        atmos_pres_a = myjson_o.get_dict_values( ['atmosphere','atm_struct_pressure'] )
        atmos_temp_a = myjson_o.get_dict_values( ['atmosphere','atm_struct_temp'] )
        atmos_height_a = myjson_o.get_dict_values( ['atmosphere','atm_struct_depth'] )
        atmos_height_a *= 1.0E-3
        indices = atmos_height_a < height_point
        atmos_height_a = atmos_height_a[indices]
        atmos_temp_a = atmos_temp_a[indices]
        ax0.plot( atmos_temp_a, atmos_height_a, '-', color=color )

        # use melt fraction to determine mixed region
        MIX = myjson_o.get_mixed_phase_boolean_array( 'basic_internal' )
        label = fig_o.get_legend_label( time )
        yy = myjson_o.get_dict_values_internal(['data','temp_b'])
        ax1.plot( yy, xx_pres, '--', color=color )
        handle, = ax1.plot( yy*MIX, xx_pres*MIX, '-', color=color, label=label )
        handle_l.append( handle )

    yticks = [10,50,100,135]
    ymax = 138
    xticks = [200,1000,2000,3000,4000,5000]

    title = 'Coupled interior-atmosphere evolution'
    ylabel = 'Height (km)'
    fig_o.set_myaxes( ax0, ylabel=ylabel, title=title, yrotation=90 )

    # titles and axes labels, legends, etc
    units = myjson_o.get_dict_units(['data','temp_b'])
    #title = '(b) Mantle temperature, {}'.format(units)

    ylabel = 'Pressure (GPa)'
    fig_o.set_myaxes( ax1, ylabel=ylabel, xlabel='$T$ (K)', xticks=xticks, ymin=0.0, ymax=ymax, yticks=yticks, yrotation=90 )
    ax1.yaxis.set_label_coords(-0.2,0.5)

    ax0.set_xlim( xticks[0], xticks[-1])
    ax0.set_ylim( 0.0, 300.0 )
    ax1.set_ylim( 0.0, ymax )

    ax1.invert_yaxis()

    fig_o.set_mylegend( ax0, handle_l, loc='upper right' )

    fig_o.savefig(2)

#====================================================================

if __name__ == "__main__":

    # arguments (run with -h to summarize)
    parser = argparse.ArgumentParser(description='SPIDER plotting script')
    parser.add_argument('-t', '--times', type=str, help='Comma-separated (no spaces) list of times');
    args = parser.parse_args()

    if not args.times:
        logger.critical( 'You must specify times in a comma-separated list (no spaces) using -t' )
        sys.exit(0)

    #plot_radius_evolution()
    plot_interior_atmosphere( args.times )
    plt.show()
