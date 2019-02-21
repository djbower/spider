#!/usr/bin/env python

import spider_utils as su
import matplotlib.pyplot as plt
import numpy as np
import argparse
from scipy.interpolate import interp1d

#====================================================================
def plot_radius_evolution( times ):

    figw = 4.7747/2.0
    figh = 4.7747/2.0

    time_l = [0,1,2]
    fig_o = su.FigureData( 1, 1, figw, figh, 'radius', times )

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

    color0 = fig_o.get_color( 0 )
    color1 = fig_o.get_color( 1 )
    color2 = fig_o.get_color( 2 )

    h1, = ax0.plot( time, radius, label='Surface', color=color0 )
    h2, = ax0.plot( time, radius_10mb, label='10 mb', color=color1 )
    h3, = ax0.plot( time, radius_1mb, label='1 mb', color=color2 )
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

    time = fig_o.time[2] # first timestep since liquidus and solidus
                         # are time-independent

    myjson_o = su.MyJSON( 'output/{}.json'.format(time) )

    TIMEYRS = myjson_o.data_d['nstep']

    xx_pres = myjson_o.get_dict_values(['data','pressure_b'])
    xx_pres *= 1.0E-9
    xx_pres_s = myjson_o.get_dict_values(['data','pressure_s'])
    xx_pres_s *= 1.0E-9
    xx_radius = myjson_o.get_dict_values_internal(['data','radius_b'])
    xx_radius *= 1.0E-3
    xx_depth = xx_radius[0] - xx_radius
    xx_radius_s = myjson_o.get_dict_values(['data','radius_s'])
    xx_radius_s *= 1.0E-3
    xx_depth_s = xx_radius_s[0] - xx_radius_s

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
        ax0.plot( temp_point, height_point, color=color, marker='_', markersize=8 )

        atmos_pres_a = myjson_o.get_dict_values( ['atmosphere','atm_struct_pressure'] )
        atmos_temp_a = myjson_o.get_dict_values( ['atmosphere','atm_struct_temp'] )
        atmos_height_a = myjson_o.get_dict_values( ['atmosphere','atm_struct_depth'] )
        atmos_height_a *= 1.0E-3
        indices = atmos_height_a < height_point
        atmos_height_a = atmos_height_a[indices]
        atmos_temp_a = atmos_temp_a[indices]
        ax0.plot( atmos_temp_a, atmos_height_a, '-', color=color )

        # use melt fraction to determine mixed region
        rho1D_o = myjson_o.get_rho_interp1d()
        temp1D_o = myjson_o.get_temp_interp1d()
        radius = su.solve_for_planetary_radius( rho1D_o )
        myargs = su.get_myargs_static_structure( rho1D_o )
        z = su.get_static_structure_for_radius( radius, *myargs )
        radius_a = su.get_radius_array_static_structure( radius, *myargs )
        pressure_a = z[:,0]
        temp_a = temp1D_o( pressure_a )
        label = fig_o.get_legend_label( time )
        depth_a = radius_a[0] - radius_a
        depth_a *= 1.0e-3 # to km
        handle, = ax1.plot( temp_a, depth_a, color=color, label=label )
        ax1.plot( temp_a[-1], depth_a[-1], color=color, marker='_', markersize=8 )

        # shade grey between liquidus and solidus
        yy_liqt = myjson_o.get_dict_values(['data','liquidus_temp_b'])
        liqt_interp1d = interp1d( xx_pres, yy_liqt, kind='linear', fill_value='extrapolate' )
        yy_solt = myjson_o.get_dict_values(['data','solidus_temp_b'])
        solt_interp1d = interp1d( xx_pres, yy_solt, kind='linear', fill_value='extrapolate' )
        liqt = liqt_interp1d( pressure_a * 1.0E-9 )
        solt = solt_interp1d( pressure_a * 1.0E-9 )
        if nn==3: 
            ax1.fill_betweenx( depth_a, liqt, solt, facecolor='gray', alpha=0.3, linewidth=0 )

        #MIX = myjson_o.get_mixed_phase_boolean_array( 'basic_internal' )
        #label = fig_o.get_legend_label( time )
        #yy = myjson_o.get_dict_values(['data','temp_s'])
        #handle, = ax1.plot( yy, xx_depth_s, '-', color=color, label=label )

        print( 'atmos_temp_a[-1]=', atmos_temp_a[-1], 'temp_a[0]=', temp_a[0] )

        #ax1.plot( yy, xx_pres, '--', color=color )
        #handle, = ax1.plot( yy*MIX, xx_pres*MIX, '-', color=color, label=label )
        handle_l.append( handle )

    #yticks = [10,50,100,135]
    #ymax = 138
    yticks = [0,1E3,2E3,3E3,3500]
    ymin = 0.0
    ymax = 3500.0
    xticks = [200,1000,2000,3000,4000,5000]

    title = 'Coupled interior-atmosphere evolution'
    ylabel = 'Atmosphere height (km)'
    fig_o.set_myaxes( ax0, ylabel=ylabel, title=title, yrotation=90 )

    # titles and axes labels, legends, etc
    units = myjson_o.get_dict_units(['data','temp_b'])
    #title = '(b) Mantle temperature, {}'.format(units)

    #ylabel = 'Pressure (GPa)'
    ylabel= 'Mantle depth (km)'
    fig_o.set_myaxes( ax1, ylabel=ylabel, xlabel='$T$ (K)', xticks=xticks, ymin=ymin, ymax=ymax, yticks=yticks, yrotation=90 )
    ax1.yaxis.set_label_coords(-0.2,0.5)
    ax1.set_ylim( ymin, ymax )
    ax1.invert_yaxis()

    ax0.set_xlim( xticks[0], xticks[-1])
    ax0.set_ylim( 0.0, 250.0 )

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

    #plot_radius_evolution( args.times )
    plot_interior_atmosphere( args.times )
    plt.show()
