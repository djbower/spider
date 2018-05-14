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
def get_mo_PT_conditions( fig_o, time_l, field1, criteria, label="", color='black' ):

    data_l = []

    time = fig_o.time[0]
    myjson_o = su.MyJSON( 'output/{}.json'.format(time) )
    xx_pres_b = myjson_o.get_scaled_field_values('pressure_b')
    xx_pres_b *= 1.0E-9

    for nn, time in enumerate( time_l ):

        myjson_o = su.MyJSON( 'output/{}.json'.format(time) )

        yy1 = myjson_o.get_scaled_field_values( field1 )
        try:
            yy2 = float(criteria)
        except TypeError:
            yy2 = myjson_o.get_scaled_field_values( criteria )

        # always need temperature
        y_temp = myjson_o.get_scaled_field_values( 'temp_b' )

        # find the pressure(s) that correspond(s) to the cross-over point
        pres = su.find_xx_for_yy( xx_pres_b, yy1, yy2 )
        ind = su.get_first_non_zero_index( pres )

        current_pres = xx_pres_b[-1]
        current_temp = y_temp[-1]
        try:
            previous_pres = data_l[-1][0]
            if ind is None and previous_pres > 120:
                data_l.append( (current_pres, current_temp, time ) )
            elif ind is None and previous_pres < 20:
                current_pres = xx_pres_b[0]
                current_temp = y_temp[0]
                data_l.append( (current_pres, current_temp, time ) )
            #if ind is None:
            #    pass
            else:
                current_pres = pres[ind]
                current_temp = y_temp[ind]
                if current_pres < previous_pres:
                    data_l.append( (current_pres, current_temp, time ) )
        except IndexError:
            #if ind is not None:
            #    current_pres = pres[ind]
            #    current_temp = y_temp[ind]
            #    data_l.append( (current_pres, current_temp, time ) )
            data_l.append( (current_pres, current_temp, time ) )

    # reshape data
    # column order: pressure, temperature, time
    data_a = np.reshape( data_l, (-1,3) )
    print(data_a)

    polydegree = 3 # hard-coded

    # just fit the part where the rheological transition is advancing
    # through the magma ocean
    # compute the start and end indices for this
    diff_a = np.diff( data_a[:,0] )
    sin = su.get_first_non_zero_index( diff_a )
    ein = 1 + len( diff_a ) - su.get_first_non_zero_index( diff_a[::-1] )


    pres_coeff = np.polyfit( data_a[:,2][sin:ein], data_a[:,0][sin:ein], polydegree ) # pressure
    pres_poly = np.poly1d( pres_coeff )
    pres_fit = pres_poly( data_a[:,2][sin:ein] )
    temp_coeff = np.polyfit( data_a[:,2][sin:ein], data_a[:,1][sin:ein], polydegree ) # temperature
    temp_poly = np.poly1d( temp_coeff )
    temp_fit = temp_poly( data_a[:,2][sin:ein] )

    tmax = np.max( data_a[:,2][sin:ein] )
    tmin = np.min( data_a[:,2][sin:ein] )

    vmax = np.max( data_a[:,2] )
    vmin = np.min( data_a[:,2] )

    print( '---------------------------------------------------' )
    print( 'rheological transition advancing through the mantle' )
    print( '---------------------------------------------------' )
    print( 'fitting coefficients for %(polydegree)sth order polynomial' % vars() )
    print( 'between t_min= %(tmin)s yrs and t_max= %(tmax)s yrs' % vars() )
    print( 'pres_coeff=', pres_coeff )
    print( 'temp_coeff=', temp_coeff )

    # now fit cooling at 135 GPa fixed pressure
    temp_coeff2 = np.polyfit( data_a[:,2][0:sin+1], data_a[:,1][0:sin+1], polydegree )
    temp_poly2 = np.poly1d( temp_coeff2 )
    temp_fit2 = temp_poly2( data_a[:,2][0:sin+1] )

    print( '------------------------------------------------' )
    print( 'cooling before rheological transition is reached' )
    print( '------------------------------------------------' )
    print( 'temp_coeff2=', temp_coeff2 )

    handle_scatter = fig_o.ax[1].scatter( data_a[:,0], data_a[:,1], c=data_a[:,2], cmap='inferno', vmin=vmin, vmax=vmax )
    handle_fit, = fig_o.ax[1].plot( pres_fit, temp_fit, color=color, label=label )

    handle_fit2 = fig_o.ax[0].plot( data_a[:,2][0:sin+1], temp_fit2, 'k-' )

    return ( handle_scatter, handle_fit )

#====================================================================
def figure8():

    width = 4.7747#/2
    height = 4.7747/2
    fig_o = su.FigureData( 1, 2, width, height )

    ax0 = fig_o.ax[0]
    ax1 = fig_o.ax[1]

    handle_l = []

    fig_o.time = su.get_all_output_times()

    hs1, hf1 = get_mo_PT_conditions( fig_o, fig_o.time, 'regime_b', 1.5, r'Polynomial fit', 'blue' )
    #hs2, hf2 = get_mo_PT_conditions( fig_o, fig_o.time, 'phi', 0.4, r'$\phi=0.4$', 'red' )

    hf_l = [hf1]#[hf1, hf2]

    #for i, txt in enumerate( liquidus_l[2] ):
    #    if txt is not None:
    #        xoff = liquidus_l[0][i]+1
    #        yoff = liquidus_l[1][i]-50
    #        ax0.annotate( txt, xy=(liquidus_l[0][i], liquidus_l[1][i]), xytext=(xoff,yoff) )


    # titles and axes labels, legends, etc
    title = r'Global magma ocean ($P$=135 GPa)'
    fig_o.set_myaxes( ax0, title=title, ylabel='$T$ (K)', xlabel='Time after $t_0$ (yrs)' )
    ax0.yaxis.set_label_coords(-0.15,0.575)

    title = r'Partially molten magma ocean'
    xticks = [0,50,100,135]
    yticks = [1000,2000,3000,4000,5000]

    fig_o.set_myaxes( ax1, title=title, ylabel='$T$ (K)', yticks=yticks, xlabel='$P$ (GPa)', xticks=xticks )
    ax1.yaxis.set_label_coords(-0.15,0.575)

    fig_o.set_mylegend( ax1, hf_l, loc='upper left', ncol=1, TITLE="" )
    cbar = fig_o.fig.colorbar( hs1 )#, ticks=[0,200,400,600,800,1000,2000] )
    cbar.set_label('Time after $t_0$ (yrs)')

    #cbar.ax.set_yticklabels(['0','0.2 kyr','0.4 kyr','0.6 kyr','0.8 kyr','1 kyr'])
    #cbar.ax.set_ylim([0,1000])

    plt.show()

    #ax0.set_xlim([0,2000])

    fig_o.savefig(8)

#====================================================================
def main():

    # arguments (run with -h to summarize)
    parser = argparse.ArgumentParser(description='SPIDER plotting script')
    parser.add_argument('-t', '--times', type=str, help='Comma-separated (no spaces) list of times');
    args = parser.parse_args()

    figure7( args.times )
    figure8()
    plt.show()

#====================================================================

if __name__ == "__main__":

    main()
