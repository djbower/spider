#!/usr/bin/env python

import spider_utils as su
import logging
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import numpy as np
import argparse
from scipy.optimize import curve_fit
import os, sys

#====================================================================
def figure9():

    #---------------
    # set parameters
    #---------------
    # flag to enable data output
    EXPORT = True
    DYNAMIC = True # dynamic criterion, otherwise critical melt fraction
    # pressure cut-offs for curve fitting and plotting
    #PMIN = 1.0 # GPa # TODO: CURRENTLY NOT USED
    #PMAX = 135.0 # GPa # TODO: CURRENTLY NOT USED
    TPRIME1_PRESSURE = 15.0 # pressure at which rheological transition is assumed to be at the surface
    XMAXFRAC = 1.05 # plot maximum x value this factor beyond TPRIME1
    # currently not used RADFIT = 1 # fit using radiative cooling model (n=4), otherwise use linear fit
    markersize = 3.0
    #---------------

    name = 'rheological_front_'
    if DYNAMIC:
        name += 'dynamic'
    else:
        name += 'phi'

    width = 4.7747
    height = 4.7747
    fig_o = su.FigureData( 2, 2, width, height, name )
    fig_o.fig.subplots_adjust(wspace=0.4,hspace=0.4)
    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]
    ax2 = fig_o.ax[1][0]
    ax3 = fig_o.ax[1][1]

    fig_o.time = su.get_all_output_times()
    time_yrs = fig_o.time
    time0 = fig_o.time[0]

    # pressure in GPa
    myjson_o = su.MyJSON( 'output/{}.json'.format(time0) )
    xx_pres_b = myjson_o.get_dict_values_internal( ('data','pressure_b') )
    xx_pres_b *= 1.0E-9 # to GPa

    data_l = []
    handle_l = []

    # by asking for all the data in one go, the code is much faster
    keys_t = ( ('rheological_front_phi','depth'),
               ('rheological_front_phi', 'pressure'),
               ('rheological_front_phi', 'temperature'),
               ('rheological_front_phi', 'below_mass_avg', 'pressure' ),
               ('rheological_front_phi', 'below_mass_avg', 'temperature' ),
               ('rheological_front_phi', 'above_mass_avg', 'pressure' ),
               ('rheological_front_phi', 'above_mass_avg', 'temperature' ),
               ('rheological_front_dynamic', 'depth'),
               ('rheological_front_dynamic', 'pressure'),
               ('rheological_front_dynamic', 'temperature'),
               ('rheological_front_dynamic', 'below_mass_avg', 'pressure' ),
               ('rheological_front_dynamic', 'below_mass_avg', 'temperature' ),
               ('rheological_front_dynamic', 'above_mass_avg', 'pressure' ),
               ('rheological_front_dynamic', 'above_mass_avg', 'temperature' ),
               ('atmosphere','temperature_surface'),
               ('data','Jtot_b') )
    data_a = su.get_dict_surface_values_for_times( keys_t, fig_o.time )

    temperature_surface = data_a[14,:]
    flux_surface = data_a[15,:]

    if DYNAMIC:
        # rheological front based on dynamic criterion
        rf_depth = data_a[7,:]
        rf_pressure = data_a[8,:] * 1.0E-9 # to GPa
        rf_temperature = data_a[9,:]
        mo_pressure = data_a[10,:] * 1.0E-9 # to GPa
        mo_temperature = data_a[11,:]
        so_pressure = data_a[12,:] * 1.0E-9 # to GPa
        so_temperature = data_a[13:,]
    else:
        # rheological front based on critical melt fraction
        rf_depth = data_a[0,:]
        rf_pressure = data_a[1,:] * 1.0E-9 # to GPa
        rf_temperature = data_a[2,:]
        mo_pressure = data_a[3,:] * 1.0E-9 # to GPa
        mo_temperature = data_a[4,:]
        so_pressure = data_a[5,:] * 1.0E-9 # to GPa
        so_temperature = data_a[6:,]

    # compute start time and end time of rheological front moving through the mantle
    (t0,n0,t1,n1) = get_t0_t1( fig_o.time, rf_pressure, xx_pres_b[-1], TPRIME1_PRESSURE )

    time_yrs_shift = fig_o.time[n0:n1] - t0 # shift time for plotting Pf and Tf
    rf_pressure_shift = rf_pressure[n0:n1]
    rf_temperature_shift = rf_temperature[n0:n1]
    mo_pressure_shift = mo_pressure[n0:n1]
    mo_temperature_shift = mo_temperature[n0:n1]
    so_pressure_shift = so_pressure[n0:n1]
    so_temperature_shift = so_temperature[n0:n1]

    print( so_temperature )

    # ---------------------------
    # ---- surface heat flux ----
    trans = transforms.blended_transform_factory(ax0.transData, ax0.transAxes)
    h1, = ax0.semilogy( fig_o.time, flux_surface, marker='o', markersize=markersize, color='0.8' )
    ax0.axvline( t0, ymin=0.1, ymax=0.8, color='0.25', linestyle=':')
    ax0.text( t0, 0.82, '$t^\prime_0$', ha='right', va='bottom', transform = trans )
    ax0.axvline( t1, ymin=0.1, ymax=0.8, color='0.25', linestyle=':')
    ax0.text( t1, 0.82, '$t^\prime_1$', ha='right', va='bottom', transform = trans )
    title = r'(a) $q_0(t)$, W/m$^2$'

    xmax = XMAXFRAC*t1
    ydata_plot = flux_surface[np.where(time_yrs<xmax)]
    yticks = get_yticks( ydata_plot, 1.0E4 )
    #yticks = get_yticks( ydata_plot, 1.0E5 )

    fig_o.set_myaxes( ax0, title=title, ylabel='$q_0$', xlabel='$t$ (yrs)' )
    ax0.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax0.yaxis.set_label_coords(-0.15,0.5)
    ax0.set_xlim( None, xmax )
    ax0.set_ylim( yticks[0], yticks[-1] )

    # --- end of surface heat flux ----
    # --------------------------------- 

    # -----------------------------
    # ---- surface temperature ----
    trans = transforms.blended_transform_factory(ax1.transData, ax1.transAxes)
    ax1.plot( fig_o.time, temperature_surface, marker='o', markersize=markersize, color='0.8' )
    ax1.axvline( t0, ymin=0.1, ymax=0.8, color='0.25', linestyle=':')
    ax1.text( t0, 0.82, '$t^\prime_0$', ha='right', va='bottom', transform = trans )
    ax1.axvline( t1, ymin=0.1, ymax=0.8, color='0.25', linestyle=':')
    ax1.text( t1, 0.82, '$t^\prime_1$', ha='right', va='bottom', transform = trans )
    title = r'(b) $Ts(t)$, K'

    xmax = XMAXFRAC*t1
    ydata_plot = temperature_surface[np.where(time_yrs<xmax)]
    #yticks = [1300,1400,1500,1600,1700]
    yticks = get_yticks( ydata_plot, 400 ) # [1500,2000,2500,3000,3500,4000]

    fig_o.set_myaxes( ax1, title=title, ylabel='$Ts$', xlabel='$t$ (yrs)', yticks=yticks )
    #ax1.set_ylim( [1300,1700])
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    #ax1.yaxis.set_label_coords(-0.2,0.575)
    ax1.yaxis.set_label_coords(-0.15,0.5)
    ax1.set_xlim( None, xmax )
    # ---- end surface temperature ----
    # ---------------------------------

    # ---------------------------
    # ---- rheological front ----
    #if EXPORT:
    #    np.savetxt( 'Pr.dat', np.column_stack((time_a,pres_a)))
    #    np.savetxt( 'Tr.dat', np.column_stack((time_a,temp_a)))
    ax2.plot( time_yrs_shift, rf_pressure_shift, marker='o', markersize=markersize, color='green' )
    ax2.plot( time_yrs_shift, mo_pressure_shift, marker='o', markersize=markersize, color='red' )
    ax2.plot( time_yrs_shift, so_pressure_shift, marker='o', markersize=markersize, color='blue' )
    title = r'(c) $P_f(t^\prime)$, GPa'
    yticks = [0,50,100,150]
    fig_o.set_myaxes( ax2, title=title, ylabel='$P_f$', xlabel='$t^\prime$ (yrs)', yticks=yticks )
    ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax2.yaxis.set_label_coords(-0.15,0.5)
    ax2.set_xlim( (0, t1-t0) )

    ax3.plot( time_yrs_shift, rf_temperature_shift, marker='o', markersize=markersize, color='green' )
    ax3.plot( time_yrs_shift, mo_temperature_shift, marker='o', markersize=markersize, color='red' )
    ax3.plot( time_yrs_shift, so_temperature_shift, marker='o', markersize=markersize, color='blue' )
    title = r'(d) $T_f(t^\prime)$, K'
    yticks = [1500,2500,3500,4500]
    fig_o.set_myaxes( ax3, title=title, ylabel='$T_f$', xlabel='$t^\prime$ (yrs)', yticks=yticks )
    ax3.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax3.set_ylim( [1500,4500] )
    #ax3.yaxis.set_label_coords(-0.15,0.575)
    ax3.yaxis.set_label_coords(-0.15,0.5)
    ax3.set_xlim( (0,t1-t0) )

    ###################
    ##### fitting #####
    ###################

    if DYNAMIC:
        deg = 2
    else:
        deg = 1

    # linear fit to pressure, using a polynomial
    # TODO: fix P_cmb, and then just fit gradient
    # TODO: 2nd order fit is better for dynamic criteria - but why?
    xx = time_yrs_shift # data_a[tprime0ii:tprime1ii,0]-tprime0
    yy_pres = mo_pressure_shift #data_a[tprime0ii:tprime1ii,2]
    rheological_pres_poly = np.polyfit( xx, yy_pres, deg )
    rheological_pres_o = np.poly1d( rheological_pres_poly )
    rheological_pres_fit = rheological_pres_o( xx )
    Rval = Rsquared( yy_pres, rheological_pres_fit )
    print('----pressure fit----')
    print( 'R^2 for pressure=', Rval )
    print( 'poly for pressure=', rheological_pres_poly )
    print( 'gradient=', rheological_pres_poly[0] )
    print( 'intercept=', rheological_pres_poly[1] )
    print( 'minimum for fit=', np.min(rheological_pres_fit) )
    print( 'maximum for fit=', np.max(rheological_pres_fit) )
    ax2.plot( xx, rheological_pres_fit )

    # linear fit to temperature, using a polynomial
    # TODO: fix Tc_cmb, and then just fit gradient
    yy_temp = mo_temperature_shift # data_a[tprime0ii:tprime1ii,3]
    rheological_temp_poly = np.polyfit( xx, yy_temp, deg )
    rheological_temp_o = np.poly1d( rheological_temp_poly )
    rheological_temp_fit = rheological_temp_o( xx )
    Rval = Rsquared( yy_temp, rheological_temp_fit )
    print('----temperature fit----')
    print( 'R^2 for temperature=', Rval )
    print( 'poly for temperature=', rheological_temp_poly )
    print( 'gradient=', rheological_temp_poly[0] )
    print( 'intercept=', rheological_pres_poly[1] )
    print( 'minimum for fit=', np.min(rheological_temp_fit) )
    print( 'maximum for fit=', np.max(rheological_temp_fit) )
    #alpha = - rheological_temp_poly[1] / rheological_pres_poly[1]
    #print( 'alpha=', alpha )
    ax3.plot( xx, rheological_temp_fit )

    fig_o.savefig(9)

#====================================================================
def get_t0_t1( time_l, pres_l, val_t0, val_t1 ):

    '''return absolute time and index when the rheological front
       leaves the CMB and arrives at the surface, as determined by
       the pressures val_t0 and val_t1, respectively'''

    FLAG_t0 = False
    FLAG_t1 = False

    for nn, time in enumerate(time_l):
        pres = pres_l[nn]
        if pres < val_t0 and not FLAG_t0:
            t0 = time
            n0 = nn
            FLAG_t0 = True
        if pres < val_t1 and not FLAG_t1:
            t1 = time
            n1 = nn
            FLAG_t1 = True

    print( 't0= ', t0 )
    print( 't1= ', t1 )

    return (t0,n0,t1,n1)

#====================================================================
def get_yticks( array_in, roundtonearest=500 ):

    # set the range of the y axis and determine the y labels based
    # on the data extent in the input array

    minimum = np.min( array_in )
    maximum = np.max( array_in )

    minimum = np.floor(minimum/roundtonearest) * roundtonearest
    maximum = np.ceil(maximum/roundtonearest) * roundtonearest

    labels = np.arange( minimum, maximum+1, roundtonearest )
    labels = labels.tolist()

    return labels

#====================================================================
def radiative_function_fit( x, *popt ):

    return (popt[0]*x+popt[1])**(-1.0/3.0) + popt[2]

#====================================================================
def Rsquared( ydata, fmodel ):

    # sum of squares of residuals
    ybar = np.mean( ydata )
    #print( 'ybar=', ybar )
    SSres = np.sum( np.square(ydata-fmodel) )
    #print( 'SSres=', SSres )
    SStot = np.sum( np.square(ydata-ybar) )
    #print( 'SStot=', SStot )
    Rsquared = 1.0 - SSres/SStot
    #print( 'Rsquared=', Rsquared )

    return Rsquared

#====================================================================
def main():

    figure9()
    plt.show()

#====================================================================

if __name__ == "__main__":

    main()
