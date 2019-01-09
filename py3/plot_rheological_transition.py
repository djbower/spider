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
    # pressure cut-offs for curve fitting and plotting
    #PMIN = 1.0 # GPa # TODO: CURRENTLY NOT USED
    #PMAX = 135.0 # GPa # TODO: CURRENTLY NOT USED
    PHI = 0.4 #-1 #0.4 # melt fraction contour, or -1 for dynamic
    TPRIME1_PRESSURE = 15.0 # pressure at which rheological transition is assumed to be at the surface
    XMAXFRAC = 1.05 # plot maximum x value this factor beyond TPRIME1
    # currently not used RADFIT = 1 # fit using radiative cooling model (n=4), otherwise use linear fit
    #---------------

    width = 4.7747
    height = 4.7747
    fig_o = su.FigureData( 2, 2, width, height, 'rheological_front' )
    fig_o.fig.subplots_adjust(wspace=0.4,hspace=0.4)
    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]
    ax2 = fig_o.ax[1][0]
    ax3 = fig_o.ax[1][1]

    fig_o.time = su.get_all_output_times()

    # pressure in GPa
    time = fig_o.time[0]
    myjson_o = su.MyJSON( 'output/{}.json'.format(time) )
    xx_pres_b = myjson_o.get_scaled_field_values('pressure_b')
    xx_pres_b *= 1.0E-9 # to GPa

    data_l = []
    handle_l = []

    for ii, time in enumerate( fig_o.time ):

        print('time=',time)
        myjson_o = su.MyJSON( 'output/{}.json'.format(time) )

        if PHI < 0.0:
            # dynamic criteria
            yy1 = myjson_o.get_scaled_field_values( 'regime_b' )
            yy2 = 1.5
            print( 'min(yy1)={}, max(yy1)={}'.format(np.min(yy1), np.max(yy1)))
        else:
            # melt fraction criteria
            yy1 = myjson_o.get_scaled_field_values( 'phi_b' )
            yy2 = PHI

        # mantle temperature (vector)
        mantle_temp = myjson_o.get_scaled_field_values( 'temp_b' )

        # mantle surface temperature (scalar)
        mantle_surface_temp = myjson_o.get_scaled_field_values( 'temperature_surface' )

        # surface flux (scalar)
        mantle_surface_flux = myjson_o.get_scaled_field_values( 'Jtot_b' )[0]

        # find the pressure(s) that correspond(s) to the cross-over point
        pres = su.find_xx_for_yy( xx_pres_b, yy1, yy2 )
        ind = su.get_first_non_zero_index( pres )
        print( 'ind=', ind )
        # if the initial condition is not great, then it can 'accidentally' trigger
        # regimes switches, even though in principle everything super-liquidus should
        # start out in the inviscid regime.  So ignore the first time.
        if int(time) is not 0 and ind is not None: # before rheological transition reaches CMB
            # save the time at which the rheological transition first
            # reaches the CMB
            print( 'ind is not None')
            try:
                tprime0
                print('try tprime0')
                #if tprime0:
                #    print('tprime0 exists')
            except NameError:
                print('in except block')
                print('setting tprime0={}'.format(time))
                tprime0 = time
                print(tprime0)
                tprime0ii = ii
            mantle_rheological_pres = xx_pres_b[ind]
            if mantle_rheological_pres < TPRIME1_PRESSURE:
                try:
                    tprime1
                except NameError:
                    print('setting tprime1={}'.format(time))
                    tprime1 = time
                    tprime1ii = ii
            mantle_rheological_temp = mantle_temp[ind]
        else:
            mantle_rheological_pres = xx_pres_b[-1]
            mantle_rheological_temp = mantle_temp[-1]

        # time is appended twice, because the second time will be corrected for the offset
        data_l.append( ( time, time, mantle_rheological_pres, mantle_rheological_temp,
            mantle_surface_temp, mantle_surface_flux ) )

    if not tprime0:
        print('ERROR: tprime0 not set during time loop')
    if not tprime1:
        print('ERROR: tprime1 not set during time loop')
    if not tprime0 or not tprime1:
        sys.exit(1)

    # reshape into numpy array
    data_a = np.reshape( data_l, (-1,6) )
    # tprime0 defines the time at which the rheological front is at the CMB
    #tprime0 = data_a[0,0]
    data_a[:,1] -= tprime0

    print('----times----')
    print('tprime0=', tprime0 )
    print('tprime1=', tprime1 )
    dt = tprime1 - tprime0
    print('delta time=', dt )

    # ---------------------------
    # ---- surface heat flux ----
    trans = transforms.blended_transform_factory(ax0.transData, ax0.transAxes)
    h1, = ax0.semilogy( data_a[:,0], data_a[:,5], marker='o', markersize=4.0, color='0.8' )
    ax0.axvline( tprime0, ymin=0.1, ymax=0.8, color='0.25', linestyle=':')
    ax0.text( tprime0, 0.82, '$t^\prime_0$', ha='right', va='bottom', transform = trans )
    ax0.axvline( tprime1, ymin=0.1, ymax=0.8, color='0.25', linestyle=':')
    ax0.text( tprime1, 0.82, '$t^\prime_1$', ha='right', va='bottom', transform = trans )
    title = r'(a) $q_0(t)$, W/m$^2$'
    fig_o.set_myaxes( ax0, title=title, ylabel='$q_0$', xlabel='$t$ (yrs)' )
    ax0.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax0.yaxis.set_label_coords(-0.15,0.5)
    xmax = XMAXFRAC*tprime1
    ax0.set_xlim( None, xmax )
    # --- end of surface heat flux ----
    # --------------------------------- 

    # -----------------------------
    # ---- surface temperature ----
    trans = transforms.blended_transform_factory(ax1.transData, ax1.transAxes)
    ax1.plot( data_a[:,0], data_a[:,4], marker='o', markersize=4.0, color='0.8' )
    ax1.axvline( tprime0, ymin=0.1, ymax=0.8, color='0.25', linestyle=':')
    ax1.text( tprime0, 0.82, '$t^\prime_0$', ha='right', va='bottom', transform = trans )
    ax1.axvline( tprime1, ymin=0.1, ymax=0.8, color='0.25', linestyle=':')
    ax1.text( tprime1, 0.82, '$t^\prime_1$', ha='right', va='bottom', transform = trans )
    title = r'(b) $Ts(t)$, K'
    #yticks = [1300,1400,1500,1600,1700]
    yticks = get_ylabels( data_a[:,4] ) # [1500,2000,2500,3000,3500,4000]
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
    ax2.plot( data_a[:,1], data_a[:,2], marker='o', markersize=4.0, color='0.8' )
    title = r'(c) $P_f(t^\prime)$, GPa'
    yticks = [0,50,100,150]
    fig_o.set_myaxes( ax2, title=title, ylabel='$P_f$', xlabel='$t^\prime$ (yrs)', yticks=yticks )
    ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax2.yaxis.set_label_coords(-0.15,0.5)
    ax2.set_xlim( (0, tprime1-tprime0) )

    ax3.plot( data_a[:,1], data_a[:,3], marker='o', markersize=4.0, color='0.8' )
    title = r'(d) $T_f(t^\prime)$, K'
    yticks = [1500,2500,3500,4500]
    fig_o.set_myaxes( ax3, title=title, ylabel='$T_f$', yticks=yticks, xlabel='$t^\prime$ (yrs)' )
    ax3.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax3.set_ylim( [1500,4500] )
    #ax3.yaxis.set_label_coords(-0.15,0.575)
    ax3.yaxis.set_label_coords(-0.15,0.5)
    ax3.set_xlim( (0, tprime1-tprime0) )

    ###################
    ##### fitting #####
    ###################

    # linear fit to pressure, using a polynomial
    # TODO: fix P_cmb, and then just fit gradient
    # TODO: 2nd order fit is better for dynamic criteria - but why?
    deg = 1
    xx = data_a[tprime0ii:tprime1ii,0]-tprime0
    yy_pres = data_a[tprime0ii:tprime1ii,2]
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
    deg = 1
    yy_temp = data_a[tprime0ii:tprime1ii,3]
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
    alpha = - rheological_temp_poly[1] / rheological_pres_poly[1]
    print( 'alpha=', alpha )
    ax3.plot( xx, rheological_temp_fit )

    fig_o.savefig(9)

#====================================================================
def get_ylabels( array_in, roundtonearest=500 ):

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
