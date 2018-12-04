#!/usr/bin/env python

import spider_utils as su
import logging
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import numpy as np
import argparse
from scipy.optimize import curve_fit
import os

#====================================================================
def figure9():

    #---------------
    # set parameters
    #---------------
    # flag to enable data output
    EXPORT = True
    # pressure cut-offs for curve fitting and plotting
    PMIN = 1.0 # GPa
    PMAX = 135.0 # GPa
    PHI = 0.4 # melt fraction contour, or -1 for dynamic
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
        else:
            # melt fraction criteria
            yy1 = myjson_o.get_scaled_field_values( 'phi_b' )
            yy2 = 0.4

        # mantle temperature (vector)
        mantle_temp = myjson_o.get_scaled_field_values( 'temp_b' )

        # mantle surface temperature (scalar)
        mantle_surface_temp = myjson_o.get_scaled_field_values( 'temperature_surface' )

        # surface flux (scalar)
        mantle_surface_flux = myjson_o.get_scaled_field_values( 'Jtot_b' )[0]

        # find the pressure(s) that correspond(s) to the cross-over point
        pres = su.find_xx_for_yy( xx_pres_b, yy1, yy2 )
        ind = su.get_first_non_zero_index( pres )
        if ind is not None:
            mantle_rheological_pres = xx_pres_b[ind]
            mantle_rheological_temp = mantle_temp[ind]
        else:
            mantle_rheological_pres = xx_pres_b[-1]
            mantle_rheological_temp = mantle_temp[-1]

        data_l.append( ( time, time, mantle_rheological_pres, mantle_rheological_temp,
            mantle_surface_temp, mantle_surface_flux ) )

        # 'None' is returned until the rheological transition first
        # reaches the CMB
        #if ind is not None:
            # time is appended twice, because the second time will be
            # corrected for the offset
        #    data_l.append( (time, time, xx_pres_b[ind], y_temp[ind]) )

    # reshape into numpy array
    # columns: time, time-offset, pressure, temperature
    data_a = np.reshape( data_l, (-1,6) )
    # tprime0 defines the time at which the rheological front is at the CMB
    tprime0 = data_a[0,0]
    data_a[:,1] -= tprime0

    # ---------------------------
    # ---- surface heat flux ----
    #surff_a = np.column_stack( (fig_o.time, fluxsurf_l) )
    #surff0_a = surff_a[:,0]
    #surff1_a = surff_a[:,1]

    trans = transforms.blended_transform_factory(ax0.transData, ax0.transAxes)
    h1, = ax0.semilogy( data_a[:,0], data_a[:,5], marker='o', markersize=4.0, color='0.8' )
    #ax2.axvline( tprime0, ymin=0.1, ymax=0.8, color='0.25', linestyle=':')
    #ax2.text( tprime0, 0.82, '$t^\prime_0$', ha='right', va='bottom', transform = trans )
    #ax2.axvline( tprime1abs, ymin=0.1, ymax=0.8, color='0.25', linestyle=':')
    #ax2.text( tprime1abs, 0.82, '$t^\prime_1$', ha='right', va='bottom', transform = trans )
    title = r'(a) $q_0(t)$, W/m$^2$'
    fig_o.set_myaxes( ax0, title=title, ylabel='$q_0$', xlabel='$t$ (yrs)' )
    ax0.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax0.yaxis.set_label_coords(-0.2,0.5)
    # --- end of surface heat flux ----
    # --------------------------------- 

    # -----------------------------
    # ---- surface temperature ----
    #surft_a = np.column_stack( (fig_o.time, tempsurf_l) )
    #ind_tcut = surft_a[:,1] > 1200.0
    #surft0_a = surft_a[:,0]#[ind_tcut] # time
    #surft1_a = surft_a[:,1]#[ind_tcut] # surface temperature

    trans = transforms.blended_transform_factory(ax1.transData, ax1.transAxes)
    ax1.plot( data_a[:,0], data_a[:,4], marker='o', markersize=4.0, color='0.8' )


    #surft_fit = np.exp( log_func_Psi_C(expon)( time_rf_a, *popt2 ) )
    #h1, = ax3.plot( tprime0 + time_rf_a, surft_fit, ':', color='blue', linewidth=2, label=r'Fit2' )
    #surft_fit2 = np.exp( log_func_Psi_C(expon)( time_rf_abs_a, *popt ) )
    #h2, = ax3.plot( time_rf_abs_a, surft_fit2, '-', color='black', linewidth=2, label=r'Fit')
    #ax3.axvline( tprime0, ymin=0.1, ymax=0.8, color='0.25', linestyle=':')
    #ax3.text( tprime0, 0.82, '$t^\prime_0$', ha='right', va='bottom', transform = trans )
    #ax3.axvline( tprime1abs, ymin=0.1, ymax=0.8, color='0.25', linestyle=':')
    #ax3.text( tprime1abs, 0.82, '$t^\prime_1$', ha='right', va='bottom', transform = trans )
    title = r'(b) $T0_m(t)$, K'
    yticks = [1500,2000,2500,3000,3500,4000]
    fig_o.set_myaxes( ax1, title=title, ylabel='$T0_m$', xlabel='$t$ (yrs)', yticks=yticks )
    ax1.set_ylim( [1000,4500])
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax1.yaxis.set_label_coords(-0.2,0.575)
    # ---- end surface temperature ----
    # ---------------------------------

    # ---------------------------
    # ---- rheological front ----
    #ind_pcut = np.where( (data_a[:,1]<PMAX) & (data_a[:,1]>PMIN) )[0]
    # time offset has already been applied, so time is 0 at t'=0
    #time_a = data_a[:,1][ind_pcut] # time with offset
    #pres_a = data_a[:,2][ind_pcut] # pressure
    #temp_a = data_a[:,3][ind_pcut] # temperature
    # save data to text files for plotting
    #if EXPORT:
    #    np.savetxt( 'Pr.dat', np.column_stack((time_a,pres_a)))
    #    np.savetxt( 'Tr.dat', np.column_stack((time_a,temp_a)))


    # indices for plotting
    #indp = np.where( data_a[:,1] > 1.0 )[0]
    #tmax = np.max( data_a[:,0][indp] )
    #tmin = np.min( data_a[:,0][indp] )
    #vmax = np.max( data_a[:,0][indp] )
    #vmin = np.min( data_a[:,0][indp] )

    ax2.plot( data_a[:,1], data_a[:,2], marker='o', markersize=4.0, color='0.8' )
    #Pf_fit = pressure_from_time_rheological_front( data_a[:,0][indp], expon, Psi, C, lambdaa, T0c )
    #h2, = ax0.plot( data_a[:,0][indp], Pf_fit, '-', color='black', linewidth=2, label=r'Fit' )
    title = r'(c) $P_f(t^\prime)$, GPa'
    yticks = [0,50,100,150]
    fig_o.set_myaxes( ax2, title=title, ylabel='$P_f$', xlabel='$t^\prime$ (yrs)', yticks=yticks )
    ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax2.yaxis.set_label_coords(-0.2,0.5)

    ax3.plot( data_a[:,1], data_a[:,3], marker='o', markersize=4.0, color='0.8' )
    #Tf_fit = temperature_from_time_rheological_front( data_a[:,0][indp], expon, Psi, C, lambdaa, T0c, alpha )
    #h2, = ax1.plot( data_a[:,0][indp], Tf_fit, '-', color='black', linewidth=2, label=r'Fit' )
    title = r'(d) $T_f(t^\prime)$, K'
    yticks = [2000,3000,4000]
    fig_o.set_myaxes( ax3, title=title, ylabel='$T_f$', yticks=yticks, xlabel='$t^\prime$ (yrs)' )
    ax3.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax3.set_ylim( [1500,4500] )
    ax3.yaxis.set_label_coords(-0.2,0.575)

    fig_o.savefig(9)

#====================================================================
def Rsquared( ydata, fmodel ):

    # sum of squares of residuals
    ybar = np.mean( ydata )
    print( 'ybar=', ybar )
    SSres = np.sum( np.square(ydata-fmodel) )
    print( 'SSres=', SSres )
    SStot = np.sum( np.square(ydata-ybar) )
    print( 'SStot=', SStot )
    Rsquared = 1.0 - SSres/SStot
    print( 'Rsquared=', Rsquared )

    return Rsquared

#====================================================================
def main():

    figure9()
    plt.show()

#====================================================================

if __name__ == "__main__":

    main()
