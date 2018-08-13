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
def log_func_Psi_C( expon ):
    '''Eq 8 in manuscript.  Used to determine Psi and C.'''
    def log_func( x, *p):
        y = p[0]*x+p[1]
        y = expon * np.log(y)
        return y
    return log_func

#====================================================================
def log_func_lambda_T0c( expon, Psi, C ):
    '''Eq 11 in manuscript.  Used to determine lambda and T0c.'''
    def log_func( x, *p):
        y = ( Psi*x + C )**expon # Eq 8 in manuscript for Psi and C
        y = np.log( p[0] ) + np.log( y - p[1] )
        return y
    return log_func

#====================================================================
def log_func_alpha( expon, Psi, C, lambdaa, T0c ):
    def log_func( x, *p ):
        Pf = pressure_from_time_rheological_front( x, expon, Psi, C, lambdaa, T0c )
        y = np.log( T0c + p[0]*Pf )
        return y
    return log_func

#====================================================================
def pressure_from_time_rheological_front( time, expon, Psi, C, lambdaa, T0c ):
    '''pressure of rheological front from time'''
    y = (Psi*time+C)**expon - T0c
    y *= lambdaa
    return y

#====================================================================
def time_from_pressure_rheological_front( pressure, expon, Psi, C, lambdaa, T0c ):
    '''Inverse calculation of pressure_from_time_rheological_front'''
    t = ( pressure / lambdaa + T0c ) ** (1.0/expon) - C
    t /= Psi
    return t

#====================================================================
def temperature_from_time_rheological_front( time, expon, Psi, C, lambdaa, T0c, alpha ):
    '''temperature of rheological front from time'''
    Pf = pressure_from_time_rheological_front( time, expon, Psi, C, lambdaa, T0c )
    Tf = T0c + alpha * Pf
    return Tf

#====================================================================
def fit_pressure_rheological_front( ind_a, data_a, expon, Psi, C ):

    time_a = data_a[:,0][ind_a] # time
    pres_a = data_a[:,1][ind_a] # pressure

    in_func = log_func_lambda_T0c( expon, Psi, C )
    popt, pcov = curve_fit( in_func, time_a, np.log(pres_a), maxfev=80000, p0=(1.0,1.0) )

    # lambda is a python keyword, hence lambdaa
    lambdaa = popt[0]
    T0c = popt[1]

    tprime1 = time_from_pressure_rheological_front( 0.0, expon, Psi, C, lambdaa, T0c )
    
    print( 'R-squared information for Pr' )
    pres_fit = pressure_from_time_rheological_front( time_a, expon, Psi, C, lambdaa, T0c )
    Rsqr = Rsquared( pres_a, pres_fit )
    Rsqr = np.round( Rsqr, 4 )
    #print ('R-squared= {}'.format( Rsqr ) )

    return popt, pcov, tprime1

#====================================================================
def fit_temperature_rheological_front( ind_a, data_a, expon, Psi, C, lambdaa, T0c ):

    time_a = data_a[:,0][ind_a] # time
    temp_a = data_a[:,2][ind_a] # temperature

    in_func = log_func_alpha( expon, Psi, C, lambdaa, T0c )
    popt, pcov = curve_fit( in_func, time_a, np.log(temp_a), maxfev=80000, p0=(1.0) )

    alpha = popt[0]

    print( 'R-squared information for Tr' )
    temp_fit = temperature_from_time_rheological_front( time_a, expon, Psi, C, lambdaa, T0c, alpha )
    Rsqr = Rsquared( temp_a, temp_fit )
    Rsqr = np.round( Rsqr, 4 )
    #print ('R-squared= {}'.format( Rsqr ) )

    return popt, pcov

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
    fig_o = su.FigureData( 1, 3, width, height, 'rheological_front' )

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

        # for dynamic criteria
        yy1 = myjson_o.get_scaled_field_values( 'regime_b' )
        yy2 = 1.5
        # for melt fraction criteria
        #yy1 = myjson_o.get_scaled_field_values( 'phi_b' )
        #yy2 = 0.4

        # always need temperature
        y_temp = myjson_o.get_scaled_field_values( 'temp_b' )

        # find the pressure(s) that correspond(s) to the cross-over point
        pres = su.find_xx_for_yy( xx_pres_b, yy1, yy2 )
        ind = su.get_first_non_zero_index( pres )

        if ind is not None:
            data_l.append( (xx_pres_b[ind], y_temp[ind], time) )

        # break loop once the rheological front has reached
        # the surface
        if ind == 1:
            break

    # reshape data
    # column order: pressure, temperature, time
    data_a = np.reshape( data_l, (-1,3) )

    # set time offset to zero, such that t0 defines the time at which the
    # rheological front is at the CMB
    data_a[:,2] -= data_a[0,2]

    # just fit the part where the rheological front is advancing
    # through the magma ocean
    #polydegree = 3 # hard-coded
    #pres_coeff = np.polyfit( data_a[:,2], data_a[:,0], polydegree ) # pressure
    #pres_poly = np.poly1d( pres_coeff )
    #pres_fit = pres_poly( data_a[:,2] )
    ind = np.where( data_a[:,0] > 5.0 )[0]
    data2fit = data_a[:,2][ind]
    data0fit = data_a[:,0][ind]
    data1fit = data_a[:,1][ind]

    popt, pcov = curve_fit( exp_func, data2fit, data0fit, maxfev=80000 )

    #temp_coeff = np.polyfit( data_a[:,2], data_a[:,1], polydegree ) # temperature
    #temp_poly = np.poly1d( temp_coeff )
    #temp_fit = temp_poly( data_a[:,2] )
    popt2, pcov2 = curve_fit( exp_func, data2fit, data1fit, maxfev=80000 )

    tmax = np.max( data_a[:,2] )
    tmin = np.min( data_a[:,2] )

    vmax = np.max( data_a[:,2] )
    vmin = np.min( data_a[:,2] )

    print( '----------------------------------------------' )
    print( 'rheological front advancing through the mantle' )
    print( '----------------------------------------------' )
    print( 'fitting coefficients' % vars() )
    print( 'between t_min= %(tmin)s yrs and t_max= %(tmax)s yrs' % vars() )
    print( 'pres_coeff=', popt )
    print( 'temp_coeff=', popt2 )

    fig_o.ax[0].scatter( data_a[:,2], data_a[:,0], c=data_a[:,2], cmap='inferno', vmin=vmin, vmax=vmax, s=3.0  )
    #h1, = fig_o.ax[0].plot( data_a[:,2], pres_fit, linestyle='-', color='green', label=r'Fit' )
    h1, = fig_o.ax[0].plot( data_a[:,2], exp_func(data_a[:,2], *popt), linestyle='-', color='green', label=r'Fit' )

    fig_o.ax[1].scatter( data_a[:,2], data_a[:,1], c=data_a[:,2], cmap='inferno', vmin=vmin, vmax=vmax, s=3.0 )
    #h2, = fig_o.ax[1].plot( data_a[:,2], temp_fit, linestyle='-', color='green', label=r'Fit' )
    h2, = fig_o.ax[1].plot( data_a[:,2], exp_func(data_a[:,2], *popt2), linestyle='-', color='green', label=r'Fit' )

    handle_scatter = fig_o.ax[2].scatter( data_a[:,0], data_a[:,1], c=data_a[:,2], cmap='inferno', vmin=vmin, vmax=vmax, s=3.0 )
    handle_fit, = fig_o.ax[2].plot( exp_func(data_a[:,2], *popt), exp_func(data_a[:,2], *popt2), color='green', linestyle='-', label=r'Fit' )

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
def figure9():

    # set parameters
    # flag to enable data output
    EXPORT = True
    # exponent of surface heat flux
    nn = 4.0
    #nn = 0.0
    # plot low pressure (linear) extension to cooling profile
    PLOT_LOWP = False
    # pressure cut-offs for curve fitting and plotting
    PMIN = 1.0 # GPa
    PMAX = 135.0 # GPa

    # this exponent appears in lots of the equations
    expon = 1.0/(1.0-nn)

    width = 4.7747
    height = 4.7747
    fig_o = su.FigureData( 2, 2, width, height, 'rheological_front' )
    fig_o.fig.subplots_adjust(wspace=0.4,hspace=0.4)

    ax0 = fig_o.ax[1][0]
    ax1 = fig_o.ax[1][1]
    ax2 = fig_o.ax[0][0]
    ax3 = fig_o.ax[0][1]

    fig_o.time = su.get_all_output_times()

    # pressure in GPa
    time = fig_o.time[0]
    myjson_o = su.MyJSON( 'output/{}.json'.format(time) )
    xx_pres_b = myjson_o.get_scaled_field_values('pressure_b')
    xx_pres_b *= 1.0E-9

    data_l = []
    tempsurf_l = []
    fluxsurf_l = []
    handle_l = []

    for ii, time in enumerate( fig_o.time ):

        print('time=',time)
        myjson_o = su.MyJSON( 'output/{}.json'.format(time) )

        # for dynamic criteria
        yy1 = myjson_o.get_scaled_field_values( 'regime_b' )
        yy2 = 1.5
        # for melt fraction criteria
        #yy1 = myjson_o.get_scaled_field_values( 'phi_b' )
        #yy2 = 0.4

        # mantle potential temperature
        y_temp = myjson_o.get_scaled_field_values( 'temp_b' )

        # mantle surface temperature
        tempsurf = myjson_o.get_scaled_field_values( 'temperature_surface' )
        tempsurf_l.append( tempsurf )

        # surface flux
        y_flux = myjson_o.get_scaled_field_values( 'Jtot_b' )
        fluxsurf_l.append( y_flux[0] )

        # find the pressure(s) that correspond(s) to the cross-over point
        pres = su.find_xx_for_yy( xx_pres_b, yy1, yy2 )
        ind = su.get_first_non_zero_index( pres )

        if ind is not None:
            data_l.append( (time, xx_pres_b[ind], y_temp[ind]) )

    data_a = np.reshape( data_l, (-1,3) )

    # set time offset to zero; tprime0 defines the time at which the
    # rheological front is at the CMB
    tprime0 = data_a[0,0]
    data_a[:,0] -= tprime0

    # -----------------------------
    # ---- surface temperature ----

    # surface temperature evolution from start to cut off temperature
    # of around 1400 K for plotting first set of fitting parameters
    # Psi and C (see manuscript)
    surft_a = np.column_stack( (fig_o.time, tempsurf_l) )
    ind_tcut = surft_a[:,1] > 1200.0
    surft0_a = surft_a[:,0][ind_tcut] # time
    surft1_a = surft_a[:,1][ind_tcut] # surface temperature
    # reasonable initial guess is quite important
    popt, pcov = curve_fit( log_func_Psi_C(expon), surft0_a, np.log(surft1_a), maxfev=80000, p0=[1.0E-13,1.0E-11]  )

    # based on the fit to the whole range of data, determine the value
    # of the surface temperature at tprime0
    logsurftC = log_func_Psi_C(expon)( tprime0, *popt )
    # temperature (K) at t'=0
    surftC = np.exp( logsurftC )
    print( 'surftC= {}'.format(surftC) )
    # the parameter C in Eq. 8 is raised to a power
    param_C = surftC**(1.0/expon)
    #print( 'param_C= {}'.format(param_C) )
    # the 1st fitting parameter (cooling rate) is the same
    # the 2nd fitting parameter (constant) needs to be updated
    #     to account for the time shift to t'=0
    popt2 = (popt[0], param_C)
    # fitting parameters for Psi and C
    print( 'Psi= {}, C= {}'.format( *popt2 ) )

    # ---- end surface temperature ----
    # ---------------------------------


    # ---------------------------
    # ---- surface heat flux ----

    surff_a = np.column_stack( (fig_o.time, fluxsurf_l) )
    #ind = tflux_a[:,0] <= tprime1abs
    surff0_a = surff_a[:,0][ind_tcut]
    surff1_a = surff_a[:,1][ind_tcut]

    # --- end of surface heat flux ----
    # --------------------------------- 


    # ---------------------------
    # ---- rheological front ----
    ind_pcut = np.where( (data_a[:,1]<PMAX) & (data_a[:,1]>PMIN) )[0]
    # time offset has already been applied, so time is 0 at t'=0
    time_a = data_a[:,0][ind_pcut] # time
    pres_a = data_a[:,1][ind_pcut] # pressure
    temp_a = data_a[:,2][ind_pcut] # temperature
    # save data to text files for plotting
    if EXPORT:
        np.savetxt( 'Pr.dat', np.column_stack((time_a,pres_a)))
        np.savetxt( 'Tr.dat', np.column_stack((time_a,temp_a)))

    # due to the change in the gradient of the melting curve around 20 GPa,
    # we should formally fit the cooling trajectory using two different
    # power laws.
    PCUT = 20.0 # GPa

    # high pressure data (P>PCUT)
    ind_a = np.where( data_a[:,1]>PCUT )[0]
    popt3, pcov3, tprime1 = fit_pressure_rheological_front( ind_a, data_a, expon, *popt2 )
    tprime1abs = tprime1 + tprime0
    print( 'start of rheological front advancing', tprime0 )
    print( 'end of rheological front advancing', tprime1abs )

    # option to plot another linear segment, since a change in the
    # gradient of the melting curve will break the constant
    # gradient assumption (lambdaa)
    #if PLOT_LOWP:
        # low pressure data (P<=PCUT)
        #indl = np.where( data_a[:,1]<PCUT )[0]
        #pop4, pcov4, tprime3 = fit_pressure_rheological_front( ind_a, data_a, expon, *popt2 )
        #tprime3abs = tprime3 + tprime0

    # clarify parameters
    # expon = expon
    # nn = nn
    Psi = popt2[0]
    C = popt2[1]
    lambdaa = popt3[0]
    T0c = popt3[1]

    # fit temperature of the rheological front
    popt4, pcov4 = fit_temperature_rheological_front( ind_a, data_a, expon, Psi, C, lambdaa, T0c )
    alpha = popt4[0]

    # only need to plot until time that rheological front reaches the
    # surface (up to tprime1)
    rel_ind = time_a <= tprime1
    abs_ind = surft0_a <= tprime1abs
    time_rf_a = time_a[rel_ind]
    time_rf_abs_a = surft0_a[abs_ind]

    trans = transforms.blended_transform_factory(
        ax3.transData, ax3.transAxes)
    ax3.plot( surft0_a, surft1_a, marker='o', markersize=4.0, color='0.8' )
    #surft_fit = np.exp( log_func_Psi_C(expon)( time_rf_a, *popt2 ) )
    #h1, = ax3.plot( tprime0 + time_rf_a, surft_fit, ':', color='blue', linewidth=2, label=r'Fit2' )
    surft_fit2 = np.exp( log_func_Psi_C(expon)( time_rf_abs_a, *popt ) )
    h2, = ax3.plot( time_rf_abs_a, surft_fit2, '-', color='black', linewidth=2, label=r'Fit')
    ax3.axvline( tprime0, ymin=0.1, ymax=0.8, color='0.25', linestyle=':')
    ax3.text( tprime0, 0.82, '$t^\prime_0$', ha='right', va='bottom', transform = trans )
    ax3.axvline( tprime1abs, ymin=0.1, ymax=0.8, color='0.25', linestyle=':')
    ax3.text( tprime1abs, 0.82, '$t^\prime_1$', ha='right', va='bottom', transform = trans )
    title = r'(b) $T0_m(t)$, K'
    yticks = [1500,2000,3000,4000]
    fig_o.set_myaxes( ax3, title=title, ylabel='$T0_m$', xlabel='$t$ (yrs)', yticks=yticks )
    ax3.set_ylim( [1000,4500])
    ax3.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax3.yaxis.set_label_coords(-0.2,0.575)

    # surface heat flux evolution
    trans = transforms.blended_transform_factory(
        ax2.transData, ax2.transAxes)
    h1, = ax2.semilogy( surff0_a[abs_ind], surff1_a[abs_ind], marker='o', markersize=4.0, color='0.8' )
    ax2.axvline( tprime0, ymin=0.1, ymax=0.8, color='0.25', linestyle=':')
    ax2.text( tprime0, 0.82, '$t^\prime_0$', ha='right', va='bottom', transform = trans )
    ax2.axvline( tprime1abs, ymin=0.1, ymax=0.8, color='0.25', linestyle=':')
    ax2.text( tprime1abs, 0.82, '$t^\prime_1$', ha='right', va='bottom', transform = trans )
    title = r'(a) $q_0(t)$, W/m$^2$'
    fig_o.set_myaxes( ax2, title=title, ylabel='$q_0$', xlabel='$t$ (yrs)' )
    ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax2.yaxis.set_label_coords(-0.2,0.5)

    # indices for plotting
    indp = np.where( data_a[:,1] > 1.0 )[0]
    tmax = np.max( data_a[:,0][indp] )
    tmin = np.min( data_a[:,0][indp] )
    vmax = np.max( data_a[:,0][indp] )
    vmin = np.min( data_a[:,0][indp] )

    ax0.plot( data_a[:,0][indp], data_a[:,1][indp], marker='o', markersize=4.0, color='0.8' )
    Pf_fit = pressure_from_time_rheological_front( data_a[:,0][indp], expon, Psi, C, lambdaa, T0c )
    h2, = ax0.plot( data_a[:,0][indp], Pf_fit, '-', color='black', linewidth=2, label=r'Fit' )
    # low pressure fit
    #if PLOT_LOWP:
    #    h3, = ax0.plot( data_a[:,0][indp], exp_func( data_a[:,0][indp], *popt_Prl ), '-', color='black', linewidth=2, label=r'Fit' )
    title = r'(c) $P_f(t^\prime)$, GPa'
    yticks = [0,50,100,150]
    fig_o.set_myaxes( ax0, title=title, ylabel='$P_f$', xlabel='$t^\prime$ (yrs)', yticks=yticks )
    ax0.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax0.yaxis.set_label_coords(-0.2,0.5)

    ax1.plot( data_a[:,0][indp], data_a[:,2][indp], marker='o', markersize=4.0, color='0.8' )
    Tf_fit = temperature_from_time_rheological_front( data_a[:,0][indp], expon, Psi, C, lambdaa, T0c, alpha )
    h2, = ax1.plot( data_a[:,0][indp], Tf_fit, '-', color='black', linewidth=2, label=r'Fit' )
    # low pressure fit
    #if PLOT_LOWP:
    #    h3, = ax1.plot( data_a[:,0][indp], exp_func( data_a[:,0][indp], *popt_Trl ), '-', color='black', linewidth=2, label=r'Fit' )
    title = r'(d) $T_f(t^\prime)$, GPa'
    yticks = [2000,3000,4000]
    fig_o.set_myaxes( ax1, title=title, ylabel='$T_f$', yticks=yticks, xlabel='$t^\prime$ (yrs)' )
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax1.set_ylim( [1500,4500] )
    ax1.yaxis.set_label_coords(-0.2,0.575)

    #tmax2 = int(np.round(tmax,0))
    #intPCUT = int(np.round(PCUT,0))
    #print('line below is for copy-paste into table')
    #print( '{:d} & {:d} & {:e} & {:e} & {:e} & {:0.4f} & {:e} & {:e} & {:e} & {:0.4f}'.format(tmax2,intPCUT,popt_Prh[0], popt_Prh[1], popt_Prh[2],RsqrPrh,popt_Trh[0],popt_Trh[1],popt_Trh[2],RsqrTrh))

    print()
    print( 'Fitting coefficients' )
    print( '--------------------' )
    print( 'n= {}'.format(nn) )
    print( 'Psi= {}'.format(Psi) )
    print( 'C= {}'.format(C) )
    print( 'lambda= {}'.format(lambdaa) )
    print( 'T0c= {}'.format(T0c) )
    print( 'alpha= {}'.format(alpha) )

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
def figure10():
 
    width = 4.7747
    height = 4.7747*3.0/2
    fig_o = su.FigureData( 3, 2, width, height, 'coefficient_surface' )
    fig_o.fig.subplots_adjust(wspace=0.4,hspace=0.4)

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[1][0]
    ax2 = fig_o.ax[2][0]
    ax3 = fig_o.ax[0][1]
    ax4 = fig_o.ax[1][1]
    ax5 = fig_o.ax[2][1]

    a_data = get_all_data()

    rad = np.unique(a_data[:,0])
    flux = np.unique(a_data[:,1])
    Tau = a_data[:,2]
    Pra = a_data[:,3]
    Prb = a_data[:,4]
    Prc = a_data[:,5]
    Tra = a_data[:,6]
    Trb = a_data[:,7]
    Trc = a_data[:,8]
 
    X, Y = np.meshgrid(rad,flux)

    num = np.size(rad)
    TAU = Tau.reshape(-1,num)
    PRA = Pra.reshape(-1,num)
    PRB = Prb.reshape(-1,num)
    PRC = Prc.reshape(-1,num)
    TRA = Tra.reshape(-1,num)
    TRB = Trb.reshape(-1,num)
    TRC = Trc.reshape(-1,num)

    dx = 0.9
    dy = 4
    aspect= dx/dy
    xticks = (0.1,0.5,1.0)
    yticks = (2,3,4,5,6)
    extent = (0.1,1.0,2,6)

    #plt.pcolormesh( X, Y, Z, cmap='RdBu', shading='gourand')
    im0 = ax0.imshow(PRA,origin='lower',interpolation='bilinear',extent=extent,aspect=aspect)
    title = r'$P_r,\ a$ coefficient'
    fig_o.set_myaxes( ax0, title=title, ylabel='$Q$', xlabel='$R$', yticks=yticks, xticks=xticks )
    plt.colorbar(im0,ax=ax0)

    #plt.pcolormesh( X, Y, Z, cmap='RdBu', shading='gourand')
    im1 = ax1.imshow(PRB,origin='lower',interpolation='bilinear',extent=extent,aspect=aspect)
    title = r'$P_r,\ b$ coefficient'
    fig_o.set_myaxes( ax1, title=title, ylabel='$Q$', xlabel='$R$', yticks=yticks, xticks=xticks )
    plt.colorbar(im1,ax=ax1)

    #plt.pcolormesh( X, Y, Z, cmap='RdBu', shading='gourand')
    im2 = ax2.imshow(PRC,origin='lower',interpolation='bilinear',extent=extent,aspect=aspect)
    title = r'$P_r,\ c$ coefficient'
    fig_o.set_myaxes( ax2, title=title, ylabel='$Q$', xlabel='$R$', yticks=yticks, xticks=xticks )
    plt.colorbar(im2,ax=ax2)

    #plt.pcolormesh( X, Y, Z, cmap='RdBu', shading='gourand')
    im3 = ax3.imshow(TRA,origin='lower',interpolation='bilinear',extent=extent,aspect=aspect)
    title = r'$T_r,\ a$ coefficient'
    fig_o.set_myaxes( ax3, title=title, ylabel='$Q$', xlabel='$R$', yticks=yticks, xticks=xticks )
    plt.colorbar(im3,ax=ax3)

    #plt.pcolormesh( X, Y, Z, cmap='RdBu', shading='gourand')
    im4 = ax4.imshow(TRB,origin='lower',interpolation='bilinear',extent=extent,aspect=aspect)
    title = r'$T_r,\ b$ coefficient'
    fig_o.set_myaxes( ax4, title=title, ylabel='$Q$', xlabel='$R$', yticks=yticks, xticks=xticks )
    plt.colorbar(im4,ax=ax4)

    #plt.pcolormesh( X, Y, Z, cmap='RdBu', shading='gourand')
    im5 = ax5.imshow(TRC,origin='lower',interpolation='bilinear',extent=extent,aspect=aspect)
    title = r'$T_r,\ c$ coefficient'
    fig_o.set_myaxes( ax5, title=title, ylabel='$Q$', xlabel='$R$', yticks=yticks, xticks=xticks )
    plt.colorbar(im5,ax=ax5)

    fig_o.savefig(10)

#====================================================================
def figure11():
 
    width = 4.7747/2
    height = 4.7747/2
    fig_o = su.FigureData( 1, 1, width, height, 'contour_plot' )
    #fig_o.fig.subplots_adjust(wspace=0.4,hspace=0.4)

    ax0 = fig_o.ax

    dir_l = ['R10Q6','R10Q5','R10Q4','R10Q3','R10Q2']

    aa = get_all_data()
    bb = (aa[14],aa[11],aa[8],aa[5],aa[2])

    for nn, directory in enumerate(dir_l):
        infile = os.path.join( directory, 'Pr.dat' )
        time, pres = np.loadtxt( infile, unpack=True )
        ax0.semilogx( time, pres )
        popt = (bb[nn][3],bb[nn][4],bb[nn][5])
        ax0.semilogx( time, exp_func(time,*popt), 'k--' )

    title = 'Cooling fits'
    fig_o.set_myaxes( ax0, title=title, ylabel='$P_f$', xlabel='$t^\prime$ (yrs)' )
    ax0.set_ylim([0.0,135])
    fig_o.savefig(11)


#====================================================================
def get_all_data():

    a_data = np.array([
        (0.1,2,22000,-7.154085e-04,1.181943e+01,8.707197e-01,-3.144027e-02,2.059165e+03,8.950754e-01),
        (0.5,2,450000,-2.491468e-01,7.020726e+01,4.258310e-01,-1.502319e+00,3.297166e+03,5.193690e-01),
        (1.0,2,1970000,-7.039320e-03,1.402555e+02,6.996766e-01,-3.920296e-02,4.381983e+03,7.774888e-01),
        (0.1,3,5000,-1.070143e-03,1.151725e+01,1.072060e+00,-3.564894e-02,2.039594e+03,1.138273e+00),
        (0.5,3,54000,-3.917908e-01,6.925599e+01,4.786142e-01,-1.847744e+00,3.259747e+03,6.167259e-01),
        (1.0,3,168000,-3.057649e-02,1.397249e+02,7.120953e-01,-1.988970e-01,4.373327e+03,7.920353e-01),
        (0.1,4,460,-8.597201e-03,1.146412e+01,1.145746e+00,-3.281326e-01,2.036894e+03,1.215488e+00),
        (0.5,4,4600,-8.011969e-01,6.833226e+01,5.287567e-01,-4.765589e+00,3.239179e+03,6.774352e-01),
        (1.0,4,15100,-1.288923e-01,1.388771e+02,7.343078e-01,-9.679921e-01,4.359172e+03,8.188002e-01),
        (0.1,5,44,-1.006908e-01,1.142039e+01,1.208162e+00,-4.450038e+00,2.034457e+03,1.282759e+00),
        (0.5,5,410,-2.062429e+00,6.748609e+01,5.799420e-01,-1.562798e+01,3.219738e+03,7.471367e-01),
        (1.0,5,1410,-5.394592e-01,1.375848e+02,7.735149e-01,-4.710575e+00,4.338459e+03,8.648539e-01),
        (0.1,6,3,-1.964233e+00,1.044135e+01,1.147085e+00,-1.102602e+02,1.984107e+03,1.187183e+00),
        (0.5,6,37,-6.564250e+00,6.659050e+01,6.414332e-01,-6.940401e+01,3.201900e+03,8.272998e-01),
        (1.0,6,125,-2.570421e+00,1.361605e+02,8.267566e-01,-2.679335e+01,4.316338e+03,9.266183e-01)])

    print(a_data)
    return a_data

#====================================================================
def main():

    # arguments (run with -h to summarize)
    #parser = argparse.ArgumentParser(description='SPIDER plotting script')
    #parser.add_argument('-t', '--times', type=str, help='Comma-separated (no spaces) list of times');
    #args = parser.parse_args()

    #figure7( args.times )
    #figure8()
    figure9()
    #figure10()
    #figure11()
    plt.show()

#====================================================================

if __name__ == "__main__":

    main()
