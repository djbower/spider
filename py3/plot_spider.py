#!/usr/bin/env python

# next is bad practice, but fastest for getting something up and running
from spider_classes import *
import json
import logging
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.transforms as transforms
import numpy as np
import os
import sys
import argparse

logger = logging.getLogger(__name__)

#====================================================================
def get_mix( melt_fraction_a ):

    # this is a boolean array
    MIX = ( melt_fraction_a < 0.999 ) & ( melt_fraction_a > 0.001 )
    # convert to float array
    MIX = MIX * 1.0
    # set single phase regions to NaNs, this prevents them from
    # plotting
    MIX[MIX == 0] = np.nan

    return MIX

#====================================================================
def bower_et_al_2018_fig3( args ):

    # article class text width is 4.7747 inches
    # http://tex.stackexchange.com/questions/39383/determine-text-width

    logger.info( 'building bower_et_al_2018_fig3' )

    fig_o = FigureData( args, 2, 2, 4.7747, 4.7747, 'bower_et_al_2018_fig3' )

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]
    ax2 = fig_o.ax[1][0]
    ax3 = fig_o.ax[1][1]

    time = fig_o.time[0] # first timestep since liquidus and solidus
                         # are time-independent

    myjson_o = MyJSON( 'output/{}.json'.format(time) )

    TIMEYRS = myjson_o.data_d['nstep']

    # hack to compute some average properties for Bower et al. (2018)
    #xx_liq, yy_liq = fig_o.get_xy_data( 'liquidus_rho', time )
    #xx_sol, yy_sol = fig_o.get_xy_data( 'solidus_rho', time )
    #diff = (yy_liq - yy_sol) / yy_sol * 100.0
    #print diff[40:]
    #print np.mean(diff[40:])
    #sys.exit(1)

    xx_pres = myjson_o.get_scaled_field_values_internal('pressure_b')
    xx_pres *= 1.0E-9

    xx_pres_s = myjson_o.get_scaled_field_values('pressure_s')
    xx_pres_s *= 1.0E-9

    # shade grey between liquidus and solidus
    yy_liq = myjson_o.get_scaled_field_values_internal('liquidus_b')
    yy_sol = myjson_o.get_scaled_field_values_internal('solidus_b')
    ax0.fill_between( xx_pres, yy_liq, yy_sol, facecolor='grey', alpha=0.35, linewidth=0 )
    yy_liqt = myjson_o.get_scaled_field_values_internal('liquidus_temp_b')
    yy_solt = myjson_o.get_scaled_field_values_internal('solidus_temp_b')
    ax1.fill_between( xx_pres, yy_liqt, yy_solt, facecolor='grey', alpha=0.35, linewidth=0 )

    # hack to compute some average properties for Bower et al. (2018)
    #print xx_sol
    #print np.mean(yy_liq[20:]-yy_sol[20:])
    #sys.exit(1)

    # dotted lines of constant melt fraction
    for xx in range( 0, 11, 2 ):
        yy_b = xx/10.0 * (yy_liq - yy_sol) + yy_sol
        if xx == 0:
            # solidus
            ax0.plot( xx_pres, yy_b, '-', linewidth=0.5, color='black' )
        elif xx == 3:
            # typically, the approximate location of the rheological transition
            ax0.plot( xx_pres, yy_b, '-', linewidth=1.0, color='white')
        elif xx == 10:
            # liquidus
            ax0.plot( xx_pres, yy_b, '-', linewidth=0.5, color='black' )
        else:
            # dashed constant melt fraction lines
            ax0.plot( xx_pres, yy_b, '--', linewidth=1.0, color='white' )

    handle_l = [] # handles for legend

    for nn, time in enumerate( fig_o.time ):
        # read json
        myjson_o = MyJSON( 'output/{}.json'.format(time) )

        color = fig_o.get_color( nn )
        # use melt fraction to determine mixed region
        yy = myjson_o.get_scaled_field_values_internal('phi_b')
        MIX = get_mix( yy )
        yy_s = myjson_o.get_scaled_field_values('phi_s')
        MIX_s = get_mix( yy_s )

        label = fig_o.get_legend_label( time )

        # entropy
        # plot staggered, since this is where entropy is defined
        # cell-wise and we can easily see the CMB boundary condition
        yy = myjson_o.get_scaled_field_values('S_s')
        #yy = myjson_o.get_scaled_field_values_internal('S_b')
        ax0.plot( xx_pres_s, yy, '-', color=color )
        handle, = ax0.plot( xx_pres_s*MIX_s, yy*MIX_s, '-', color=color, label=label )
        handle_l.append( handle )
        # temperature
        yy = myjson_o.get_scaled_field_values_internal('temp_b')
        ax1.plot( xx_pres, yy, '--', color=color )
        ax1.plot( xx_pres*MIX, yy*MIX, '-', color=color )
        # melt fraction
        yy = myjson_o.get_scaled_field_values_internal('phi_b')
        ax2.plot( xx_pres, yy, '-', color=color )
        # viscosity
        visc_const = 1 # this is used for the arcsinh scaling
        visc_fmt = MyFuncFormatter( visc_const )
        yy = myjson_o.get_scaled_field_values_internal('visc_b', visc_fmt)
        ax3.plot( xx_pres, yy, '-', color=color )

    xticks = [0,50,100,135]
    xmax = 138

    # titles and axes labels, legends, etc
    units = myjson_o.get_field_units('S_b')
    title = '(a) Entropy, {}'.format(units)
    yticks = [700,1600,2400,3200]
    # Bower et al. (2018)
    #yticks = [1600,2000,2400,2800,3200]
    # DJB used this next range for work with Bayreuth
    #yticks = [300,1000,1600,2000,2400,2800,3200]
    fig_o.set_myaxes( ax0, title=title, ylabel='$S$', xticks=xticks, xmax=xmax, yticks=yticks )
    ax0.yaxis.set_label_coords(-0.075,0.59)
    units = myjson_o.get_field_units('temp_b')
    title = '(b) Temperature, {}'.format(units)
    yticks = [300,1000,2000,3000,4000,5000]
    # Bower et al. (2018)
    #yticks= [1000,2000,3000,4000,5000]
    # DJB used this next range for work with Bayreuth
    #yticks= [300,1000,2000,3000,4000,5000]
    fig_o.set_myaxes( ax1, title=title, ylabel='$T$', xticks=xticks, xmax=xmax, yticks=yticks )
    ax1.set_xlim( xticks[0], 138 )
    ax1.yaxis.set_label_coords(-0.075,0.59)
    fig_o.set_mylegend( ax1, handle_l, loc=4, ncol=2 )
    #fig_o.set_mylegend( ax0, handle_l, loc=2, ncol=2 )
    title = '(c) Melt fraction'
    yticks = [0,0.2,0.4,0.6,0.8,1.0]
    fig_o.set_myaxes( ax2, title=title, xlabel='$P$ (GPa)',
        ylabel='$\phi$', xticks=xticks, xmax=xmax, yticks=yticks )
    ax2.yaxis.set_label_coords(-0.075,0.475)
    ax2.set_ylim( [0, 1] )
    units = myjson_o.get_field_units('visc_b')
    title = '(d) Viscosity, ' + units
    yticks = [1.0E-2, 1.0E6, 1.0E12, 1.0E18, 1.0E21]
    fig_o.set_myaxes( ax3, title=title, xlabel='$P$ (GPa)',
        ylabel='$\eta$', xticks=xticks, xmax=xmax, yticks=yticks, fmt=visc_fmt )
    ax3.yaxis.set_label_coords(-0.075,0.67)

    if 0:
        # add sol and liq text boxes to mark solidus and liquidus
        # these have to manually adjusted according to the solidus
        # and liquidus used for the model
        prop_d = {'boxstyle':'round', 'facecolor':'white', 'linewidth': 0}
        ax0.text(120, 2925, r'liq', fontsize=8, bbox=prop_d )
        ax0.text(120, 2075, r'sol', fontsize=8, bbox=prop_d )
        # DJB commented out for now, below labels 0.2 melt fraction
        # so the reader knows that dashed white lines denote melt
        # fraction contours
        ax0.text(110, 2275, r'$\phi=0.2$', fontsize=6 )

    fig_o.savefig(1)

#====================================================================
def bower_et_al_2018_fig4( args ):

    logger.info( 'building bower_et_al_2018_fig4' )

    fig_o = FigureData( args, 2, 2, 4.7747, 4.7747, 'bower_et_al_2018_fig4' )

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]
    ax2 = fig_o.ax[1][0]
    ax3 = fig_o.ax[1][1]

    handle_l = []

    # for arcsinh scaling
    flux_const = 1.0E6
    flux_fmt = MyFuncFormatter( flux_const )

    for nn, time in enumerate( fig_o.time ):
        # uncomment below to plot every other line
        if (nn-1) % 2:
            continue

        # read json
        myjson_o = MyJSON( 'output/{}.json'.format(time) )

        color = fig_o.get_color( nn )
        # use melt fraction to determine mixed region
        yy = myjson_o.get_scaled_field_values_internal('phi_b')
        MIX = get_mix( yy )

        label = fig_o.get_legend_label( time )

        # pressure for x-axis
        xx_pres = myjson_o.get_scaled_field_values_internal('pressure_b')
        xx_pres *= 1.0E-9

        # Jconv_b
        yy = myjson_o.get_scaled_field_values_internal('Jconv_b', flux_fmt)
        ax0.plot( xx_pres, yy, '--', color=color )
        handle, = ax0.plot( xx_pres*MIX, yy*MIX, '-', label=label,
            color=color )
        handle_l.append( handle )
        # Jmix_b
        yy = myjson_o.get_scaled_field_values_internal('Jmix_b', flux_fmt)
        ax2.plot( xx_pres, yy, '--', color=color )
        ax2.plot( xx_pres*MIX, yy*MIX, '-', color=color )
        # Jgrav_b
        yy = myjson_o.get_scaled_field_values_internal('Jgrav_b', flux_fmt)
        ax1.plot( xx_pres, yy, '--', color=color )
        ax1.plot( xx_pres*MIX, yy*MIX, '-', color=color )
        # Jtot_b
        yy = myjson_o.get_scaled_field_values_internal('Jtot_b', flux_fmt)
        ax3.plot( xx_pres, yy, '--', color=color )
        ax3.plot( xx_pres*MIX, yy*MIX, '-', color=color )

    # titles and axes labels, legends, etc.
    xticks = [0,50,100,135]
    xmax = 138
    yticks = [-1E15,-1E12, -1E6, -1E0, 0, 1E0, 1E6, 1E12, 1E15]
    units = myjson_o.get_field_units('Jconv_b')
    title = '(a) Convective flux, {}'.format(units)
    fig_o.set_myaxes( ax0, title=title, ylabel='$F_\mathrm{conv}$',
        yticks=yticks, xticks=xticks, xmax=xmax, fmt=flux_fmt )
    ax0.yaxis.set_label_coords(-0.16,0.54)
    fig_o.set_mylegend( ax0, handle_l, ncol=2 )
    units = myjson_o.get_field_units('Jmix_b')
    title = '(c) Mixing flux, {}'.format(units)
    fig_o.set_myaxes( ax2, title=title, xlabel='$P$ (GPa)',
        ylabel='$F_\mathrm{mix}$', yticks=yticks, xticks=xticks, xmax=xmax, fmt=flux_fmt )
    ax2.yaxis.set_label_coords(-0.16,0.54)
    units = myjson_o.get_field_units('Jgrav_b')
    title = '(b) Separation flux, {}'.format(units)
    fig_o.set_myaxes( ax1, title=title,
        ylabel='$F_\mathrm{grav}$', yticks=yticks, xticks=xticks, xmax=xmax, fmt=flux_fmt )
    ax1.yaxis.set_label_coords(-0.16,0.54)
    units = myjson_o.get_field_units('Jtot_b')
    title = '(d) Total flux, {}'.format(units)
    fig_o.set_myaxes( ax3, title=title, xlabel='$P$ (GPa)',
        ylabel='$F_\mathrm{tot}$', yticks=yticks, xticks=xticks, xmax=xmax, fmt=flux_fmt )
    ax3.yaxis.set_label_coords(-0.16,0.54)

    fig_o.savefig(2)

#====================================================================
def bower_et_al_2018_fig5( args ):

    logger.info( 'building bower_et_al_2018_fig5' )

    # keep y the same by scaling (7.1621)
    fig_o = FigureData( args, 3, 2, 4.7747, 7.1621, 'bower_et_al_2018_fig5' )

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]
    ax2 = fig_o.ax[1][0]
    ax3 = fig_o.ax[1][1]
    ax4 = fig_o.ax[2][0]
    ax5 = fig_o.ax[2][1]

    handle_l = []

    eddy_const = 1.0
    eddy_fmt = MyFuncFormatter( eddy_const )
    alpha_const = 1.0E6
    alpha_fmt = MyFuncFormatter( alpha_const )
    dSdr_const = 1.0E15
    dSdr_fmt = MyFuncFormatter( dSdr_const )

    for nn, time in enumerate( fig_o.time ):
        # uncomment below to plot every other line
        if (nn-1) % 2:
            continue

        # read json
        myjson_o = MyJSON( 'output/{}.json'.format(time) )

        color = fig_o.get_color( nn )
        # use melt fraction to determine mixed region
        yy = myjson_o.get_scaled_field_values_internal('phi_b')
        MIX = get_mix( yy )

        label = fig_o.get_legend_label( time )

        # pressure for x-axis
        xx_pres = myjson_o.get_scaled_field_values_internal('pressure_b')
        xx_pres *= 1.0E-9

        # eddy diffusivity
        yy = myjson_o.get_scaled_field_values_internal('kappah_b', eddy_fmt)
        ax0.plot( xx_pres, yy, '--', color=color )
        handle, = ax0.plot( xx_pres*MIX, yy*MIX, '-', label=label,
            color=color )
        handle_l.append( handle )
        # density
        yy = myjson_o.get_scaled_field_values_internal('rho_b')
        ax1.plot( xx_pres, yy, '--', color=color )
        ax1.plot( xx_pres*MIX, yy*MIX, '-', color=color )
        # thermal expansion coefficient
        yy = myjson_o.get_scaled_field_values_internal('alpha_b', alpha_fmt)
        ax3.plot( xx_pres, yy, '--', color=color )
        ax3.plot( xx_pres*MIX, yy*MIX, '-', color=color )
        # heat capacity
        yy = myjson_o.get_scaled_field_values_internal('cp_b')
        ax4.plot( xx_pres, yy, '--', color=color )
        ax4.plot( xx_pres*MIX, yy*MIX, '-', color=color )
        # dTdrs
        yy = myjson_o.get_scaled_field_values_internal('dTdrs_b')
        ax5.plot( xx_pres, yy, '--', color=color )
        ax5.plot( xx_pres*MIX, yy*MIX, '-', color=color )
        # entropy gradient
        yy = myjson_o.get_scaled_field_values_internal('dSdr_b', dSdr_fmt )
        yy *= -1.0
        ax2.plot( xx_pres, yy, '--', color=color )
        ax2.plot( xx_pres*MIX, yy*MIX, '-', color=color )

    xticks = [0,50,100,135]
    xmax = 138

    # titles and axes labels, legends, etc.
    units = myjson_o.get_field_units('kappah_b')
    title = '(a) Eddy diffusivity, {}'.format(units)
    yticks = [1.0, 1.0E3, 1.0E6, 1.0E9]
    fig_o.set_myaxes( ax0, title=title, ylabel='$\kappa_h$',
        yticks=yticks, fmt=eddy_fmt, xticks=xticks, xmax=xmax )
    ax0.yaxis.set_label_coords(-0.1,0.475)

    units = myjson_o.get_field_units('rho_b')
    title = '(b) Density, {}'.format(units)
    yticks = [2000,3000,4000,5000,6000]
    fig_o.set_myaxes( ax1, title=title, ylabel='$\\rho$',
        yticks=yticks, xticks=xticks, xmax=xmax )
    ax1.yaxis.set_label_coords(-0.1,0.6)
    fig_o.set_mylegend( ax1, handle_l, loc=4, ncol=2 )

    units = myjson_o.get_field_units('dSdr_b')
    title = '(c) Entropy grad, {}'.format(units)
    #yticks= [-1.0E-3, -1.0E-6, -1.0E-9, -1.0E-12, -1.0E-15]
    yticks = [1E-15, 1E-12, 1E-9, 1E-6, 1E-3]
    fig_o.set_myaxes( ax2, title=title,
        yticks=yticks, xticks=xticks, xmax=xmax,
        ylabel='$-\\frac{\partial S}{\partial r}$', fmt=dSdr_fmt)
    ax2.yaxis.set_label_coords(-0.1,0.565)

    fig_o.savefig(3)


    units = myjson_o.get_field_units('alpha_b')
    title = '(d) Thermal expansion, {}'.format(units)
    yticks = [1.0E-5, 1.0E-4, 1.0E-3, 1.0E-2]
    fig_o.set_myaxes( ax3, title=title, ylabel='$\\alpha$',
        yticks=yticks, fmt=alpha_fmt, xticks=xticks, xmax=xmax )
    ax3.yaxis.set_label_coords(-0.1,0.475)

    units = myjson_o.get_field_units('cp_b')
    title = '(e) Heat capacity, {}'.format(units)
    yticks = [0,5000,10000,15000]
    fig_o.set_myaxes( ax4, title=title, xlabel='$P$ (GPa)', ylabel='$c$',
        yticks=yticks, xticks=xticks, xmax=xmax )
    ax4.yaxis.set_label_coords(-0.1,0.475)

    units = myjson_o.get_field_units('dTdrs_b')
    title = '(f) Adiabatic grad, {}'.format(units)
    yticks = [-3E-3, -2E-3, -1E-3, 0]
    fig_o.set_myaxes( ax5, title=title, xlabel='$P$ (GPa)',
        yticks=yticks, xticks=xticks, xmax=xmax,
        ylabel='$\\left(\\frac{\partial T}{\partial r}\\right)_S$')
    ax5.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
    ax5.yaxis.set_label_coords(-0.15,0.43)

    fig_o.savefig(3)

#====================================================================
def dep2pres( dep, R0 ):

    dep *= R0
    pres = 4078.95095544*10 / 1.1115348931000002e-07
    pres *= np.exp( 1.1115348931000002e-07 * dep) - 1.0
    pres *= 1.0E-9
    return pres

#====================================================================
def bower_et_al_2018_fig2( args ):

    logger.info( 'building bower_et_al_2018_fig2' )

    # define scalings (hard-coded)
    # these were used for the original Bower et al. (2018) work
    # but could in principle be different in future models
    S0 = 2993.025100070677
    R0 = 6371000.0
    T0 = 4033.6070755893948
    D0 = 4613.109568155063
    FLUX0 = D0 * (S0*T0)**(3.0/2.0)

    width = 4.7747 # * 0.5
    height = 4.7747 # * 0.5
    fig_o = FigureData( args, 2, 2, width, height, 'bower_et_al_2018_fig2' )

    dd = fig_o.data_d

    # below, only basic internal nodes

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]
    ax2 = fig_o.ax[1][0]
    ax3 = fig_o.ax[1][1]

    # directory containing simplified model data
    prefix = 'examples/bower_et_al_2018/simplified_model/fig2'

    ###############
    # plot liquidus
    Sliq_file = os.path.join(prefix,'Sliq_processed.dat')
    rad, Sliq = np.loadtxt( Sliq_file, unpack=True )
    depth = 1.0 - rad
    pres = dep2pres( depth, R0 )
    const = S0 # scaling for entropy but could change
    Sliq *= const
    # titles and axes labels, legends, etc
    title = '(d) Liquidus, J kg$^{-1}$ K$^{-1}$'
    yticks = [1600,2000,2400,2800,3200]
    xticks = [0,50,100,135]
    fig_o.set_myaxes( ax3, title=title, xlabel='$P$ (GPa)', ylabel='$S_{\\rm liq}$', xticks=xticks, yticks=yticks )
    ax3.yaxis.set_label_coords(-0.1,0.59)
    handle, = ax3.plot( pres, Sliq, 'k--' )

    ###############################
    # plot semi-analytical solution
    const = S0/R0

    dSdr_file = os.path.join(prefix,'dSdr_processed.dat')
    rad, dSdr = np.loadtxt( dSdr_file, unpack=True )
    depth = 1.0 - rad
    pres = dep2pres( depth, R0 )
    dSdr *= const * -1.0

    # dSliqdr was not plotted in the end
    #dSliqdr_file = os.path.join(prefix,'dSliqdr_processed.dat')
    #rad2, dSliqdr = np.loadtxt( dSliqdr_file, unpack=True )
    #depth2 = 1.0 - rad2
    #pres2 = dep2pres( depth2, R0 )
    #dSliqdr *= const
    # must scale by 1/2
    #dSliqdr *= 0.5

    dSdr_const = 1.0E15
    dSdr_fmt = MyFuncFormatter( dSdr_const )

    dSdr = ascale( dSdr, dSdr_fmt.const )
    handle2, = ax1.plot( pres, dSdr, 'k-' )

    title = '(b) Entropy grad, J kg$^{-1}$ K$^{-1}$ m$^{-1}$'
    #yticks= [-1.0E-3, -1.0E-6, -1.0E-9, -1.0E-12, -1.0E-15]#, 1.0E-3]
    yticks = [1.0E-15, 1.0E-12, 1.0E-9,1.0E-6,1.0E-3]
    fig_o.set_myaxes( ax1, yticks=yticks, xticks=xticks,
        ylabel='$-\\frac{\partial S}{\partial r}$', fmt=dSdr_fmt, title=title)
    ax1.yaxis.set_label_coords(-0.11,0.565)

    ############
    # plot Fconv
    Fconv_file = os.path.join(prefix,'FluxConv_processed.dat')
    rad, Fconv = np.loadtxt( Fconv_file, unpack=True )
    depth = 1.0 - rad
    pres = dep2pres( depth, R0 )

    const = FLUX0
    Fconv *= const

    # for arcsinh scaling
    flux_const = 1.0E6
    flux_fmt = MyFuncFormatter( flux_const )
    Fconv = ascale( Fconv, flux_fmt.const )
    handle3, = ax0.plot( pres, Fconv, 'k-' )

    # titles and axes labels, legends, etc.
    #yticks = [-1.0E12, -1.0E6, -1.0E0, -1.0E-3, 1.0E-3, 1.0E0, 1.0E6, 1E12]
    yticks = [-1.0E15, -1.0E12, -1.0E6, -1.0E0, 0, 1.0E0, 1.0E6, 1.0E12, 1.0E15]

    title = '(a) Convective flux, W m$^{-2}$'
    fig_o.set_myaxes( ax0, title=title, ylabel='$F_\mathrm{conv}$',
        yticks=yticks, fmt=flux_fmt, xticks=xticks )
    ax0.yaxis.set_label_coords(-0.16,0.54)

    ###########
    # plot Fmix
    Fmix_file = os.path.join(prefix,'FluxMix_processed.dat')
    rad, Fmix = np.loadtxt( Fmix_file, unpack=True )
    depth = 1.0 - rad
    pres = dep2pres( depth, R0 )

    const = FLUX0
    Fmix *= const

    Fmix = ascale( Fmix, flux_fmt.const )
    handle4, = ax2.plot( pres, Fmix, 'k-' )

    #fig_o.set_mylegend( ax3, handle_l, ncol=2 )
    title = '(c) Mixing flux, W m$^{-2}$'
    fig_o.set_myaxes( ax2, title=title, xlabel='$P$ (GPa)', ylabel='$F_\mathrm{mix}$',
        yticks=yticks, fmt=flux_fmt, xticks=xticks )
    ax2.yaxis.set_label_coords(-0.16,0.54)

    ###########
    # plot Ftot
    pres = np.array([0.0, 140.0])
    Ftot = np.array([10.0**6, 10.0**6])
    Ftot = ascale( Ftot, flux_fmt.const )
    ax0.plot( pres, Ftot, 'k--' )
    ax2.plot( pres, Ftot, 'k--' )

    fig_o.savefig(4)

#====================================================================
def bower_et_al_2018_fig1( args ):

    logger.info( 'building bower_et_al_2018_fig1' )

    prefix = 'examples/bower_et_al_2018/simplified_model/fig1'

    # define scalings (hard-coded)
    # these were used for the original Bower et al. (2018) work
    # but could in principle be different in future models
    S0 = 2993.025100070677
    R0 = 6371000.0
    T0 = 4033.6070755893948
    D0 = 4613.109568155063

    width = 4.7747 * 0.5
    height = 4.7747 * 0.5
    fig_o = FigureData( args, 1, 1, width, height, 'bower_et_al_2018_fig1' )

    ax0 = fig_o.ax

    const = S0/R0

    dSdr_const = 1.0E15
    dSdr_fmt = MyFuncFormatter( dSdr_const )
    yticks = [1.0E-15, 1.0E-12, 1.0E-9,1.0E-6,1.0E-3]
    #yticks = [-1E-3,-1E-6,-1E-9,-1E-12,-1E-15]

    dSliqdr_const = 1.0E6 #1.0E15
    dSliqdr_fmt = MyFuncFormatter( dSliqdr_const )
    xticks = [-1E-2,-1E-4,0,1E-4,1E-2]#,1E-3,1]

    #fig_o.set_mylegend( ax0, handle_l, ncol=2 )
    title = 'Flux solution space'# for negative fluxes'
    xlabel = '$\\frac{dS_{\\rm liq}}{dr}$'
    ylabel = '$-\\frac{\partial S}{\partial r}$'
    fig_o.set_myaxes( ax0, title=title, xlabel=xlabel, ylabel=ylabel,
        yticks=yticks, fmt=dSdr_fmt, xticks=xticks, xfmt=dSliqdr_fmt )
    ax0.yaxis.set_label_coords(-0.1,0.565)
    #ax0.set_xlim( -5.2, 0.4)

    # plot asymptote
    FLAG = 1
    if FLAG:
        # not sure why I need a negative number less than 5.2?
        # probably to account for the x shift?
        xx2 = np.linspace( -20.4, 0.8, 1000000 )
        xx2 *= const
        yy2 = 0.5 * xx2 * -1.0
        xx2 = ascale( xx2, dSliqdr_fmt.const )
        yy2 = ascale( yy2, dSdr_fmt.const )
        ax0.plot( xx2, yy2, 'k--' )
        ax0.fill_between( xx2, yy2, 0, facecolor='0.75' )

    # plot flux minimum
    if 1:
        #xx2 = np.linspace( -20.4, 0.8, 1000000 )
        xx2 = np.linspace( -20.4, -0.0001, 1000000 )
        xx2 *= const
        yy2 = 1.0/6.0 * xx2 * -1.0
        xx2 = ascale( xx2, dSliqdr_fmt.const )
        yy2 = ascale( yy2, dSdr_fmt.const )
        xx2 -= 0.4
        yy2 -= 0.4
        ax0.plot( xx2, yy2, color='0.25', linestyle=':' )

    # positive heat fluxes
    for cc, nn in enumerate([12,9,6]):
        filein = os.path.join(prefix,'dSdr_10p%(nn)d_processed.dat' % vars())
        xx, yy = np.loadtxt( filein, unpack=True )
        xx *= const
        xx = ascale( xx, dSliqdr_fmt.const )
        yy *= const * -1.0
        yy = ascale( yy, dSdr_fmt.const )
        # offset y for visual clarity
        xx += 0.4
        yy += 0.4
        ax0.plot( xx, yy, '-', color=fig_o.get_color(cc) ) #, 'k-' )

    # negative heat fluxes
    for cc, nn in enumerate([12,9,6]):
        filein = os.path.join(prefix,'dSdr_n10p%(nn)d_processed.dat' % vars())
        xx, yy = np.loadtxt( filein, unpack=True )
        xx *= const
        xx = ascale( xx, dSliqdr_fmt.const )
        yy *= const * -1.0
        yy = ascale( yy, dSdr_fmt.const )
        # offset y for visual clarity
        xx -= 0.4
        yy -= 0.4
        ax0.plot( xx, yy, '-', color=fig_o.get_color(7-cc) )

    # domains of F>0 and F<0
    if 1:
        prop_d = {'boxstyle':'round', 'facecolor':'white', 'linewidth': 0}
        ax0.text(0.23, 0.4, '$\widetilde{F}<0$', fontsize=8, bbox=prop_d, transform=ax0.transAxes, horizontalalignment='left' )
        ax0.text(0.77, 0.4, '$\widetilde{F}>0$', fontsize=8, transform=ax0.transAxes, horizontalalignment='right' )

    # negative heat flux labels
    if 1:
        ax0.text(0.0, 0.86, '-10$^{12}$', fontsize=8, rotation=50, transform=ax0.transAxes )
        ax0.text(0.01, 0.36, '-10$^9$', fontsize=8, rotation=50, transform=ax0.transAxes )
        ax0.text(0.18, 0.09, '-10$^6$', fontsize=8, rotation=50, transform=ax0.transAxes )

    # positive heat flux labels
    if 1:
        ax0.text(0.75, 0.85, '10$^{12}$', fontsize=8, rotation=0, transform=ax0.transAxes )
        ax0.text(0.75, 0.57, '10$^9$', fontsize=8, rotation=-50, transform=ax0.transAxes )
        ax0.text(0.75, 0.07, '10$^6$', fontsize=8, rotation=-50, transform=ax0.transAxes )

    # flux minimum and zero label
    if 1:
        ax0.text(0.13, 0.77, '$\widetilde{F}_{\\rm min}$', fontsize=8, rotation=-33, transform=ax0.transAxes )
        ax0.text(0.54, 0.16, '$\widetilde{F}=0$', fontsize=8, rotation=-90, transform=ax0.transAxes, horizontalalignment='center' )


    fig_o.savefig(5)

#====================================================================
def figure6( args ):

    logger.info( 'building atmosphere' )

    width = 4.7747
    height = 4.7747
    fig_o = FigureData( args, 2, 2, width, height, 'atmosphere' )

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]
    ax2 = fig_o.ax[1][0]
    ax3 = fig_o.ax[1][1]

    fig_o.time = get_all_output_times()
    timeMyr_a = fig_o.time * 1.0E-6 # Myrs

    mass_liquid_a = get_single_values_for_times( 'mass_liquid', fig_o.time )
    mass_solid_a = get_single_values_for_times( 'mass_solid', fig_o.time )
    mass_mantle_a = get_single_values_for_times( 'mass_mantle', fig_o.time )
    mass_mantle = mass_mantle_a[0] # time independent

    CO2_initial_a = get_single_values_for_times( 'CO2_initial', fig_o.time )
    CO2_liquid_a = get_single_values_for_times( 'CO2_liquid', fig_o.time )
    CO2_solid_a = get_single_values_for_times( 'CO2_solid', fig_o.time )
    CO2_atmos_a = get_single_values_for_times( 'CO2_atmosphere', fig_o.time )

    H2O_initial_a = get_single_values_for_times( 'H2O_initial', fig_o.time )
    H2O_liquid_a = get_single_values_for_times( 'H2O_liquid', fig_o.time )
    H2O_solid_a = get_single_values_for_times( 'H2O_solid', fig_o.time )
    H2O_atmos_a = get_single_values_for_times( 'H2O_atmosphere', fig_o.time )

    # compute total mass (kg) in each reservoir
    CO2_liquid_kg_a = ppm_to_mass_fraction( CO2_liquid_a ) * mass_liquid_a
    CO2_solid_kg_a = ppm_to_mass_fraction( CO2_solid_a ) * mass_solid_a
    CO2_total_kg_a = ppm_to_mass_fraction( CO2_initial_a ) * mass_mantle
    CO2_total_kg = CO2_total_kg_a[0] # time-independent
    # TODO: below mass is conserved by definition, but can also
    # compute directly from partial pressure
    CO2_atmos_kg_a = CO2_total_kg - CO2_liquid_kg_a - CO2_solid_kg_a

    H2O_liquid_kg_a = ppm_to_mass_fraction( H2O_liquid_a ) * mass_liquid_a
    H2O_solid_kg_a = ppm_to_mass_fraction( H2O_solid_a ) * mass_solid_a
    H2O_total_kg_a = ppm_to_mass_fraction( H2O_initial_a ) * mass_mantle
    H2O_total_kg = H2O_total_kg_a[0] # time-independent
    # TODO: below mass is conserved by definition, but can also
    # compute directly from partial pressure
    H2O_atmos_kg_a = H2O_total_kg - H2O_liquid_kg_a - H2O_solid_kg_a

    temperature_surface_a = get_single_values_for_times( 'temperature_surface', fig_o.time )
    emissivity_a = get_single_values_for_times( 'emissivity', fig_o.time )

    xticks = [1E-2,1E-1,1E0,1E1,1E2]
    xlabel = 'Time (Myr)'

    red = fig_o.get_color(3)
    blue = fig_o.get_color(-3)

    # to plot melt fraction contours on figure (a)
    # compute time at which desired melt fraction is reached
    phi_a = mass_liquid_a / mass_mantle
    phi_cont_l = [0.2, 0.8] # 20% and 80% melt fraction
    phi_time_l = [] # contains the times for each contour
    for cont in phi_cont_l:
        time_temp_l = find_value( timeMyr_a, phi_a, cont )
        index = get_first_non_zero_index( time_temp_l )
        phi_time_l.append( timeMyr_a[index] )

    ##########
    # figure a
    ##########
    title = '(a) Partial pressure (bar)'
    ylabel = '$p$'
    trans = transforms.blended_transform_factory(
        ax0.transData, ax0.transAxes)
    for cc, cont in enumerate(phi_cont_l):
        ax0.axvline( phi_time_l[cc], ymin=0.1, ymax=0.9, color='0.25', linestyle=':' )
        label = int(cont*100) # as percent
        ax0.text( phi_time_l[cc], 0.9, '{:2d}'.format(label), va='bottom', ha='center', transform=trans )
    ax0.text( 0.1, 0.9, '$\phi (\%)$', ha='center', va='bottom', transform=ax0.transAxes )
    h1, = ax0.semilogx( timeMyr_a, CO2_atmos_a, color=red, linestyle='-', label='CO$_2$')
    h2, = ax0.semilogx( timeMyr_a, H2O_atmos_a, color=blue, linestyle='-', label='H$_2$O')
    handle_l = [h1,h2]
    fig_o.set_myaxes( ax0, title=title, ylabel=ylabel, xticks=xticks )
    fig_o.set_mylegend( ax0, handle_l, loc='center left', ncol=1, TITLE="" )
    ax0.yaxis.set_label_coords(-0.1,0.575)

    ##########
    # figure b
    ##########
    title = '(b) Volatile reservoirs'
    #h5, = ax1.semilogx( timeMyr_a, mass_liquid_a / mass_mantle, 'k--', label='melt' )
    h1, = ax1.semilogx( timeMyr_a, (CO2_liquid_kg_a+CO2_solid_kg_a) / CO2_total_kg, color=red, linestyle='-', label='Magma' )
    h2, = ax1.semilogx( timeMyr_a, CO2_atmos_kg_a / CO2_total_kg, color=red, linestyle='--', label='Atmos' )
    h3, = ax1.semilogx( timeMyr_a, (H2O_liquid_kg_a+H2O_solid_kg_a) / H2O_total_kg, color=blue, linestyle='-', label='Magma' )
    h4, = ax1.semilogx( timeMyr_a, H2O_atmos_kg_a / H2O_total_kg, color=blue, linestyle='--', label='Atmos')
    fig_o.set_myaxes( ax1, title=title, ylabel='$x$', xticks=xticks )
    handle_l = [h1,h2]
    fig_o.set_mylegend( ax1, handle_l, loc='center left', ncol=1, TITLE="" )
    ax1.yaxis.set_label_coords(-0.1,0.46)

    ##########
    # figure c
    ##########
    title = '(c) Surface temperature'
    ylabel = '$T$'
    ax2.semilogx( timeMyr_a, temperature_surface_a, 'k-' )
    fig_o.set_myaxes( ax2, title=title, xlabel=xlabel, ylabel=ylabel, xticks=xticks )
    ax2.yaxis.set_label_coords(-0.1,0.5)

    ##########
    # figure d
    ##########
    title = '(d) Emissivity'
    ylabel = '$\epsilon$'
    ax3.semilogx( timeMyr_a, emissivity_a, 'k-' )
    fig_o.set_myaxes( ax3, title=title, xlabel=xlabel, ylabel=ylabel, xticks=xticks )
    ax3.yaxis.set_label_coords(-0.1,0.55)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    fig_o.savefig(6)

#====================================================================
def figure7( args ):

    '''Extract magma ocean depth for partitioning of moderately
       siderophile elements'''

    # TODO: this is an intermediary plot that has been superseded
    # by fig. 8.  It can probably be removed soon, once I have the
    # desired plotting for fig. 8.

    width = 4.7747
    height = 4.7747
    fig_o = FigureData( args, 2, 2, width, height )

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]

    handle_l = []

    time = fig_o.time[0]
    myjson_o = MyJSON( 'output/{}.json'.format(time) )
    xx_pres_b = myjson_o.get_scaled_field_values_internal('pressure_b')
    xx_pres_b *= 1.0E-9

    for nn, time in enumerate( fig_o.time ):
        # read json
        myjson_o = MyJSON( 'output/{}.json'.format(time) )

        color = fig_o.get_color( nn )
        # use melt fraction to determine mixed region
        yy = myjson_o.get_scaled_field_values_internal( 'phi_b' )
        MIX = get_mix( yy )

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
def ppm_to_mass_fraction( ppm ):

    return ppm / 1.0E6

#====================================================================
def get_mo_PT_conditions( fig_o, time_l, field1, criteria, label="", color='black' ):

    data_l = []

    time = fig_o.time[0]
    myjson_o = MyJSON( 'output/{}.json'.format(time) )
    xx_pres_b = myjson_o.get_scaled_field_values('pressure_b')
    xx_pres_b *= 1.0E-9

    for nn, time in enumerate( time_l ):

        myjson_o = MyJSON( 'output/{}.json'.format(time) )

        yy1 = myjson_o.get_scaled_field_values( field1 )
        try:
            yy2 = float(criteria)
        except TypeError:
            yy2 = myjson_o.get_scaled_field_values( criteria )

        # always need temperature
        y_temp = myjson_o.get_scaled_field_values( 'temp_b' )

        # find the pressure(s) that correspond(s) to the cross-over point
        pres = find_value( xx_pres_b, yy1, yy2 )
        ind = get_first_non_zero_index( pres )

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
    sin = get_first_non_zero_index( diff_a )
    ein = 1 + len( diff_a ) - get_first_non_zero_index( diff_a[::-1] )

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
def figure8( args ):

    width = 4.7747#/2
    height = 4.7747/2
    fig_o = FigureData( args, 1, 2, width, height )

    ax0 = fig_o.ax[0]
    ax1 = fig_o.ax[1]

    handle_l = []

    fig_o.time = get_all_output_times()

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
def find_value( xx, yy, yywant ):

    a = yy - yywant

    s = sign_change( a )

    # for ease, just add zero at the beginning to enable us to
    # have the same length array.  Could equally add to the end, or
    # interpolate

    s = np.insert(s,0,0)

    result = xx * s

    return result

#====================================================================
def get_first_non_zero_index( myList ):

    # https://stackoverflow.com/questions/19502378/python-find-first-instance-of-non-zero-number-in-list

    index = next((i for i, x in enumerate(myList) if x), None)

    return index

#====================================================================
def sign_change( a ):

    s = (np.diff(np.sign(a)) != 0)*1

    return s

#====================================================================
def get_single_values_for_times( field, time_l ):

    data_l = []

    for nn, time in enumerate( time_l ):

        # read json
        myjson_o = MyJSON( 'output/{}.json'.format(time) )

        yy = myjson_o.get_scaled_field_values( field )
        data_l.append( yy )

    data_a = np.array( data_l )

    return data_a

#====================================================================
def get_all_output_times():

    '''get all times (in Myrs) from the json files located in the
       output directory'''

    # locate times to process based on files located in output/
    odir = 'output'
    file_l = [f for f in os.listdir(odir) if os.path.isfile(os.path.join(odir,f))]
    if not file_l:
        print('ERROR: output directory contains no files')
        sys.exit(1)

    time_l = [fname for fname in file_l]
    time_l = list(filter(lambda a: a.endswith('json'), time_l))
    time_l = [int(time.split('.json')[0]) for time in time_l]
    # ascending order
    time_l = sorted( time_l, key=int)
    time_a = np.array( time_l )

    return time_a

#====================================================================
def main():

    # arguments (run with -h to summarize)
    parser = argparse.ArgumentParser(description='SPIDER plotting script')
    parser.add_argument('times',type=str, help='Comma-separated (no spaces) list of times');
    parser.add_argument('-f1', '--fig1', help='Plot figure 1', action="store_true")
    parser.add_argument('-f2', '--fig2', help='Plot figure 2', action="store_true")
    parser.add_argument('-f3', '--fig3', help='Plot figure 3', action="store_true")
    parser.add_argument('-f4', '--fig4', help='Plot figure 4', action="store_true")
    parser.add_argument('-f5', '--fig5', help='Plot figure 5', action="store_true")
    parser.add_argument('-f6', '--fig6', help='Plot figure 6', action="store_true")
    parser.add_argument('-f7', '--fig7', help='Plot figure 7', action="store_true")
    parser.add_argument('-f8', '--fig8', help='Plot figure 8', action="store_true")
    args = parser.parse_args()

    # If nothing specified, choose a default set
    if not (args.fig1 or args.fig2 or args.fig3 or args.fig4 or args.fig5 or args.fig6 or args.fig7 or args.fig8) :
        args.fig3 = True;
        args.fig4 = True;
        args.fig5 = True;

    # Old-style arguments as expected by the plotting functions
    old_args = [None,args.times]

    # setup logger

    # create logger with 'main'
    #logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    fh = logging.FileHandler('main.log')
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

    # simplified model in Bower et al. (2018)
    # i.e., figs 1,2
    if args.fig1 :
        bower_et_al_2018_fig1( old_args )
    if args.fig2 :
        bower_et_al_2018_fig2( old_args )

    # reproduce staple figures in Bower et al. (2018)
    # i.e., figs 3,4,5,6,7,8
    if args.fig3 :
        bower_et_al_2018_fig3( old_args )
    if args.fig4 :
        bower_et_al_2018_fig4( old_args )
    if args.fig5 :
        bower_et_al_2018_fig5( old_args )

    # additional figures
    if args.fig6 :
        figure6( old_args )
    if args.fig7 :
        figure7( old_args )
    if args.fig8 :
        figure8( old_args )

    plt.show()

#====================================================================

if __name__ == "__main__":

    main()
