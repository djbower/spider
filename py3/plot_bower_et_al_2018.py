#!/usr/bin/env python

import spider_utils as su
import argparse
import logging
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import os
import sys

logger = su.get_my_logger(__name__)

# global constants are the scalings used in Bower et al. (2018)
# note that these could change for other models
radius0 = 6371000.0
entropy0 = 2993.025100070677
temperature0 = 4033.6070755893948
density0 = 4613.109568155063
dSdr0 = entropy0 / radius0
flux0 = density0 * (entropy0 * temperature0)**(3.0/2.0)

#====================================================================
def dep2pres( dep ):

    '''Adams-Williamson EOS used to relate pressure and depth'''

    rho_surface = 4078.95095544
    gravity = 10.0
    beta = 1.1115348931000002e-07

    dep *= radius0
    pres = rho_surface*gravity / beta
    pres *= np.exp( beta * dep) - 1.0 
    pres *= 1.0E-9 # to GPa
    return pres

#====================================================================
def bower_et_al_2018_fig1():

    logger.info( 'building bower_et_al_2018_fig1' )

    # Note: this function expects that the user has set SPIDERPATH
    # so that it can find data files
    spiderpath = os.environ.get('SPIDERPATH')
    if not spiderpath :
        raise RuntimeError('You must define SPIDERPATH in your environment to point to the root directory of the SPIDER installation');

    prefix = os.path.join(spiderpath,'py3/bower_et_al_2018/simplified_model/fig1')

    width = 4.7747 * 0.5 
    height = 4.7747 * 0.5 
    fig_o = su.FigureData( 1, 1, width, height, 'bower_et_al_2018_fig1')

    ax0 = fig_o.ax

    const = dSdr0

    dSdr_const = 1.0E15
    dSdr_fmt = su.MyFuncFormatter( dSdr_const )
    yticks = [1.0E-15, 1.0E-12, 1.0E-9,1.0E-6,1.0E-3]
    #yticks = [-1E-3,-1E-6,-1E-9,-1E-12,-1E-15]

    dSliqdr_const = 1.0E6 #1.0E15
    dSliqdr_fmt = su.MyFuncFormatter( dSliqdr_const )
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
        xx2 = dSliqdr_fmt.ascale( xx2 )
        yy2 = dSdr_fmt.ascale( yy2 )
        ax0.plot( xx2, yy2, 'k--' )
        ax0.fill_between( xx2, yy2, 0, facecolor='0.75' )

    # plot flux minimum
    if 1:
        #xx2 = np.linspace( -20.4, 0.8, 1000000 )
        xx2 = np.linspace( -20.4, -0.0001, 1000000 )
        xx2 *= const
        yy2 = 1.0/6.0 * xx2 * -1.0
        xx2 = dSliqdr_fmt.ascale( xx2 )
        yy2 = dSdr_fmt.ascale( yy2 )
        xx2 -= 0.4
        yy2 -= 0.4
        ax0.plot( xx2, yy2, color='0.25', linestyle=':' )

    # positive heat fluxes
    for cc, nn in enumerate([12,9,6]):
        filein = os.path.join(prefix,'dSdr_10p{:d}_processed.dat'.format(nn))
        xx, yy = np.loadtxt( filein, unpack=True )
        xx *= const
        xx = dSliqdr_fmt.ascale( xx )
        yy *= const * -1.0
        yy = dSdr_fmt.ascale( yy )
        # offset y for visual clarity
        xx += 0.4
        yy += 0.4
        ax0.plot( xx, yy, '-', color=fig_o.get_color(cc) ) #, 'k-' )

    # negative heat fluxes
    for cc, nn in enumerate([12,9,6]):
        filein = os.path.join(prefix,'dSdr_n10p{:d}_processed.dat'.format(nn))
        xx, yy = np.loadtxt( filein, unpack=True )
        xx *= const
        xx = dSliqdr_fmt.ascale( xx )
        yy *= const * -1.0
        yy = dSdr_fmt.ascale( yy )
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
def bower_et_al_2018_fig2():

    logger.info( 'building bower_et_al_2018_fig2' )

    # Note: this function expects that the user has set SPIDERPATH
    # so that it can find data files
    spiderpath = os.environ.get('SPIDERPATH')
    if not spiderpath :
        raise RuntimeError('You must define SPIDERPATH in your environment to point to the root directory of the SPIDER installation');

    prefix = os.path.join(spiderpath,'py3/bower_et_al_2018/simplified_model/fig2')

    width = 4.7747 # * 0.5
    height = 4.7747 # * 0.5
    fig_o = su.FigureData( 2, 2, width, height, 'bower_et_al_2018_fig2' )

    dd = fig_o.data_d

    # below, only basic internal nodes

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]
    ax2 = fig_o.ax[1][0]
    ax3 = fig_o.ax[1][1]

    ###############
    # plot liquidus
    Sliq_file = os.path.join(prefix,'Sliq_processed.dat')
    rad, Sliq = np.loadtxt( Sliq_file, unpack=True )
    depth = 1.0 - rad
    pres = dep2pres( depth )
    const = entropy0 # scaling for entropy but could change
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
    const = dSdr0

    dSdr_file = os.path.join(prefix,'dSdr_processed.dat')
    rad, dSdr = np.loadtxt( dSdr_file, unpack=True )
    depth = 1.0 - rad
    pres = dep2pres( depth )
    dSdr *= const * -1.0

    # dSliqdr was not plotted in the end
    #dSliqdr_file = os.path.join(prefix,'dSliqdr_processed.dat')
    #rad2, dSliqdr = np.loadtxt( dSliqdr_file, unpack=True )
    #depth2 = 1.0 - rad2
    #pres2 = dep2pres( depth2 )
    #dSliqdr *= const
    # must scale by 1/2
    #dSliqdr *= 0.5

    dSdr_const = 1.0E15
    dSdr_fmt = su.MyFuncFormatter( dSdr_const )

    dSdr = dSdr_fmt.ascale( dSdr )
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
    pres = dep2pres( depth )

    const = flux0
    Fconv *= const

    # for arcsinh scaling
    flux_const = 1.0E6
    flux_fmt = su.MyFuncFormatter( flux_const )
    Fconv = flux_fmt.ascale( Fconv )
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
    pres = dep2pres( depth )

    const = flux0
    Fmix *= const

    Fmix = flux_fmt.ascale( Fmix )
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
    Ftot = flux_fmt.ascale( Ftot )
    ax0.plot( pres, Ftot, 'k--' )
    ax2.plot( pres, Ftot, 'k--' )

    fig_o.savefig(4)

#====================================================================
def bower_et_al_2018_fig3( times ):

    # article class text width is 4.7747 inches
    # http://tex.stackexchange.com/questions/39383/determine-text-width

    logger.info( 'building bower_et_al_2018_fig3' )

    fig_o = su.FigureData( 2, 2, 4.7747, 4.7747, 'bower_et_al_2018_fig3', times )

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]
    ax2 = fig_o.ax[1][0]
    ax3 = fig_o.ax[1][1]

    time = fig_o.time[0] # first timestep since liquidus and solidus
                         # are time-independent

    myjson_o = su.MyJSON( 'output/{}.json'.format(time) )

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
        myjson_o = su.MyJSON( 'output/{}.json'.format(time) )

        color = fig_o.get_color( nn )
        # use melt fraction to determine mixed region
        MIX = myjson_o.get_mixed_phase_boolean_array( 'basic_internal' )
        MIX_s = myjson_o.get_mixed_phase_boolean_array( 'staggered' )

        label = fig_o.get_legend_label( time )

        # entropy
        # plot staggered, since this is where entropy is defined
        # cell-wise and we can easily see the CMB boundary condition
        yy = myjson_o.get_scaled_field_values('S_s')
        #yy = myjson_o.get_scaled_field_values_internal('S_b')
        ax0.plot( xx_pres_s, yy, '--', color=color )
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
        visc_fmt = su.MyFuncFormatter( visc_const )
        yy = myjson_o.get_scaled_field_values_internal('visc_b', visc_fmt)
        ax3.plot( xx_pres, yy, '-', color=color )

    xticks = [0,50,100,135]
    xmax = 138

    # titles and axes labels, legends, etc
    units = myjson_o.get_field_units('S_b')
    title = '(a) Entropy, {}'.format(units)
    #yticks = [1000,2000,2400,2800,3200]
    # Bower et al. (2018)
    yticks = [1600,2000,2400,2800,3200]
    # DJB used this next range for work with Bayreuth
    #yticks = [300,1000,1600,2000,2400,2800,3200]
    fig_o.set_myaxes( ax0, title=title, ylabel='$S$', xticks=xticks, xmax=xmax, yticks=yticks )
    ax0.yaxis.set_label_coords(-0.075,0.59)
    units = myjson_o.get_field_units('temp_b')
    title = '(b) Temperature, {}'.format(units)
    # Bower et al. (2018)
    yticks= [1000,2000,3000,4000,5000]
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
    yticks = [1.0E2, 1.0E6, 1.0E12, 1.0E18, 1.0E21]
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
def bower_et_al_2018_fig4( times ):

    logger.info( 'building bower_et_al_2018_fig4' )

    fig_o = su.FigureData( 2, 2, 4.7747, 4.7747, 'bower_et_al_2018_fig4', times )

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]
    ax2 = fig_o.ax[1][0]
    ax3 = fig_o.ax[1][1]

    handle_l = []

    # for arcsinh scaling
    flux_const = 1.0E6
    flux_fmt = su.MyFuncFormatter( flux_const )

    for nn, time in enumerate( fig_o.time ):
        # In Bower et al. (2018), just every other line
        # is plotted for visual clarity (uncomment below)
        if (nn-1) % 2:
            continue

        # read json
        myjson_o = su.MyJSON( 'output/{}.json'.format(time) )

        color = fig_o.get_color( nn )
        # use melt fraction to determine mixed region
        MIX = myjson_o.get_mixed_phase_boolean_array( 'basic_internal')

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
def bower_et_al_2018_fig5( times ):

    logger.info( 'building bower_et_al_2018_fig5' )

    # keep y the same by scaling (7.1621)
    fig_o = su.FigureData( 3, 2, 4.7747, 7.1621, 'bower_et_al_2018_fig5', times )

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]
    ax2 = fig_o.ax[1][0]
    ax3 = fig_o.ax[1][1]
    ax4 = fig_o.ax[2][0]
    ax5 = fig_o.ax[2][1]

    handle_l = []

    eddy_const = 1.0
    eddy_fmt = su.MyFuncFormatter( eddy_const )
    alpha_const = 1.0E6
    alpha_fmt = su.MyFuncFormatter( alpha_const )
    dSdr_const = 1.0E15
    dSdr_fmt = su.MyFuncFormatter( dSdr_const )

    for nn, time in enumerate( fig_o.time ):
        # In Bower et al. (2018), just every other line
        # is plotted for visual clarity (uncomment below)
        if (nn-1) % 2:
            continue

        # read json
        myjson_o = su.MyJSON( 'output/{}.json'.format(time) )

        color = fig_o.get_color( nn )
        # use melt fraction to determine mixed region
        MIX = myjson_o.get_mixed_phase_boolean_array( 'basic_internal' )

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
def main():

    # arguments (run with -h to summarize)
    parser = argparse.ArgumentParser(description='SPIDER plotting script')
    parser.add_argument('-t', '--times', type=str, help='Comma-separated (no spaces) list of times');
    parser.add_argument('-f1', '--fig1', help='Plot figure 1', action="store_true")
    parser.add_argument('-f2', '--fig2', help='Plot figure 2', action="store_true")
    parser.add_argument('-f3', '--fig3', help='Plot figure 3', action="store_true")
    parser.add_argument('-f4', '--fig4', help='Plot figure 4', action="store_true")
    parser.add_argument('-f5', '--fig5', help='Plot figure 5', action="store_true")
    args = parser.parse_args()

    # if nothing specified, choose a default set
    if not (args.fig1 or args.fig2 or args.fig3 or args.fig4 or args.fig5):
        args.fig3 = True;
        args.fig4 = True;
        args.fig5 = True;

    if args.fig3 or args.fig4 or args.fig5:
        if not args.times:
            logger.critical( 'You must specify times in a comma-separated list (no spaces) using -t' )
            sys.exit(0)

    # simplified model in Bower et al. (2018)
    # i.e., figs 1,2
    if args.fig1 :
        bower_et_al_2018_fig1()
    if args.fig2 :
        bower_et_al_2018_fig2()

    # reproduce staple figures in Bower et al. (2018)
    # i.e., figs 3,4,5,6,7,8
    if args.fig3 :
        bower_et_al_2018_fig3( times=args.times )
    if args.fig4 :
        bower_et_al_2018_fig4( times=args.times )
    if args.fig5 :
        bower_et_al_2018_fig5( times=args.times )

    plt.show()

#====================================================================

if __name__ == "__main__":

    main()
