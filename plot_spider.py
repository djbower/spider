#!/usr/bin/env python

import json
import logging
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import os, sys 

logger = logging.getLogger(__name__)

#===================================================================
class MyFuncFormatter( object ):

    '''the default function formatter from
       matplotlib.ticker.FuncFormatter(func) only accepts two
       arguments, which is not enough to scale an arcsinh function.
       But by creating our own class here we can attach the scaling
       to the object which can then be accessed in __call__'''

    def __init__( self, arcsinh_scale ):
        self.const = arcsinh_scale

    def __call__( self, x, pos ):
        y = invascale( x, self.const )
        fmt = sci_notation( y, 0 )
        return fmt

#====================================================================
class MyJSON( object ):

    '''load and access json data'''

    def __init__( self, filename ):
        self.filename = filename
        self._load()

    def _load( self ):
        '''load and store json data from file'''
        try:
            json_data  = open( self.filename )
        except FileNotFoundError:
            logger.critical('cannot find file: %s', self.filename )
            logger.critical('please specify times for which data exists')
            sys.exit(1)
        self.data_d = json.load( json_data )
        json_data.close()

    def get_field_data( self, field ):
        '''get all data relating to a particular field'''
        field_l = self.data_d['data']
        for ii in range( len(field_l) ):
            if field_l[ii]['name'] == field:
                fdata_d = field_l[ii]
        try:
            return fdata_d
        except NameError:
            logger.critical('data for %s does not exist', field )
            sys.exit(1)

    def get_field_units( self, field ):
        '''get the units (SI) of a particular field'''
        fdata_d = self.get_field_data( field )
        units = fdata_d['units']
        units = None if units == 'None' else units
        return units

    def get_scaled_field_values( self, field, fmt_o='' ):
        '''get the scaled values for a particular field'''
        fdata_d = self.get_field_data( field )
        scaling = fdata_d['scaling']
        values_a = np.array( [float(value) for value in fdata_d['values']] )
        scaled_values_a = scaling * values_a
        if fmt_o:
            scaled_values_a = ascale( scaled_values_a, fmt_o.const )
        return scaled_values_a

    def get_scaled_field_values_internal( self, field, fmt_o='' ):
        '''get the scaled values for the internal nodes (ignore top
           and bottom nodes)'''
        scaled_values_a = self.get_scaled_field_values( field, fmt_o )
        return scaled_values_a[1:-1]
 
#====================================================================
class FigureData( object ):

    def __init__( self, args, nrows, ncols, width, height, outname='fig' ):
        dd = {}
        self.data_d = dd
        dd['time_l'] = args[1]
        dd['time_units'] = 'kyr' # hard-coded
        dd['time_decimal_places'] = 2 # hard-coded
        dd['outname'] = outname
        self.process_time_list()
        self.set_properties( nrows, ncols, width, height )

    def get_color( self, num ):
        dd = self.data_d
        return dd['colors_l'][num]

    def get_legend_label( self, time ):
        dd = self.data_d
        units = dd['time_units']
        dp = dd['time_decimal_places']
        age = float(time)
        if units == 'yr':
            age = round( age, 0 )
            label = '%d'
        elif units == 'kyr':
            age /= 1.0E3
            label = '%0.1f'
        elif units == 'Myr':
            age /= 1.0E6
            label = '%0.2f'
        elif units == 'Byr':
            age /= 1.0E9
            label = '%0.2f'
        #label = '%0.{}e'.format( dp )
        #label = '%0.{}f'.format( dp )
        label = label % age
        return label

    def process_time_list( self ):
        dd = self.data_d
        time_l = dd['time_l']
        try:
            time_l = [int(time_l)]
        except ValueError:
            time_l = [int(time) for time in time_l.split(',')]
        self.time = time_l

    def make_figure( self ):
        dd = self.data_d
        nrows = dd['nrows']
        ncols = dd['ncols']
        fig, ax = plt.subplots( nrows, ncols )
        fig.subplots_adjust(wspace=0.3,hspace=0.3)
        fig.set_size_inches( dd['width'], dd['height'] )
        self.fig = fig
        self.ax = ax

    def savefig( self, num ):
        dd = self.data_d
        if dd['outname']:
            outname = dd['outname'] + '.pdf'
        else:
            outname = 'fig{}.pdf'.format( num)
        self.fig.savefig(outname, transparent=True, bbox_inches='tight',
            pad_inches=0.05, dpi=dd['dpi'])

    def set_colors( self, num=6 ):
        dd = self.data_d
        cmap = plt.get_cmap('jet')
        #colors_l = [cmap(i) for i in np.linspace(0, 1, num)]
        # color scheme from Tim.  Nice reds and blues
        #colors_l = ['#2364A4',
        #            '#1695F9',
        #            '#95D5FD',
        #            '#8B0000',
        #            '#CD5C5C',
        #            '#FA141B',
        #            '#FFA07A']
        # color scheme 'bkr8' for light background from Crameri
        # see f_Colours.m at http://www.fabiocrameri.ch/visualisation.php
        # this is actually very similar (same?) as Tim's scheme above
        # used in Bower et al. (2018)
        colors_l = [(0.0,0.0,0.3),
                    (0.1,0.1,0.5),
                    (0.2,0.2,0.7),
                    (0.4,0.4,0.8),
                    (0.8,0.4,0.4),
                    (0.7,0.2,0.2),
                    (0.5,0.1,0.1),
                    (0.3,0.0,0.0)]
        colors_l.reverse()
        dd['colors_l'] = colors_l

    def set_properties( self, nrows, ncols, width, height ):
        dd = self.data_d
        dd['nrows'] = nrows
        dd['ncols'] = ncols
        dd['width'] = width # inches
        dd['height'] = height # inches
        # TODO: breaks for MacOSX, since I don't think Mac comes
        # with serif font.  But whatever it decides to switch to
        # also looks OK and LaTeX-like.
        font_d = {'family' : 'serif',
                  #'style': 'normal',
                  #'weight' : 'bold'
                  'serif': ['Computer Modern Roman'],
                  'size'   : '8'}
        mpl.rc('font', **font_d)
        # Use TeX font for labels etc.
        plt.rc( 'text', usetex=True )
        dd['dpi'] = 300
        dd['extension'] = 'png'
        dd['fontsize_legend'] = 6
        dd['fontsize_title'] = 8
        dd['fontsize_xlabel'] = 8
        dd['fontsize_ylabel'] = 8
        self.set_colors( len(self.time) )
        self.make_figure()

    def set_myaxes( self, ax, title='', xlabel='', xticks='',
        ylabel='', yticks='', fmt='', xfmt='' ):
        if title:
            self.set_mytitle( ax, title )
        if xlabel:
            self.set_myxlabel( ax, xlabel )
        if xticks:
            self.set_myxticks( ax, xticks, xfmt )
        if ylabel:
            self.set_myylabel( ax, ylabel )
        if yticks:
            self.set_myyticks( ax, yticks, fmt )

    def set_mylegend( self, ax, handles, loc=4, ncol=1, TITLE=1 ):
        dd = self.data_d
        units = dd['time_units']
        if TITLE==1:
            title = r'Time ({0})'.format( units )
        else:
            title = TITLE
        fontsize = self.data_d['fontsize_legend']
        legend = ax.legend(title=title, handles=handles, loc=loc,
            ncol=ncol, fontsize=fontsize)
        plt.setp(legend.get_title(),fontsize=fontsize)

    def set_mytitle( self, ax, title ):
        dd = self.data_d
        fontsize = dd['fontsize_title']
        title = r'{}'.format( title )
        ax.set_title( title, fontsize=fontsize )

    def set_myxlabel( self, ax, label ):
        dd = self.data_d
        fontsize = dd['fontsize_xlabel']
        label = r'{}'.format( label )
        ax.set_xlabel( label, fontsize=fontsize )

    def set_myylabel( self, ax, label ):
        dd = self.data_d
        fontsize = dd['fontsize_ylabel']
        rotation = 'horizontal'
        label = r'{}'.format( label )
        ax.set_ylabel( label, fontsize=fontsize, rotation=rotation )

    def set_myxticks( self, ax, xticks, fmt ):
        dd = self.data_d
        if fmt:
            xticks = ascale( np.array(xticks), fmt.const )
            ax.xaxis.set_major_formatter(
                mpl.ticker.FuncFormatter(fmt))
        ax.set_xticks( xticks)
        # set x limits to match extent of ticks
        ax.set_xlim( xticks[0], xticks[-1] )

    def set_myyticks( self, ax, yticks, fmt ):
        dd = self.data_d
        if fmt:
            yticks = ascale( np.array(yticks), fmt.const )
            ax.yaxis.set_major_formatter(
                mpl.ticker.FuncFormatter(fmt))
        ax.set_yticks( yticks)
        # set y limits to match extent of ticks
        ax.set_ylim( yticks[0], yticks[-1] )

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

        label = fig_o.get_legend_label( time )

        # entropy
        yy = myjson_o.get_scaled_field_values_internal('S_b')
        ax0.plot( xx_pres, yy, '--', color=color )
        handle, = ax0.plot( xx_pres*MIX, yy*MIX, '-', color=color, label=label )
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

    # titles and axes labels, legends, etc
    units = myjson_o.get_field_units('S_b')
    title = '(a) Entropy, {}'.format(units)
    yticks = [1600,2000,2400,2800,3200]
    # DJB used this next range for work with Bayreuth
    #yticks = [300,1000,1600,2000,2400,2800,3200]
    fig_o.set_myaxes( ax0, title=title, ylabel='$S$', xticks=xticks, yticks=yticks )
    ax0.yaxis.set_label_coords(-0.075,0.59)
    units = myjson_o.get_field_units('temp_b')
    title = '(b) Temperature, {}'.format(units)
    yticks= [1000,2000,3000,4000,5000]
    # DJB used this next range for work with Bayreuth
    #yticks= [300,1000,2000,3000,4000,5000]
    fig_o.set_myaxes( ax1, title=title, ylabel='$T$', xticks=xticks, yticks=yticks )
    ax1.yaxis.set_label_coords(-0.075,0.59)
    fig_o.set_mylegend( ax1, handle_l, loc=4, ncol=2 )
    #fig_o.set_mylegend( ax0, handle_l, loc=2, ncol=2 )
    title = '(c) Melt fraction'
    yticks = [0,0.2,0.4,0.6,0.8,1.0]
    fig_o.set_myaxes( ax2, title=title, xlabel='$P$ (GPa)',
        ylabel='$\phi$', xticks=xticks, yticks=yticks )
    ax2.yaxis.set_label_coords(-0.075,0.475)
    ax2.set_ylim( [0, 1] )
    units = myjson_o.get_field_units('visc_b')
    title = '(d) Viscosity, ' + units
    yticks = [1.0E2, 1.0E6, 1.0E12, 1.0E18, 1.0E21]
    fig_o.set_myaxes( ax3, title=title, xlabel='$P$ (GPa)',
        ylabel='$\eta$', xticks=xticks, yticks=yticks, fmt=visc_fmt )
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
    yticks = [-1E15,-1E12, -1E6, -1E0, 0, 1E0, 1E6, 1E12, 1E15]
    units = myjson_o.get_field_units('Jconv_b')
    title = '(a) Convective flux, {}'.format(units)
    fig_o.set_myaxes( ax0, title=title, ylabel='$F_\mathrm{conv}$',
        yticks=yticks, xticks=xticks, fmt=flux_fmt )
    ax0.yaxis.set_label_coords(-0.16,0.54)
    fig_o.set_mylegend( ax0, handle_l, ncol=2 )
    units = myjson_o.get_field_units('Jmix_b')
    title = '(c) Mixing flux, {}'.format(units)
    fig_o.set_myaxes( ax2, title=title, xlabel='$P$ (GPa)',
        ylabel='$F_\mathrm{mix}$', yticks=yticks, xticks=xticks, fmt=flux_fmt )
    ax2.yaxis.set_label_coords(-0.16,0.54)
    units = myjson_o.get_field_units('Jgrav_b')
    title = '(b) Separation flux, {}'.format(units)
    fig_o.set_myaxes( ax1, title=title,
        ylabel='$F_\mathrm{grav}$', yticks=yticks, xticks=xticks, fmt=flux_fmt )
    ax1.yaxis.set_label_coords(-0.16,0.54)
    units = myjson_o.get_field_units('Jtot_b')
    title = '(d) Total flux, {}'.format(units)
    fig_o.set_myaxes( ax3, title=title, xlabel='$P$ (GPa)',
        ylabel='$F_\mathrm{tot}$', yticks=yticks, xticks=xticks, fmt=flux_fmt )
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

    # titles and axes labels, legends, etc.
    units = myjson_o.get_field_units('kappah_b')
    title = '(a) Eddy diffusivity, {}'.format(units)
    yticks = [1.0, 1.0E3, 1.0E6, 1.0E9]
    fig_o.set_myaxes( ax0, title=title, ylabel='$\kappa_h$',
        yticks=yticks, fmt=eddy_fmt, xticks=xticks )
    ax0.yaxis.set_label_coords(-0.1,0.475)

    units = myjson_o.get_field_units('rho_b')
    title = '(b) Density, {}'.format(units)
    yticks = [2000,3000,4000,5000,6000]
    fig_o.set_myaxes( ax1, title=title, ylabel='$\\rho$',
        yticks=yticks, xticks=xticks )
    ax1.yaxis.set_label_coords(-0.1,0.6)
    fig_o.set_mylegend( ax1, handle_l, loc=4, ncol=2 )

    units = myjson_o.get_field_units('dSdr_b')
    title = '(c) Entropy grad, {}'.format(units)
    #yticks= [-1.0E-3, -1.0E-6, -1.0E-9, -1.0E-12, -1.0E-15]
    yticks = [1E-15, 1E-12, 1E-9, 1E-6, 1E-3]
    fig_o.set_myaxes( ax2, title=title,
        yticks=yticks, xticks=xticks,
        ylabel='$-\\frac{\partial S}{\partial r}$', fmt=dSdr_fmt)
    ax2.yaxis.set_label_coords(-0.1,0.565)

    fig_o.savefig(3)


    units = myjson_o.get_field_units('alpha_b')
    title = '(d) Thermal expansion, {}'.format(units)
    yticks = [1.0E-5, 1.0E-4, 1.0E-3, 1.0E-2]
    fig_o.set_myaxes( ax3, title=title, ylabel='$\\alpha$',
        yticks=yticks, fmt=alpha_fmt, xticks=xticks )
    ax3.yaxis.set_label_coords(-0.1,0.475)

    units = myjson_o.get_field_units('cp_b')
    title = '(e) Heat capacity, {}'.format(units)
    yticks = [0,5000,10000,15000]
    fig_o.set_myaxes( ax4, title=title, xlabel='$P$ (GPa)', ylabel='$c$',
        yticks=yticks, xticks=xticks )
    ax4.yaxis.set_label_coords(-0.1,0.475)

    units = myjson_o.get_field_units('dTdrs_b')
    title = '(f) Adiabatic grad, {}'.format(units)
    yticks = [-3E-3, -2E-3, -1E-3, 0]
    fig_o.set_myaxes( ax5, title=title, xlabel='$P$ (GPa)',
        yticks=yticks, xticks=xticks,
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

    width = 4.7747
    height = 4.7747
    fig_o = FigureData( args, 2, 2, width, height )

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]
    ax2 = fig_o.ax[1][0]
    ax3 = fig_o.ax[1][1]

    data = atmosphere_data_to_array( fig_o )

    # bit clunky, but works for now
    time = data[:,0] * 1.0E-6 # in Myrs

    Mliq = data[:,1]
    Msol = data[:,2]
    Mtot = data[:,3]
    tsurf = data[:,4]
    tau = data[:,5]
    emissivity = data[:,6]
    CO2liq = data[:,7]
    CO2sol = data[:,8]
    CO2atm = data[:,9]
    CO2tot = data[:,10]
    CO2p = data[:,11]
    CO2dpdx = data[:,12]
    CO2tau = data[:,13]
    H2Oliq = data[:,14]
    H2Osol = data[:,15]
    H2Oatm = data[:,16]
    H2Otot = data[:,17]
    H2Op = data[:,18]
    H2Odpdx = data[:,19]
    H2Otau = data[:,20]

    # mass fractions
    CO2liq /= CO2tot
    CO2sol /= CO2tot
    CO2atm /= CO2tot

    # mass fractions
    H2Oliq /= H2Otot
    H2Osol /= H2Otot
    H2Oatm /= H2Otot

    xticks = [1E-2,1E-1,1E0,1E1,1E2]
    xlabel = 'Time (Myr)'

    red = fig_o.get_color(3)
    blue = fig_o.get_color(-3)

    ### new to make regime plot
    yy = Mliq/Mtot
    xx = time
    #cont_l = np.linspace(0.2,0.8,3)
    cont_l = [0.2,0.8]#,0.8]
    val_l = []
    for cont in cont_l:
        # this is probbaly broken now
        item, yval = find_value( xx, yy, cont )
        val_l.append( yval )

    # figure a
    title = '(a) Partial pressure (bar)'
    ylabel = '$p$'
    # for testing
    #h3, = ax0.semilogx( time, 400*Mliq/Mtot, 'k--' )
    for cc, cont in enumerate(cont_l):
        h3, = ax0.semilogx( [val_l[cc],val_l[cc]],[0,400], color='0.25', linestyle=':')
        label = 100*cont
        ax0.text( val_l[cc], 400, '%(label)0d' % vars(), ha='center' )#, rotation='vertical' )
        # for massive atmosphere
        #h3, = ax0.semilogx( [val_l[cc],val_l[cc]],[0,4000], color='0.25', linestyle=':')
        #ax0.text( val_l[cc], 4000, '%(label)0d' % vars(), ha='center' )#, rotation='vertical' )
    #phishift = val_l[cc] - 1.8 #0.135
    #ax0.text( phishift, 4000, '$\phi (\%)$' )
    phishift = val_l[cc] - 0.135
    ax0.text( phishift, 400, '$\phi (\%)$' )
    h1, = ax0.semilogx( time, CO2p*1.0E-5, color=red, linestyle='-', label='CO$_2$')
    h2, = ax0.semilogx( time, H2Op*1.0E-5, color=blue, linestyle='-', label='H$_2$O')
    handle_l = [h1,h2]
    fig_o.set_myaxes( ax0, title=title, ylabel=ylabel, xticks=xticks )
    fig_o.set_mylegend( ax0, handle_l, loc='center right', ncol=1, TITLE="" )
    ax0.yaxis.set_label_coords(-0.1,0.575)

    # figure b
    title = '(b) Volatile reservoirs'
    #h5, = ax0.semilogx( time, Mliq/Mtot, 'k--', label='melt' )
    h1, = ax1.semilogx( time, CO2liq+CO2sol, color=red, linestyle='-', label='Magma' )
    #h2, = ax0.semilogx( time, CO2sol, 'b-', label='CO2 sol' )
    h2, = ax1.semilogx( time, CO2atm, color=red, linestyle='--', label='Atmos' )
    h3, = ax1.semilogx( time, H2Oliq+H2Osol, color=blue, linestyle='-')#, label='H2O magma' )
    #h5, = ax0.semilogx( time, H2Osol, 'b--', label='H2O sol' )
    h4, = ax1.semilogx( time, H2Oatm, color=blue, linestyle='--')#, label='H2O atmos' )
    #ax0.semilogx( time, CO2liq+CO2sol+CO2atm, 'k-' )
    fig_o.set_myaxes( ax1, title=title, ylabel='$x$', xticks=xticks )
    handle_l = [h1,h2]#,h3,h4,h5]
    fig_o.set_mylegend( ax1, handle_l, loc='center left', ncol=1, TITLE="" )
    ax1.yaxis.set_label_coords(-0.1,0.46)

    # figure c
    #title = '(c) H$_2$O reservoirs'
    #h4, = ax2.semilogx( time, Mliq/Mtot, 'k--', label='Melt frac' )
    #h1, = ax2.semilogx( time, H2Oliq, 'r-', label='Liquid' )
    #h2, = ax2.semilogx( time, H2Osol, 'b-', label='Solid' )
    #h3, = ax2.semilogx( time, H2Oatm, 'g-', label='Atmos' )
    #ax2.semilogx( time, H2Oliq+H2Osol+H2Oatm, 'k-' )
    #fig_o.set_myaxes( ax2, title=title, xlabel=xlabel, ylabel='$x$', xticks=xticks )
    #handle_l = [h1,h3,h2,h4]
    #fig_o.set_mylegend( ax2, handle_l, loc='center left', ncol=1, TITLE=0 )
    #ax2.yaxis.set_label_coords(-0.1,0.46)

    # figure c
    title = '(c) Surface temperature'
    ylabel = '$T$'
    ax2.semilogx( time, tsurf, 'k-' )
    fig_o.set_myaxes( ax2, title=title, xlabel=xlabel, ylabel=ylabel, xticks=xticks )
    ax2.yaxis.set_label_coords(-0.1,0.5)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    # figure d
    title = '(d) Emissivity'
    ylabel = '$\epsilon$'
    ax3.semilogx( time, emissivity, 'k-' )
    fig_o.set_myaxes( ax3, title=title, xlabel=xlabel, ylabel=ylabel, xticks=xticks )
    ax3.yaxis.set_label_coords(-0.1,0.55)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    fig_o.savefig(6)

#====================================================================
def figure7( args ):

    '''Extract magma ocean depth for partitioning of moderately
       siderophile elements'''

    width = 4.7747
    height = 4.7747
    fig_o = FigureData( args, 2, 2, width, height )

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]

    handle_l = []

    for nn, time in enumerate( fig_o.time ):
        color = fig_o.get_color( nn )
        # use melt fraction to determine mixed region
        xx, yy = fig_o.get_xy_data( 'phi', time )
        MIX = get_mix( yy )

        label = fig_o.get_legend_label( time )

        # regime
        xx, yy = fig_o.get_xy_data( 'regime', time )
        ax0.plot( xx, yy, '--', color=color )
        handle, = ax0.plot( xx*MIX, yy*MIX, '-', color=color, label=label )
        handle_l.append( handle )

        # regime based on melt fraction
        xx, yy = fig_o.get_xy_data( 'phi', time )
        regime = np.copy(yy)
        regime[:] = 0
        regime[yy>=0.4] = 1.0
        regime[yy<0.4] = 2.0 

        ax1.plot( xx, regime, '--', color=color )
        handle, = ax1.plot( xx*MIX, regime*MIX, '-', color=color, label=label )
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

    for nn, time in enumerate( time_l ):

        xx1, yy1 = fig_o.get_xy_data( field1, time )
        try:
            xx2, yy2 = fig_o.get_xy_data( criteria, time )
        except ValueError:
            yy2 = criteria

        # always need temperature
        x_temp, y_temp = fig_o.get_xy_data( 'temp', time )

        # find the pressure(s) that correspond(s) to the cross-over point
        pres = find_value( xx1, yy1, yy2 )
        ind = get_first_non_zero_index( pres )

        current_pres = x_temp[-1]
        current_temp = y_temp[-1]
        try:
            previous_pres = data_l[-1][0]
            if ind is None and previous_pres > 120:
                data_l.append( (current_pres, current_temp, time ) )
            elif ind is None and previous_pres < 20:
                current_pres = x_temp[0]
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

    polydegree = 5 # hard-coded

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

    fig_o.time = range(0,3000,10)
    #fig_o.time = [40,100,200,400,600,800,1000,1400,1800]
    #fig_o.time = [40,100,200,400,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200]

    hs1, hf1 = get_mo_PT_conditions( fig_o, fig_o.time, 'regime', 1.5, r'Polynomial fit', 'blue' )
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
def atmosphere_data_to_array( fig_o ):

    # locate times to process based on files located in output/
    odir = 'output'
    file_l = [f for f in os.listdir(odir) if os.path.isfile(os.path.join(odir,f))]
    if not file_l:
        print('ERROR: output directory contains no files')
        sys.exit(1)

    time_l = [fname.split('.')[-1] for fname in file_l]

    # remove '.info' suffix
    time_l = list(filter(lambda a: a !='info', time_l))
    # remove '.m' suffix
    time_l = list(filter(lambda a: a !='m', time_l))
    # remove gitignore
    time_l = list(filter(lambda a: a !='gitignore', time_l))
    # convert to int
    time_l = [int(val) for val in time_l]
    # ascending order
    time_l = sorted( time_l, key=int)

    data = fig_o.get_data( 'atmosphere_data', time_l[0] )

    # must be some data to process
    if data is None:
        print('atmosphere_data_to_array: no data')
        sys.exit(0)

    ncols = np.size( data ) + 1
    nrows = len( time_l )

    data_a = np.zeros( (nrows, ncols+2) )

    for nn, time in enumerate( time_l ):

        data = fig_o.get_data( 'atmosphere_data', time )
        data_a[nn,0] = time
        data_a[nn,1:-2] = data

        dSdr_b_aug = fig_o.get_data( 'dSdr_b_aug', time )
        data_a[nn,-2] = dSdr_b_aug[0] # x0
        data_a[nn,-1] = dSdr_b_aug[1] # x1

    return data_a

#====================================================================
def ascale( yy, constant ):
    '''convert data to log-like values (scaled arcsinh)'''
    return np.arcsinh( yy*constant )

#====================================================================
def invascale( yy, constant ):
    '''recover original data from log-like values (scaled arcsinh)'''
    return np.sinh(yy) / constant

#====================================================================
# Define function for string formatting of scientific notation
def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """

    # plotting zero is useful to emphasize that we are plotting both
    # positive and negative values, e.g. for the heat fluxes
    if num==0:
        fmt = r"$0$"
        return fmt

    if not exponent:
        exponent = abs(num)
        exponent = np.log10( exponent )
        exponent = np.floor( exponent )
        exponent = int( exponent )

    coeff = round(num / float(10**exponent), decimal_digits)
    # sometimes, probably due to floating point precision? the coeff
    # is not less than ten.  Correct for that here
    if np.abs(coeff) >= 10.0:
        coeff /= 10.0
        exponent += 1
    if not precision:
        precision = decimal_digits

    if coeff < 0.0:
        fmt = r"$-10^{{{0}}}$".format(exponent)
        #fmt= r"${{{0}}}$".format(exponent)
    else:
        fmt = r"$10^{{{0}}}$".format(exponent)

    return fmt
    #return r"${0:.{2}f}\cdot10^{{{1:d}}}$".format(coeff, exponent, precision)

#====================================================================

# REMOVE

#def get_vector_index( instring ):
#    """
#    the order of fields in data_l must match exactly the order in
#    which PetscVecs are written to the output Petsc binary.
#    See monitor.c
#    """

#    data_l = [
        # extraneous
#        'time_data',
#        'constant_data',
#        'atmosphere_data',
#        'dSdr_b_aug',
#        # stored in M->meshVecs_b
#        'area_b',
#        'dPdr_b',
#        'pressure_b',
#        'radius_b',
#        'mix_b',
#        # stored in M->meshVecs_s
#        'pressure_s',
#        'radius_s',
#        'volume_s',
#        'dPdr_s',
#        'area_s',
#        'rho_s',
#        'mass_s',
#        # stored in S->solutionVecs_b
#        'alpha',
#        'alpha_mix',
#        'cond',
#        'cp',
#        'cp_mix',
#        'dfusdr',
#        'dfusdr_temp',
#        'dSdr',
#        'dSliqdr',
#        'dSsoldr',
#        'dTdrs',
#        'dTdrs_mix',
#        'Etot',
#        'fusion',
#        'fusion_curve',
#        'fusion_curve_temp',
#        'fusion_rho',
#        'fusion_temp',
#        'fwtl',
#        'fwts',
#        'gphi',
#        'gsuper',
#        'Jcond',
#        'Jconv',
#        'Jgrav',
#        'Jmix',
#        'Jtot',
#        'kappah',
#        'liquidus',
#        'liquidus_rho',
#        'liquidus_temp',
#        'nu',
#        'phi',
#        'regime',
#        'rho',
#        'S',
#        'solidus',
#        'solidus_rho',
#        'solidus_temp',
#        'temp',
#        'visc',
         # stored in S->solutionVecs_s
#        'cp_s',
#        'cp_mix_s',
#        'dSdt_s',
#        'fusion_s',
#        'fusion_curve_s',
#        'fusion_curve_temp_s',
#        'fusion_temp_s',
#        'fwtl_s',
#        'fwts_s',
#        'gphi_s',
#        'Hradio_s',
#        'Htidal_s',
#        'Htot_s',
#        'lhs_s',
#        'liquidus_rho_s',
#        'liquidus_s',
#        'liquidus_temp_s',
#        'phi_s',
#        'rho_s',
#        'S_s',
#        'solidus_s',
#        'solidus_rho_s',
#        'solidus_temp_s',
#        'temp_s'
#    ]

#    index = data_l.index( instring )

#    return index

#====================================================================
def main( args ):

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

    # reproduce staple figures in Bower et al. (2018)
    # i.e., figs 3,4,5,6,7,8
    if 0:
        if len(args) < 2 :
            raise Exception('You must provide an argument consisting of comma-separated times.')
        bower_et_al_2018_fig3( args )
        bower_et_al_2018_fig4( args )
        bower_et_al_2018_fig5( args )

    # simplified model in Bower et al. (2018)
    # i.e., figs 1,2
    if 0:
        bower_et_al_2018_fig1( args )
        bower_et_al_2018_fig2( args )

    #figure6( args )
    #figure7( args )
    #figure8( args )
    plt.show()

#====================================================================

if __name__ == "__main__":

    main( sys.argv )

    sys.exit(1)
