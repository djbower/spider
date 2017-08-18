#!/usr/bin/env python

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import os
import sys 
PETSC_DIR = os.getenv('PETSC_DIR')
if not PETSC_DIR :
    print('You must define PETSC_DIR in your environment')
    sys.exit(1)
sys.path.append(os.path.join(PETSC_DIR,'bin'))
import PetscBinaryIO

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
class FigureData( object ):

    def __init__( self, args, nrows, ncols, width, height ):
        dd = {}
        self.data_d = dd
        dd['time_l'] = args[1]
        dd['time_units'] = 'Myr' # hard-coded
        dd['time_decimal_places'] = 2 # hard-coded
        self.process_time_list()
        self.set_properties( nrows, ncols, width, height )

    def petsc_bin_filename( self, time ):
        '''filename for petsc binary file output'''
        dd = self.data_d
        out = 'output'
        filename = 'petscbin.{0}'.format( time )
        tname = os.path.join( out, filename )
        if not os.path.isfile( tname ):
            print 'petsc_bin_filename: ERROR', tname, 'does not exist'
            print '    please specify times for which data exists in output/'
            sys.exit(1)
        return tname

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

    def get_data( self, field, time ):

        infile = self.petsc_bin_filename( time )
        io = PetscBinaryIO.PetscBinaryIO()
        objects_l = io.readBinaryFile( infile )
        index = get_vector_index( field )
        yy = objects_l[index]
        return yy

    def get_xy_data( self, field, time, fmt_o='' ):
        yy_b = self.get_data( field, time )

        # next index (4) should correlate to pressure_b
        xx_b = self.get_data( 'pressure_b', time ) * 1.0E-9 # to GPa
        xx_b = xx_b[1:-1] # basic internal nodes only
        # DJB HACKY
        # want to flip the sign of the entropy gradient to make it
        # cleaner to plot
        if field == 'dSdr':
            yy_b *= -1.0
        yy = yy_b[1:-1] # only basic internal nodes
        if fmt_o:
            yy = ascale( yy, fmt_o.const )
        return ( xx_b, yy )

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
        outname = 'fig{}.pdf'.format( num)
        self.fig.savefig(outname, transparent=True, bbox_inches='tight',
            pad_inches=0.05, dpi=dd['dpi'])

    def set_colors( self, num=6 ):
        dd = self.data_d
        cmap = plt.get_cmap('jet')
        colors_l = [cmap(i) for i in np.linspace(0, 1, num)]
        # color scheme from Tim.  Nice reds and blues
        #colors = ['#2364A4',
        #          '#1695F9',
        #          '#95D5FD',
        #          '#8B0000',
        #          '#CD5C5C',
        #          '#FA141B',
        #          '#FFA07A']
        dd['colors_l'] = colors_l

    def set_properties( self, nrows, ncols, width, height ):
        dd = self.data_d
        dd['nrows'] = nrows
        dd['ncols'] = ncols
        dd['width'] = width # inches
        dd['height'] = height # inches
        font_d = {'family' : 'serif',
                  #'style': 'normal',
                  #'weight' : 'bold',
                  'serif': ['computer modern roman'],
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
        if TITLE:
            title = r'Time ({0})'.format( units )
        else:
            title = ''
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
def figure1( args ):

    # article class text width is 4.7747 inches
    # http://tex.stackexchange.com/questions/39383/determine-text-width

    fig_o = FigureData( args, 2, 2, 4.7747, 4.7747 )

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]
    ax2 = fig_o.ax[1][0]
    ax3 = fig_o.ax[1][1]
 
    time = fig_o.time[0] # first timestep since liquidus and solidus
                          # are time-independent

    # shade grey between liquidus and solidus
    xx_liq, yy_liq = fig_o.get_xy_data( 'liquidus', time )
    xx_sol, yy_sol = fig_o.get_xy_data( 'solidus', time )
    ax0.fill_between( xx_liq, yy_liq, yy_sol, facecolor='grey', alpha=0.35, linewidth=0 )
    xx_liqt, yy_liqt = fig_o.get_xy_data( 'liquidus_temp', time )
    xx_solt, yy_solt = fig_o.get_xy_data( 'solidus_temp', time )
    ax1.fill_between( xx_liqt, yy_liqt, yy_solt, facecolor='grey', alpha=0.35, linewidth=0 )

    # dotted lines of constant melt fraction
    for xx in range( 0, 11, 2 ):
        yy_b = xx/10.0 * (yy_liq - yy_sol) + yy_sol
        if xx == 0:
            # solidus
            ax0.plot( xx_liq, yy_b, '-', linewidth=0.5, color='black' )
        elif xx == 3:
            # typically, the approximate location of the rheological transition
            ax0.plot( xx_liq, yy_b, '-', linewidth=1.0, color='white')
        elif xx == 10:
            # liquidus
            ax0.plot( xx_liq, yy_b, '-', linewidth=0.5, color='black' )
        else:
            # dashed constant melt fraction lines
            ax0.plot( xx_liq, yy_b, ':', linewidth=1.0, color='white' )

    handle_l = [] # handles for legend

    for nn, time in enumerate( fig_o.time ):
        color = fig_o.get_color( nn )
        # use melt fraction to determine mixed region
        xx, yy = fig_o.get_xy_data( 'phi', time )
        MIX = get_mix( yy )

        label = fig_o.get_legend_label( time )

        # entropy
        xx, yy = fig_o.get_xy_data( 'S', time )
        ax0.plot( xx, yy, '--', color=color )
        handle, = ax0.plot( xx*MIX, yy*MIX, '-', color=color, label=label )
        handle_l.append( handle )
        # temperature
        xx, yy = fig_o.get_xy_data( 'temp', time )
        ax1.plot( xx, yy, '--', color=color )
        ax1.plot( xx*MIX, yy*MIX, '-', color=color )
        # melt fraction
        xx, yy = fig_o.get_xy_data( 'phi', time )
        ax2.plot( xx, yy, '-', color=color )
        # viscosity
        visc_const = 1 # this is used for the arcsinh scaling
        visc_fmt = MyFuncFormatter( visc_const )
        xx, yy = fig_o.get_xy_data( 'visc', time, visc_fmt)
        ax3.plot( xx, yy, '-', color=color )

    xticks = [0,50,100,135]

    # titles and axes labels, legends, etc
    title = '(a) Entropy, J kg$^{-1}$ K$^{-1}$'
    yticks = [1600,2000,2400,2800,3200]
    fig_o.set_myaxes( ax0, title=title, ylabel='$S$', xticks=xticks, yticks=yticks )
    ax0.yaxis.set_label_coords(-0.075,0.59)
    title = '(b) Temperature, K'
    yticks= [1000,2000,3000,4000,5000]
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
    title = '(d) Viscosity, Pa$\cdot$s'
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
def figure2( args ):

    fig_o = FigureData( args, 2, 2, 4.7747, 4.7747 )

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
        #if (nn-1) % 2:
        #    continue

        color = fig_o.get_color( nn )
        # use melt fraction to determine mixed region
        xx, yy = fig_o.get_xy_data( 'phi', time )
        MIX = get_mix( yy )

        label = fig_o.get_legend_label( time )

        # Jconv_b
        xx, yy = fig_o.get_xy_data( 'Jconv', time, flux_fmt )
        ax0.plot( xx, yy, '--', color=color )
        handle, = ax0.plot( xx*MIX, yy*MIX, '-', label=label, 
            color=color )
        handle_l.append( handle )
        # Jmix_b
        xx, yy = fig_o.get_xy_data( 'Jmix', time, flux_fmt )
        ax2.plot( xx, yy, '--', color=color )
        ax2.plot( xx*MIX, yy*MIX, '-', color=color )
        # Jgrav_b
        xx, yy = fig_o.get_xy_data( 'Jgrav', time, flux_fmt )
        ax1.plot( xx, yy, '--', color=color )
        ax1.plot( xx*MIX, yy*MIX, '-', color=color )
        # Jtot_b
        xx, yy = fig_o.get_xy_data( 'Jtot', time, flux_fmt )
        ax3.plot( xx, yy, '--', color=color )
        ax3.plot( xx*MIX, yy*MIX, '-', color=color )

    # titles and axes labels, legends, etc.
    xticks = [0,50,100,135]
    yticks = [-1E15,-1E12, -1E6, -1E0, 0, 1E0, 1E6, 1E12, 1E15]
    title = '(a) Convective flux, W m$^{-2}$'
    fig_o.set_myaxes( ax0, title=title, ylabel='$F_\mathrm{conv}$',
        yticks=yticks, xticks=xticks, fmt=flux_fmt )
    ax0.yaxis.set_label_coords(-0.16,0.54)
    fig_o.set_mylegend( ax0, handle_l, ncol=2 )
    title = '(c) Mixing flux, W m$^{-2}$'
    fig_o.set_myaxes( ax2, title=title, xlabel='$P$ (GPa)',
        ylabel='$F_\mathrm{mix}$', yticks=yticks, xticks=xticks, fmt=flux_fmt )
    ax2.yaxis.set_label_coords(-0.16,0.54)
    title = '(b) Separation flux, W m$^{-2}$'
    fig_o.set_myaxes( ax1, title=title,
        ylabel='$F_\mathrm{grav}$', yticks=yticks, xticks=xticks, fmt=flux_fmt )
    ax1.yaxis.set_label_coords(-0.16,0.54)
    title = '(d) Total flux, W m$^{-2}$'
    fig_o.set_myaxes( ax3, title=title, xlabel='$P$ (GPa)',
        ylabel='$F_\mathrm{tot}$', yticks=yticks, xticks=xticks, fmt=flux_fmt )
    ax3.yaxis.set_label_coords(-0.16,0.54)

    fig_o.savefig(2)

#====================================================================
def figure3( args ):

    # keep y the same by scaling (7.1621)
    fig_o = FigureData( args, 3, 2, 4.7747, 7.1621 )

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
        #if (nn-1) % 2:
        #    continue
        color = fig_o.get_color( nn )
        # use melt fraction to determine mixed region
        xx, yy = fig_o.get_xy_data( 'phi', time )
        MIX = get_mix( yy )

        label = fig_o.get_legend_label( time )

        # eddy diffusivity
        xx, yy = fig_o.get_xy_data( 'kappah', time, eddy_fmt )
        ax0.plot( xx, yy, '--', color=color )
        handle, = ax0.plot( xx*MIX, yy*MIX, '-', label=label, 
            color=color )
        handle_l.append( handle )
        # density
        xx, yy = fig_o.get_xy_data( 'rho', time )
        ax1.plot( xx, yy, '--', color=color )
        ax1.plot( xx*MIX, yy*MIX, '-', color=color )
        # thermal expansion coefficient
        xx, yy = fig_o.get_xy_data( 'alpha', time, alpha_fmt )
        ax3.plot( xx, yy, '--', color=color )
        ax3.plot( xx*MIX, yy*MIX, '-', color=color )
        # heat capacity
        xx, yy = fig_o.get_xy_data( 'cp', time )
        ax4.plot( xx, yy, '--', color=color )
        ax4.plot( xx*MIX, yy*MIX, '-', color=color )
        # dTdrs
        xx, yy = fig_o.get_xy_data( 'dTdrs', time )
        ax5.plot( xx, yy, '--', color=color )
        ax5.plot( xx*MIX, yy*MIX, '-', color=color )
        # entropy gradient
        xx, yy = fig_o.get_xy_data( 'dSdr', time, dSdr_fmt )
        ax2.plot( xx, yy, '--', color=color )
        ax2.plot( xx*MIX, yy*MIX, '-', color=color )

    xticks = [0,50,100,135]

    # titles and axes labels, legends, etc.
    title = '(a) Eddy diffusivity, W m$^{-1}$ K$^{-1}$'
    yticks = [1.0, 1.0E3, 1.0E6, 1.0E9]
    fig_o.set_myaxes( ax0, title=title, ylabel='$\kappa_h$',
        yticks=yticks, fmt=eddy_fmt, xticks=xticks )
    ax0.yaxis.set_label_coords(-0.1,0.475)

    title = '(b) Density, kg m$^{-3}$'
    yticks = [2000,3000,4000,5000,6000]
    fig_o.set_myaxes( ax1, title=title, ylabel='$\\rho$',
        yticks=yticks, xticks=xticks )
    ax1.yaxis.set_label_coords(-0.1,0.6)
    fig_o.set_mylegend( ax1, handle_l, loc=4, ncol=2 )

    title = '(d) Thermal expansion, K$^{-1}$'
    yticks = [1.0E-5, 1.0E-4, 1.0E-3, 1.0E-2]
    fig_o.set_myaxes( ax3, title=title, ylabel='$\\alpha$',
        yticks=yticks, fmt=alpha_fmt, xticks=xticks )
    ax3.yaxis.set_label_coords(-0.1,0.475)

    title = '(e) Heat capacity, J kg$^{-1}$ K$^{-1}$'
    yticks = [0,5000,10000,15000]
    fig_o.set_myaxes( ax4, title=title, xlabel='$P$ (GPa)', ylabel='$c$',
        yticks=yticks, xticks=xticks )
    ax4.yaxis.set_label_coords(-0.1,0.475)

    title = '(f) Adiabatic grad, K m$^{-1}$'
    yticks = [-3E-3, -2E-3, -1E-3, 0]
    fig_o.set_myaxes( ax5, title=title, xlabel='$P$ (GPa)',
        yticks=yticks, xticks=xticks,
        ylabel='$\\left(\\frac{dT}{dr}\\right)_S$')
    ax5.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
    ax5.yaxis.set_label_coords(-0.15,0.43)

    title = '(c) Entropy grad, K m$^{-1}$'
    #yticks= [-1.0E-3, -1.0E-6, -1.0E-9, -1.0E-12, -1.0E-15]
    yticks = [1E-15, 1E-12, 1E-9, 1E-6, 1E-3]
    fig_o.set_myaxes( ax2, title=title,
        yticks=yticks, xticks=xticks,
        ylabel='$-\\frac{dS}{dr}$', fmt=dSdr_fmt)
    ax2.yaxis.set_label_coords(-0.1,0.565)

    fig_o.savefig(3)

#====================================================================
def figure4( args ):

    width = 4.7747 # * 0.5
    height = 4.7747 # * 0.5
    fig_o = FigureData( args, 2, 2, width, height )

    dd = fig_o.data_d

    # FIXME: below will now break because I am not loading in the 
    # python model
    dep2pres = fig_o.data_d['solver_o'].sol_o.mesh_o._calc_pressure
    mesh_o = dd['solver_o'].sol_o.mesh_o
    # below, only basic internal nodes
    xconst = 55692628932.82327
    xconst *= 1.0E-9 # to GPa

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]
    ax2 = fig_o.ax[1][0]
    ax3 = fig_o.ax[1][1]

    ###############
    # plot liquidus
    Sliq_file = '/Users/dan/Documents/research/my_papers/in_prep/magma_ocean_method/toymodel/Sliq_processed.dat'
    rad, Sliq = np.loadtxt( Sliq_file, unpack=True )
    depth = 1.0 - rad
    # FIXME
    pres = dep2pres( depth )
    pres *= xconst
    const = 2993.025100070677
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
    dSdr_file = '/Users/dan/Documents/research/my_papers/in_prep/magma_ocean_method/toymodel/flux1E6dim_processed.dat'
    dSliqdr_file = '/Users/dan/Documents/research/my_papers/in_prep/magma_ocean_method/toymodel/dSliqdr_processed.dat'

    rad, dSdr = np.loadtxt( dSdr_file, unpack=True )
    rad2, dSliqdr = np.loadtxt( dSliqdr_file, unpack=True )
    depth = 1.0 - rad
    depth2 = 1.0 - rad2
    # FIXME
    pres = dep2pres( depth )
    pres2 = dep2pres( depth2 )

    pres *= xconst
    pres2 *= xconst
    const = 0.0004697889028520918
    dSdr *= const * -1.0 # DJB HACK
    dSliqdr *= const

    # must scale by 1/2
    dSliqdr *= 0.5

    dSdr_const = 1.0E15
    dSdr_fmt = MyFuncFormatter( dSdr_const )

    dSdr = ascale( dSdr, dSdr_fmt.const )
    handle2, = ax1.plot( pres, dSdr, 'k-' )

    title = '(b) Entropy grad, K m$^{-1}$'
    #yticks= [-1.0E-3, -1.0E-6, -1.0E-9, -1.0E-12, -1.0E-15]#, 1.0E-3]
    yticks = [1.0E-15, 1.0E-12, 1.0E-9,1.0E-6,1.0E-3]
    fig_o.set_myaxes( ax1, yticks=yticks, xticks=xticks,
        ylabel='$-\\frac{dS}{dr}$', fmt=dSdr_fmt, title=title)
    ax1.yaxis.set_label_coords(-0.11,0.565)

    ############
    # plot Fconv
    Fconv_file = '/Users/dan/Documents/research/my_papers/in_prep/magma_ocean_method/toymodel/FluxConv_processed.dat'
    rad, Fconv = np.loadtxt( Fconv_file, unpack=True )
    depth = 1.0 - rad
    # FIXME
    pres = dep2pres( depth )
    pres *= xconst

    const = 193508342723647.66
    Fconv *= const

    # for arcsinh scaling
    #flux_const = 1.0E6
    flux_const = 1.0E6
    flux_fmt = MyFuncFormatter( flux_const )
    Fconv = ascale( Fconv, flux_fmt.const )
    handle3, = ax0.plot( pres, Fconv, 'k-' )

    # titles and axes labels, legends, etc.
    #yticks = [-1.0E12, -1.0E6, -1.0E0, -1.0E-3, 1.0E-3, 1.0E0, 1.0E6, 1E12]
    yticks = [-1.0E15, -1.0E12, -1.0E6, -1.0E0,0, 1.0E0, 1.0E6, 1.0E12, 1.0E15]

    title = '(a) Convective flux, W m$^{-2}$'
    fig_o.set_myaxes( ax0, title=title, ylabel='$F_\mathrm{conv}$',
        yticks=yticks, fmt=flux_fmt, xticks=xticks )
    ax0.yaxis.set_label_coords(-0.16,0.54)

    ###########
    # plot Fmix
    Fmix_file = '/Users/dan/Documents/research/my_papers/in_prep/magma_ocean_method/toymodel/FluxMix_processed.dat'
    rad, Fmix = np.loadtxt( Fmix_file, unpack=True )
    depth = 1.0 - rad
    # FIXME
    pres = dep2pres( depth )
    pres *= xconst

    const = 193508342723647.66
    Fmix *= const

    Fmix = ascale( Fmix, flux_fmt.const )
    handle4, = ax2.plot( pres, Fmix, 'k-' )
    #handle4, = ax2.plot( pres, Fmix, 'k-' )

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
def figure5( args ):

    width = 4.7747 * 0.5
    height = 4.7747 * 0.5
    fig_o = FigureData( args, 1, 1, width, height )

    ax0 = fig_o.ax

    const = 0.0004697889028520918

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
    ylabel = '$-\\frac{dS}{dr}$'
    fig_o.set_myaxes( ax0, title=title, xlabel=xlabel, ylabel=ylabel,
        yticks=yticks, fmt=dSdr_fmt, xticks=xticks, xfmt=dSliqdr_fmt )
    ax0.yaxis.set_label_coords(-0.1,0.565)
    #ax0.set_xlim( -5.2, 0.4)

    prefix = '/Users/dan/Documents/research/my_papers/in_prep/magma_ocean_method/toymodel/regime_diagram/'

    # plot asymptote
    FLAG = 1
    if FLAG:
        # not sure why I need a negative number less than 5.2?
        # probably to account for the x shift?
        xx2 = np.linspace( -20.4, 0.8, 1000000 )
        xx2 *= const
        yy2 = 0.5 * xx2 * -1.0 # DJB HACK TO FLIP Y
        xx2 = ascale( xx2, dSliqdr_fmt.const )
        yy2 = ascale( yy2, dSdr_fmt.const )
        ax0.plot( xx2, yy2, 'k--' )
        ax0.fill_between( xx2, yy2, 0, facecolor='0.75' )

    # plot flux minimum
    if 1:
        #xx2 = np.linspace( -20.4, 0.8, 1000000 )
        xx2 = np.linspace( -20.4, -0.0001, 1000000 )
        xx2 *= const
        yy2 = 1.0/6.0 * xx2 * -1.0 # DJB HACK TO FLIP Y
        xx2 = ascale( xx2, dSliqdr_fmt.const )
        yy2 = ascale( yy2, dSdr_fmt.const )
        xx2 -= 0.4
        yy2 -= 0.4
        ax0.plot( xx2, yy2, color='0.25', linestyle=':' )

    # positive heat fluxes
    for nn in [6,9,12]:
        filein = prefix + 'dSdr_10p%(nn)d_processed.dat' % vars()
        xx, yy = np.loadtxt( filein, unpack=True )
        xx *= const
        xx = ascale( xx, dSliqdr_fmt.const )
        yy *= const * -1.0 # DJB HACK TO FLIP Y
        yy = ascale( yy, dSdr_fmt.const )
        # offset y for visual clarity
        xx += 0.4
        yy += 0.4
        ax0.plot( xx, yy, '-' ) #, 'k-' )

    # negative heat fluxes
    for nn in [6,9,12]:
        filein = prefix + 'dSdr_n10p%(nn)d_processed.dat' % vars()
        xx, yy = np.loadtxt( filein, unpack=True )
        xx *= const
        xx = ascale( xx, dSliqdr_fmt.const )
        yy *= const * -1.0 # DJB HACK TO FLIP Y
        yy = ascale( yy, dSdr_fmt.const )
        # offset y for visual clarity
        xx -= 0.4
        yy -= 0.4
        ax0.plot( xx, yy, '-' )

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
    M0 = data[:,1]
    Mliq = data[:,2]
    Msol = data[:,3]
    dMliqdt = data[:,4]
    tau = data[:,5]
    emissivity = data[:,6]
    x0init = data[:,7]
    x0 = data[:,8]
    dx0dt = data[:,9]
    p0 = data[:,10]
    dp0dx = data[:,11]
    m0 = data[:,12]
    tau0 = data[:,13]
    x1init = data[:,14]
    x1 = data[:,15]
    dx1dt = data[:,16]
    p1 = data[:,17]
    dp1dx = data[:,18]
    m1 = data[:,19]
    tau1 = data[:,20]

    # THIS PROVES THAT THE INITIAL CONDITION IS APPLIED
    # CORRECTLY BOTH CO2 AND H2O
    #HENRY = 4.4E-12
    #HENRY_POW = 1.0
    #alpha = 4.0*np.pi*6371000.0**2.0
    #alpha /= 10.0 * M0
    #alpha *= 100.0**(1.0-HENRY_POW)
    #alpha /= HENRY**HENRY_POW
    #result = x0 + alpha*x0**HENRY_POW - x0init
    #print result
    # GOOD TO HERE
    #sys.exit(1)

    # based on above, what should partial pressure of the atmosphere be?
    #p = ((x0/100.0)/HENRY)**HENRY_POW
    #print p, p0
    # GOOD TO HERE
    #sys.exit(1)

    # what about mass of the atmosphere?
    #Matm = 4.0*np.pi*6371000.0**2.0 / 10.0 * p
    #print Matm, m0
    # GOOD TO HERE
    #sys.exit(1)

    # what about mass of the volatiles in the liquid?
    #Mvolliq = (x0/100.0) * M0
    #print Mvolliq
    #sys.exit(1)

    #Minit = (x0init/100.0) * M0
    #print Minit
    #print Mvolliq + Matm
    #sys.exit(1)


    #Minit = (x0init/100.0) * M0
    #Mliq = (x0/100.0) * M0
    #Matm = m0
    #Matm2 = 4.0*np.pi*6371000.0**2.0 * p0 / 10.0
    #print Minit, Mvolliq, Matm, Matm2
    #sys.exit(1)

    # mass fractions
    CO2init = x0init * M0
    CO2liq = x0 * Mliq
    CO2liq /= CO2init # normalised
    CO2atm = m0 * 100.0 / CO2init # normalised
    # FIXME: kdist is hard-coded here
    CO2sol = 5.0E-4 * x0 * Msol
    CO2sol /= CO2init # normalised
    CO2tot = CO2liq + CO2atm + CO2sol # normalised
    #CO2atm2 = 4.0*np.pi*6371000.0**2.0 / 10.0 * p0
    #CO2atm2 /= CO2init
    #print CO2atm2
    #print np.min(CO2tot), np.max(CO2tot)
    #print 'CO2init=', CO2init
    #print 'CO2sol=', CO2sol
    #print 'CO2liq=', CO2liq
    #print 'CO2atm=', CO2atm
    #print 'CO2tot=', CO2tot

    # more testing
    #Matm = 4.0*np.pi*6371000.0**2 / 10.0 * p0
    #print 'Matm=', Matm
    #print 'Mtot2=', Mliq2 + Matm
    #print 'Mtotinit=', (x0init/100.0) * M0

    # mass fractions
    H2Oinit = x1init * M0
    H2Oliq = x1 * Mliq
    H2Oliq /= H2Oinit # normalised
    H2Oatm = m1 * 100.0 / H2Oinit # normalised
    # FIXME: kdist is hard-coded here
    H2Osol = 1.0E-4 * x1 * Msol
    H2Osol /= H2Oinit # normalised
    H2Otot = H2Oliq + H2Oatm + H2Osol # normalised
    #print np.min(H2Otot), np.max(H2Otot)
    #sys.exit(1)

    print np.min(x0), np.max(x0)
    print np.min(x1), np.max(x1)

    xticks = [1E-2,1E-1,1E0,1E1,1E2,1E3,4.55E3]
    xlabel = 'Time (Myr)'

    # figure a
    title = '(a) CO$_2$ reservoirs'
    h4, = ax0.semilogx( time, Msol/M0, 'k--', label='Mantle' )
    h1, = ax0.semilogx( time, CO2liq, 'r-', label='Liquid' )
    h2, = ax0.semilogx( time, CO2sol, 'b-', label='Solid' )
    h3, = ax0.semilogx( time, CO2atm, 'g-', label='Atmos' )
    #ax0.semilogx( time, CO2tot, 'k-' )
    fig_o.set_myaxes( ax0, title=title, ylabel='$x$', xticks=xticks )
    handle_l = [h1,h3,h2,h4]
    fig_o.set_mylegend( ax0, handle_l, loc='center right', ncol=1, TITLE=0 )
    ax0.yaxis.set_label_coords(-0.1,0.46)

    # solid mass fraction is now plotted on figures a and b along with
    # the mass fraction of the volatile in each reservoir (liq, sol, atm)
    # figure b
    #title = '(b) Mass fraction'
    #h1, = ax1.semilogx( time, Mliq/M0, label='Liquid' )
    #h2, = ax1.semilogx( time, Msol/M0, label='Solid' )
    #h3, = ax1.semilogx( time, (Mliq+Msol)/M0 )
    #fig_o.set_myaxes( ax1, title=title, ylabel='Frac', xticks=xticks )
    #handle_l = [h1,h2]#,h3]
    #fig_o.set_mylegend( ax1, handle_l, loc='center right', ncol=1, TITLE=0 )
    #ax1.yaxis.set_label_coords(-0.1,0.46)

    # figure c
    title = '(c) H$_2$O reservoirs'
    h4, = ax2.semilogx( time, Msol/M0, 'k--', label='Mantle' )
    h1, = ax2.semilogx( time, H2Oliq, 'r-', label='Liquid' )
    h2, = ax2.semilogx( time, H2Osol, 'b-', label='Solid' )
    h3, = ax2.semilogx( time, H2Oatm, 'g-', label='Atmos' )
    #ax1.semilogx( time, H2Otot, 'k-' )
    fig_o.set_myaxes( ax2, title=title, xlabel=xlabel, ylabel='$x$', xticks=xticks )
    handle_l = [h1,h3,h2,h4]
    fig_o.set_mylegend( ax2, handle_l, loc='center right', ncol=1, TITLE=0 )
    ax2.yaxis.set_label_coords(-0.1,0.46)

    # figure e
    title = '(b) Partial pressure (bar)'
    ylabel = '$p$'
    h1, = ax1.semilogx( time, p0*1.0E-5, 'r-', label='CO$_2$')
    h2, = ax1.semilogx( time, p1*1.0E-5, 'b-', label='H$_2$O')
    handle_l = [h1,h2]
    fig_o.set_myaxes( ax1, title=title, ylabel=ylabel, xticks=xticks )
    fig_o.set_mylegend( ax1, handle_l, loc='center right', ncol=1, TITLE=0 )
    ax1.yaxis.set_label_coords(-0.1,0.575)

    # figure f
    title = '(d) Emissivity'
    ylabel = '$\epsilon$'
    ax3.semilogx( time, emissivity )
    fig_o.set_myaxes( ax3, title=title, xlabel=xlabel, ylabel=ylabel, xticks=xticks )
    ax3.yaxis.set_label_coords(-0.1,0.575)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    fig_o.savefig(6)

#====================================================================
def atmosphere_data_to_array( fig_o ):

    # locate times to process based on files located in output/
    odir = 'output'
    file_l = [f for f in os.listdir(odir) if os.path.isfile(os.path.join(odir,f))]
    if not file_l:
        print 'ERROR: output directory contains no files'
        sys.exit(1)

    time_l = [fname.split('.')[-1] for fname in file_l]

    # remove '.info' suffix
    time_l = list(filter(lambda a: a !='info', time_l))
    # remove '.m' suffix
    time_l = list(filter(lambda a: a !='m', time_l))
    # convert to int
    time_l = [int(val) for val in time_l]
    # ascending order
    time_l = sorted( time_l, key=int)

    data = fig_o.get_data( 'atmosphere_data', time_l[0] )

    # must be some data to process
    if data is None:
        print 'atmosphere_data_to_array: no data'
        sys.exit(0)

    ncols = np.size( data ) + 1
    nrows = len( time_l )

    data_a = np.zeros( (nrows, ncols) )

    for nn, time in enumerate( time_l ):

        data = fig_o.get_data( 'atmosphere_data', time )
        data_a[nn,0] = time
        data_a[nn,1:] = data

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

    # Bit hacky, but flag allows us to plot 0 on the y axis, which is
    # useful to emphasize that we are plotting positive and negative
    # values, e.g. for the heat fluxes
    FLAG = 0

    if not exponent:
        try:
            exponent = int(np.floor(np.log10(abs(num))))
        except OverflowError:
            exponent = 0
            FLAG=1
    coeff = round(num / float(10**exponent), decimal_digits)
    # sometimes, probably due to floating point precision? the coeff
    # is not less than ten.  Correct for that here
    if np.abs(coeff) >= 10.0:
        coeff /= 10.0
        exponent += 1
    if not precision:
        precision = decimal_digits

    if coeff < 0.0 and not FLAG:
        fmt = r"$-10^{{{0}}}$".format(exponent)
        #fmt= r"${{{0}}}$".format(exponent)
    elif coeff > 0.0 and not FLAG:
        fmt = r"$10^{{{0}}}$".format(exponent)
    else:
        fmt = r"$0$"

    return fmt
    #return r"${0:.{2}f}\cdot10^{{{1:d}}}$".format(coeff, exponent, precision)

#====================================================================
def get_vector_index( instring ):
    """
    the order of fields in data_l must match exactly the order in
    which PetscVecs are written to the output Petsc binary.
    See monitor.c
    """

    data_l = [
        # extraneous
        'time_data',
        'atmosphere_data',
        'dSdr_b_aug',
        # stored in M->meshVecs_b
        'area_b',
        'dPdr_b',
        'pressure_b',
        'radius_b',
        'mix_b',
        # stored in M->meshVecs_s
        'pressure_s',
        'radius_s',
        'volume_s',
        'dPdr_s',
        'area_s',
        'rho_s',
        'mass_s',
        # stored in S->solutionVecs_b
        'alpha',
        'alpha_mix',
        'cond',
        'cp',
        'cp_mix',
        'dfusdr',
        'dfusdr_temp',
        'dSdr',
        'dSliqdr',
        'dSsoldr',
        'dTdrs',
        'dTdrs_mix',
        'Etot',
        'fusion',
        'fusion_curve',
        'fusion_curve_temp',
        'fusion_rho',
        'fusion_temp',
        'fwtl',
        'fwts',
        'gphi',
        'gsuper',
        'Jcond',
        'Jconv',
        'Jgrav',
        'Jmix',
        'Jtot',
        'kappah',
        'liquidus',
        'liquidus_rho',
        'liquidus_temp',
        'nu',
        'phi',
        'rho',
        'S',
        'solidus',
        'solidus_rho',
        'solidus_temp',
        'temp',
        'visc',
         # stored in S->solutionVecs_s
        'cp_s',
        'cp_mix_s',
        'dSdt_s',
        'fusion_s',
        'fusion_curve_s',
        'fusion_curve_temp_s',
        'fusion_temp_s',
        'fwtl_s',
        'fwts_s',
        'gphi_s',
        'Hradio_s',
        'Htidal_s',
        'Htot_s',
        'lhs_s',
        'liquidus_rho_s',
        'liquidus_s',
        'liquidus_temp_s',
        'phi_s',
        'rho_s',
        'S_s',
        'solidus_s',
        'solidus_rho_s',
        'solidus_temp_s',
        'temp_s'
    ]

    index = data_l.index( instring )

    return index

#====================================================================
def main( args ):

    #figure1( args )
    #figure2( args )
    #figure3( args )
    #figure4( args )
    #figure5( args )
    figure6( args )
    plt.show()

#====================================================================

if __name__ == "__main__":

    main( sys.argv )

    sys.exit(1)