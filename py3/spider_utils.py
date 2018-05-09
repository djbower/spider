#!/usr/bin/env python

import json
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.transforms as transforms
import numpy as np
import logging
import os

#===================================================================
# CLASSES
#===================================================================

#===================================================================
class MyFuncFormatter( object ):

    '''the default function formatter from
       matplotlib.ticker.FuncFormatter(func) only accepts two
       arguments, which is not enough to scale an arcsinh function.
       But by creating our own class here we can attach the scaling
       to the object which can then be accessed in __call__'''

    def __init__( self, arcsinh_scale ):
        self.const = arcsinh_scale

    def ascale( self, xx ):
        '''map input to log-like values (scaled arcsinh)'''
        yy = np.arcsinh( xx*self.const )
        return yy

    def _invascale( self, yy ):
        '''map input from log-like values (inverse transform)'''
        xx = np.sinh(yy) / self.const
        return xx

    def _sci_notation( self, num, decimal_digits=1, precision=None, exponent=None):
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

    def __call__( self, x, pos ):
        y = self._invascale( x )
        fmt = self._sci_notation( y, 0 )
        return fmt

#===================================================================
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
        scaling = float(fdata_d['scaling'])
        if len( fdata_d['values'] ) == 1:
            values_a = float( fdata_d['values'][0] )
        else:
            values_a = np.array( [float(value) for value in fdata_d['values']] )
        scaled_values_a = scaling * values_a
        if fmt_o:
            scaled_values_a = fmt_o.ascale( scaled_values_a )
        return scaled_values_a

    def get_scaled_field_values_internal( self, field, fmt_o='' ):
        '''get the scaled values for the internal nodes (ignore top
           and bottom nodes)'''
        scaled_values_a = self.get_scaled_field_values( field, fmt_o )
        return scaled_values_a[1:-1]

    def get_mixed_phase_boolean_array( self, nodes='basic' ):
        '''this array enables us to plot different linestyles for 
           mixed phase versus single phase quantities'''
        if nodes == 'basic':
            phi = self.get_scaled_field_values( 'phi_b' )
        elif nodes == 'basic_internal':
            phi = self.get_scaled_field_values_internal( 'phi_b' )
        elif nodes == 'staggered':
            phi = self.get_scaled_field_values( 'phi_s' )
        # define mixed phase by these threshold values
        MIX = (phi<0.999) & (phi>0.001)
        MIX = MIX * 1.0 # convert to float array
        # set single phase region to nan to prevent plotting
        MIX[MIX==0] = np.nan

        return MIX

#===================================================================
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
        ylabel='', yticks='', fmt='', xfmt='', xmax='' ):
        if title:
            self.set_mytitle( ax, title )
        if xlabel:
            self.set_myxlabel( ax, xlabel )
        if xticks:
            self.set_myxticks( ax, xticks, xmax, xfmt )
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

    def set_myxticks( self, ax, xticks, xmax, fmt ):
        dd = self.data_d
        if fmt:
            xticks = fmt.ascale( np.array(xticks) )
            ax.xaxis.set_major_formatter(
                mpl.ticker.FuncFormatter(fmt))
        ax.set_xticks( xticks)
        # set x limits to match extent of ticks
        if not xmax: xmax=xticks[-1]
        ax.set_xlim( xticks[0], xmax )

    def set_myyticks( self, ax, yticks, fmt ):
        dd = self.data_d
        if fmt:
            yticks = fmt.ascale( np.array(yticks) )
            ax.yaxis.set_major_formatter(
                mpl.ticker.FuncFormatter(fmt))
        ax.set_yticks( yticks)
        # set y limits to match extent of ticks
        ax.set_ylim( yticks[0], yticks[-1] )

#====================================================================

#====================================================================
# FUNCTIONS
#====================================================================
def get_all_output_times( odir='output' ):

    '''get all times (in Myrs) from the json files located in the
       output directory'''

    # locate times to process based on files located in odir/
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
def ppm_to_mass_fraction( ppm ):

    return ppm * 1.0E-6

#====================================================================
