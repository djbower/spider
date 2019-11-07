#!/usr/bin/env python

import logging
import spider_utils as su
import matplotlib.transforms as transforms
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.ticker as ticker

logger = su.get_my_logger(__name__)

#====================================================================
def plot_atmosphere():

    logger.info( 'building atmosphere' )
    wdir = 'output'

    width = 10 * 3.0/2.0
    height = 10 / 2.0
    fig_o = su.FigureData( 1, 3, width, height, wdir, units='Myr' )
    fig_o.fig.subplots_adjust(wspace=0.4,hspace=0.5)

    ax0 = fig_o.ax[0]
    ax1 = fig_o.ax[1]
    ax2 = fig_o.ax[2]
    #ax3 = fig_o.ax[1][1]

    fig_o.time = su.get_all_output_times(odir=wdir)
    timeMyr_a = fig_o.time * 1.0E-6 # Myrs

    keys_t = ( ('atmosphere','mass_liquid'),
               ('atmosphere','mass_solid'),
               ('atmosphere','mass_mantle'),
           #    ('atmosphere','CO2','liquid_kg'),
           #    ('atmosphere','CO2','solid_kg'),
           #    ('atmosphere','CO2','initial_kg'),
          #     ('atmosphere','CO2','atmosphere_kg'),
          #     ('atmosphere','CO2','atmosphere_bar'),
               ('atmosphere','H2O','liquid_kg'),
               ('atmosphere','H2O','solid_kg'),
               ('atmosphere','H2O','initial_kg'),
               ('atmosphere','H2O','atmosphere_kg'),
               ('atmosphere','H2O','atmosphere_bar'),
               ('atmosphere','H2','liquid_kg'),
               ('atmosphere','H2','solid_kg'),
               ('atmosphere','H2','initial_kg'),
               ('atmosphere','H2','atmosphere_kg'),
               ('atmosphere','H2','atmosphere_bar'),
               #('atmosphere','O2','liquid_kg'),
               #('atmosphere','O2','solid_kg'),
               #'atmosphere','O2','initial_kg'),
               #('atmosphere','O2','atmosphere_kg'),
               #('atmosphere','O2','atmosphere_bar'),
               ('atmosphere','temperature_surface'),
               ('atmosphere','emissivity'),
               ('rheological_front_phi','phi_global'),
               ('atmosphere','Fatm'),
               ('atmosphere','H2O','reaction_kg'),
               ('atmosphere','H2','reaction_kg'))

    data_a = su.get_dict_surface_values_for_times( keys_t, fig_o.time, indir=wdir )

    out_a = np.column_stack( (timeMyr_a,data_a[15,:],data_a[12,:]) )
    np.savetxt( 'out.dat', out_a )

    mass_liquid_a = data_a[0,:]
    mass_solid_a = data_a[1,:]
    mass_mantle_a = data_a[2,:]
    mass_mantle = mass_mantle_a[0] # time independent

    # compute total mass (kg) in each reservoir
    #CO2_liquid_kg_a = data_a[3,:]
    #CO2_solid_kg_a = data_a[4,:]
    #CO2_total_kg_a = data_a[5,:]
    #CO2_total_kg = CO2_total_kg_a[0] # time-independent
    #CO2_atmos_kg_a = data_a[6,:]
    #CO2_atmos_a = data_a[7,:]
    #CO2_escape_kg_a = CO2_total_kg - CO2_liquid_kg_a - CO2_solid_kg_a - CO2_atmos_kg_a

    H2O_liquid_kg_a = data_a[3,:]
    H2O_solid_kg_a = data_a[4,:]
    H2O_total_kg_a = data_a[5,:]
    H2O_total_kg = H2O_total_kg_a[0] # time-independent
    H2O_atmos_kg_a = data_a[6,:]
    H2O_atmos_a = data_a[7,:]
    H2O_reaction_kg_a = data_a[17,:]
    H2O_total2_kg_a = H2O_liquid_kg_a + H2O_solid_kg_a + H2O_atmos_kg_a# + H2O_reaction_kg_a
    H2O_escape_kg_a = H2O_total_kg - H2O_liquid_kg_a - H2O_solid_kg_a - H2O_atmos_kg_a - H2O_reaction_kg_a

    H2_liquid_kg_a = data_a[8,:]
    H2_solid_kg_a = data_a[9,:]
    H2_total_kg_a = data_a[10,:]
    H2_total_kg = H2_total_kg_a[0] # time-independent
    H2_atmos_kg_a = data_a[11,:]
    H2_atmos_a = data_a[12,:]
    H2_reaction_kg_a = data_a[18,:]
    H2_total2_kg_a =  H2_liquid_kg_a + H2_solid_kg_a  + H2_atmos_kg_a# + H2_reaction_kg_a
    H2_escape_kg_a = H2_total_kg - H2_liquid_kg_a - H2_solid_kg_a - H2_atmos_kg_a - H2_reaction_kg_a

    #O2_liquid_kg_a = data_a[18,:]
    #O2_solid_kg_a = data_a[19,:]
    #O2_total_kg_a = data_a[10,:]
    #O2_total_kg = O2_total_kg_a[0] # time-independent
    #O2_atmos_kg_a = data_a[21,:]
    #O2_atmos_a = data_a[22,:]
    #O2_escape_kg_a = O2_total_kg - O2_liquid_kg_a - O2_solid_kg_a - O2_atmos_kg_a

    temperature_surface_a = data_a[13,:]
    emissivity_a = data_a[14,:]
    phi_global = data_a[15,:]
    Fatm = data_a[16,:]

    H2O_err = 100.0 * (H2O_total2_kg_a-H2O_total_kg) / H2O_total_kg
    H2_err =  100.0 * (H2_total2_kg_a-H2_total_kg) / H2_total_kg

    print('H2O_err MAX=', np.max(H2O_err))
    print('H2_err MAX=', np.max(H2_err))

    #xticks = [1E-5,1E-4,1E-3,1E-2,1E-1]#,1]
    #xticks = [1.0E-2, 1.0E-1, 1.0E0, 1.0E1, 1.0E2,1.0E3] #[1E-6,1E-4,1E-2,1E0,1E2,1E4,1E6]#,1]
    #xticks = (1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1E0)
    xlabel = 'Time (Myr)'
    #xlim = (1.0E-2, 4550)
    xlim = (1e-6,1)

    red = (0.5,0.1,0.1)
    blue = (0.1,0.1,0.5)
    yellow = (1,1,0)
    green = (0.1,0.5,0.1)
    black = 'black'

    # to plot melt fraction contours on figure (a)
    # compute time at which desired melt fraction is reached
    #phi_a = mass_liquid_a / mass_mantle
    #phi_cont_l = [0.75, 0.5,0.25,0.1,0.01]
    #phi_time_l = [] # contains the times for each contour
    #for cont in phi_cont_l:
    #    time_temp_l = su.find_xx_for_yy( timeMyr_a, phi_global, cont )
    #    index = su.get_first_non_zero_index( time_temp_l )
    #    if index is None:
    #        out = 0.0
    #    else:
    #        out = timeMyr_a[index]
    #    phi_time_l.append( out )


    ##########
    if 1:
        title = r'\textbf{(a) Pressure and global melt fraction}'
        ylabel = '$p$ (bar)'
        trans = transforms.blended_transform_factory(
            ax0.transData, ax0.transAxes)
        #for cc, cont in enumerate(phi_cont_l):
        #    ax0.axvline( phi_time_l[cc], ymin=0.05, ymax=0.7, color='0.25', linestyle=':' )
        #    label = cont #int(cont*100) # as percent
        #    ax0.text( phi_time_l[cc], 0.40, '{:2d}'.format(label), va='bottom', ha='center', rotation=90, bbox=dict(facecolor='white'), transform=trans )
        #ax0.text( 0.1, 0.9, '$\phi (\%)$', ha='center', va='bottom', transform=ax0.transAxes )
        #h1, = ax0.semilogx( timeMyr_a, CO2_atmos_a, color=red, linestyle='-', label=r'CO$_2$')
        h2, = ax0.loglog( timeMyr_a, H2O_atmos_a, color=blue, linestyle='-', label=r'H$_2$O')
        h3, = ax0.loglog( timeMyr_a, H2_atmos_a, color=green, linestyle='-', label=r'H$_2$')
        #h4, = ax0.semilogx( timeMyr_a, O2_atmos_a, color=green, linestyle='-', label=r'O$_2$')
        ax0b = ax0.twinx()
        h5, = ax0b.semilogx( timeMyr_a, phi_global, color=black, linestyle='--', label=r'Melt, $\phi_g$')

        fig_o.set_myaxes(ax0, title=title, ylabel=ylabel, xlabel=xlabel )#, xticks=xticks )
        #ax0.legend()
        #ax0b.legend()
        #handle_l = [h1,h2,h3,h5]
        fig_o.set_myaxes(ax0, title=title, ylabel=ylabel, xlabel=xlabel )#, xticks=xticks )
        #fig_o.set_mylegend(ax0, handle_l, loc='upper center', ncol=1 )

        ax0.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20) )
        ax0.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=(0.2,0.4,0.6,0.8), numticks=20))
        ax0.xaxis.set_minor_formatter(ticker.NullFormatter())
        ax0.set_xlim( *xlim )
        #ax0.set_ylim( -10, 500)
        ax0.yaxis.set_label_coords(-0.15,0.5)
        ax0b.set_ylabel( r'$\phi_g$', rotation=0 )
        ax0b.yaxis.set_label_coords(1.1,0.525)

    ##########
    # figure b
    ##########
    #print((CO2_liquid_kg_a+CO2_solid_kg_a) / CO2_total_kg)
    #print(CO2_total_kg, H2O_total_kg, H2_total_kg)
    if 1:
        title = r'\textbf{(b) Reservoir mass fraction}'
        #h5, = ax1.semilogx( timeMyr_a, mass_liquid_a / mass_mantle, 'k--', label='melt' )
        #h1, = ax1.semilogx( timeMyr_a, (CO2_liquid_kg_a+CO2_solid_kg_a) / CO2_total_kg, color=red, linestyle='-', label=r'CO$_2$ interior' )
        #h2, = ax1.semilogx( timeMyr_a, CO2_atmos_kg_a / CO2_total_kg, color=red, linestyle='--', label=r'CO$_2$ atmos' )
        #h2b, = ax1.semilogx( timeMyr_a, CO2_escape_kg_a / CO2_total_kg, color=red, linestyle=':', label='Escape' )
        h3, = ax1.loglog( timeMyr_a, H2O_liquid_kg_a / H2O_total_kg, color=blue, linestyle='-', label=r'H$_2$O liquid' )
        h4, = ax1.loglog( timeMyr_a, H2O_solid_kg_a / H2O_total_kg, color=blue, linestyle=':', label=r'H$_2$O solid' )
        h5, = ax1.loglog( timeMyr_a, H2O_atmos_kg_a / H2O_total_kg, color=blue, linestyle='--', label=r'H$_2$O atmos')
        #h6, = ax1.loglog( timeMyr_a, H2O_reaction_kg_a / H2O_total_kg, color=blue, linestyle='-.', label=r'H$_2$O (+ve) reaction')
        h7, = ax1.loglog( timeMyr_a, H2O_total2_kg_a / H2O_total_kg, color=black, linestyle='-', label=r'H$_2$O total')

        #h4b, = ax1.semilogx( timeMyr_a, H2O_escape_kg_a / H2O_total_kg, color=blue, linestyle=':', label='Atmos' )
        h8, = ax1.loglog( timeMyr_a, H2_liquid_kg_a / H2_total_kg, color=green, linestyle='-', label=r'H$_2$ liquid' )
        h9, = ax1.loglog( timeMyr_a, H2_solid_kg_a / H2_total_kg, color=green, linestyle=':', label=r'H$_2$ solid')
        h10, = ax1.loglog( timeMyr_a, H2_atmos_kg_a / H2_total_kg, color=green, linestyle='--', label=r'H$_2$ atmos')
        #h11, = ax1.loglog( timeMyr_a, -H2_reaction_kg_a / H2_total_kg, color=green, linestyle='-.', label=r'H$_2$ (-ve) reaction')
        h12, = ax1.loglog( timeMyr_a, H2_total2_kg_a / H2_total_kg, color=yellow, linestyle='-', label=r'H$_2$ total')

        fig_o.set_myaxes( ax1, title=title, ylabel='$x$', xlabel=xlabel) #,xticks=xticks )
        ax1.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20) )
        ax1.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=(0.2,0.4,0.6,0.8), numticks=20))
        ax1.xaxis.set_minor_formatter(ticker.NullFormatter())
        ax1.set_xlim( *xlim )
        ax1.set_ylim(1.0E-8, 2.0)
        ax1.legend()
        #handle_l = [h3,h4,h5,h6,h7,h8,h9,h10,h11,h12] #[h1,h2,h3,h4]#,h2b]
        handle_l = [h3,h4,h5,h7,h8,h9,h10,h12] #[h1,h2,h3,h4]#,h2b]
        #fig_o.set_mylegend( ax1, handle_l, loc='center left', ncol=1 )
        ax1.yaxis.set_label_coords(-0.1,0.47)

    ##########
    # figure c
    ##########
    title = r'\textbf{(c) Surface temp and emissivity}'
    ylabel = '$T_s$ (K)'
    yticks = range(500,4001,500)
    h1, = ax2.semilogx( timeMyr_a, temperature_surface_a, 'k-', label=r'Surface temp, $T_s$' )
    ax2b = ax2.twinx()
    h2, = ax2b.loglog( timeMyr_a, emissivity_a, 'k--', label=r'Emissivity, $\epsilon$' )
    fig_o.set_myaxes( ax2, title=title, xlabel=xlabel, ylabel=ylabel, yticks=yticks )#, xticks=xticks )
    ax2.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20) )
    ax2.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=(0.2,0.4,0.6,0.8), numticks=20))
    ax2.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax2.set_xlim( *xlim )
    ax2.yaxis.set_label_coords(-0.175,0.48)
    #ax2.set_ylim( 1050, 1850 )
    #ax2.set_xlim( 1E-5 , 1 )
    ax2.set_ylim(200,4100)
    #ax2b.set_ylim( 4E-4, 2E-3 )
    handle_l = [h1,h2]
    #fig_o.set_mylegend( ax2, handle_l, loc='upper right', ncol=1 )
    ax2b.set_ylabel( r'$\epsilon$', rotation=0)
    ax2b.yaxis.set_label_coords(1.1,0.52)

    ##########
    # figure d
    ##########
    #if 1:
    #    title = r'\textbf{(d) Fatm}'
    #    ylabel = '$F_{atm}$'
    #    yticks = (1.0E3,1.0E4,1.0E5)
    #    h1, = ax3.loglog( timeMyr_a, Fatm,'k', )
    #   fig_o.set_myaxes( ax3, title=title, xlabel=xlabel, ylabel=ylabel, yticks=yticks)#, xticks=xticks )
    #    ax3.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20) )
    #    ax3.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=(0.2,0.4,0.6,0.8), numticks=20))
    #   ax3.set_xlim( *xlim )

    #title = '(d) Emissivity'
    #ylabel = '$\epsilon$'
    #ax3.loglog( timeMyr_a, emissivity_a, 'k-' )
    #fig_o.set_myaxes( ax3, title=title, xlabel=xlabel, ylabel=ylabel, xticks=xticks )
    #ax3.yaxis.set_label_coords(-0.1,0.55)
    #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    #ax3.set_ylim( 1E-4, 1E-2 )
    #ax3.set_xlim( 1E-5, 1 )

    # output (effective) emissivity for Bonati et al. (2018), a&a paper
    #out_a = np.column_stack( (timeMyr_a, temperature_surface_a, emissivity_a ) )
    #np.savetxt( 'out.dat', out_a )

    #fig_o.savefig(6, jpg=True)
    fig_o.savefig(6)

#====================================================================
def main():

    plot_atmosphere()
    plt.show()

#====================================================================

if __name__ == "__main__":

    main()