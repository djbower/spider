#!/usr/bin/env python

import logging
import spider_utils as su
import matplotlib.transforms as transforms
import matplotlib.pyplot as plt
import numpy as np
import os

logger = su.get_my_logger(__name__)

#====================================================================
def plot_atmosphere():

    logger.info( 'building atmosphere' )

    width = 4.7747
    height = 4.7747
    fig_o = su.FigureData( 2, 2, width, height, 'atmosphere', units='Myr' )

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]
    ax2 = fig_o.ax[1][0]
    ax3 = fig_o.ax[1][1]

    fig_o.time = su.get_all_output_times()
    timeMyr_a = fig_o.time * 1.0E-6 # Myrs

    keys_t = ( ('atmosphere','mass_liquid'),
               ('atmosphere','mass_solid'),
               ('atmosphere','mass_mantle'),
               ('atmosphere','CO2','liquid_kg'),
               ('atmosphere','CO2','solid_kg'),
               ('atmosphere','CO2','initial_kg'),
               ('atmosphere','CO2','atmosphere_kg'),
               ('atmosphere','CO2','atmosphere_bar'),
               ('atmosphere','H2O','liquid_kg'),
               ('atmosphere','H2O','solid_kg'),
               ('atmosphere','H2O','initial_kg'),
               ('atmosphere','H2O','atmosphere_kg'),
               ('atmosphere','H2O','atmosphere_bar'),
               ('atmosphere','temperature_surface'),
               ('atmosphere','emissivity') )

    data_a = su.get_dict_surface_values_for_times( keys_t, fig_o.time )

    mass_liquid_a = data_a[0,:]
    mass_solid_a = data_a[1,:]
    mass_mantle_a = data_a[2,:]
    mass_mantle = mass_mantle_a[0] # time independent

    # compute total mass (kg) in each reservoir
    CO2_liquid_kg_a = data_a[3,:]
    CO2_solid_kg_a = data_a[4,:]
    CO2_total_kg_a = data_a[5,:]
    CO2_total_kg = CO2_total_kg_a[0] # time-independent
    CO2_atmos_kg_a = data_a[6,:]
    CO2_atmos_a = data_a[7,:]
    CO2_escape_kg_a = CO2_total_kg - CO2_liquid_kg_a - CO2_solid_kg_a - CO2_atmos_kg_a

    H2O_liquid_kg_a = data_a[8,:]
    H2O_solid_kg_a = data_a[9,:]
    H2O_total_kg_a = data_a[10,:]
    H2O_total_kg = H2O_total_kg_a[0] # time-independent
    H2O_atmos_kg_a = data_a[11,:]
    H2O_atmos_a = data_a[12,:]
    H2O_escape_kg_a = H2O_total_kg - H2O_liquid_kg_a - H2O_solid_kg_a - H2O_atmos_kg_a

    temperature_surface_a = data_a[13,:] #su.get_dict_values_for_times( ['atmosphere','temperature_surface'], fig_o.time )
    emissivity_a = data_a[14,:] #su.get_dict_values_for_times( ['atmosphere','emissivity'], fig_o.time )

    #xticks = [1E-5,1E-4,1E-3,1E-2,1E-1]#,1]
    xticks = [1.0E-2, 1.0E-1, 1.0E0, 1.0E1, 1.0E2,1.0E3,5.0E3] #[1E-6,1E-4,1E-2,1E0,1E2,1E4,1E6]#,1]
    xlabel = 'Time (Myr)'

    red = fig_o.get_color(3)
    blue = fig_o.get_color(-3)

    # to plot melt fraction contours on figure (a)
    # compute time at which desired melt fraction is reached
    phi_a = mass_liquid_a / mass_mantle
    phi_cont_l = [0.2, 0.8] # 20% and 80% melt fraction
    phi_time_l = [] # contains the times for each contour
    for cont in phi_cont_l:
        time_temp_l = su.find_xx_for_yy( timeMyr_a, phi_a, cont )
        index = su.get_first_non_zero_index( time_temp_l )
        if index is None:
            out = 0.0
        else:
            out = timeMyr_a[index]
        phi_time_l.append( out )

    ##########
    # figure a
    ##########
    if 1:
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
    if 1:
        title = '(b) Volatile reservoirs'
        #h5, = ax1.semilogx( timeMyr_a, mass_liquid_a / mass_mantle, 'k--', label='melt' )
        h1, = ax1.semilogx( timeMyr_a, (CO2_liquid_kg_a+CO2_solid_kg_a) / CO2_total_kg, color=red, linestyle='-', label='Magma' )
        h2, = ax1.semilogx( timeMyr_a, CO2_atmos_kg_a / CO2_total_kg, color=red, linestyle='--', label='Atmos' )
        h2b, = ax1.semilogx( timeMyr_a, CO2_escape_kg_a / CO2_total_kg, color=red, linestyle=':', label='Escape' )
        h3, = ax1.semilogx( timeMyr_a, (H2O_liquid_kg_a+H2O_solid_kg_a) / H2O_total_kg, color=blue, linestyle='-', label='Magma' )
        h4, = ax1.semilogx( timeMyr_a, H2O_atmos_kg_a / H2O_total_kg, color=blue, linestyle='--', label='Atmos')
        h4b, = ax1.semilogx( timeMyr_a, H2O_escape_kg_a / H2O_total_kg, color=blue, linestyle=':', label='Atmos' )
        fig_o.set_myaxes( ax1, title=title, ylabel='$x$', xticks=xticks )
        handle_l = [h1,h2,h2b]
        fig_o.set_mylegend( ax1, handle_l, loc='center left', ncol=1, TITLE="" )
        ax1.yaxis.set_label_coords(-0.1,0.46)

    ##########
    # figure c
    ##########
    title = '(c) Surface temperature'
    ylabel = '$T$'
    yticks = range(200,2601,400)
    ax2.semilogx( timeMyr_a, temperature_surface_a, 'k-' )
    fig_o.set_myaxes( ax2, title=title, xlabel=xlabel, ylabel=ylabel, xticks=xticks, yticks=yticks )
    ax2.yaxis.set_label_coords(-0.1,0.5)
    #ax2.set_ylim( 1050, 1850 )
    #ax2.set_xlim( 1E-5 , 1 )
    ax2.set_ylim(260,2600)

    ##########
    # figure d
    ##########
    title = '(d) Emissivity'
    ylabel = '$\epsilon$'
    ax3.loglog( timeMyr_a, emissivity_a, 'k-' )
    fig_o.set_myaxes( ax3, title=title, xlabel=xlabel, ylabel=ylabel, xticks=xticks )
    ax3.yaxis.set_label_coords(-0.1,0.55)
    #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax3.set_ylim( 1E-4, 1E-2 )
    #ax3.set_xlim( 1E-5, 1 )

    # output (effective) emissivity for Bonati et al. (2018), a&a paper
    #out_a = np.column_stack( (timeMyr_a, temperature_surface_a, emissivity_a ) )
    #np.savetxt( 'out.dat', out_a )

    fig_o.savefig(6)

#====================================================================
def main():

    plot_atmosphere()
    plt.show()

#====================================================================

if __name__ == "__main__":

    main()
