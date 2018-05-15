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
    fig_o = su.FigureData( 2, 2, width, height, 'atmosphere' )

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[0][1]
    ax2 = fig_o.ax[1][0]
    ax3 = fig_o.ax[1][1]

    fig_o.time = su.get_all_output_times()
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
    CO2_liquid_kg_a = su.ppm_to_mass_fraction( CO2_liquid_a ) * mass_liquid_a
    CO2_solid_kg_a = su.ppm_to_mass_fraction( CO2_solid_a ) * mass_solid_a
    CO2_total_kg_a = su.ppm_to_mass_fraction( CO2_initial_a ) * mass_mantle
    CO2_total_kg = CO2_total_kg_a[0] # time-independent
    # TODO: below mass is conserved by definition, but can also
    # compute directly from partial pressure
    CO2_atmos_kg_a = CO2_total_kg - CO2_liquid_kg_a - CO2_solid_kg_a

    H2O_liquid_kg_a = su.ppm_to_mass_fraction( H2O_liquid_a ) * mass_liquid_a
    H2O_solid_kg_a = su.ppm_to_mass_fraction( H2O_solid_a ) * mass_solid_a
    H2O_total_kg_a = su.ppm_to_mass_fraction( H2O_initial_a ) * mass_mantle
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
        time_temp_l = su.find_xx_for_yy( timeMyr_a, phi_a, cont )
        index = su.get_first_non_zero_index( time_temp_l )
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
def get_single_values_for_times( field, time_l ):

    data_l = []

    for nn, time in enumerate( time_l ):

        # read json
        myjson_o = su.MyJSON( 'output/{}.json'.format(time) )

        yy = myjson_o.get_scaled_field_values( field )
        data_l.append( yy )

    data_a = np.array( data_l )

    return data_a

#====================================================================
def main():

    plot_atmosphere()
    plt.show()

#====================================================================

if __name__ == "__main__":

    main()