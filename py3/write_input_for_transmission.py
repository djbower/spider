#!/usr/bin/env python

import numpy as np
import spider_utils as su

#====================================================================
def main():

    time_l = su.get_all_output_times()
    # FIXME: hack for testing
    #time_l = time_l[0:1]
    #time_l = [0,150000,250000,450000,100000000]

    # times for reference model in paper
    time_l = [0,150000,500000,575000000,4550000000]

    # planetary radii estimates
    outradiusfile = open('radius.dat','w')
    outradiusfile.write( '# time (yr), radius (m), 10mB height (m), 1mB height (m)\n' )

    for time in time_l:

        outfilename = 'input_transmission_{}yr.dat'.format(time)
        print( outfilename )
        outfile = open(outfilename,'w')

        myjson_o = su.MyJSON( 'output/{}.json'.format(time) )
        rho_interp1d = myjson_o.get_rho_interp1d()
        radius = su.solve_for_planetary_radius( rho_interp1d )

        outfile.write('planetary radius (m)= {}\n'.format(radius))
        pres_interp1d = myjson_o.get_atm_struct_depth_interp1d()
        height1mb = pres_interp1d( 1.0E-3 )
        height10mb = pres_interp1d( 1.0E-2 )

        outfile.write('height of 1mb contour (m)= {}\n'.format( height1mb ) )
        outfile.write('height of 10mb contour (m)= {}\n'.format( height10mb ) )
        write_value_to_file( myjson_o, outfile, ['atmosphere','molar_mass'], time, 'atmosphere mean molar mass' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','CO2','molar_mass'], time, 'CO2 molar mass' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','H2O','molar_mass'], time, 'H2O molar mass' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','CO2','mixing_ratio'], time, 'CO2 mixing ratio' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','H2O','mixing_ratio'], time, 'H2O mixing ratio' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','CO2','atmosphere_kg'], time, 'CO2 mass (kg)' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','H2O','atmosphere_kg'], time, 'H2O mass (kg)' )

        outfile.write('Temp (K), Pressure (bar), Optical depth, Height above surface (m)\n')

        temperature_a = myjson_o.get_dict_values( ['atmosphere','atm_struct_temp'] )
        pressure_a = myjson_o.get_dict_values( ['atmosphere','atm_struct_pressure'] )
        opticaldepth_a = myjson_o.get_dict_values( ['atmosphere','atm_struct_tau'] )
        height_a = myjson_o.get_dict_values( ['atmosphere','atm_struct_depth'] )
        npoints = np.size( temperature_a ) # number of points in the atmosphere

        outfile.write('npoints= {}\n'.format(npoints))
        outfile.close()

        atmos_struct_a = np.column_stack( (temperature_a, pressure_a, opticaldepth_a, height_a ) )
        f = open( outfilename, 'ab' ) # must open in binary append mode
        np.savetxt( f, atmos_struct_a )
        f.close()

        # write radius data to file
        outradiusfile.write( '{} {} {} {}\n'.format( time, radius, radius+height10mb, radius+height1mb ) )

    outradiusfile.close()

#====================================================================
def write_value_to_file( myjson_o, outfile, keys, time, label ):

    qty = myjson_o.get_dict_values( keys )

    outfile.write( '{}= {}\n'.format(label, qty) )

#====================================================================

if __name__ == "__main__":

    main()
