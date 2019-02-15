#!/usr/bin/env python

import numpy as np
import spider_utils as su

#====================================================================
def main():

    time_l = su.get_all_output_times()
    # FIXME: hack for testing
    time_l = time_l[0:1]
    time_l = [0,150000,250000,450000,100000000]

    for time in time_l:

        outfilename = 'input_transmission_{}yr.dat'.format(time)
        outfile = open(outfilename,'w')

        myjson_o = su.MyJSON( 'output/{}.json'.format(time) )
        rho_interp1d = myjson_o.get_rho_interp1d()
        radius = su.solve_for_planetary_radius( rho_interp1d )

        outfile.write('planetary radius (m)= {}\n'.format(radius))
        write_value_to_file( myjson_o, outfile, ['atmosphere','molecular_mass'], time, 'atmosphere mean molecular mass' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','CO2','molecular_mass'], time, 'CO2 molecular mass' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','H2O','molecular_mass'], time, 'H2O molecular mass' )
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

#====================================================================
def write_value_to_file( myjson_o, outfile, keys, time, label ):

    qty = myjson_o.get_dict_values( keys )

    outfile.write( '{}= {}\n'.format(label, qty) )

#====================================================================

if __name__ == "__main__":

    main()
