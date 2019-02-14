#!/usr/bin/env python

import numpy as np
import spider_utils as su

#====================================================================
def main():

    #time = 0
    #time = 150000
    #time = 250000
    #time = 450000
    time = 100000000

    myjson_o = su.MyJSON( 'output/{}.json'.format(time) )

    atmos_d = myjson_o.data_d['atmosphere']

    outfile = open('input_transmission_at_{}yr.dat'.format(time),'w')

    # need to manually add planetary radii from static structure calculation

    radius = 6371000.0
    outfile.write('planetary radius (m)= {}\n'.format(radius))

    write_value_to_file( outfile, ['molecular_mass'], time, 'atmosphere mean molecular mass' )
    write_value_to_file( outfile, ['CO2','molecular_mass'], time, 'CO2 molecular mass' )
    write_value_to_file( outfile, ['H2O','molecular_mass'], time, 'H2O molecular mass' )
    write_value_to_file( outfile, ['CO2','atmosphere_kg'], time, 'CO2 mass (kg)' )
    write_value_to_file( outfile, ['H2O','atmosphere_kg'], time, 'H2O mass (kg)' )

    outfile.write('Temp (K), Pressure (bar), Optical depth, Height above surface (m)\n')

    # TODO: must manually update number of atmosphere points here
    outfile.write('npoints= 500\n')

    temp_a = get_single_value_at_time( ['atm_struct_temp'], time )
    pres_a = get_single_value_at_time( ['atm_struct_pressure'], time )
    tau_a = get_single_value_at_time( ['atm_struct_tau'], time )
    depth_a = get_single_value_at_time( ['atm_struct_depth'], time )

    data_a = np.column_stack( (temp_a, pres_a, tau_a, depth_a ) ) 

    for entry in data_a:
        col0 = entry[0]
        col1 = entry[1]
        col2 = entry[2]
        col3 = entry[3]
        outfile.write( '{} {} {} {}\n'.format(col0,col1,col2,col3) )

    #print( data_a)

    #write_value_to_file( outfile, ['atm_struct_temp'], time, 'Temperature (K)' )
    #write_value_to_file( outfile, ['atm_struct_pressure'], time, 'Pressure (bar)' )
    #write_value_to_file( outfile, ['atm_struct_tau'], time, 'Optical depth' )
    #write_value_to_file( outfile, ['atm_struct_depth'], time, 'Height (m)' )

    outfile.close()

#====================================================================
def get_single_values_for_times( keys, time_l ):

    data_l = []

    for nn, time in enumerate( time_l ):

        data_l.append( get_single_value_at_time( keys, time ) )

    data_a = np.array( data_l )

    return data_a

#====================================================================
def write_value_to_file( outfile, keys, time, label ):

    qty = get_single_value_at_time( keys, time )

    outfile.write( '{}= {}\n'.format(label, qty) )

#====================================================================
def get_single_value_at_time( keys, time ):

    # read json
    myjson_o = su.MyJSON( 'output/{}.json'.format(time) )

    # all atmosphere-related output is stored here
    atmos_d = myjson_o.data_d['atmosphere']

    fdata_d = recursive_get(atmos_d,keys)
    scaling = float(fdata_d['scaling'])

    if len( fdata_d['values'] ) == 1:
        values_a = float( fdata_d['values'][0] )
    else:
        values_a = np.array( [float(value) for value in fdata_d['values']] )

    scaled_values_a = scaling * values_a

    return scaled_values_a

#====================================================================

if __name__ == "__main__":

    main()
