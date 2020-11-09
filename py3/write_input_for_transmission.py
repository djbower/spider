#!/usr/bin/env python

import numpy as np
import spider_utils as su

#====================================================================
def main():

    #casenum = '1m'
    #casenum = '9m'
    # for Mercury models
    casenum = 'SV'
    #casenum = 'LV'
    strcase = 'case' + str(casenum)

    # this is for short-cutting to generate data for transmission
    # and emission spectra
    if 1:
        if casenum == '1m':
            # Bower et al. (2019), case 1
            #time_l = [0,94400,5450000,78400000,560000000,4550000000]
            time_l = [0,94300,250000,500000,2250000]
        elif casenum == '9m':
            # Bower et al. (2019), case 9
            #time_l = [0,3850000,12050000,141000000,802000000,4550000000]
            time_l = [0,3850000,11150000,16300000,21250000]
        elif casenum == 'S3':
            time_l = [9437]
        # for Mercury models
        # 2400 K, 2000 K, 1600 K
        elif casenum == 'SV':
            time_l = [30,1302,3550]
        elif casenum == 'LV':
            time_l = [50,2609,7100]

    else:
        # this is for looping over all date to compute the static structure
        # as a function of time
        time_l = su.get_all_output_times()

    # planetary radii estimates
    outradiusfile = open(strcase+'_evolving_radius.dat','w')
    outradiusfile.write( '# time (yr), radius (m), 10mB height (m), 1mB height (m), surface temp (K), phi_g\n' )

    nn_total = len(time_l)

    for nn, time in enumerate(time_l):

        print( 'processing output {} of {}'.format( nn, nn_total ) )
        percent = np.around( float(nn)/float(nn_total)*100.0, decimals=2 )
        print( '{} percent complete'.format(percent) )
        outfilename = strcase + '_input_transmission_{}yr.dat'.format(time)
        print( outfilename )
        outfile = open(outfilename,'w')

        myjson_o = su.MyJSON( 'output/{}.json'.format(time) )

        # commented out some sections for Mercury models, since radius not compatible
        #rho_interp1d = myjson_o.get_rho_interp1d()
        #radius = su.solve_for_planetary_radius( rho_interp1d )

        #outfile.write('planetary radius (m)= {}\n'.format(radius))
        #pres_interp1d = myjson_o.get_atm_struct_depth_interp1d()
        #height1mb = pres_interp1d( 1.0E-3 )
        #height10mb = pres_interp1d( 1.0E-2 )

        #outfile.write('height of 1mb contour (m)= {}\n'.format( height1mb ) )
        #outfile.write('height of 10mb contour (m)= {}\n'.format( height10mb ) )
        write_value_to_file( myjson_o, outfile, ['atmosphere','molar_mass'], time, 'atmosphere mean molar mass' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','CO2','molar_mass'], time, 'CO2 molar mass' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','H2O','molar_mass'], time, 'H2O molar mass' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','CO','molar_mass'], time, 'CO molar mass' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','H2','molar_mass'], time, 'H2 molar mass' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','CO2','mixing_ratio'], time, 'CO2 mixing ratio' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','H2O','mixing_ratio'], time, 'H2O mixing ratio' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','CO','mixing_ratio'], time, 'CO mixing ratio' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','H2','mixing_ratio'], time, 'H2 mixing ratio' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','CO2','atmosphere_kg'], time, 'CO2 mass (kg)' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','H2O','atmosphere_kg'], time, 'H2O mass (kg)' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','CO','atmosphere_kg'], time, 'CO mass (kg)' )
        write_value_to_file( myjson_o, outfile, ['atmosphere','H2','atmosphere_kg'], time, 'H2 mass (kg)' )

        outfile.write('Temp (K), Pressure (bar), Optical depth, Height above surface (m)\n')

        temperature_a = myjson_o.get_dict_values( ['atmosphere','atm_struct_temp'] )
        pressure_a = myjson_o.get_dict_values( ['atmosphere','atm_struct_pressure'] )
        opticaldepth_a = myjson_o.get_dict_values( ['atmosphere','atm_struct_tau'] )
        height_a = myjson_o.get_dict_values( ['atmosphere','atm_struct_depth'] )
        npoints = np.size( temperature_a ) # number of points in the atmosphere
        surface_temperature = temperature_a[-1]

        phi_g = myjson_o.get_dict_values( ['rheological_front_phi','phi_global'] )

        outfile.write('npoints= {}\n'.format(npoints))
        outfile.close()

        atmos_struct_a = np.column_stack( (temperature_a, pressure_a, opticaldepth_a, height_a ) )
        f = open( outfilename, 'ab' ) # must open in binary append mode
        np.savetxt( f, atmos_struct_a )
        f.close()

        # write radius data to file
        #outradiusfile.write( '{} {} {} {} {} {}\n'.format( time, radius, radius+height10mb, radius+height1mb, surface_temperature, phi_g ) )

    #outradiusfile.close()

#====================================================================
def write_value_to_file( myjson_o, outfile, keys, time, label ):

    qty = myjson_o.get_dict_values( keys )

    outfile.write( '{}= {}\n'.format(label, qty) )

#====================================================================

if __name__ == "__main__":

    main()
