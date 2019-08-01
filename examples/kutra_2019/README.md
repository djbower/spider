bu_input.opts
default 

CO2H2O_varML.opts
mixing length: 1(variable)
volatiles tracked: CO2,H2O

CO2H2O_constML.opts 
mixing length: 2(constant)
volatiles tracked: CO2,H2O

constML.opts
mixing length: 2(constant)
volatiles tracked: CO2,H2O,H2

varML.opts
mixing length: 1(variable)
volatiles tracked: CO2,H2O,H2

constML_1kabs.opts 
mixing length: 2(constant)
volatiles tracked: CO2,H2O,H2
H2_kabs = 1.0

constML_0kabs.opts 
mixing length: 2(constant)
volatiles tracked: CO2,H2O,H2
H2_kabs = 0.0

constML_coupled.opts
mixing length: 2(constant)
volatiles tracked: CO2, H2O, H2
coefficients for coupling: added 