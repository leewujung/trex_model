function run_field(mode_fname)

runfield = which( 'field.exe' );
flp_filename = ['field_',mode_fname];
eval( [ '! "' runfield '" ' mode_fname ' < ' flp_filename '.flp > ' flp_filename '.prt' ] );
