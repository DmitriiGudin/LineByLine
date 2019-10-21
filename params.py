# USER

program_folder = '/afs/crc.nd.edu/user/d/dgudin/LineByLine_xcor/RAVE/' # The location of the folder containing the program (including this file).
file_list = '/scratch365/dgudin/RAVE_spectra/RAVE_file_list.txt' # The list of the spectral files to process.
resampling_min_wavelength = 3700 # Minimum wavelength in the resampled spectra.
resampling_max_wavelength = 5000 # Maximum wavelength in the resampled spectra.
resampling_del_wavelength = 0.5 # The distance between the wavelengths in the resampled spectra.




# IRAF

wave_start_header_varname = 'CRVAL1' # The variable in the spectra headers corresponding to the starting wavelength.
wave_delta_header_varname = 'CD1_1' # The variable in the spectra headers corresponding to the difference between neighboring wavelengths.




# DEV

lambda_file = 'lambda.cl' # Lambda script.
name_list = 'name.list' # Generic star name list.
generic_spectra_list = 'S.list' # Generic spectra list.
star_name_prefix = 'S' # Beginning of the generic star name.
ablines_file = 'ablines.dat' # ablines.dat location.
brad_fortrane_file = 'hbrad2.for' # hbrad2.for location.
brad_exec_file = 'hbrad2.run' # Executable hbrad2 file.
run_script = 'run.cl' # IRAF run script.
dad_clean = 'dad_clean.sh' # Cleaner of the dad-files.


