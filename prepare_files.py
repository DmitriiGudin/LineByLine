import numpy as np
import csv
import os
import sys
import params


if __name__=='__main__':

    program_folder = params.program_folder
    if not(program_folder[-1]=='/'):
        programfolder+='/'
    padding_length = 9-len(params.star_name_prefix)
     
    spectra_list = np.genfromtxt (params.file_list, dtype=str, delimiter='\t')
    
    print "Retrieving the spectra files and performing generic renaming procedure..."
    for i, f in enumerate(spectra_list):
        os.system("cp "+f+" "+program_folder+params.star_name_prefix+str(i+1).zfill(padding_length)+".fits")
        print str(i+1) +" / "+str(len(spectra_list))+'\r',
        sys.stdout.flush()
    print "\n"
    print "..."
    
    print "Generating the lambda script "+params.lambda_file+"..."
    print "..."
    with open (program_folder+params.lambda_file, 'w') as f:
        writer=csv.writer(f, delimiter='\t')
        f.write('ctio')
        for i in range(len(spectra_list)):
            f.write('\n')
            f.write('lambda input='+params.star_name_prefix+str(i+1).zfill(padding_length)+'.fits start='+params.wave_start_header_varname+' delta='+params.wave_delta_header_varname+' number=no > '+params.star_name_prefix+str(i+1).zfill(padding_length)+'.dad')
    
    print "Generating the list of the generic star names names.list..."
    print "..."
    with open (program_folder+params.name_list, 'w') as f:
        writer=csv.writer(f, delimiter='\t')
        for i in range(len(spectra_list)):
            f.write(params.star_name_prefix+str(i+1).zfill(padding_length))
            f.write('\n')
    
    print "Compiling "+params.brad_fortrane_file+" to produce executable "+params.brad_exec_file+"..."
    print "..."
    os.system("gfortran -w "+params.brad_fortrane_file+" -o "+params.brad_exec_file)
    
    print "Building the IRAF run script "+params.run_script+"..."
    print "..."
    with open (params.run_script, 'w') as f:
        writer=csv.writer(f, delimiter='\t')
        f.write ('print "Generating the spectra list..."')
        f.write ('\n')
        f.write ('ls -1 '+params.star_name_prefix+'*.fits > '+params.generic_spectra_list)
        f.write ('\n')
        f.write ('print "Resampling the spectra..."')
        f.write ('\n')
        f.write ('noao')
        f.write ('\n')
        f.write ('onedspec')
        f.write ('\n')
        f.write ('dispcor input=@'+params.generic_spectra_list+' output=@'+params.generic_spectra_list+' w1='+str(params.resampling_min_wavelength)+' w2='+str(params.resampling_max_wavelength)+' dw='+str(params.resampling_del_wavelength))
        f.write ('\n')
        f.write ('print "Generating the spectral text files..."')
        f.write ('\n')
        f.write ('cl < '+ params.lambda_file)
        f.write ('\n')
        f.write ('print "Ready. Exit IRAF and run ./'+params.dad_clean+' (ONCE!!!) to continue."')
    
    print "Generating the fake RV-correction file rvc.out..."
    print "..."
    with open ('rvc.out', 'w') as f:
        writer=csv.writer(f, delimiter='\t')
        for i in range(0,5):
            f.write ('#')
            f.write ('\n')
        for i in range(len(spectra_list)):
            f.write ('0000000.00000 0.00 0.00 0.00 0.00 0.00 0.00 0.00')
            f.write ('\n')           
 
 
    print "Generating the *.dad file cleaner "+params.dad_clean+"..."
    print "..."
    with open (params.dad_clean, 'w') as f:
        writer=csv.writer(f, delimiter='\t')
        f.write ('for file in *.dad')
        f.write ('\n')    
        f.write ('do')
        f.write ('\n')
        f.write ("    sed -i -e '1,2d' $file")
        f.write ('\n')
        f.write ('done')
        f.write ('\n')
        f.write ('echo Ready. Run ./'+params.brad_exec_file+' and use '+params.name_list+' as the input.')
    os.system('chmod 777 '+params.dad_clean)

    print "Ready. Start IRAF and run 'cl < "+params.run_script+"' to proceed."
            
