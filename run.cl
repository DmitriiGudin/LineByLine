print Generating the spectra list...
ls -1 S*.fits > S.list
print Resampling the spectra...
noao
onedspec
dispcor input=@S.list output=@S.list w1=3700 w2=5000 dw=0.5
print Generating the spectral text files...
cl < lambda.cl
print Ready. Run ./hbrad2.run and use name.list as the input.
            
