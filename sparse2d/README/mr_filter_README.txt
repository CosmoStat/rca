The progam read/write FITS files. THe format is supported by IDL/Matlab and 
the fitsio library is available here: http://heasarc.gsfc.nasa.gov/fitsio

To run the filtering, the command is
$ISAP/sparse2d/bin/mr_filter in_map out_map
An online is help is obtained just by running mr_filter without options.

 
Example of command:

- Isotropic wavelet filtering with 5 wavelet scales
mr_filter -n5  input.fits  output.fits 


- filtering with the 7/9 filter bank, 5 wavelets and 3sigma hard thresholding 
mr_filter -n5 t24 input.fits  output.fits 



- filtering with the 7/9 filter bank, 5 wavelets and multiscale Wiener 
mr_filter -n5 t24 -f3 input.fits  output.fits 



Wavelet filtering is described in the book:
 J.-L. Starck , F. Murtagh and J. Fadili, 
 "Sparse Image & Signal Processing: wavelets, curvelets, morphological diversity ",  
 Cambridge University Press, Cambridge (GB),  2010.
   
