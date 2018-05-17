#CC="-L/usr/lib/gcc/i386-redhat-linux/4.0.1 -lstdc++"
#export FMKMF_F90="ifort -I/opt/intel/mkl/8.0/include -no-ipo -no-prec-div -xP"
#export FMKMF_F90="ifort  -I$MKL_ROOT $galahad_modules -O4 -ipo"
#export FMKMF_F90="ifort  -I$MKL_ROOT $galahad_modules -O4 -parallel -no-ipo"
# seems -O0 causes trouble in case of RNase:
export FMKMF_F90="ifort  $lapack95_modules $galahad_modules -openmp -g -C -traceback -no-ipo"
#export FMKMF_LINKOPTS="-no-ipo"
fmkmf -tag f90 \
-l "-lstdc++ $CCP4_LIBRARIES -L/usr/lib64 -L/usr/local/lib -lgrafpack -lfftw3_threads -lfftw3 $galahad_mkllib -L/usr/local/lib -larpack_Linux -L. -lro" \
fmlsq.f90
