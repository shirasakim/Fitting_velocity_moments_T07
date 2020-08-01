#!/bin/sh

dir=current_directory

cd ${dir}/gen_pdf
make clean

cd ${dir}/save_mass_variance
make clean

cd ${dir}/fit
make clean