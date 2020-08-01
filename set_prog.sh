#!/bin/sh

dir=current_directory

cd ${dir}/gen_pdf
make

cd ${dir}/save_mass_variance
make

cd ${dir}/fit
make fit_mean
make fit_variance