#!/bin/bash

mkdir -p /home/bart/DiagHam_latest/run/SFCI_data_2/plots/finite_size/script_output/q128/n6
cd /home/bart/DiagHam_latest/run/SFCI_data_2/plots/finite_size/script_output/q128/n6
/home/bart/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel  -p 6 -X 16 -Y 8 -x 3 -y 6 --t3hop 1.31 --t6hop -0.25 --t9hop -0.25 -m 96000 -S --processors 6 -n 5 --lanczos-precision 1e-10 --u-potential 1 --v-potential 0 --v2-potential 0 --v3-potential 0 --auto-addprojector;
/home/bart/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel  -p 6 -X 16 -Y 8 -x 3 -y 6 --t3hop -1.81 --t6hop 0.25 --t9hop 0.25 -m 96000 -S --processors 6 -n 5 --lanczos-precision 1e-10 --u-potential 1 --v-potential 0 --v2-potential 0 --v3-potential 0 --auto-addprojector;
mkdir -p /home/bart/DiagHam_latest/run/SFCI_data_2/plots/finite_size/script_output/q128/n7
cd /home/bart/DiagHam_latest/run/SFCI_data_2/plots/finite_size/script_output/q128/n7
/home/bart/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel  -p 7 -X 16 -Y 8 -x 3 -y 7 --t3hop 1.31 --t6hop -0.25 --t9hop -0.25 -m 96000 -S --processors 6 -n 5 --lanczos-precision 1e-10 --u-potential 1 --v-potential 0 --v2-potential 0 --v3-potential 0 --auto-addprojector;
/home/bart/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel  -p 7 -X 16 -Y 8 -x 3 -y 7 --t3hop -1.81 --t6hop 0.25 --t9hop 0.25 -m 96000 -S --processors 6 -n 5 --lanczos-precision 1e-10 --u-potential 1 --v-potential 0 --v2-potential 0 --v3-potential 0 --auto-addprojector;
mkdir -p /home/bart/DiagHam_latest/run/SFCI_data_2/plots/finite_size/script_output/q128/n8
cd /home/bart/DiagHam_latest/run/SFCI_data_2/plots/finite_size/script_output/q128/n8
/home/bart/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel  -p 8 -X 16 -Y 8 -x 3 -y 8 --t3hop 1.31 --t6hop -0.25 --t9hop -0.25 -m 96000 -S --processors 6 -n 5 --lanczos-precision 1e-10 --u-potential 1 --v-potential 0 --v2-potential 0 --v3-potential 0 --auto-addprojector;
/home/bart/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel  -p 8 -X 16 -Y 8 -x 3 -y 8 --t3hop -1.81 --t6hop 0.25 --t9hop 0.25 -m 96000 -S --processors 6 -n 5 --lanczos-precision 1e-10 --u-potential 1 --v-potential 0 --v2-potential 0 --v3-potential 0 --auto-addprojector;
mkdir -p /home/bart/DiagHam_latest/run/SFCI_data_2/plots/finite_size/script_output/q128/n9
cd /home/bart/DiagHam_latest/run/SFCI_data_2/plots/finite_size/script_output/q128/n9
/home/bart/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel  -p 9 -X 64 -Y 2 -x 1 -y 27 --t3hop 1.31 --t6hop -0.25 --t9hop -0.25 -m 96000 -S --processors 6 -n 5 --lanczos-precision 1e-10 --u-potential 1 --v-potential 0 --v2-potential 0 --v3-potential 0 --auto-addprojector;
/home/bart/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel  -p 9 -X 64 -Y 2 -x 1 -y 27 --t3hop -1.81 --t6hop 0.25 --t9hop 0.25 -m 96000 -S --processors 6 -n 5 --lanczos-precision 1e-10 --u-potential 1 --v-potential 0 --v2-potential 0 --v3-potential 0 --auto-addprojector;
