#!/bin/bash

runs() {
echo ~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 20 --lanczos-precision 1e-10 --atan_0-potential 1
echo ~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 20 --lanczos-precision 1e-10 --atan_1_4-potential 1
echo ~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 20 --lanczos-precision 1e-10 --atan_1_3-potential 1
echo ~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 20 --lanczos-precision 1e-10 --atan_1_2-potential 1
echo ~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 20 --lanczos-precision 1e-10 --atan_2_3-potential 1
echo ~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 20 --lanczos-precision 1e-10 --atan_3_4-potential 1
echo ~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 20 --lanczos-precision 1e-10 --atan_1-potential 1
}
export -f runs
runs | nohup nice parallel -j 4 > /dev/null &
