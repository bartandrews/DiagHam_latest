#!/bin/bash

runs() {
echo ~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --u-potential 0.36787944117144233 --v-potential 0.0 --v2-potential 0.0 --v3-potential 0.0
echo ~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --u-potential 0.36787944117144233 --v-potential 0.0049787068367863905 --v2-potential 3.059023205018258e-08 --v3-potential 3.7751345442790845e-12
echo ~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --u-potential 0.36787944117144233 --v-potential 0.009957413673572781 --v2-potential 6.118046410036516e-08 --v3-potential 7.550269088558169e-12
echo ~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --u-potential 0.36787944117144233 --v-potential 0.014936120510359172 --v2-potential 9.177069615054775e-08 --v3-potential 1.1325403632837253e-11
echo ~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --u-potential 0.36787944117144233 --v-potential 0.019914827347145562 --v2-potential 1.2236092820073031e-07 --v3-potential 1.5100538177116338e-11
echo ~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --u-potential 0.36787944117144233 --v-potential 0.024893534183931948 --v2-potential 1.529511602509129e-07 --v3-potential 1.887567272139542e-11
echo ~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --u-potential 0.36787944117144233 --v-potential 0.029872241020718344 --v2-potential 1.835413923010955e-07 --v3-potential 2.2650807265674507e-11
echo ~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --u-potential 0.36787944117144233 --v-potential 0.034850947857504734 --v2-potential 2.1413162435127808e-07 --v3-potential 2.642594180995359e-11
echo ~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --u-potential 0.36787944117144233 --v-potential 0.039829654694291124 --v2-potential 2.4472185640146063e-07 --v3-potential 3.0201076354232676e-11
echo ~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --u-potential 0.36787944117144233 --v-potential 0.04480836153107751 --v2-potential 2.753120884516432e-07 --v3-potential 3.3976210898511755e-11
echo ~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -X 18 -Y 9 -x 3 -y 6 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --u-potential 0.36787944117144233 --v-potential 0.049787068367863896 --v2-potential 3.059023205018258e-07 --v3-potential 3.775134544279084e-11
}
export -f runs
runs | nohup nice parallel -j 4 &
