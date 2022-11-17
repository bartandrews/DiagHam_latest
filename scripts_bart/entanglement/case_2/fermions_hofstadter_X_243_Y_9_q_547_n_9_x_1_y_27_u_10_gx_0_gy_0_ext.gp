set xrange [-1:28]
set yrange [-12.7155235759771:-12.7154465478575]
set xlabel "k_x*L_y + k_y" font "Helvetica,28"
set ylabel "E(k)" font "Helvetica,28"
#set size 1.0, 0.6
set terminal postscript eps size 3.5,3 enhanced colour "Helvetica" 16
set output "fermions_hofstadter_X_243_Y_9_q_547_n_9_x_1_y_27_u_10_gx_0_gy_0_ext.eps"
plot "fermions_hofstadter_X_243_Y_9_q_547_n_9_x_1_y_27_u_10_gx_0_gy_0_ext.dat" using ($1*27+$2):3 title "L_x=1, L_y=27, MUC: p=547, q=243 x 9" with points lc rgbcolor "red" lt 6 pt 1 ps 2
