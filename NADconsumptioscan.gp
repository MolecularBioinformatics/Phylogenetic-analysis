set terminal postscript eps color enhanced 'Arial' 12
set output 'NADconsumptionscan_wdiff.eps'
set multiplot

set xlabel 'cell devision rate per hour'

set log x
#set format x "%g"

set key right
set size 0.5,0.5
set origin 0.0,0.5

set xrange [0.000001*3600:0.0001*3600]



set yrange [0:0.3]
set title 'NAD consumption'
set ylabel 'NAD consumption flux ({/Symbol m}M/s)'
set label 'A' font 'Arial,16' at screen 0.0, screen 0.98

plot 'NADconsumptionscan_0_10.txt' using ($6*3600):($8*1000) title '10 U NC NamPRT only' w lines lc rgb "#00AA00" lw 3 dt 1, 'NADconsumptionscan_10_10.txt' using ($6*3600):($8*1000) title '10 U NC  NamPRT + NNMT' w lines  lc  rgb "#3333FF" lw 3 dt 1,  'enyzme_contribution_0_0_200.txt' using ($6*3600):($8*1000) title '50 U NC NamPRT only' w lines lc rgb "#00AA00" lw 3 dt 3, 'enyzme_contribution_0_10_200.txt' using ($6*3600):($8*1000) title '50 U NC NamPRT + NNMT' w lines  lc  rgb "#3333FF" lw 3 dt 3,  'NADconsumptionscan_0_200.txt' using ($6*3600):($8*1000) title '200 U NC NamPRT only' w lines lc rgb "#00AA00" lw 3 dt 5, 'NADconsumptionscan_10_200.txt' using ($6*3600):($8*1000) title '200 U NC NamPRT + NNMT' w lines  lc  rgb "#3333FF" lw 3 dt 5;




set size 0.5,0.5
set origin 0.0,0.0

set yrange [0:0.14]
set title 'NAD consumption benefit with NNMT'
set ylabel 'NAD consumption flux benefit ({/Symbol m}M/s)'
set label 'C' font 'Arial,16' at screen 0.0, screen 0.48

plot  'NADconsumptionscan_0_10-NADconsumptionscan_10_10.txt'  using ($6*3600):(-$8*1000) title '10 U NC' w lines lc black lw 3 dt 1, 'enyzme_contribution_0_0_200-enyzme_contribution_0_10_200.txt'using ($6*3600):(-$8*1000) title '50 U NC' w lines lc rgb "#000000" lw 3 dt 3,'NADconsumptionscan_0_200-NADconsumptionscan_10_200.txt' using ($6*3600):(-$8*1000) title '200 U NC' w lines lc rgb "#000000" lw 3 dt 5, ;

unset title
unset ylabel

unset key
unset ylabel
set label 'B' font 'Arial, 16' at screen 0.52, screen 0.98
set lmargin at screen 0.55
set rmargin at screen 0.98
set bmargin at screen 0.572
set tmargin at screen 0.715
set border 1+2+8
set xtics nomirror
set yrange [0:8]
set label 'concentration of free NAD  in {/Symbol m}M' at screen 0.5,screen 0.57 rotate by 90
set ytics 2



plot 'NADconsumptionscan_0_10.txt' using ($6*3600):($7*1000) title '10 U NC NamPRT only' w lines lc rgb "#00AA00" lw 3 dt 1, 'NADconsumptionscan_10_10.txt' using ($6*3600):($7*1000) title '10 U NC  NamPRT + NNMT' w lines  lc  rgb "#3333FF" lw 3 dt 1,  'enyzme_contribution_0_0_200.txt' using ($6*3600):($7*1000) title '50 U NC NamPRT only' w lines lc rgb "#00AA00" lw 3 dt 3, 'enyzme_contribution_0_10_200.txt' using ($6*3600):($7*1000) title '50 U NC NamPRT + NNMT' w lines  lc  rgb "#3333FF" lw 3 dt 3,  'NADconsumptionscan_0_200.txt' using ($6*3600):($7*1000) title '200 U NC NamPRT only' w lines lc rgb "#00AA00" lw 3 dt 5, 'NADconsumptionscan_10_200.txt' using ($6*3600):($7*1000) title '200 U NC NamPRT + NNMT' w lines  lc  rgb "#3333FF" lw 3 dt 5;





set border 2+4+8
set bmargin at screen 0.725
set tmargin at screen 0.925
unset xtics

unset xlabel
set yrange [10:180]
set ytics 20

set arrow from screen 0.545,screen 0.71 to screen 0.555,screen 0.72 nohead
set arrow from screen 0.545,screen 0.72 to screen 0.555,screen 0.73 nohead
set arrow from screen 0.975,screen 0.71 to screen 0.985,screen 0.72 nohead
set arrow from screen 0.975,screen 0.72 to screen 0.985,screen 0.73 nohead
set title 'NAD concentration'

plot 'NADconsumptionscan_0_10.txt' using ($6*3600):($7*1000) title '10 U NC NamPRT only' w lines lc rgb "#00AA00" lw 3 dt 1, 'NADconsumptionscan_10_10.txt' using ($6*3600):($7*1000) title '10 U NC  NamPRT + NNMT' w lines  lc  rgb "#3333FF" lw 3 dt 1,  'enyzme_contribution_0_0_200.txt' using ($6*3600):($7*1000) title '50 U NC NamPRT only' w lines lc rgb "#00AA00" lw 3 dt 3, 'enyzme_contribution_0_10_200.txt' using ($6*3600):($7*1000) title '50 U NC NamPRT + NNMT' w lines  lc  rgb "#3333FF" lw 3 dt 3,  'NADconsumptionscan_0_200.txt' using ($6*3600):($7*1000) title '200 U NC NamPRT only' w lines lc rgb "#00AA00" lw 3 dt 5, 'NADconsumptionscan_10_200.txt' using ($6*3600):($7*1000) title '200 U NC NamPRT + NNMT' w lines  lc  rgb "#3333FF" lw 3 dt 5;


set lmargin at screen 0.55
set rmargin at screen 0.98
set bmargin at screen 0.072
set tmargin at screen 0.215
set border 1+2+8
set xtics nomirror


set label 'concentration reduction in free NAD  in {/Symbol m}M' at screen 0.5,screen 0.02 rotate by 90
unset ylabel
set size 0.5,0.5
set origin 0.5,0.0
unset key
unset title
set label 'D' font 'Arial, 16' at screen 0.52, screen 0.48
set ytics 2
set yrange [0:10]
set xlabel 'cell devision rate per hour'
plot  'NADconsumptionscan_0_10-NADconsumptionscan_10_10.txt'  using ($6*3600):($7*1000) title '10 U NC' w lines lc black lw 3 dt 1,  'NADconsumptionscan_0_105-NADconsumptionscan_10_105.txt' using ($6*3600):($7*1000) title '105 U NC' w lines lc rgb "#000000" lw 3 dt 3,  'NADconsumptionscan_0_200-NADconsumptionscan_10_200.txt' using ($6*3600):($7*1000) title '200 U NC' w lines lc rgb "#000000" lw 3 dt 5;


unset xtics
unset xlabel
set ytics 20


set border 2+4+8
set bmargin at screen 0.225
set tmargin at screen 0.425

set title 'NAD concentration reduction with NNMT'

set yrange [14:180]

set arrow from screen 0.545,screen 0.21 to screen 0.555,screen 0.22 nohead
set arrow from screen 0.545,screen 0.22 to screen 0.555,screen 0.23 nohead
set arrow from screen 0.975,screen 0.21 to screen 0.985,screen 0.22 nohead
set arrow from screen 0.975,screen 0.22 to screen 0.985,screen 0.23 nohead

plot  'NADconsumptionscan_0_10-NADconsumptionscan_10_10.txt'  using ($6*3600):($7*1000) title '10 U NC' w lines lc black lw 3 dt 1,  'NADconsumptionscan_0_105-NADconsumptionscan_10_105.txt' using ($6*3600):($7*1000) title '105 U NC' w lines lc rgb "#000000" lw 3 dt 3,  'NADconsumptionscan_0_200-NADconsumptionscan_10_200.txt' using ($6*3600):($7*1000) title '200 U NC' w lines lc rgb "#000000" lw 3 dt 5;