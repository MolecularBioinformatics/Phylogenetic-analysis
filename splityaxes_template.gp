
set xlabel 'K_M NamPRT (mM)' 
set size 0.3,0.3
set origin 0.0, 0.0
set ylabel 'NAD-consumption flux (nM/s)' 
set label  'E' font 'Arial,9' at screen 0.0, screen 0.3
set yrange [0:0.16]
unset grid
set key left at screen 0.06,screen 0.22
plot  'NamPRT_Km_0_1e-05.txt' using 4:($8*1000) title 'NamPRT only' w lines lt rgb "#00AA00" lw 3, 'NamPRT_Km_10_1e-05.txt' using 4:($8*1000) title 'NamPRT + NNMT' w lines  lt  rgb "#3333FF" lw 3 ;


set lmargin at screen 0.35
set rmargin at screen 0.6
set bmargin at screen 0.072
set tmargin at screen 0.165
set border 1+2+8
set xtics nomirror
set ytics 2 

unset ylabel
set label  'F' font 'Arial,9' at screen 0.3, screen 0.3


set label 'Concentration of free NAD  in {/Symbol m}M' at screen 0.32,screen 0.085 rotate by 90

set yrange [1:10]
unset key



plot  'NamPRT_Km_0_1e-05.txt' using 4:($7*1000) title 'NamPRT only' w lines lt rgb "#00AA00" lw 3, 'NamPRT_Km_10_1e-05.txt' using 4:($7*1000) title 'NamPRT + NNMT' w lines  lt  rgb "#3333FF" lw 3 ;

set yrange [113:124]

unset xtics
set ytics 4 
unset xlabel
set key left at screen 0.37,screen 0.22
set border 2+4+8
set bmargin at screen 0.175
set tmargin at screen 0.285

set arrow from screen 0.345,screen 0.16 to screen 0.355,screen 0.17 nohead
set arrow from screen 0.345,screen 0.17 to screen 0.355,screen 0.18 nohead
set arrow from screen 0.595,screen 0.16 to screen 0.605,screen 0.17 nohead
set arrow from screen 0.595,screen 0.17 to screen 0.605,screen 0.18 nohead

plot  'NamPRT_Km_0_1e-05.txt' using 4:($7*1000) title 'NamPRT only' w lines lt rgb "#00AA00" lw 3, 'NamPRT_Km_10_1e-05.txt' using 4:($7*1000) title 'NamPRT + NNMT' w lines  lt  rgb "#3333FF" lw 3 ;

