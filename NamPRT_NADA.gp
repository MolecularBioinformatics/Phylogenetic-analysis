set terminal postscript eps color enhanced 'Arial' 12
set output 'NamPRT_NADA.eps'
set multiplot

set xlabel 'cell division rate (per h)'

set log x
set format x "%g
set size 0.5,0.5
set origin 0.0,0.5
set xrange [0.000001*3600:0.0001*3600]

set yrange [0:0.6]

set label 'A' font 'Arial, 16' at screen 0.0, screen 0.95
set title 'NAD consumption flux ({/Symbol m}M/s)'
set ylabel 'Flux NAD consumption(nM/s)'
plot 'enyzme_contribution_10_0_0.txt' using ($6*3600):($8*1000) title 'NADA only' w lines lc rgb "#AA0000" lt 1 , 'enyzme_contribution_10_10_0.txt' using ($6*3600):($8*1000) title 'NADA+NNMT' w lines lc rgb "#AA00AA" lt 1, 'enyzme_contribution_10_10_200.txt' using ($6*3600):($8*1000) title 'NADA+NamPRT+NNMT' w lines lc rgb "#000000" lt 1, 'enyzme_contribution_10_0_200.txt' using ($6*3600):($8*1000) title 'NADA+NamPRT' w lines lc rgb "#DDDD00" lt 1, 'enyzme_contribution_0_0_200.txt' using ($6*3600):($8*1000) title 'NamPRT only' w lines lc rgb "#00AA00" lt 1, 'enyzme_contribution_0_10_200.txt' using ($6*3600):($8*1000) title 'NamPRT+NNMT' w lines lc  "#3333FF" lt 1;
set title 'Effect on NAD concentration'


set yrange [-0.001:350]


set ylabel offset 2,0 ' concentration of free NAD ({/Symbol m}M)'


set size 0.5,0.5
set origin 0.5,0.5
set label 'B' font 'Arial,16' at screen 0.5, screen 0.95
#set yrange [-0.001:0.08]



plot 'enyzme_contribution_10_0_0.txt' using ($6*3600):($7*1000) title 'NADA only' w lines lc rgb "#AA0000" lt 1 , 'enyzme_contribution_10_10_0.txt' using ($6*3600):($7*1000) title 'NADA+NNMT' w lines lc rgb "#AA00AA" lt 1, 'enyzme_contribution_10_10_200.txt' using ($6*3600):($7*1000) title 'NADA+NamPRT+NNMT' w lines lc rgb "#000000" lt 1, 'enyzme_contribution_10_0_200.txt' using ($6*3600):($7*1000) title 'NADA+NamPRT' w lines lc rgb "#DDDD00" lt 1, 'enyzme_contribution_0_0_200.txt' using ($6*3600):($7*1000) title 'NamPRT only' w lines lc rgb "#00AA00" lt 1, 'enyzme_contribution_0_10_200.txt' using ($6*3600):($7*1000) title 'NamPRT+NNMT' w lines lc  "#3333FF" lt 1;


set size 0.5,0.5
set origin 0.0,0.0
set xrange [0.000001*3600:0.0001*3600]

set yrange [0:0.22]

set label 'A' font 'Arial, 16' at screen 0.0, screen 0.55
set title 'NAD consumption flux ({/Symbol m}M/s)'
set ylabel 'Flux NAD consumption(nM/s)'
plot 'enyzme_contribution_10_10_0.txt' using ($6*3600):($8*1000) title 'NADA+NNMT' w lines lc rgb "#AA00AA" lt 1, 'enyzme_contribution_10_10_200.txt' using ($6*3600):($8*1000) title 'NADA+NamPRT+NNMT' w lines lc rgb "#000000" lt 1, 'enyzme_contribution_0_0_200.txt' using ($6*3600):($8*1000) title 'NamPRT only' w lines lc rgb "#00AA00" lt 1, 'enyzme_contribution_0_10_200.txt' using ($6*3600):($8*1000) title 'NamPRT+NNMT' w lines lc  "#3333FF" lt 1;
set title 'Effect on NAD concentration'


set yrange [-0.001:40]


set ylabel offset 2,0 ' concentration of free NAD ({/Symbol m}M)'


set size 0.5,0.5
set origin 0.5,0.0
set label 'B' font 'Arial,16' at screen 0.5, screen 0.55
#set yrange [-0.001:0.08]



plot 'enyzme_contribution_10_10_0.txt' using ($6*3600):($7*1000) title 'NADA+NNMT' w lines lc rgb "#AA00AA" lt 1, 'enyzme_contribution_10_10_200.txt' using ($6*3600):($7*1000) title 'NADA+NamPRT+NNMT' w lines lc rgb "#000000" lt 1, 'enyzme_contribution_0_0_200.txt' using ($6*3600):($7*1000) title 'NamPRT only' w lines lc rgb "#00AA00" lt 1, 'enyzme_contribution_0_10_200.txt' using ($6*3600):($7*1000) title 'NamPRT+NNMT' w lines lc  "#3333FF" lt 1;
