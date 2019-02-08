set terminal postscript eps color enhanced 'Arial' 10
set output 'NamPRTNNMTkm_woANC.eps'
set multiplot
set format x "%g"
set format y "%g"
set view 60,15
set parametric

set xlabel offset 0,-2  'K_m NamPRT (mM)'
set ylabel offset 1,-0.75 'K_m NNMT (mM)'



#set size 1, 1
set key font "Arial, 10"
set pm3d
set yrange[0.001:1]
set xrange[1e-07:0.001]
#set logscale xy
set xyplan 0
unset key
#unset clabel
set logscale xy
set grid x, y, z
set size 0.5,0.5
set origin 0.0,0.5

set label "*" at 0.000003,0.3,0.15 tc rgb "white" font ",20" front
set title "NAD consumption ({/Symbol m}M/s)"
set label 'A' font 'Arial, 16' at screen 0.0, screen 0.9

splot 'NamPRTNNMT_Km_woANC.txt' using 4:5:($8*1000)  with pm3d;

unset label

set label "*" at 0.000003,0.4,8.4 tc rgb "white" font ",20" front

set size 0.5,0.5
set origin 0.5,0.5
set label 'B' font 'Arial,16' at screen 0.5, screen 0.9
set title "NAD concentration ({/Symbol m}M)"

splot 'NamPRTNNMT_Km_woANC.txt' using 4:5:($7*1000)  with pm3d

