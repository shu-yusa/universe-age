 set term postscript eps enhanced color
 set output "universe.eps"
 unset key
 set size 0.6,0.6
 set xlabel "t [Gyr]"
 set ylabel "a / a_0"
 set xrange [-20:50]
 set yrange [0:4]
 set label "(K>0,{/Symbol L}=0)" at -2,0.5
 set label "(K=0,{/Symbol L}=0)" at 32,2.5
 set label "(K<0,{/Symbol L}=0)" at 22,3.6
 set label "(K=0,{/Symbol L}>0)" at  4,3.5
 set arrow from 25,3.4 to 27,3
 pl "model1.dat" w l, "model2.dat" w l, "model3.dat" w l, "model4.dat" w l
 set term x11
