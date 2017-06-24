set terminal epslatex
set out 'magnvstime.eps'

set   autoscale # scale axes automatically
set xr [0.0:2500]
set yr [-0.5:1.05] 
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set xlabel "MCS"
set ylabel "$M [J/(TN)]$"
set title "Rilassamento verso l'equilibrio termodinamico"
plot "magnvstime.dat"  with points notitle


set out
