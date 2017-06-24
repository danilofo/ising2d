set terminal epslatex
set out 'chivsT.eps'

set   autoscale                        # scale axes automatically
set xr [1.0:5.0]
set yr [0.0:1]
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Suscettibilità magnetica in funzione della temperatura (unità naturali) "

set xlabel "$T[J/k_B]$"
set ylabel "$\chi[\mu/k_B]$"
set arrow from 2.26, graph 0 to 2.26, graph 1 nohead lc rgb "red"
plot "chivsT.dat" using 1:2 with points pointtype 7 notitle

set out
