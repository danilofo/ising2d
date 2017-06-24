set terminal epslatex
set out 'magnetization.eps'

#set   autoscale                        # scale axes automatically
set xr [1.0:4.0]
set yr [0.0:1.05]
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set xlabel "$T [J/k_b]$"
set ylabel "$|M| [J/(TN)]$"
set title "Magnetizzazione (valore assoluto) in funzione della temperatura (unità natuali)"
set arrow from 2.26, graph 0 to 2.26, graph 1 nohead lc rgb "red"
plot "magnvsT.dat"  with points pointtype 7 notitle

set out
