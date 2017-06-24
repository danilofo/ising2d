set terminal epslatex
set out 'heatcap.eps'

set   autoscale                        # scale axes automatically
set xr [1.0:5.0]
set yr [0.0:0.005]
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Capacità termica in funzione della temperatura (unità naturali) "

#power_law(x,a,b, left, right) = (x < left || x> right ? 1/0 : abs(x-2.26)**(a)+b )

set xlabel "$T[J/k_B]$"
set ylabel "$C_V[J/k_B^2]$"
set arrow from 2.26, graph 0 to 2.26, graph 1 nohead lc rgb "red"
plot "heatcap.dat" with points pointtype 7 notitle
#power_law(x,-2*10**(-22),0.005,1,2.26),\
#power_law(x,-2*10**(-10),1*10**(-10),2.26,5)
   

set out
