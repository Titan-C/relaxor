reset;
set size 1.0, 1.0;
set origin 0.0, 0.0;

set multiplot;
set size 1.0, 0.5;
set origin 0.0, 0.5;
set xlabel "Pasos de la simulación";
set ylabel "Energía";
set title "Evolución de la energía del sistema";
plot "energy_log.dat";

set size 0.5, 0.5;
set origin 0.0, 0.0;
set xlabel "Temperatura";
set ylabel "Proporción";
set title "Proporción de polarización congelada";
plot "Congelamiento.dat" using 5:1 w l, \
      "Congelamiento.dat" using 5:2 w l,\
      "Congelamiento.dat" using 5:3 w l,\
      "Congelamiento.dat" using 5:4 w l;

set origin 0.5, 0.0;
set xlabel "Temperatura";
set ylabel "Susceptibilidad";
set title "Susceptibilidad del sistema"
plot "Susceptibilidad.dat" using 5:1 w l,\
      "Susceptibilidad.dat" using 5:2 w l,\
      "Susceptibilidad.dat" using 5:3 w l,\
      "Susceptibilidad.dat" using 5:4 w l;
unset multiplot
pause 10