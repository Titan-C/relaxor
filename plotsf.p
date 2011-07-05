reset;
set terminal png

set output '2_Dip_cong.png'
set xrange [0:2.5]
set xlabel "Temperatura";
set ylabel "Proporción";
set title "Proporción de polarización congelada";
plot "Dip_cong_Efixcool.dat" using 1:2 w linespoints, \
      "Dip_cong_Efixcool.dat" using 1:3 w linespoints;

set output '3_Sus_E0.png'
set xlabel "Temperatura";
set ylabel "Susceptibilidad";
set title "Susceptibilidad del sistema"
set xrange [0:2.7]
set yrange [0:1]
plot "Sus_Efix_cool.dat" using 1:2 w linespoints,\
      "Sus_Efix_cool.dat" using 1:4,\
      "Sus_Efix_cool.dat" using 1:5,\
      "Sus_Efix_cool.dat" using 1:6;

set output '4_Dip_cong_E.png'
set xlabel "Temperatura";
set ylabel "Proporción";
set title "Proporción de polarización congelada E=0;0.8;2";
set xrange [0:2.2]
set yrange [0:1]
plot "Dip_cong_90_cool.dat" using 1:2 w linespoints,\
      "Dip_cong_90_cool.dat" using 1:8 w linespoints,\
      "Dip_cong_90_cool.dat" using 1:13 w linespoints;

set output '5_Sus_E.png'
set xlabel "Temperatura";
set ylabel "Susceptibilidad";
set title "Susceptibilidad del sistema E=0;0.8;2"
set yrange [0:0.8]
plot "Sus_90_cool.dat" using 1:2 w linespoints,\
      "Sus_90_cool.dat" using 1:8 w linespoints,\
      "Sus_90_cool.dat" using 1:13 w linespoints;

set output '6a_Xmax.png'
set xlabel "Campo";
set ylabel "Xmax";
set title "Susceptibilidad maxima"
set yrange [0.1:0.8]
plot "Xmax_Tmaxcool.dat" using 1:2 w linespoints;

set output '6b_Tmax.png'
set xlabel "Campo";
set ylabel "Tmax";
set title "Temperatura del maximo"
set yrange [0:5]
plot "Xmax_Tmaxcool.dat" using 1:3 w linespoints;

set output '7_Polarizacion.png'
set xlabel "Temperatura";
set ylabel "polarizacion";
set title "polarización material calentando E=0";
set xrange [0:3]
set yrange [0:0.35]
plot "polarizacion_heatpol.dat" using 1:2,\
      "polarizacion_heatunpol.dat" using 1:2;

set output '8a_Polarizacion.png'
set xlabel "Temperatura";
set ylabel "polarizacion";
set title "polarización material calentando E=0.2";
set xrange [0:3]
set yrange [0:0.35]
plot "polarizacion_heatpol.dat" using 1:3,\
      "polarizacion_heatunpol.dat" using 1:3;

set output '8b_Polarizacion.png'
set xlabel "Temperatura";
set ylabel "polarizacion";
set title "polarización material calentando E=0.6";
set xrange [0:3]
set yrange [0.2:0.55]
plot "polarizacion_heatpol.dat" using 1:4,\
      "polarizacion_heatunpol.dat" using 1:4;

set output '9_Polarizacion.png'
set xlabel "Temperatura";
set ylabel "polarizacion";
set title "polarización material enfriando E=0;0.2;0.6";
set xrange [0:3]
set yrange [0:0.55]
plot "polarizacion_cool.dat" using 1:2,\
      "polarizacion_cool.dat" using 1:4,\
      "polarizacion_cool.dat" using 1:7;

set output '10_Polarizacion.png'
set xlabel "Temperatura";
set ylabel "polarizacion";
set title "polarización material enfriando E=0";
set xrange [0:1.2]
set yrange [0:0.1]
plot "polarizacion_coolpol.dat",\
      "polarizacion_coolpol_t.dat",\
      "polarizacion_coolpol_t_t.dat",\
      "polarizacion_coolpol_t_t_t.dat";
