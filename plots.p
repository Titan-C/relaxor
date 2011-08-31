binwidth=0.005
bin(x,width)=width*floor(x/width)
plot "polinfo.dat" using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause 10