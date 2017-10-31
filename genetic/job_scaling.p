set term postscript eps enhanced font 'Helvetica,20' linewidth 4
set output 'job_scaling.eps'

set title "Increasing Number of Jobs Workflow Time"

set key left

set auto x

plot 'new.dat' using ($1):($2) with lines lt 3 title "New", 'orig.dat' using ($1):($2) with lines lt 1 title "Orig"
