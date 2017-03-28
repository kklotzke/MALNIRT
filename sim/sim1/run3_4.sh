#!/bin/bash
until Rscript sim_run3.R; do
echo Retry3
done

until Rscript sim_run4.R; do 
echo Retry4
done  
