#!/bin/bash
until Rscript sim_run1.R; do
echo Retry1
done

until Rscript sim_run2.R; do 
echo Retry2
done  
