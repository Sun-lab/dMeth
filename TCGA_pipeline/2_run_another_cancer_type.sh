#!/bin/bash


for type in $( cat 2_cancer_types.txt )
do
     R CMD BATCH --quiet --no-save \
        "--args cancer_type='${type}'" 2_run_another_cancer_type.R 2_run_another_cancer_type_${type}.Rout &
done
