#!/bin/bash

fn=$1;
cfn=$2;
if [ -z "$1" ]; then
    echo scan-igme.sh TPM_file_to_handle
    exit
elif [ ! -e $fn ]; then
    echo scan-igme.sh : error : cannot open $fn
    exit
fi 
if [ ! -e $cfn ]; then
    echo scan-igme.sh : error : cannot open $cfn
    exit
fi

if [ -z "$2" ]||[ -z "$cfn" ]; then
    cfn=$fn
fi

len=`cat $fn|wc -l`; length=`cat $cfn|wc -l`;
dim=`matrics $fn dim|head -n1`

peq_items=`seq 1 $dim|tr '\n' ' '`
peq=`tail -n1 $fn | matrics $dim con vl | matrics $dim con t | matrics $dim con print --"$peq_items" | awk 'NF>1{sum=0;for(i=1;i<=NF;i++){sum+=$i};for(i=1;i<=NF;i++)printf(" %.15g",$i/sum)}'`

printf "#%4s %4s " "tauK" "tauG"
printf "%18s " "rmse"
for ((i=1;i<dim;i++)); do printf "%12s " "its$i"; done
printf "%10s %10s" "ln(T_hat)" "lnA"
printf "\n"

for ((taug=2;taug<=len;taug++)); do for ((tauk=1;tauk<=taug-2;tauk++)); do
    printf " %4s %4s " $tauk $taug
    igme_result=`bash ./do-igme.sh $fn -b $tauk -e $taug 2>&1`;
    lnT=`echo "$igme_result" | grep ln_T_hat | sed s/ln_T_hat=//g`
    lnA=`echo "$igme_result" | grep ln_A | sed s/ln_A=//g`

    # rmse
    for ((i=1;i<=length;i++)); do echo $i; done | matrics --"$lnT" . con | matrics --"$lnA" + con | matrics $dim 1 exp con | matrics $cfn - con | matrics -D"$peq" . con | matrics $dim con fn | matrics 1 con ^ 2 | matrics 1 con sum | matrics 1 con / $((dim*dim*length)) | matrics 1 1 sqrt con | awk '{printf("%18.12g ",$1)}';

    # ITS
    lnThat=`matrics $dim 1 exp --"$lnT" -dd | awk 'NF>1{sum=0;for(i=1;i<=NF;i++){if($i<0)$i=0;sum+=$i};for(i=1;i<=NF;i++)printf(" %.15g",$i/sum)}NF<=1{printf("\n")}'`;
    matrics --"$lnThat" ev | matrics $dim con diag | awk '{for(i=2;i<=NF;i++)printf("%12.9g ",-2/log($i^2))}';

    # matrix
    matrics $dim --"$lnT" -%.15g, | awk '{printf("%s ",substr($0,1,length($0)-1))}' 
    matrics $dim --"$lnA" -%.15g, | awk '{printf("%s ",substr($0,1,length($0)-1))}'

    printf "\n";
done; done;
