#!/bin/bash

argc=0; show_help=0; success=1;
for args in $*; do argc=$((argc+1)); argv[$argc]=$args; done

begin=0;
end=-1;
fn="";

for ((i=1;i<=argc;i++)); do key=${argv[i]}; key_prefix_char=`echo ${argv[i]}|awk '{print substr($1,1,1)}'`;
  if [ "$key" == "-h" ]||[ "$key" == "-help" ]||[ "$key" == "--help" ]; then
    show_help=1;
  elif [ "$key" == "-f" ]||[ "$key" == "--f" ]; then
    if [ $i -lt $argc ]&&[ "`echo ${argv[i+1]}|awk '{print substr($1,1,1)}'`" != "-" ]; then
        i=$((i+1)); fn=${argv[i]};
    fi
  elif [ "$key" == "-b" ]||[ "$key" == "--b" ]||[ "$key" == "-begin" ]||[ "$key" == "--begin" ]||[ "$key" == "-from" ]||[ "$key" == "--from" ]; then
    if [ $i -lt $argc ]&&[ "`echo ${argv[i+1]}|awk '{print substr($1,1,1)}'`" != "-" ]; then
        i=$((i+1)); begin=${argv[i]};
    fi
  elif [ "$key" == "-e" ]||[ "$key" == "--e" ]||[ "$key" == "-end" ]||[ "$key" == "--end" ]||[ "$key" == "-to" ]||[ "$key" == "--to" ]; then
    if [ $i -lt $argc ]&&[ "`echo ${argv[i+1]}|awk '{print substr($1,1,1)}'`" != "-" ]; then
        i=$((i+1)); end=${argv[i]};
    fi
  else
    if [ "`echo ${argv[i]}|awk '{print substr($1,1,1)}'`" != "-" ]; then
        if [ -z "$fn" ]; then
            fn=${argv[i]};
        elif [ $begin -eq 0 ]; then
            begin=${argv[i]};
        elif [ $end -le 0 ]; then
            end=${argv[i]};
        else
            printf "do-igme.sh : warning : extra parameters ignored : \""${argv[i]}"\"\n";
        fi;
    fi
    #echo "$0 : argv[$i] : unrecognizable parameter ${argv[i]}"
    #success=0;
  fi
done

if [ $show_help -eq 1 ]||[ -z "$fn" ]; then
    echo "do-igme.sh : build Memoryless non-Markov Model"
    echo "            (c) Cao Siqin, Mar 16, 2021"
    echo "   [-f]       TPM file to handle"
    echo "   -b, -begin begin frame (included, default 1)"
    echo "   -e, -end   end frame (included, default not set)"
    exit;
fi
if [ ! -e "$fn" ]; then
    echo do-igme.sh : error : cannot open $fn
    exit;
fi 

dim=`matrics $fn dim|head -n1`;
max_lines=`cat $fn | awk 'END{print NR}'`

if [ $end -lt $begin ]; then end=$max_lines; fi
if [ $begin -lt 1 ]; then begin=1; fi

ts_line=`matrics $dim 1 ln $fn -%.20g\  | awk '{printf("%3d  %s\n",NR,$0)}' | awk -v begin=$begin -v end=$end '{if(NR>=begin&&(end<begin||NR<=end))print $0}'`;
fitting_lines=`for ((i=2;i<=$((dim*dim+1));i++)); do echo "$ts_line" | lsa-linear -x 0 1 -y $i -%28.20g | sed s/+\ //g | sed s/-\ /-/g; done`;
ln_T_hat=`echo "$fitting_lines" | awk '{printf("%s ",$4)}'`;
ln_A=`echo "$fitting_lines" | awk '{printf("%s ",$3)}'`;
T_hat=`echo "$fitting_lines" | awk '{printf("%s ",$4)}' | matrics $dim 1 exp con -%.20g\ `;
A_matrix=`echo "$fitting_lines" | awk '{printf("%s ",$3)}' | matrics $dim 1 exp con -%.20g\ `;


begin_text=$begin; if [ $begin -lt 1 ]; then begin_text="begin"; fi
end_text=$end; if [ $end -lt $begin ]||[ $end -ge $max_lines ]; then end_text="end"; fi
if [ $begin_text == "begin" ]&&[ $end_text == "end" ]; then
    printf "MNM of $fn \n"
else
    printf "MNM of $fn :: $begin_text ~ $end_text \n"
fi
printf "    <R>= %g ± %g\n" `echo "$fitting_lines" | awk '{printf("%s\n",$NF)}' | matrics 1 con average -%f` `echo "$fitting_lines" -%.20g, | awk '{printf("%s\n",$NF)}' | matrics 1 con stdev -%f -%.20g,`;
printf "    RMSE= "; for ((i=1;i<=$end;i++)); do echo $i; done | matrics --"$T_hat" ^ con -%.20g, | matrics --"$A_matrix" . con -%.20g, | matrics $dim con - $fn -%.20g, | awk -v begin=$begin -v end=$end '{for(i=1;i<=NF;i++){if(NR>=begin&&(end<begin||NR<=end)){rmse+=$i^2;n++}; rmse2+=$i^2;n2++; }}END{printf("%g (totally: %g)\n",n>0?sqrt(rmse/n):0,n2>0?sqrt(rmse2/n2):0)}' 
printf "    ITS_bound="; matrics --"$T_hat" ev -%.20g, | matrics $dim con print-diag -%.20g\ | awk '{for(i=2;i<=NF;i++)printf($i<0?" %g+iπ":" %g",-2/log($i^2))}'; printf "\n";
printf "    y^2=ev(A)= "; matrics --"$A_matrix" ev -%.20g, | matrics $dim con print-diag -%.20g\ | awk '{for(i=1;i<=NF;i++)printf(i<NF?"%g ":"%g\n",$i)}' 

printf "  Conclusions:\n"
printf "    ln_T_hat= "; echo $ln_T_hat;
printf "    ln_A= "; echo $ln_A;
#printf "  Approximatedly:\n"
#printf "    T_hat= "; echo $T_hat;
#printf "    A_matrix= "; echo $A_matrix;


