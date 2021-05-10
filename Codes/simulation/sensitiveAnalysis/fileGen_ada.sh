
#############################################################################
################## File generation for PhyloBCG simulation      #############
#############################################################################


declare -a par=("lambda" "h" "IGsig2" "IGv0")
declare -a val=(1 2500 1e-3 1e-3) # Default hyper parameter values (lambda, h, IGsig2 (for the tree scale), IGv0 (for the spike variance) )
declare -a lam=(0.5 1 2 3 4 5)
declare -a h=(2000 2500 3000 3500 4000 5000)
declare -a abp=(1e-6 1e-5 1e-4 1e-3 1e-2 1e-1)
declare -a abv=(1e-6 1e-5 1e-4 1e-3 1e-2 1e-1)

for ip in ${par[@]};
do
echo "$ip"
mkdir -p "$ip"
done

sdir="/Users/heecheolchung/Dropbox/Research/TAMU/graphical/p50cdf/sensitiveAnalysis/subs.sh"
rdir="/Users/heecheolchung/Dropbox/Research/TAMU/graphical/p50cdf/sensitiveAnalysis/tree.R"


for im in {1..4}
    do
    where="/Users/heecheolchung/Dropbox/Research/TAMU/graphical/p50cdf/sensitiveAnalysis/${par[im-1]}/"
    cd $where
    for iset in {1..6}
        do
        cp $(ls -1 $sdir) "subs"_"$iset.sh"
        cp $(ls -1 $rdir) "${par[$im-1]}"_"$iset.R"
        rfile="${par[$im-1]}"_"$iset.R"
        sfile="subs"_"$iset.sh"
        sed -i.bak "s/m1.R/${par[$im-1]}"_"$iset.R/" $sfile
        if [[ $im -eq 1 ]]
        then
        sed -i.bak "s/${par[$im-1]} = ${val[$im-1]}/${par[$im-1]} = ${lam[$iset-1]}/" $rfile
        sed -i.bak "s/Gibbs_cdfT/${par[$im-1]}_$iset/" $rfile
        sed -i.bak "s/R_fix/sens_${par[$im-1]}"_"$iset/" $sfile
        sed -i.bak "s/fix./s_${par[$im-1]}"_"$iset./" $sfile
        rm $(ls -1 *.bak)
        elif [[ $im -eq 2 ]]
        then
        sed -i.bak "s/${par[$im-1]} = ${val[$im-1]}/${par[$im-1]} = ${h[$iset-1]}/" $rfile
        sed -i.bak "s/Gibbs_cdfT/${par[$im-1]}_$iset/" $rfile
        sed -i.bak "s/R_fix/sens_${par[$im-1]}"_"$iset/" $sfile
        sed -i.bak "s/fix./s_${par[$im-1]}"_"$iset./" $sfile
        rm $(ls -1 *.bak)
        elif  [[ $im -eq 3 ]]
        then
        sed -i.bak "s/${par[$im-1]} = rep(1e-3,2)/${par[$im-1]} = rep(${abp[$iset-1]},2)/" $rfile
        sed -i.bak "s/Gibbs_cdfT/${par[$im-1]}_$iset/" $rfile
        sed -i.bak "s/R_fix/sens_${par[$im-1]}"_"$iset/" $sfile
        sed -i.bak "s/fix./s_${par[$im-1]}"_"$iset./" $sfile
        rm $(ls -1 *.bak)
        else
        sed -i.bak "s/${par[$im-1]} = rep(1e-3,2)/${par[$im-1]} = rep(${abv[$iset-1]},2)/" $rfile
        sed -i.bak "s/Gibbs_cdfT/${par[$im-1]}_$iset/" $rfile
        sed -i.bak "s/R_fix/sens_${par[$im-1]}"_"$iset/" $sfile
        sed -i.bak "s/fix./s_${par[$im-1]}"_"$iset./" $sfile
        rm $(ls -1 *.bak)
        fi
    done
    cd ..
done


###################
# Job submission  #
###################
declare -a par=("lambda" "h" "IGsig2" "IGv0")
for im in {1..4}
    do
    where="/general/home/hcchung/graphical/p50cdf/sensitiveAnalysis/${par[im-1]}/"
    cd $where
    for h in $(ls -1 *.sh)
        do
        bsub < $h
    done
    cd ..
done



