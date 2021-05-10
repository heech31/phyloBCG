#cd /Users/heecheolchung/Dropbox/Research/TAMU/graphical/p50cdf/
##########################################################################################
################## File generation for Bayesian copula graphical model       #############
##########################################################################################
####################################################################################
# Ada cluster
declare -a mod=("m1" "m2" "m3" "m4" "m5")

for im in ${mod[@]};
do
echo "$im"
mkdir -p "$im"
done

declare -a subfiles=("subjob1" "subjob2" "subjob3" "subjob4" "subjob5")

for im in {1..5}
do
sdir="/general/home/hcchung/graphical/p50cdf/${subfiles[im-1]}.sh"
rdir="/general/home/hcchung/graphical/p50cdf/${mod[im-1]}.R"
where="/general/home/hcchung/graphical/p50cdf/${mod[im-1]}/"
cd $where
for idata in {1..10}
do
cp $(ls -1 $rdir) "${mod[im-1]}"_"$idata.R"
rfile="${mod[im-1]}"_"$idata.R"
sed -i.bak "s/iseed=1/iseed=$idata/" $rfile
rm $(ls -1 *.bak)

cp $(ls -1 $sdir) "${subfiles[im-1]}"_"$idata.sh"
shfile="${subfiles[im-1]}"_"$idata.sh"
sed -i.bak "s/${mod[im-1]}.R/${mod[im-1]}"_"$idata.R/" $shfile
rm $(ls -1 *.bak)
done
cd ..
done



for h in $(ls -1 *.sh)
do
echo $h
done


for h in $(ls -1 *.sh)
do
bsub < $h
done




bjobs -u hcchung

#bkill -u hcchung

for im in {1..1}
do
echo ${mod[im-1]}
done




##########################################################################################
## Sapelo2
declare -a d1=("d1_1" "d1_2" "d1_3")
declare -a d2=("d2_1" "d2_2" "d2_3" "d2_4")

for id1 in ${d1[@]};
do for id2 in ${d2[@]};
do
echo "$id1"_"$id2"
mkdir -p "$id1"_"$id2"
done
done



sdir="/home/hc32681/copulaLDA/sim_censored/jnt_all/sub.sh"
rdir="/home/hc32681/copulaLDA/sim_censored/jnt_all/sim_JNT.R"

for i1 in {1..3}
    do for i2 in {1..3}
        do where="/home/hc32681/copulaLDA/sim_censored/jnt_all/${d1[i1-1]}"_"${d2[i2-1]}/"
            cd $where
        for irep in {1..50}
        do
        cp $(ls -1 $rdir) "sim_JNT"_"$i1"_"$i2"_"$irep.R"
        rfile="sim_JNT"_"$i1"_"$i2"_"$irep.R"
        sed -i.bak "s/idelta1 = 1/idelta1 = $i1/" $rfile
        sed -i.bak "s/idelta2 = 1/idelta2 = $i2/" $rfile
        sed -i.bak "s/irep = 1/irep = $irep/" $rfile
        rm $(ls -1 *.bak)

        cp $(ls -1 $sdir) "sub"_"$i1"_"$i2"_"$irep.sh"
        shfile="sub"_"$i1"_"$i2"_"$irep.sh"
        sed -i.bak "s/sim_JNT.R/sim_JNT"_"$i1"_"$i2"_"$irep.R/" $shfile
        rm $(ls -1 *.bak)
        done
    cd ..
    done
done




for i1 in {1..3}
    do for i2 in {1..3}
        do where="/home/hc32681/copulaLDA/sim_censored/jnt_all/${d1[i1-1]}"_"${d2[i2-1]}/"
            cd $where
            echo $PWD
            for h in $(ls -1 *.sh)
                do
                qsub $h
            done
            cd ..
    done
done







#----------------------------------------------------------------
for h in $( seq 3238031 3238916)
do
qdel $h
done







#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

## Sapelo2
declare -a s=("s_1" "s_2" "s_3")
declare -a d2=("d2_1" "d2_2" "d2_3" "d2_4")

for is in ${s[@]};
do for id2 in ${d2[@]};
do
echo "$is"_"$id2"
mkdir -p "$is"_"$id2"
done
done



sdir="/home/hc32681/copulaLDA/sim_censored/lda_all/sub.sh"
rdir="/home/hc32681/copulaLDA/sim_censored/lda_all/sim_LDA.R"

for is in {1..3}
do for i2 in {1..4}
do where="/home/hc32681/copulaLDA/sim_censored/lda_all/${s[is-1]}"_"${d2[i2-1]}/"
cd $where
for irep in {1..50}
do
cp $(ls -1 $rdir) "sim_LDA"_"$is"_"$i2"_"$irep.R"
rfile="sim_LDA"_"$is"_"$i2"_"$irep.R"
sed -i.bak "s/icsize = 1/icsize = $is/" $rfile
sed -i.bak "s/idelta2 = 1/idelta2 = $i2/" $rfile
sed -i.bak "s/irep = 1/irep = $irep/" $rfile
rm $(ls -1 *.bak)

cp $(ls -1 $sdir) "sub"_"$is"_"$i2"_"$irep.sh"
shfile="sub"_"$is"_"$i2"_"$irep.sh"
sed -i.bak "s/sim_LDA.R/sim_LDA"_"$is"_"$i2"_"$irep.R/" $shfile
rm $(ls -1 *.bak)
done
cd ..
done
done




for is in {1..3}
    do for i2 in {1..2}
        do where="/home/hc32681/copulaLDA/sim_censored/lda_all/${s[is-1]}"_"${d2[i2-1]}/"
        cd $where
        echo $PWD
        for h in $(ls -1 *.sh)
            do
            qsub $h
        done
        cd ..
    done
done



for is in {1..3}
do for i2 in {3..4}
do where="/home/hc32681/copulaLDA/sim_censored/lda_all/${s[is-1]}"_"${d2[i2-1]}/"
cd $where
echo $PWD
for h in $(ls -1 *.sh)
do
qsub $h
done
cd ..
done
done









#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

## Sapelo2
declare -a s=("d1_1" "d1_2" "d1_3")
declare -a d2=("d2_1" "d2_2" "d2_3" "d2_4")

for is in ${s[@]};
do for id2 in ${d2[@]};
do
echo "$is"_"$id2"
mkdir -p "$is"_"$id2"
done
done



sdir="/home/hc32681/copulaLDA/sim_censored/jnt_all/sub.sh"
rdir="/home/hc32681/copulaLDA/sim_censored/jnt_all/sim_JNT.R"

for is in {1..3}
do for i2 in {1..4}
do where="/home/hc32681/copulaLDA/sim_censored/jnt_all/${s[is-1]}"_"${d2[i2-1]}/"
cd $where
for irep in {1..50}
do
cp $(ls -1 $rdir) "sim_JNT"_"$is"_"$i2"_"$irep.R"
rfile="sim_JNT"_"$is"_"$i2"_"$irep.R"
sed -i.bak "s/idelta1 = 1/idelta1 = $is/" $rfile
sed -i.bak "s/idelta2 = 1/idelta2 = $i2/" $rfile
sed -i.bak "s/irep = 1/irep = $irep/" $rfile
rm $(ls -1 *.bak)

cp $(ls -1 $sdir) "sub"_"$is"_"$i2"_"$irep.sh"
shfile="sub"_"$is"_"$i2"_"$irep.sh"
sed -i.bak "s/sim_JNT.R/sim_JNT"_"$is"_"$i2"_"$irep.R/" $shfile
rm $(ls -1 *.bak)
done
cd ..
done
done




for is in {1..3}
do for i2 in {1..2}
do where="/home/hc32681/copulaLDA/sim_censored/jnt_all/${s[is-1]}"_"${d2[i2-1]}/"
cd $where
echo $PWD
for h in $(ls -1 *.sh)
do
qsub $h
done
cd ..
done
done



for is in {1..3}
do for i2 in {3..4}
do where="/home/hc32681/copulaLDA/sim_censored/jnt_all/${s[is-1]}"_"${d2[i2-1]}/"
cd $where
echo $PWD
for h in $(ls -1 *.sh)
do
qsub $h
done
cd ..
done
done




########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################












## Sapelo2
declare -a s=("d1_1" "d1_2" "d1_3")
declare -a d2=("d2_1" "d2_2" "d2_3" "d2_4")

for is in ${s[@]};
do for id2 in ${d2[@]};
do
echo "$is"_"$id2"
mkdir -p "$is"_"$id2"
done
done



sdir="/home/hc32681/copulaLDA/sim_censored/jnt_par/subp.sh"
rdir="/home/hc32681/copulaLDA/sim_censored/jnt_par/sim_JNTp.R"

for is in {1..3}
do for i2 in {1..4}
do where="/home/hc32681/copulaLDA/sim_censored/jnt_par/${s[is-1]}"_"${d2[i2-1]}/"
cd $where

cp $(ls -1 $rdir) "sim_JNTp"_"$is"_"$i2"_".R"
rfile="sim_JNTp"_"$is"_"$i2"_".R"
sed -i.bak "s/idelta1 = 1/idelta1 = $is/" $rfile
sed -i.bak "s/idelta2 = 1/idelta2 = $i2/" $rfile
rm $(ls -1 *.bak)

cp $(ls -1 $sdir) "subp"_"$is"_"$i2"_".sh"
shfile="subp"_"$is"_"$i2"_".sh"
sed -i.bak "s/sim_JNTp.R/sim_JNTp"_"$is"_"$i2"_".R/" $shfile
rm $(ls -1 *.bak)
cd ..
done
done




for is in {1..3}
do for i2 in {1..4}
do where="/home/hc32681/copulaLDA/sim_censored/jnt_par/${s[is-1]}"_"${d2[i2-1]}/"
cd $where
echo $PWD
for h in $(ls -1 *.sh)
do
qsub $h
done
cd ..
done
done




