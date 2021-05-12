#cd /Users/heecheolchung/Dropbox/Research/TAMU/graphical/p50cdf/
##########################################################################################
################## File generation for Bayesian copula graphical model       #############
##########################################################################################
####################################################################################
# Ada cluster
declare -a mod=("m1" "m2" "m3" "m4" "m5") 
#m1: Oracle
#m2: Dist
#m3: PhyloBCG
#m4: Spiec-Easi
#m5: SPRING

for im in ${mod[@]};
do
echo "$im"
mkdir -p "$im"
done

# job submission file names
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


# Job submission code 
for h in $(ls -1 *.sh)
do
bsub < $h
done

