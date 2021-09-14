#!/bin/bash
# run-cmp30-17.sh -- tab:cmp-on-OncoNEM-simul

source init.sh
DIR=result/TEST-CMP30-17P
OUTPUT=run-cmp30-17
LEN=10

rm -rf $DIR
mkdir -p $DIR

SAMPLE_OPT="50 100 50"
LOC_OPT="20 50 30"
MV_OPT="0.07 0.071 1"
ZO_OPT="0.015 0.0151 0.2"
OZ_OPT="0.10 0.101 0.1"
RUN2=,onco,steiner,bitphyl,scite,sasc,sciphi,sifit,siclonefit,

#$SCRIPT/test-gen-eval-all.sh $DIR -rep $LEN -mut 2 2.1 1 -ozr $OZ_OPT -zor $ZO_OPT -mvr $MV_OPT -locus $LOC_OPT -sample $SAMPLE_OPT -clone 5 10 5 -step 5 5 1 -gen onco -run2 $RUN2 -keep &>$DIR/$OUTPUT.eo
#for REP in `seq 1 $LEN`; do 
#for RUN2_ELEMENT in onco steiner bitphyl scite sasc sciphi sifit siclonefit ; do
#echo RUNNIG $RUN2_ELEMENT $REP
#sem -j 20 $SCRIPT/test-gen-eval-all.sh $DIR -rep 1 -mut 2 2.1 1 -ozr $OZ_OPT -zor $ZO_OPT -mvr $MV_OPT -locus $LOC_OPT -sample $SAMPLE_OPT -clone 5 10 5 -step 5 5 1 -gen onco -run2 ,$RUN2_ELEMENT, -keep -seed $((77777*REP)) &>$DIR/$OUTPUT-$RUN2_ELEMENT-$REP.eo
#done
#done
#sem --wait
#
#for RUN2_ELEMENT in onco steiner bitphyl scite sasc sciphi sifit siclonefit ; do
#for REP in `seq 1 $LEN`; do 
#cat $DIR/$OUTPUT-$RUN2_ELEMENT-$REP.eo
#done
#done > $DIR/$OUTPUT.eo

for REP in `seq 1 $LEN`; do 
#for REP in 2 3 4 6 7 8 9 10 ; do 
echo $DIR/$REP folder created
echo sem -j $LEN $SCRIPT/test-gen-eval-all.sh $DIR/$REP -rep 1 -mut 2 2.1 1 -ozr $OZ_OPT -zor $ZO_OPT -mvr $MV_OPT -locus $LOC_OPT -sample $SAMPLE_OPT -clone 5 10 5 -step 5 5 1 -gen onco -run2 $RUN2 -keep -seed $((77777*REP)) 
mkdir $DIR/$REP
sem -j $LEN $SCRIPT/test-gen-eval-all.sh $DIR/$REP -rep 1 -mut 2 2.1 1 -ozr $OZ_OPT -zor $ZO_OPT -mvr $MV_OPT -locus $LOC_OPT -sample $SAMPLE_OPT -clone 5 10 5 -step 5 5 1 -gen onco -run2 $RUN2 -keep -seed $((77777*REP)) -timeout 10h 36000 &>$DIR/$REP/$OUTPUT-$REP.eo
done
sem --wait

for REP in `seq 1 $LEN`; do 
cat $DIR/$REP/$OUTPUT-$REP.eo
done > $DIR/$OUTPUT.eo

cat $DIR/$OUTPUT.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s O?" "\\s I?" "\\s B?" "\\s S?" "\\s sasc?" "\\s sciphi?" "\\s sifit?" "\\s siclonefit?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pO?" "\\s pI?" "\\s pB?" "\\s pS?" "\\s psasc?" "\\s psciphi?" "\\s psifit?" "\\s psiclonefit?" --sg1 '^([0-9.]*)user[ \t]*(\S+)' donco? dimpt? dbphl? dscte? dsasc? dsciphi? dsifit? dsiclonefit?  --sg1-fault RUN_END | python $SCRIPT/eq-to-csv-2.py -mean-na nan --len $LEN -a O? onco vonco -a pO? ponco vponco -a B? bphl vbphl -a pB? pbphl vpbphl -a S? scte vscte -a pS? pscte vpscte -a I? impt vimpt -a pI? pimpt vpimpt -a sasc? sasc vsasc -a psasc? psasc vpsasc -a sciphi? sciphi vsciphi -a psciphi? psciphi vpsciphi -a sifit? sifit vsifit -a psifit? psifit vpsifit -a siclonefit? siclonefit vsiclonefit -a psiclonefit? psiclonefit vpsiclonefit -a donco? donco vdonco -a dimpt? dimpt vdimpt -a dbphl? dbphl vdbphl -a dscte? dscte vdscte -a dsasc? dsasc vdsasc -a dsciphi? dsciphi vdsciphi -a dsifit? dsifit vdsifit -a dsiclonefit? dsiclonefit vdsiclonefit > $DIR/$OUTPUT.txt

cat $DIR/$OUTPUT.txt | python $SCRIPT/calculate-average.py > $DIR/$OUTPUT-avg.txt

echo "Done"
