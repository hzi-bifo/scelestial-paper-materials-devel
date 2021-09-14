#!/bin/bash
## Mimincs old run-cmp6-all12-3-5

source init.sh


DIR=result/TEST-CMP30-1-IMPUT
OUTPUT=run-cmp30-1-imput
DIR_OLD=result/TEST-CMP30-1-OS

rm -rf $DIR
mkdir -p $DIR

SAMPLE_OPT="20 100 20"
LOC_OPT="200 200 1"
MV_OPT="0.07 0.071 1"
ZO_OPT="0.015 0.0151 0.2"
OZ_OPT="0.10 0.101 0.1"
#RUN2=,onco,steiner,bitphyl,scite,sasc,sciphi,sifit,siclonefit,
RUN2=,steiner,
#$SCRIPT/test-gen-eval-all.sh $DIR -rep $LEN -mut 10 10.01 4 -ozr $OZ_OPT -zor $ZO_OPT -mvr $MV_OPT -locus $LOC_OPT -sample $SAMPLE_OPT -clone 0 9 1 -step 5 5 1 -gen synth -run2 $RUN2 -keep -exit-on-error-no -cpu-limit &>$DIR/$OUTPUT.eo

SDIR=$DIR
GEN=synth
SEED=0

DATA_FILES=NONE
REPEAT_COUNT=$LEN
STEPMUTATIONRATE_MIN=10
STEPMUTATIONRATE_MAX=10.01
STEPMUTATIONRATE_INC=4
STEPMUTATIONRATE_RANGE=`python -c "import numpy; print(' '.join([str(f) for f in numpy.arange($STEPMUTATIONRATE_MIN, $STEPMUTATIONRATE_MAX, $STEPMUTATIONRATE_INC)]));"`
OZRATE_MIN=0.10
OZRATE_MAX=0.101
OZRATE_INC=0.1
OZRATE_RANGE=`python -c "import numpy; print(' '.join([str(f) for f in numpy.arange($OZRATE_MIN, $OZRATE_MAX, $OZRATE_INC)]));"`
ZORATE_MIN=0.015
ZORATE_MAX=0.0151
ZORATE_INC=0.2
ZORATE_RANGE=`python -c "import numpy; print(' '.join([str(f) for f in numpy.arange($ZORATE_MIN, $ZORATE_MAX, $ZORATE_INC)]));"`
MISSING_VALUE_RATE_MIN=0.07
MISSING_VALUE_RATE_MAX=0.071
MISSING_VALUE_RATE_INC=1
MISSING_VALUE_RATE_RANGE=`python -c "import numpy; print(' '.join([str(f) for f in numpy.arange($MISSING_VALUE_RATE_MIN, $MISSING_VALUE_RATE_MAX, $MISSING_VALUE_RATE_INC)]));"`
LOCUS_MIN=200
LOCUS_MAX=200
LOCUS_INC=1
SAMPLE_MIN=20
SAMPLE_MAX=100
SAMPLE_INC=20
CLONE_MIN=0
CLONE_MAX=9
CLONE_INC=1
STEP_MIN=5
STEP_MAX=5
STEP_INC=1

ALLCNT=0
for DATA_FILE in ${DATA_FILES//,/ } ; do
for ((STEP=STEP_MIN; STEP<=STEP_MAX; STEP+=STEP_INC)); do
for ((CLONE=CLONE_MIN; CLONE<=CLONE_MAX; CLONE+=CLONE_INC)); do
for ((STEP=STEP_MIN; STEP<=STEP_MAX; STEP+=STEP_INC)); do
for ((SAMPLE=SAMPLE_MIN; SAMPLE<=SAMPLE_MAX; SAMPLE+=SAMPLE_INC)); do
for ((LOCUS=LOCUS_MIN; LOCUS<=LOCUS_MAX; LOCUS+=LOCUS_INC)); do
for MISSING_VALUE_RATE in $MISSING_VALUE_RATE_RANGE; do
for ZORATE in $ZORATE_RANGE; do
for OZRATE in $OZRATE_RANGE; do
for STEPMUTATIONRATE in $STEPMUTATIONRATE_RANGE; do
	ALLCNT=$((ALLCNT+1))
	echo "DIR=$SDIR SAMPLE=$SAMPLE LOCUS=$LOCUS STEP=$STEP MVR=$MISSING_VALUE_RATE 0->1=$ZORATE 1->0:$OZRATE STEPMUT:$STEPMUTATIONRATE CLONE:$CLONE ALLCNT:$ALLCNT DIV_RATE=$DIV_RATE DEATH_RATE=$DEATH_RATE "
	for ((CNT=0; CNT<REPEAT_COUNT; CNT++)); do
		DIR_OLDS=$DIR_OLD/$ALLCNT-$CNT/
		echo python bin/calc-imputation-accuracy.py $DIR_OLDS/input-scite.txt $DIR_OLDS/true-seq.txt $DIR_OLDS/impute-seq.txt I
		python bin/calc-imputation-accuracy.py $DIR_OLDS/input-scite.txt $DIR_OLDS/true-seq.txt $DIR_OLDS/impute-seq.txt I
	done
done 
done 
done 
done 
done 
done 
done 
done 
done 
done > $DIR/$OUTPUT.eo

echo "RUN_END" >> $DIR/$OUTPUT.eo

#cat $DIR/$OUTPUT.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^name:I" "\\s I?" "\\s T?" --pg2 "^GHIR? [0-9]+" "\\s pT?" "\\s pI?" --sg1 '^GHIR' dimpt?  --sg1-fault RUN_END | python $SCRIPT/eq-to-csv-2.py --len $LEN -a I? impt vimpt  -a T? true vtrue > $DIR/$OUTPUT.txt
cat $DIR/$OUTPUT.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^name:I" "\\s I?" "\\s T?" --pg2 "^GHIR? [0-9]+" "\\s pT?" "\\s pI?" --sg1 '^GHIR' dimpt? --print RUN_END DIR | python $SCRIPT/eq-to-csv-2.py --len $LEN -a I? impt vimpt  -a T? true vtrue > $DIR/$OUTPUT.txt
#cat $DIR/$OUTPUT.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^name:I" "\\s I?" "\\s T?" --pg2 "^GHIR? [0-9]+" "\\s pT?" "\\s pI?" --sg1 '^GHIR' dimpt?  --sg1-fault RUN_END | python $SCRIPT/eq-merge.py --on clone | python $SCRIPT/eq-to-csv-2.py --len $LEN -a I? impt vimpt  -a T? true vtrue -ignore-missing > $DIR/$OUTPUT-avg.txt
cat $DIR/$OUTPUT.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^name:I" "\\s I?" "\\s T?" --pg2 "^GHIR? [0-9]+" "\\s pT?" "\\s pI?" --sg1 '^GHIR' dimpt? --print RUN_END DIR | python $SCRIPT/eq-merge.py --on clone | python $SCRIPT/eq-to-csv-2.py --len $LEN -a I? impt vimpt  -a T? true vtrue > $DIR/$OUTPUT-avg.txt


##cat $DIR/$OUTPUT.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s O?" "\\s I?" "\\s B?" "\\s S?" "\\s sasc?" "\\s sciphi?" "\\s sifit?" "\\s siclonefit?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pO?" "\\s pI?" "\\s pB?" "\\s pS?" "\\s psasc?" "\\s psciphi?" "\\s psifit?" "\\s psiclonefit?" --sg1 '^([0-9.]*)user[ \t]*(\S+)' donco? dimpt? dbphl? dscte? dsasc? dsciphi? dsifit? dsiclonefit?  --sg1-fault RUN_END | python $SCRIPT/eq-to-csv-2.py -mean-na nan --len $LEN -a O? onco vonco -a pO? ponco vponco -a B? bphl vbphl -a pB? pbphl vpbphl -a S? scte vscte -a pS? pscte vpscte -a I? impt vimpt -a pI? pimpt vpimpt -a sasc? sasc vsasc -a psasc? psasc vpsasc -a sciphi? sciphi vsciphi -a psciphi? psciphi vpsciphi -a sifit? sifit vsifit -a psifit? psifit vpsifit -a siclonefit? siclonefit vsiclonefit -a psiclonefit? psiclonefit vpsiclonefit -a donco? donco vdonco -a dimpt? dimpt vdimpt -a dbphl? dbphl vdbphl -a dscte? dscte vdscte -a dsasc? dsasc vdsasc -a dsciphi? dsciphi vdsciphi -a dsifit? dsifit vdsifit -a dsiclonefit? dsiclonefit vdsiclonefit > $DIR/$OUTPUT.txt
##
##cat $DIR/$OUTPUT.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s O?" "\\s I?" "\\s B?" "\\s S?" "\\s sasc?" "\\s sciphi?" "\\s sifit?" "\\s siclonefit?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pO?" "\\s pI?" "\\s pB?" "\\s pS?" "\\s psasc?" "\\s psciphi?" "\\s psifit?" "\\s psiclonefit?" --sg1 '^([0-9.]*)user[ \t]*(\S+)' donco? dimpt? dbphl? dscte? dsasc? dsciphi? dsifit? dsiclonefit?  --sg1-fault RUN_END | python $SCRIPT/eq-merge.py --on clone | python $SCRIPT/eq-to-csv-2.py -mean-na nan --len $LEN -a O? onco vonco -a pO? ponco vponco -a B? bphl vbphl -a pB? pbphl vpbphl -a S? scte vscte -a pS? pscte vpscte -a I? impt vimpt -a pI? pimpt vpimpt -a sasc? sasc vsasc -a psasc? psasc vpsasc -a sciphi? sciphi vsciphi -a psciphi? psciphi vpsciphi -a sifit? sifit vsifit -a psifit? psifit vpsifit -a siclonefit? siclonefit vsiclonefit -a psiclonefit? psiclonefit vpsiclonefit -a donco? donco vdonco -a dimpt? dimpt vdimpt -a dbphl? dbphl vdbphl -a dscte? dscte vdscte -a dsasc? dsasc vdsasc -a dsciphi? dsciphi vdsciphi -a dsifit? dsifit vdsifit -a dsiclonefit? dsiclonefit vdsiclonefit > $DIR/$OUTPUT-avg.txt
##
#
#cat $DIR/$OUTPUT.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s I?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pI?" --sg1 '^([0-9.]*)user[ \t]*(\S+)' dimpt?  --sg1-fault RUN_END | python $SCRIPT/eq-to-csv-2.py -mean-na nan --len $LEN -a I? impt vimpt -a pI? pimpt vpimpt -a dimpt? dimpt vdimpt | sed 's/\t$//' > $DIR/$OUTPUT-1.txt
#python bin/remove-column.py I0,T0,clone,dimpt,dimpt0,impt,mvr,ozr,pI0,pT0,pimpt,sample,site,step,stepmut,vdimpt,vimpt,vpimpt,zor < $DIRM/$OUTPUTM.txt > $DIR/$OUTPUT-rm.txt
#paste <(sed '$ d' <$DIR/$OUTPUT-1.txt) $DIR/$OUTPUT-rm.txt > $DIR/$OUTPUT.txt
#
#cat $DIR/$OUTPUT.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s I?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pI?" --sg1 '^([0-9.]*)user[ \t]*(\S+)' dimpt?  --sg1-fault RUN_END > $DIR/$OUTPUT-pat.txt
#cat $DIRM/$OUTPUTM.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s O?" "\\s I?" "\\s B?" "\\s S?" "\\s sasc?" "\\s sciphi?" "\\s sifit?" "\\s siclonefit?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pO?" "\\s pI?" "\\s pB?" "\\s pS?" "\\s psasc?" "\\s psciphi?" "\\s psifit?" "\\s psiclonefit?" --sg1 '^([0-9.]*)user[ \t]*(\S+)' donco? dimpt? dbphl? dscte? dsasc? dsciphi? dsifit? dsiclonefit?  --sg1-fault RUN_END > $DIR/$OUTPUT-mpat.txt
#
#LEN=10
#paste -d' ' <(sed '$ d' <$DIR/$OUTPUT-pat.txt) $DIR/$OUTPUT-mpat.txt | python $SCRIPT/eq-merge.py --on clone | python $SCRIPT/eq-to-csv-2.py -mean-na nan --len $LEN -a O? onco vonco -a pO? ponco vponco -a B? bphl vbphl -a pB? pbphl vpbphl -a S? scte vscte -a pS? pscte vpscte -a I? impt vimpt -a pI? pimpt vpimpt -a sasc? sasc vsasc -a psasc? psasc vpsasc -a sciphi? sciphi vsciphi -a psciphi? psciphi vpsciphi -a sifit? sifit vsifit -a psifit? psifit vpsifit -a siclonefit? siclonefit vsiclonefit -a psiclonefit? psiclonefit vpsiclonefit -a donco? donco vdonco -a dimpt? dimpt vdimpt -a dbphl? dbphl vdbphl -a dscte? dscte vdscte -a dsasc? dsasc vdsasc -a dsciphi? dsciphi vdsciphi -a dsifit? dsifit vdsifit -a dsiclonefit? dsiclonefit vdsiclonefit -ignore-missing > $DIR/$OUTPUT-avg.txt




echo "Done"
