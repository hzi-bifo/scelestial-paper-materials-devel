#!/bin/bash
## Mimincs old run-cmp6-all11.sh
# run-cmp30-5e3.sh -- fig:cmp-on-tumor-simul

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

source init.sh

DIR=result/TEST-CMP30-5E3
OUTPUT=run-cmp30-5e3
LEN=10

#rm -rf $DIR
mkdir -p $DIR

SAMPLE_OPT="50 50 1"
LOC_OPT="200 200 1"
MV_OPT="0.07 0.071 1"
ZO_OPT="0.015 0.0151 0.2"
OZ_OPT="0.10 0.101 0.1"
RUN2=,onco,steiner,bitphyl,scite,sasc,sciphi,sifit,siclonefit,

$SCRIPT/test-gen-eval-all.sh $DIR -rep $LEN -mut 0 10.01 1 -ozr $OZ_OPT -zor $ZO_OPT -mvr $MV_OPT -locus $LOC_OPT -sample $SAMPLE_OPT -clone 5 5 5 -step 13 13 1 -gen synth -run2 $RUN2 -keep &>$DIR/$OUTPUT.eo

#script/extract-log.sh "$DIR" "$OUTPUT" $LEN

cat $DIR/$OUTPUT.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s O?" "\\s I?" "\\s B?" "\\s S?" "\\s sasc?" "\\s sciphi?" "\\s sifit?" "\\s siclonefit?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pO?" "\\s pI?" "\\s pB?" "\\s pS?" "\\s psasc?" "\\s psciphi?" "\\s psifit?" "\\s psiclonefit?" --sg1 '^([0-9.]*)user[ \t]*(\S+)' donco? dimpt? dbphl? dscte? dsasc? dsciphi? dsifit? dsiclonefit?  --sg1-fault RUN_END | python $SCRIPT/eq-to-csv-2.py -mean-na nan --len $LEN -a O? onco vonco -a pO? ponco vponco -a B? bphl vbphl -a pB? pbphl vpbphl -a S? scte vscte -a pS? pscte vpscte -a I? impt vimpt -a pI? pimpt vpimpt -a phiscs? phiscs vphiscs -a pphiscs? pphiscs vpphiscs -a sasc? sasc vsasc -a psasc? psasc vpsasc -a sciphi? sciphi vsciphi -a psciphi? psciphi vpsciphi -a sifit? sifit vsifit -a psifit? psifit vpsifit -a siclonefit? siclonefit vsiclonefit -a psiclonefit? psiclonefit vpsiclonefit > $DIR/$OUTPUT.txt

echo "Done"
