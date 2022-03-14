#!/bin/bash
## Mimincs old run-cmp6-all12-4-4
# DO NOT RUN IT DIRECTLY, IT SHOULD BE CALLED THROUGH run-cmp30-5b.sh

CLONECNT=$1
DIR=result/TEST-CMP30-5bs/$CLONECNT
OUTPUT=run-cmp30-5bs$CLONECNT
LEN=1

#rm -rf $DIR
mkdir -p $DIR

SAMPLE_OPT="500 900 100"
#SAMPLE_OPT="50 50 100"
LOC_OPT="100 100 100"
MV_OPT="0.07 0.071 1"
ZO_OPT="0.015 0.0151 0.2"
OZ_OPT="0.10 0.101 0.1"
RUN2=,onco,steiner,bitphyl,scite,sasc,sciphi,sifit,siclonefit,

$SCRIPT/test-gen-eval-all.sh $DIR -rep $LEN -mut 2 2.01 1 -ozr $OZ_OPT -zor $ZO_OPT -mvr $MV_OPT -locus $LOC_OPT -sample $SAMPLE_OPT -clone $CLONECNT $CLONECNT 1 -step 5 5 1 -gen synth -run2 $RUN2 -keep -cpu-limit -timeout 24h 86400 -sifit-iter 0 -siclonefit-iter 0 -exit-on-error-no &>$DIR/$OUTPUT.eo

cat $DIR/$OUTPUT.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s O?" "\\s I?" "\\s B?" "\\s S?" "\\s sasc?" "\\s sciphi?" "\\s sifit?" "\\s siclonefit?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pO?" "\\s pI?" "\\s pB?" "\\s pS?" "\\s psasc?" "\\s psciphi?" "\\s psifit?" "\\s psiclonefit?" --sg1 '^([0-9.]*)user[ \t]*(\S+)' donco? dimpt? dbphl? dscte? dsasc? dsciphi? dsifit? dsiclonefit?  --sg1-fault RUN_END | python $SCRIPT/eq-to-csv-2.py -mean-na nan --len $LEN -a O? onco vonco -a pO? ponco vponco -a B? bphl vbphl -a pB? pbphl vpbphl -a S? scte vscte -a pS? pscte vpscte -a I? impt vimpt -a pI? pimpt vpimpt -a sasc? sasc vsasc -a psasc? psasc vpsasc -a sciphi? sciphi vsciphi -a psciphi? psciphi vpsciphi -a sifit? sifit vsifit -a psifit? psifit vpsifit -a siclonefit? siclonefit vsiclonefit -a psiclonefit? psiclonefit vpsiclonefit -a donco? donco vdonco -a dimpt? dimpt vdimpt -a dbphl? dbphl vdbphl -a dscte? dscte vdscte -a dsasc? dsasc vdsasc -a dsciphi? dsciphi vdsciphi -a dsifit? dsifit vdsifit -a dsiclonefit? dsiclonefit vdsiclonefit > $DIR/$OUTPUT.txt

cat $DIR/$OUTPUT.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s O?" "\\s I?" "\\s B?" "\\s S?" "\\s sasc?" "\\s sciphi?" "\\s sifit?" "\\s siclonefit?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pO?" "\\s pI?" "\\s pB?" "\\s pS?" "\\s psasc?" "\\s psciphi?" "\\s psifit?" "\\s psiclonefit?" --sg1 '^([0-9.]*)user[ \t]*(\S+)' donco? dimpt? dbphl? dscte? dsasc? dsciphi? dsifit? dsiclonefit?  --sg1-fault RUN_END | python $SCRIPT/eq-merge.py --on clone | python $SCRIPT/eq-to-csv-2.py -mean-na nan --len $LEN -a O? onco vonco -a pO? ponco vponco -a B? bphl vbphl -a pB? pbphl vpbphl -a S? scte vscte -a pS? pscte vpscte -a I? impt vimpt -a pI? pimpt vpimpt -a sasc? sasc vsasc -a psasc? psasc vpsasc -a sciphi? sciphi vsciphi -a psciphi? psciphi vpsciphi -a sifit? sifit vsifit -a psifit? psifit vpsifit -a siclonefit? siclonefit vsiclonefit -a psiclonefit? psiclonefit vpsiclonefit -a donco? donco vdonco -a dimpt? dimpt vdimpt -a dbphl? dbphl vdbphl -a dscte? dscte vdscte -a dsasc? dsasc vdsasc -a dsciphi? dsciphi vdsciphi -a dsifit? dsifit vdsifit -a dsiclonefit? dsiclonefit vdsiclonefit > $DIR/$OUTPUT-avg.txt


#cat $DIR/$OUTPUT.eo | python extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s O?" "\\s I?" "\\s B?" "\\s S?" "\\s phiscs?" "\\s sasc?" "\\s sciphi?" "\\s sifit?" "\\s siclonefit?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pO?" "\\s pI?" "\\s pB?" "\\s pS?" "\\s pphiscs?" "\\s psasc?" "\\s psciphi?" "\\s psifit?" "\\s psiclonefit?" --sg1 '^([0-9.]*)user[ \t]*(\S+)' donco? dimpt? dbphl? dscte? dphiscs? dsasc? dsciphi? dsifit? dsiclonefit?  --sg1-fault RUN_END | python eq-to-csv-2.py -mean-na nan --len $LEN -a O? onco vonco -a pO? ponco vponco -a B? bphl vbphl -a pB? pbphl vpbphl -a S? scte vscte -a pS? pscte vpscte -a I? impt vimpt -a pI? pimpt vpimpt -a phiscs? phiscs vphiscs -a pphiscs? pphiscs vpphiscs -a sasc? sasc vsasc -a psasc? psasc vpsasc -a sciphi? sciphi vsciphi -a psciphi? psciphi vpsciphi -a sifit? sifit vsifit -a psifit? psifit vpsifit -a siclonefit? siclonefit vsiclonefit -a psiclonefit? psiclonefit vpsiclonefit -a donco? donco vdonco -a dimpt? dimpt vdimpt -a dbphl? dbphl vdbphl -a dscte? dscte vdscte -a dphiscs? dphiscs vdphiscs -a dsasc? dsasc vdsasc -a dsciphi? dsciphi vdsciphi -a dsifit? dsifit vdsifit -a dsiclonefit? dsiclonefit vdsiclonefit > $DIR/$OUTPUT.txt
#
#cat $DIR/$OUTPUT.eo | python extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s O?" "\\s I?" "\\s B?" "\\s S?" "\\s phiscs?" "\\s sasc?" "\\s sciphi?" "\\s sifit?" "\\s siclonefit?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pO?" "\\s pI?" "\\s pB?" "\\s pS?" "\\s pphiscs?" "\\s psasc?" "\\s psciphi?" "\\s psifit?" "\\s psiclonefit?" --sg1 '^([0-9.]*)user[ \t]*(\S+)' donco? dimpt? dbphl? dscte? dphiscs? dsasc? dsciphi? dsifit? dsiclonefit?  --sg1-fault RUN_END | python script/eq-merge.py --on clone | python eq-to-csv-2.py -mean-na nan --len $LEN -a O? onco vonco -a pO? ponco vponco -a B? bphl vbphl -a pB? pbphl vpbphl -a S? scte vscte -a pS? pscte vpscte -a I? impt vimpt -a pI? pimpt vpimpt -a phiscs? phiscs vphiscs -a pphiscs? pphiscs vpphiscs -a sasc? sasc vsasc -a psasc? psasc vpsasc -a sciphi? sciphi vsciphi -a psciphi? psciphi vpsciphi -a sifit? sifit vsifit -a psifit? psifit vpsifit -a siclonefit? siclonefit vsiclonefit -a psiclonefit? psiclonefit vpsiclonefit -a donco? donco vdonco -a dimpt? dimpt vdimpt -a dbphl? dbphl vdbphl -a dscte? dscte vdscte -a dphiscs? dphiscs vdphiscs -a dsasc? dsasc vdsasc -a dsciphi? dsciphi vdsciphi -a dsifit? dsifit vdsifit -a dsiclonefit? dsiclonefit vdsiclonefit > $DIR/$OUTPUT-avg.txt

#python eq-to-csv.py 		--avg imput `python repeat.py 'I?' 0 $LEN` 		--avg2 onco `python repeat.py 'O?' 0 $LEN` 		--avg3 bphl `python repeat.py 'B?' 0 $LEN` 		--avg4 scte `python repeat.py 'S?' 0 $LEN` 		--var vimput `python repeat.py 'I?' 0 $LEN` 		--var2 vonco `python repeat.py 'O?' 0 $LEN` 		--var3 vbphl `python repeat.py 'B?' 0 $LEN` 		--var4 vscte `python repeat.py 'S?' 0 $LEN` 		--avg5 pimput `python repeat.py 'pI?' 0 $LEN` 		--avg6 ponco `python repeat.py 'pO?' 0 $LEN` 		--avg7 pbphl `python repeat.py 'pB?' 0 $LEN` 		--avg8 pscte `python repeat.py 'pS?' 0 $LEN` 		--var5 vpimput `python repeat.py 'pI?' 0 $LEN` 		--var6 vponco `python repeat.py 'pO?' 0 $LEN` 		--var7 vpbphl `python repeat.py 'pB?' 0 $LEN` 		--var8 vpscte `python repeat.py 'pS?' 0 $LEN` 		> $DIR/$OUTPUT.txt

#script/extract-log.sh "$DIR" "$OUTPUT" $LEN
#cat $DIR/$OUTPUT.eo | 	python extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s O?" "\\s I?" "\\s B?" "\\s S?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pO?" "\\s pI?" "\\s pB?" "\\s pS?" --sg1 user donco? dimpt? dbphl? dscte? | python eq-to-csv.py 		--avg impt `python repeat.py 'I?' 0 $LEN` 		--avg2 onco `python repeat.py 'O?' 0 $LEN` 		--avg3 bphl `python repeat.py 'B?' 0 $LEN` 		--avg4 scte `python repeat.py 'S?' 0 $LEN` 		--var vimpt `python repeat.py 'I?' 0 $LEN` 		--var2 vonco `python repeat.py 'O?' 0 $LEN` 		--var3 vbphl `python repeat.py 'B?' 0 $LEN` 		--var4 vscte `python repeat.py 'S?' 0 $LEN` 		--avg5 pimpt `python repeat.py 'pI?' 0 $LEN` 		--avg6 ponco `python repeat.py 'pO?' 0 $LEN` 		--avg7 pbphl `python repeat.py 'pB?' 0 $LEN` 		--avg8 pscte `python repeat.py 'pS?' 0 $LEN` 		--var5 vpimpt `python repeat.py 'pI?' 0 $LEN` 		--var6 vponco `python repeat.py 'pO?' 0 $LEN` 		--var7 vpbphl `python repeat.py 'pB?' 0 $LEN` 		--var8 vpscte `python repeat.py 'pS?' 0 $LEN` --avg9 donco `python repeat.py 'donco?' 0 $LEN` --avg10 dimpt `python repeat.py 'dimpt?' 0 $LEN` --avg11 dbphl `python repeat.py 'dbphl?' 0 $LEN` --avg12 dscte `python repeat.py 'dscte?' 0 $LEN`	> $DIR/$OUTPUT.txt
#
#LEN=10
#cat $DIR/$OUTPUT.eo | 	python extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s O?" "\\s I?" "\\s B?" "\\s S?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pO?" "\\s pI?" "\\s pB?" "\\s pS?" --sg1 user donco? dimpt? dbphl? dscte? | python script/eq-merge.py --on clone | python eq-to-csv.py 		--avg impt `python repeat.py 'I?' 0 $LEN` 		--avg2 onco `python repeat.py 'O?' 0 $LEN` 		--avg3 bphl `python repeat.py 'B?' 0 $LEN` 		--avg4 scte `python repeat.py 'S?' 0 $LEN` 		--var vimpt `python repeat.py 'I?' 0 $LEN` 		--var2 vonco `python repeat.py 'O?' 0 $LEN` 		--var3 vbphl `python repeat.py 'B?' 0 $LEN` 		--var4 vscte `python repeat.py 'S?' 0 $LEN` 		--avg5 pimpt `python repeat.py 'pI?' 0 $LEN` 		--avg6 ponco `python repeat.py 'pO?' 0 $LEN` 		--avg7 pbphl `python repeat.py 'pB?' 0 $LEN` 		--avg8 pscte `python repeat.py 'pS?' 0 $LEN` 		--var5 vpimpt `python repeat.py 'pI?' 0 $LEN` 		--var6 vponco `python repeat.py 'pO?' 0 $LEN` 		--var7 vpbphl `python repeat.py 'pB?' 0 $LEN` 		--var8 vpscte `python repeat.py 'pS?' 0 $LEN` --avg9 donco `python repeat.py 'donco?' 0 $LEN` --avg10 dimpt `python repeat.py 'dimpt?' 0 $LEN` --avg11 dbphl `python repeat.py 'dbphl?' 0 $LEN` --avg12 dscte `python repeat.py 'dscte?' 0 $LEN`	> $DIR/$OUTPUT-avg.txt
echo "Done"

exit

#for ((i=0; i<=9; i++)); do cat script/run-cmp30-5bsXX.sh | sed "s/XX/$i/g" > script/run-cmp30-5bs$i.sh; done
#for ((i=0; i<=9; i++)); do qsub -l mem_free=64G -l h_vmem=64G /home/hforoughmand/sc-phyl-infer-steiner/test/test-setup-1/script/run-cmp30-5bs$i.sh ; done

OUTPUT=run-cmp30-5b
(for ((i=0; i<=9; i++)); do cat result/TEST-CMP30-5bs$i/run-cmp30-5bs$i.eo | grep -v "Done TEST-GEN-AL"; done; echo "Done TEST-GEN-ALL"; ) > $DIR/$OUTPUT.eo

cat $DIR/$OUTPUT.eo | python extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s O?" "\\s I?" "\\s B?" "\\s S?" "\\s sasc?" "\\s sciphi?" "\\s sifit?" "\\s siclonefit?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pO?" "\\s pI?" "\\s pB?" "\\s pS?" "\\s psasc?" "\\s psciphi?" "\\s psifit?" "\\s psiclonefit?" --sg1 '^([0-9.]*)user[ \t]*(\S+)' donco? dimpt? dbphl? dscte? dsasc? dsciphi? dsifit? dsiclonefit?  --sg1-fault RUN_END | python eq-to-csv-2.py -mean-na nan --len $LEN -a O? onco vonco -a pO? ponco vponco -a B? bphl vbphl -a pB? pbphl vpbphl -a S? scte vscte -a pS? pscte vpscte -a I? impt vimpt -a pI? pimpt vpimpt -a sasc? sasc vsasc -a psasc? psasc vpsasc -a sciphi? sciphi vsciphi -a psciphi? psciphi vpsciphi -a sifit? sifit vsifit -a psifit? psifit vpsifit -a siclonefit? siclonefit vsiclonefit -a psiclonefit? psiclonefit vpsiclonefit -a donco? donco vdonco -a dimpt? dimpt vdimpt -a dbphl? dbphl vdbphl -a dscte? dscte vdscte -a dsasc? dsasc vdsasc -a dsciphi? dsciphi vdsciphi -a dsifit? dsifit vdsifit -a dsiclonefit? dsiclonefit vdsiclonefit > $DIR/$OUTPUT.txt

cat $DIR/$OUTPUT.eo | python extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s O?" "\\s I?" "\\s B?" "\\s S?" "\\s sasc?" "\\s sciphi?" "\\s sifit?" "\\s siclonefit?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pO?" "\\s pI?" "\\s pB?" "\\s pS?" "\\s psasc?" "\\s psciphi?" "\\s psifit?" "\\s psiclonefit?" --sg1 '^([0-9.]*)user[ \t]*(\S+)' donco? dimpt? dbphl? dscte? dsasc? dsciphi? dsifit? dsiclonefit?  --sg1-fault RUN_END | python script/eq-merge.py --on clone | python eq-to-csv-2.py -mean-na nan --len $LEN -a O? onco vonco -a pO? ponco vponco -a B? bphl vbphl -a pB? pbphl vpbphl -a S? scte vscte -a pS? pscte vpscte -a I? impt vimpt -a pI? pimpt vpimpt -a sasc? sasc vsasc -a psasc? psasc vpsasc -a sciphi? sciphi vsciphi -a psciphi? psciphi vpsciphi -a sifit? sifit vsifit -a psifit? psifit vpsifit -a siclonefit? siclonefit vsiclonefit -a psiclonefit? psiclonefit vpsiclonefit -a donco? donco vdonco -a dimpt? dimpt vdimpt -a dbphl? dbphl vdbphl -a dscte? dscte vdscte -a dsasc? dsasc vdsasc -a dsciphi? dsciphi vdsciphi -a dsifit? dsifit vdsifit -a dsiclonefit? dsiclonefit vdsiclonefit > $DIR/$OUTPUT-avg.txt
