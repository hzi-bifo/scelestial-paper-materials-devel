#!/bin/bash
source init.sh
DIR=result/TEST-CMP30-3bs-OS/
LEN=1
OUTPUT=run-cmp30-3bs-os
DIRM=result/TEST-CMP30-3bs/
OUTPUTM=run-cmp30-3b

for ((i=0; i<=9; i++)); do sem -j 20 --id cmp30-3b ./run-cmp30-3bsXX-only-scelestial.sh $i; done
sem --wait --id cmp30-3b

OUTPUT=run-cmp30-3b-os
(for ((i=0; i<=9; i++)); do cat result/TEST-CMP30-3bs-OS/$i/run-cmp30-3bs-os$i.eo | grep -v "Done TEST-GEN-AL"; done; echo "Done TEST-GEN-ALL"; ) > $DIR/$OUTPUT.eo

LEN=1
cat $DIR/$OUTPUT.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s I?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pI?" --sg1 '^([0-9.]*)user[ \t]*(\S+)' dimpt?  --sg1-fault RUN_END > $DIR/$OUTPUT-pat.txt
cat $DIRM/$OUTPUTM.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s O?" "\\s I?" "\\s B?" "\\s S?" "\\s sasc?" "\\s sciphi?" "\\s sifit?" "\\s siclonefit?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pO?" "\\s pI?" "\\s pB?" "\\s pS?" "\\s psasc?" "\\s psciphi?" "\\s psifit?" "\\s psiclonefit?" --sg1 '^([0-9.]*)user[ \t]*(\S+)' donco? dimpt? dbphl? dscte? dsasc? dsciphi? dsifit? dsiclonefit?  --sg1-fault RUN_END > $DIR/$OUTPUT-mpat.txt  

cat $DIR/$OUTPUT-pat.txt | grep -v 'sample:50 site:170 step:5 mvr:0.07 zor:0.015 ozr:0.1 stepmut:2.0 clone:4' | grep -v 'sample:50 site:190 step:5 mvr:0.07 zor:0.015 ozr:0.1 stepmut:2.0 clone:4' | grep -v 'sample:50 site:170 step:5 mvr:0.07 zor:0.015 ozr:0.1 stepmut:2.0 clone:8' | grep -v 'sample:50 site:190 step:5 mvr:0.07 zor:0.015 ozr:0.1 stepmut:2.0 clone:8' | grep -v 'sample:50 site:150 step:5 mvr:0.07 zor:0.015 ozr:0.1 stepmut:2.0 clone:8' | grep -v 'sample:50 site:150 step:5 mvr:0.07 zor:0.015 ozr:0.1 stepmut:2.0 clone:4' | sed 's/ *$//' > $DIR/$OUTPUT-pat2.txt
cat $DIR/$OUTPUT-mpat.txt | grep -v 'sample:50 site:150 step:5 mvr:0.07 zor:0.015 ozr:0.1 stepmut:2.0 clone:8' | grep -v 'sample:50 site:150 step:5 mvr:0.07 zor:0.015 ozr:0.1 stepmut:2.0 clone:4' | sed 's/ *$//' > $DIR/$OUTPUT-mpat2.txt


LEN=10
paste $DIR/$OUTPUT-pat2.txt $DIR/$OUTPUT-mpat2.txt -d' ' | python $SCRIPT/eq-to-csv-2.py -mean-na nan --len $LEN -a O? onco vonco -a pO? ponco vponco -a B? bphl vbphl -a pB? pbphl vpbphl -a S? scte vscte -a pS? pscte vpscte -a I? impt vimpt -a pI? pimpt vpimpt -a sasc? sasc vsasc -a psasc? psasc vpsasc -a sciphi? sciphi vsciphi -a psciphi? psciphi vpsciphi -a sifit? sifit vsifit -a psifit? psifit vpsifit -a siclonefit? siclonefit vsiclonefit -a psiclonefit? psiclonefit vpsiclonefit -a donco? donco vdonco -a dimpt? dimpt vdimpt -a dbphl? dbphl vdbphl -a dscte? dscte vdscte -a dsasc? dsasc vdsasc -a dsciphi? dsciphi vdsciphi -a dsifit? dsifit vdsifit -a dsiclonefit? dsiclonefit vdsiclonefit -ignore-missing > $DIR/$OUTPUT.txt

#cat $DIR/$OUTPUT.eo | python $SCRIPT/extract-log-results.py --rep $LEN --print "^DIR=" "Done TEST-GEN-ALL" --pat1 "DIR=" "SAMPLE=" sample "LOCUS=" site "STEP=" step "MVR=" mvr "0->1=" zor "1->0:" ozr "STEPMUT:" stepmut "CLONE:" clone --pg1 "^True? [0-9]+" "\\s T?" "\\s O?" "\\s I?" "\\s B?" "\\s S?" "\\s sasc?" "\\s sciphi?" "\\s sifit?" "\\s siclonefit?" --pg2 "^pTrue? [0-9]+" "\\s pT?" "\\s pO?" "\\s pI?" "\\s pB?" "\\s pS?" "\\s psasc?" "\\s psciphi?" "\\s psifit?" "\\s psiclonefit?" --sg1 '^([0-9.]*)user[ \t]*(\S+)' donco? dimpt? dbphl? dscte? dsasc? dsciphi? dsifit? dsiclonefit?  --sg1-fault RUN_END | 
paste $DIR/$OUTPUT-pat2.txt $DIR/$OUTPUT-mpat2.txt -d' ' | python $SCRIPT/eq-merge.py --on clone | python $SCRIPT/eq-to-csv-2.py -mean-na nan --len $LEN -a O? onco vonco -a pO? ponco vponco -a B? bphl vbphl -a pB? pbphl vpbphl -a S? scte vscte -a pS? pscte vpscte -a I? impt vimpt -a pI? pimpt vpimpt -a sasc? sasc vsasc -a psasc? psasc vpsasc -a sciphi? sciphi vsciphi -a psciphi? psciphi vpsciphi -a sifit? sifit vsifit -a psifit? psifit vpsifit -a siclonefit? siclonefit vsiclonefit -a psiclonefit? psiclonefit vpsiclonefit -a donco? donco vdonco -a dimpt? dimpt vdimpt -a dbphl? dbphl vdbphl -a dscte? dscte vdscte -a dsasc? dsasc vdsasc -a dsciphi? dsciphi vdsciphi -a dsifit? dsifit vdsifit -a dsiclonefit? dsiclonefit vdsiclonefit -ignore-missing > $DIR/$OUTPUT-avg.txt
