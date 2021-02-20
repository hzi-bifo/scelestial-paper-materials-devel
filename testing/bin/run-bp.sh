#!/bin/bash
source `conda info --base`/etc/profile.d/conda.sh

DIR=$1
shift 1

MATRIX_DISTANCE_NORMALIZATION="-no-normal"
while :
do
	case "$1" in
	-m-distance-norm-pair)
		MATRIX_DISTANCE_NORMALIZATION="-normal-mean"
		shift
		;;
	-m-distance-norm-no)
		MATRIX_DISTANCE_NORMALIZATION="-no-normal"
		shift
		;;
	--) # End of all options
		shift
		break;
		;;
	-*)
		echo "Error: Unknown option: $1" >&2
		exit 1
		;;
	*)  # No more options
		break
		;;
	esac
done

conda activate py2
time python $SCRIPT/BitPhylogeny analyse $DIR/input-bp.txt $DIR/temp/bit-output -n 200 -b 10 -t 5 -mode "mutation" -seed 1234
conda deactivate
python $SCRIPT/bitphylogeny-2-tree-clone.py $DIR/temp/bit-output/input-bp.txt/treescripts/nodes-`ls $DIR/temp/bit-output/input-bp.txt/treescripts/nodes-*.graphml | sed 's/^.*-\([0-9]*\)\.graphml$/\1/' | sort -nr | head -n 1`.graphml $DIR/bp-tree.txt $DIR/bp-clones.txt
python $SCRIPT/tree-clone-2-distance-matrix.py $DIR/bp-tree.txt $DIR/bp-clones.txt $MATRIX_DISTANCE_NORMALIZATION > $DIR/bp-distance-matrix
python $SCRIPT/tree-part-sim.py $DIR/bp-tree.txt $DIR/bp-clones.txt > $DIR/bp-part
##python visualize-tree.py BPhyl $DIR/bp-tree.txt $DIR/bp-clones.txt $DIR/output-bp

