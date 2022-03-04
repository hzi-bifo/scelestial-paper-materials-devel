DIR=$1
ITER=$2
shift 2

OPT=

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
	-opt)
		OPT=$2
		shift 2
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



awk '{printf("%d %s\n", NR-1, $0)}' $DIR/input-scite.txt > $DIR/input-sifit.txt

if [ "$ITER" == "0" ]; then
OPT_ITER=""
else
OPT_ITER="-iter $ITER"
fi
time java $OPT -jar $SCRIPT/SiFit/SiFit.jar -n `cat $DIR/input-scite.txt | wc -l` -m `head -1 $DIR/input-scite.txt | grep -oE [[:digit:]] | wc -l` $OPT_ITER -df 0 -ipMat $DIR/input-sifit.txt > $DIR/sifit.txt

python $SCRIPT/sifit-2-tree-clone.py $DIR/sifit.txt $DIR/input-scite.txt $DIR/sifit-tree.txt $DIR/sifit-clones.txt
python $SCRIPT/tree-clone-2-distance-matrix.py $DIR/sifit-tree.txt $DIR/sifit-clones.txt $MATRIX_DISTANCE_NORMALIZATION > $DIR/sifit-distance-matrix
python $SCRIPT/tree-part-sim.py $DIR/sifit-tree.txt $DIR/sifit-clones.txt > $DIR/sifit-part
#python visualize-tree.py SCITE $DIR/sifit-tree.txt $DIR/sifit-clones.txt $DIR/output-sifit

