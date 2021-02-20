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


awk '{printf("%d %s\n", NR-1, $0)}' $DIR/input-scite.txt > $DIR/input-siclonefit.txt
head -1 $DIR/input-scite.txt | grep -oE [[:digit:]] | awk '{printf("%d ", NR, $0)}' > $DIR/siclonefit-cellNames.txt

mkdir -p $DIR/siclonefit
echo "time java $OPT -jar $SCRIPT/SiCloneFit/SiCloneFiTComplete.jar -n `cat $DIR/input-scite.txt | wc -l` -m `head -1 $DIR/input-scite.txt | grep -oE [[:digit:]] | wc -l` -iter $ITER -ipMat $DIR/input-siclonefit.txt -cellNames $DIR/siclonefit-cellNames.txt -outDir $DIR/siclonefit"
time java $OPT -jar $SCRIPT/SiCloneFit/SiCloneFiTComplete.jar -n `cat $DIR/input-scite.txt | wc -l` -m `head -1 $DIR/input-scite.txt | grep -oE [[:digit:]] | wc -l` -iter $ITER -ipMat $DIR/input-siclonefit.txt -cellNames $DIR/siclonefit-cellNames.txt -outDir $DIR/siclonefit

python $SCRIPT/siclonefit-2-tree-clone.py $DIR/siclonefit/samples/best/best_MAP_tree.txt $DIR/input-scite.txt $DIR/siclonefit-tree.txt $DIR/siclonefit-clones.txt
python $SCRIPT/tree-clone-2-distance-matrix.py $DIR/siclonefit-tree.txt $DIR/siclonefit-clones.txt $MATRIX_DISTANCE_NORMALIZATION > $DIR/siclonefit-distance-matrix
python $SCRIPT/tree-part-sim.py $DIR/siclonefit-tree.txt $DIR/siclonefit-clones.txt > $DIR/siclonefit-part
#python visualize-tree.py SCITE $DIR/siclonefit-tree.txt $DIR/siclonefit-clones.txt $DIR/output-siclonefit

