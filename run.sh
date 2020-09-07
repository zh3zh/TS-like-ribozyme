DNAFILE1=$1
DNAFILE2=$2
RNAFILE1=$3
RNAFILE2=$4
FASTAFILE=$5
OUTPUTPATH=$6
Weighting_factor=$7




#REQUIRED PROGRAMS
#SeqPrep
#samtools
#unpigz
#gsplit
#pigz


fastaSeq=`tail -1 $FASTAFILE`


echo "merge reads..."
pear -f $2 -r $1 -o LD
pear -f $4 -r $3 -o LR

echo "generate barcode-variant map..."
python GenerateBarcodeMapMP.py --runPath $OUTPUTPATH --refSeq $FASTAFILE
rm $OUTPUTPATH/tmp*
rm $OUTPUTPATH/pseudo*

echo "grep rna reads..."
python ChopRNAReadsMP.py --runPath $OUTPUTPATH --refSeq $FASTAFILE

echo "count variants..."
python CountVariantsMP.py --runPath $OUTPUTPATH --refSeq $FASTAFILE

echo "write multiple sequence alignment..."
python WriteMSA.py --runPath $OUTPUTPATH --refSeq $FASTAFILE --RA 0.5

echo "read relative activity..."
python AnalysisRA.py --runPath $OUTPUTPATH --refSeq $FASTAFILE --count $OUTPUTPATH/var.count >$OUTPUTPATH/var.ra

echo "change activity file format..."
python raInfoToPosRAInfo.py $OUTPUTPATH/var.ra $fastaSeq $OUTPUTPATH/var.pos.ra

echo "calculate ps score..."
python PredictContact.py $OUTPUTPATH/var.ra $OUTPUTPATH/var.pos.ra $fastaSeq $OUTPUTPATH/var.features $OUTPUTPATH/pred.mtx

echo "run monte carlo..."
java -Xmx2000m codaMC/RunPredict $fastaSeq $OUTPUTPATH/pred.mtx $Weighting_factor >$OUTPUTPATH/pred.ss
