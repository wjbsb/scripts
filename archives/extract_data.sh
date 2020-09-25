for i in `ls ~/TRAVAIL/DATA/INMG_SingleCell/data/*.tar.gz`
do
	echo $i
	tar -zxvf $i
done