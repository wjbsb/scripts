$b = "/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/Chipseq/GSM1090108/GSM1090108_2010-4428.ip.CT_density.scaled.wig.bedGraph"
cat $b | awk '{print $1 "\t" $2 "\t" $3}' > "/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/Chipseq/GSM1090108/GSM1090108CT.bed"
sort -k 3n,3 "/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/Chipseq/GSM1090108/GSM1090108CT.bed" > "/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/Chipseq/GSM1090108/GSM1090108CT_vsort.bed"





cat "/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/Chipseq/GSE44756/GSM1090108_2010-4428.ip.CT_density.scaled.wig.bedGraph" | awk '{ if ($4 -gt 1500) print $1 "\t" $2 "\t" $3 "\t" $4;}' > "/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/Chipseq/GSE52939/GSM1090108CT.bed"
cat "/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/Chipseq/GSE44756/GSM1090107_2011-947.ip.IP_density.scaled.wig.bedGraph" | sed '1,3d' | awk '{ if ($4 -gt 1500) print $1 "\t" $2 "\t" $3 "\t" $4;}' > "/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/Chipseq/GSE44756/GSM1090108CT.bed"






a="/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/Chipseq/GSE44756/GSM1090108_2010-4428.ip.CT_density.scaled.wig.bedGraph" 
outputa="/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/Chipseq/GSE44756/GSM1090108CT.bed"

b="/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/Chipseq/GSE44756/GSM1090107_2011-947.ip.IP_density.scaled.wig.bedGraph"
outputb="/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/Chipseq/GSE44756/GSM1090107IP.bed"

awk '{ if ($4 -ge 60) { print $1 "\t" $2 "\t" $3 "\t" $4 }; };' $a | sed '1,3d' | head > $outputa



awk -F "\t" '{ if ($4 >= 1500) { print $1 "\t" $2 "\t" $3 "\t" $4} };' $a | sed '1,3d' | head > $outputa

 > $outputb
 
 
 
head $a -n 100
cat $a  | sed '1,3d' | awk '($4 >= 60) {print}' | head


cat $a  | sed '1,3d' | awk '{print $4}'  > scoring_peak.txt