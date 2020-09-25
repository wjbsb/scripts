########DO NOT WORK !!!!!!!!!!


####PARAMETRES
# chemin racine FUTUR $1
#p=/Users/williamjarassier/TRAVAIL
p=$1
# chemin index aligner FUTUR $2
index=/Users/williamjarassier/TRAVAIL/references/Mus_musculus/Ensembl/GRCm38/Sequence/HISAT2Index/genome
# dossier input
folder=$p/$2
# fastq.gz input
listeread=`ls $folder/DATA/*.fastq.gz`
folderoutput='/Volumes/LaCie/InstitutNeuroMyoGene_DATA/'$2

echo "RACINE " $p
echo "INDEX  " $index
echo "INPUT  " $folder
echo $listeread
echo "OUTPUT " $folderoutput

mkdir -p $folderoutput
mkdir -p $folderoutput/results_exon
mkdir -p $folderoutput/results_transcript
echo "WORKING ON EXONS"

#traitement sur les exons
nextflow run nf-core/rnaseq \
-profile docker \
--genome GRCm38 \
--singleEnd \
--aligner hisat2 \
--hisat2_index $index \
--reads $listeread \
-c $folder/scripts/nfcore_RNA.config \
-w $folderoutput/results_exon \
--email w.jarassier@gmail.com \
--email_on_fail w.jarassier@gmail.com \
--fc_count_type exon

echo "WORKING ON TRANSCRIPTS"

#traitement sur les transcripts
nextflow run nf-core/rnaseq \
-profile docker \
--genome GRCm38 \
--singleEnd \
--aligner hisat2 \
--hisat2_index $index \
--reads '$listeread' \
-c $folder/scripts/config_nfcorepipeline.config \
-w $folderoutput/results_transcript \
--email w.jarassier@gmail.com \
--email_on_fail w.jarassier@gmail.com \
--fc_count_type transcript


