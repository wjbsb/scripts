#######download SRR

accession_number=$1
pathoutput="/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/Chipseq/"
#GSE52939
fastq-dump SRR1041910 -O $pathoutput
fastq-dump SRR1041911 -O $pathoutput
#GSE44756
fastq-dump SRR771227 -O $pathoutput
fastq-dump SRR771228 -O $pathoutput
#GSE44755 == GSE44756
#GSE11062 = chipchip

nextflow run nf-core/chipseq \
--genome GRCm38 \
-profile docker \
--single_end \
--fasta \
--gtf \
--bwa_index \
--input /Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/Chipseq/GSE52939/design_chipseq.csv



