#!/bin/bash

cd ~/private/DSC180BFinalProject/
mkdir gwas_data

curl https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010774/GCST010774_buildGRCh37.tsv > gwas_data/essential_hypertension.tsv
<<com
curl https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010773/GCST010773_buildGRCh37.tsv > gwas_data/abdominal_hernia
curl https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010772/GCST010772_buildGRCh37.tsv > gwas_data/hyperlipidemia
curl https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010771/GCST010771_buildGRCh37.tsv > gwas_data/osteoarthrosis
curl https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010770/GCST010770_buildGRCh37.tsv > gwas_data/cardiac_dysrhythmias
curl https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010769/GCST010769_buildGRCh37.tsv > gwas_data/asthma
curl https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010768/GCST010768_buildGRCh37.tsv > gwas_data/cataract
curl https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010767/GCST010767_buildGRCh37.tsv > gwas_data/coronary_atherosclerosis
curl https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010766/GCST010766_buildGRCh37.tsv > gwas_data/type_2_diabetes
curl https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010765/GCST010765_buildGRCh37.tsv > gwas_data/parkinsons_disease
curl https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010764/GCST010764_buildGRCh37.tsv > gwas_data/alzheimers_disease
curl https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010763/GCST010763_buildGRCh37.tsv > gwas_data/schizophrenia
com
