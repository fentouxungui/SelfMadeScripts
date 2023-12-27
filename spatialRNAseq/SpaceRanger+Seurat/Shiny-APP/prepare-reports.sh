#!/bin/bash
mkdir -p ../Analysis-Pipeline-APCmin-SI-LI-DSS-Runx1-QianNN/APCmin_LI-DSS_1M_D 
mkdir -p ../Analysis-Pipeline-APCmin-SI-LI-DSS-Runx1-QianNN/APCmin_LI-DSS_3M_C 
mkdir -p ../Analysis-Pipeline-APCmin-SI-LI-DSS-Runx1-QianNN/APCmin_SI_A 
mkdir -p ../Analysis-Pipeline-APCmin-SI-LI-DSS-Runx1-QianNN/APCmin_SI_B 
ln -s /data2/shared_data_backup/qiannn/SpatialRNAseq/SpatialRNAseq_Mouse/Downstream-Analysis/Single-Sample/APCmin_LI-DSS_1M_D/APCmin_LI-DSS_1M_D.cluster.markers.csv ../Analysis-Pipeline-APCmin-SI-LI-DSS-Runx1-QianNN/APCmin_LI-DSS_1M_D/APCmin_LI-DSS_1M_D.cluster.markers.csv 
ln -s /data2/shared_data_backup/qiannn/SpatialRNAseq/SpatialRNAseq_Mouse/Downstream-Analysis/Single-Sample/APCmin_LI-DSS_1M_D/seurat-spatial-pipeline.html ../Analysis-Pipeline-APCmin-SI-LI-DSS-Runx1-QianNN/APCmin_LI-DSS_1M_D/seurat-spatial-pipeline.html 
ln -s /data2/shared_data_backup/qiannn/SpatialRNAseq/SpatialRNAseq_Mouse/Downstream-Analysis/Single-Sample/APCmin_LI-DSS_3M_C/APCmin_LI-DSS_3M_C.cluster.markers.csv ../Analysis-Pipeline-APCmin-SI-LI-DSS-Runx1-QianNN/APCmin_LI-DSS_3M_C/APCmin_LI-DSS_3M_C.cluster.markers.csv 
ln -s /data2/shared_data_backup/qiannn/SpatialRNAseq/SpatialRNAseq_Mouse/Downstream-Analysis/Single-Sample/APCmin_LI-DSS_3M_C/seurat-spatial-pipeline.html ../Analysis-Pipeline-APCmin-SI-LI-DSS-Runx1-QianNN/APCmin_LI-DSS_3M_C/seurat-spatial-pipeline.html 
ln -s /data2/shared_data_backup/qiannn/SpatialRNAseq/SpatialRNAseq_Mouse/Downstream-Analysis/Single-Sample/APCmin_SI_A/APCmin_SI_A.cluster.markers.csv ../Analysis-Pipeline-APCmin-SI-LI-DSS-Runx1-QianNN/APCmin_SI_A/APCmin_SI_A.cluster.markers.csv 
ln -s /data2/shared_data_backup/qiannn/SpatialRNAseq/SpatialRNAseq_Mouse/Downstream-Analysis/Single-Sample/APCmin_SI_A/seurat-spatial-pipeline.html ../Analysis-Pipeline-APCmin-SI-LI-DSS-Runx1-QianNN/APCmin_SI_A/seurat-spatial-pipeline.html 
ln -s /data2/shared_data_backup/qiannn/SpatialRNAseq/SpatialRNAseq_Mouse/Downstream-Analysis/Single-Sample/APCmin_SI_B/APCmin_SI_B.cluster.markers.csv ../Analysis-Pipeline-APCmin-SI-LI-DSS-Runx1-QianNN/APCmin_SI_B/APCmin_SI_B.cluster.markers.csv 
ln -s /data2/shared_data_backup/qiannn/SpatialRNAseq/SpatialRNAseq_Mouse/Downstream-Analysis/Single-Sample/APCmin_SI_B/seurat-spatial-pipeline.html ../Analysis-Pipeline-APCmin-SI-LI-DSS-Runx1-QianNN/APCmin_SI_B/seurat-spatial-pipeline.html 
for file in $(find ../Analysis-Pipeline-APCmin-SI-LI-DSS-Runx1-QianNN -type l); do if [ ! -e $file ]; then echo "rm $file"; rm -f $file; fi; done