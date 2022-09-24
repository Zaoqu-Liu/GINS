# Classfy new data

load('data.rda') ## here is an example data, you shold input your data here.

load('Classifier.rda') ## Classifier

GINS = predict(rf, data)

table(GINS)

##########################################################################################################################################################

# GINS
Gene interaction perturbation network deciphers a high-resolution taxonomy in colorectal cancer

## Background network
The gene interaction-perturbation pipeline applies the discovery cohort as the tumor sample input and the GTEx cohort as the normal sample input. A pathway-derived analysis requires constructing the protein interaction functional networks projected by candidate genes(Chen et al., 2021). Data from the Reactome Pathway Database (https://reactome.org) focus on the reaction, which were utilized to establish the biological interaction networks(Jassal et al., 2020). Network nodes that absented in our cohorts were removed and existing nodes were integrated into the large background network with 148942 interactions in total.

## Interaction-perturbation-based program
Individual lesion extent, reflecting the biological state of patients, can be assessed through the perturbation degree of the background network that inherently derived from altered gene interactions(Chen et al., 2021). The interaction perturbation within the network can quantify the interaction change for each gene pair. Thus, the overall perturbation of all gene pairs in the background network is reasonably and effectively utilized to characterize the pathological condition at the individual level. Notably, gene interactions are highly conservative within normal samples, and interaction perturbations are rare(Sahni et al., 2015). Hence, we hypothesized that the background network is stable across all normal samples, and then the interactions within normal samples served as the benchmark network. Gene interactions in each patient should be compared with the benchmark network, and the corresponding difference represents the gene interaction perturbation for that patient. Collectively, the interaction-perturbation-based program includes three steps (Figure): (1) convert the gene expression matrix into the gene expression rank matrix; (2) calculate the delta rank matrix; (3) further generate the interaction-perturbation matrix that can assess the gene interaction perturbations for each patient. As an effective measurement of the sample-specific gene interaction perturbation, the interaction-perturbation matrix is used for subsequent clustering analysis.

<img width="477" alt="image" src="https://user-images.githubusercontent.com/68080738/159143569-c9b6b7d9-2634-4d02-8f6d-6020fc72d92e.png">

## Discovery of the gene interaction network-based subtypes
The gene interaction network-based subtypes (GINS) were deciphered in the discovery cohort using the consensus clustering approach that required two elements, including features in rows and samples in columns. The clustering features needed to meet two aspects: (1) being able to significantly distinguish tumor from normal samples; (2) maintaining strong heterogeneity within all tumor samples with high variability. First, the differential analysis between tumor and normal samples for each interaction was performed. The top 30,000 interactions (approximately 20%) that differed significantly between tumor and normal samples and the top 30,000 interactions with high standard deviation (SD) among all CRC samples were selected. Subsequently, the intersection of two sets with 30,000 gene interactions over all CRC samples were retained for clustering work. Furthermore, the perturbation matrix with the selected interactions was subjected to the consensus clustering procedure(Wilkerson and Hayes, 2010) using the partitioning around medoids approach and 1-Spearman correlation distance. The procedure was performed 5000 iterations on 80% of interactions and 80% of samples selected randomly. All derived partitions for a cluster K (2-10) were summarized by clustering the (samples × samples) co-classification matrix. The optimal cluster K was synthetically determined by the consensus score matrix, cumulative distribution function (CDF) curve, proportion of ambiguous clustering (PAC) score, and Nbclust(Senbabaoglu et al., 2014).

## Generation of the GINS classifier
To identify GINS subtypes in novel datasets using a small list of genes, we developed a gene expression-based classifier in three steps. First, our analysis included only samples that statistically belonged to the core of each subtype. Excluding samples with negative silhouette width has been shown to minimize the impact of sample outliers on the identification of subtype markers, as described in TCGA glioblastoma classification(Verhaak et al., 2010). Second, the significance analysis of microarrays (SAM)(Tusher et al., 2001) was employed to identify 762 genes significantly differentially expressed across the GINS subtypes. The threshold was set to Benjamini-Hochberg-corrected false discovery rate (FDR) <0.0001. Third, the differentially expressed genes were further trained by prediction analysis for microarrays (PAM)(Tibshirani et al., 2002) to build a classifier. The leave-one-out cross-validation (LOOCV) was performed to ultimately identify 289 gene that had the lowest misclassification error (1.8%). Using this strategy, a centroid-based classifier was built, and six centroids were computed on the gene-median-centered discovery cohort. For each validation dataset, the distance to the six centroids of each sample was computed and samples were assigned to the closest centroid subtype.

## Statistical analysis
All data processing, statistical analysis, and plotting were conducted in R 4.0.5 software. All statistical tests were two-sided. P < 0.05 was regarded as statistically significant.

## References
Chen, Y., Gu, Y., Hu, Z., and Sun, X. (2021). Sample-specific perturbation of gene interactions identifies breast cancer subtypes. Brief Bioinform 22(4). doi: 10.1093/bib/bbaa268.

Jassal, B., Matthews, L., Viteri, G., Gong, C., Lorente, P., Fabregat, A., et al. (2020). The reactome pathway knowledgebase. Nucleic Acids Res 48(D1), D498-D503. doi: 10.1093/nar/gkz1031.

Sahni, N., Yi, S., Taipale, M., Fuxman Bass, J.I., Coulombe-Huntington, J., Yang, F., et al. (2015). Widespread macromolecular interaction perturbations in human genetic disorders. Cell 161(3), 647-660. doi: 10.1016/j.cell.2015.04.013.

Senbabaoglu, Y., Michailidis, G., and Li, J.Z. (2014). Critical limitations of consensus clustering in class discovery. Sci Rep 4, 6207. doi: 10.1038/srep06207.

Tibshirani, R., Hastie, T., Narasimhan, B., and Chu, G. (2002). Diagnosis of multiple cancer types by shrunken centroids of gene expression. Proc Natl Acad Sci U S A 99(10), 6567-6572. doi: 10.1073/pnas.082099299.

Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance analysis of microarrays applied to the ionizing radiation response. Proc Natl Acad Sci U S A 98(9), 5116-5121. doi: 10.1073/pnas.091062498.

Verhaak, R.G., Hoadley, K.A., Purdom, E., Wang, V., Qi, Y., Wilkerson, M.D., et al. (2010). Integrated genomic analysis identifies clinically relevant subtypes of glioblastoma characterized by abnormalities in PDGFRA, IDH1, EGFR, and NF1. Cancer Cell 17(1), 98-110. doi: 10.1016/j.ccr.2009.12.020.

Wilkerson, M.D., and Hayes, D.N. (2010). ConsensusClusterPlus: a class discovery tool with confidence assessments and item tracking. Bioinformatics 26(12), 1572-1573. doi: 10.1093/bioinformatics/btq170.


R scipt "Perturbation.R" is from https://github.com/Marscolono/SSPGI, which was published in Brief in bioinformatics (PMID: 33126248).

#########################################################################################################################################################
