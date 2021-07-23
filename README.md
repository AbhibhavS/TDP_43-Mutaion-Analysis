# TDP_43-Mutaion-Analysis
# Unravelling the effects of disease associated mutations in TDP-43 protein via molecular dynamics simulation and machine learning.
Autohor: Abhibhav Sharma; Pinki Dey

The molecular mechanisms leading to the accumulation of TAR DNA-binding protein 43 (TDP-43) in central nervous system is a key feature of several common neurological disorders in ageing societies, such as frontotemporal   (FTD), Alzheimer’s disease (AD), amyotrophic lateral sclerosis (ALS) and limbic predominant age-related TDP-43 encephalopathy (LATE)1. The TDP-43 protein is highly conserved and plays a significant role in RNA regulation such as splicing, transcriptional regulation, mRNA stabilization2, 3 etc. Moreover, the TDP-43 is a ubiquitously expressed member of the large heterogeneous nuclear ribonucleoprotein (hnRNP) family that shows specific RNA/DNA binding ability by the highly conserved RNA recognition motifs (RRMs) of the proteins4. But during pathological conditions, several post-translational modifications occur in the protein that leads to their cytoplasmic aggregation causing TDP-43 proteinopathies5-7. Infact, ~97% of all the cases of ALS, ~75% of patients with severe AD and ~45% of all the cases of FTDL involve the aggregation of TDP-438, 9. And all the four diseases which are together known as TDP-proteinopathies2, 10 constitute the major cause of dementia in the world and are expected to rise notoriously in the coming years11.
The generation of C-terminal fragments, TDP-25 (25kD fragment) and TDP-35 (35kD fragment) by caspase-3 and caspase-7 is the most prominent step for the clearance of TDP-4323. Several experimental studies have also reported significant delay in cell death on blockage of the caspase digestion23-25. Moreover, it is shown that out of the four prominent caspase cleavage consensus sites, three sites lie in the RNA recognition motifs (RRM) of the TDP-43 protein26.  The cleavage at D89-A90 of N-terminal domain (NTD) of TDP-43 generates TDP-35 that is still capable of folding correctly27. But, both the cleavage sites at D169-G170 and D174-C175 of the RRM of TDP-43 generates TDP-25 which lacks the NTD, nuclear localization signal (NLS) and most of RRM1, trapping the protein in the cytoplasm28 and thus enhancing its cytoplasmic aggregation29. Interestingly, certain mutations, particularly D169G in the RRM1 domain is reported to increase the thermal stability of the protein which becomes more accessible to cleavage by caspase-3 resulting in the early onset of diseases such as ALS26. They also showed that the neighbouring I168 residue is also very crucial for protein folding. Recently, mutations in the RRMs are also shown to influence the DNA or RNA binding specificities30 indicating the role of nucleic acid binding in TDP-43 aggregation. However, studies on the effect of disease-causing mutations on the RRM domains of the TDP-43 is still very limited. Given the crucial role played by the disease-causing mutations in the RRM domain in inducing conformational stability, the nucleic acid binding and their role in affecting the cleavage sites remains largely unexplored.
In this paper, we perform an extensive structural analysis on the impact of disease-causing mutations and nucleic acid binding in the RRM1 domain of the TDP-43 protein. We address this question by integrating molecular dynamics and machine learning approaches in mainly two ways (i) structural analysis of the wild type and disease-causing mutant proteins bound to nucleic acid by molecular dynamics simulations and (ii) identifying the functionally important regions and biological descriptors of TDP-43 in different mutated states that are crucial in explaining the effect of nucleic acid binding and disease causing mutations in the caspase cleavage sites in RRM motifs of TDP-43 using Machine Learning models.


The Pipeline

(I) The Simulation dataset is first feed to generate the inter residual distance

(II) The empirical analysis is carried out

(III) The machine learning models were then implemented on the filtered inverse dataset to etract features

(IV) The importance profile is generated

(V) The important intercatons were identified

The individual pair wise analysis could be performed by the sepraret machine learning modules 

