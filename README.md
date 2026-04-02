# Bioinformatics Portfolio
This portfolio showcases my bioinformatics research projects focused on evolutionary genomics and plant adaptation.
My research work investigates adaptive evolution in chloroplast genes of Begonia species, integrating sequence analysis, phylogenetics, and SNP–environment association to identify genetic mechanisms underlying shade tolerance.

## Begonia Research (February 2024 - September 2025)
This project investigates adaptive evolution in chloroplast genes of Begonia species under different light environments (deep shade vs semi-shade). Using comparative genomics and SNP–environment association analysis, this study identifies genetic variants in rpoC1 that are significantly associated with shade adaptation.

### Key finding
- SNP 671 (A/C) and SNP 1679 (A/G) in rpoC1 show the strongest association with shade adaptation based on adjusted residual analysis.  
  - Allele A is associated with semi-shade environments  
  - Allele C (site 671) and G (site 1679) are associated with deep-shade environments  
- Amino acid substitutions:
  - Glu224 and His560 are associated with semi-shade conditions  
  - Arg560 is associated with deep-shade environments  
![SNP Data Analysis](https://github.com/RYD8803/Bioinformatics_Portfolio_Rafael8803/blob/Begonia-Research/SNP%20data%20analysis.png?raw=true)

Mutations at codons 224 and 560 in rpoC1 may influence protein stability and function, potentially affecting plastid-encoded RNA polymerase (PEP) activity.

These changes could alter transcription efficiency under low-light conditions, contributing to adaptive evolution in deep-shade Begonia species.

## Methods
- Sequence retrieval: NCBI GenBank
- Gene parsing: Biopython
- Alignment: MUSCLE (MEGA11)
- Phylogenetic analysis: Maximum Likelihood (Tamura-Nei model)
- Positive selection: PAML (M7 vs M8)
- SNP association: Fisher Exact Test & Adjusted Residuals (SPSS)

### List of Begonia Species
- 59 Begonia chloroplast genomes from NCBI; [See the Full Accessions](https://github.com/RYD8803/Bioinformatics_Portfolio_Rafael8803/blob/Begonia-Research/List%20of%20Accessions)
- Indonesian samples from Bogor Botanic Gardens
- Species span Asian, African, and Neotropical lineages
- 2 Outgroups (only used for phylogenetic tree comparison, not in PAML Analysis)

  ![Phylogeny Tree](https://github.com/RYD8803/Bioinformatics_Portfolio_Rafael8803/blob/Begonia-Research/Phylogeny%20Tree.png?raw=true)


## References
Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, et al. 2009. Biopython: freely available python tools for computational molecular biology and bioinformatics. Bioinformatics. 25(11):1422–1423. DOI:10.1093/bioinformatics/btp163.

Tamura K, Stecher G, Kumar S. 2021. MEGA11: molecular evolutionary genetics analysis version 11. Molecular Biology and Evolution. 38(7):3022–3027. DOI:10.1093/molbev/msab120.

Xiong C, Huang Y, Li Z, Wu L, Liu Z, Zhu W, Li J, Xu R, Hong X. 2023. Comparative chloroplast genomics reveals the phylogeny and the adaptive evolution of Begonia in China. BMC Genomics. 24(1). DOI:10.1186/s12864-023-09563-3. https://pubmed.ncbi.nlm.nih.gov/37891463/ 

Yang Z. 2007. PAML 4: Phylogenetic analysis by maximum likelihood. Molecular Biology and Evolution. 24(8):1586–1591. DOI:10.1093/molbev/msm088. https://academic.oup.com/mbe/article/24/8/1586/1103731 

Yang Z, Wong WSW, Nielsen R. 2005. Bayes empirical Bayes inference of amino acid sites under positive selection. Molecular Biology and Evolution. 22(4):1107–1118. DOI:10.1093/molbev/msi097. https://academic.oup.com/mbe/article/22/4/1107/1083468 
