# Bioinformatics Portfolio
This portfolio showcases my bioinformatics research projects focused on evolutionary genomics and plant adaptation.
My research work investigates adaptive evolution in chloroplast genes of Begonia species, integrating sequence analysis, phylogenetics, and SNP–environment association to identify genetic mechanisms underlying shade tolerance.

## Project Overview Begonia Research (February 2024 - September 2025)
For my final Bachelor Biotechnology Project, I was tasked to do a positive selection analysis on Begonia DNA Chloroplast. Based on Xiong et al. (2023) study, there are eight genes that are detected positively selected which are: _rpoC1, rpoB, psbE, psbK, petA, rpl22, rpl2_, and _rps12_. I was tasked to do another positive selection analysis on Indonesian Begonia samples taken from Bogor Botanic Gardens (Kebun Raya Bogor) genes (partial _rpoC1_ and _petA_) in comparison with Asian Begonia's that are available in NCBI. 

### Key finding
- SNP 671 A/C and SNP 1679 A/G in rpoC1 is the most significant associated with the shade adaptation on the Adjusted Residuals result. Allele A contributes to semi-shade adaptation, while allele C at site 671 or G at site 1679 contributes to deep shade adaptation.
- Glu 224 and His 560 are significantly linked to Semi-shade, while Arg 560 is significantly linked to deep shade.

## Methods
- Sequence retrieval: NCBI GenBank
- Gene parsing: Biopython
- Alignment: MUSCLE (MEGA11)
- Phylogenetic analysis: Maximum Likelihood (Tamura-Nei model)
- Positive selection: PAML (M7 vs M8)
- SNP association: Fisher Exact Test & Adjusted Residuals (SPSS)

### List of Begonia Species
- 59 Begonia chloroplast genomes from NCBI
  ![List of Accessions](https://github.com/RYD8803/Bioinformatics_Portfolio_Rafael8803/blob/Begonia-Research/List%20of%20Accessions)
- Indonesian samples from Bogor Botanic Gardens
- Species span Asian, African, and Neotropical lineages

![SNP Data Analysis](https://github.com/RYD8803/Bioinformatics_Portfolio_Rafael8803/blob/Begonia-Research/SNP%20data%20analysis.png?raw=true)

## References
Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, et al. 2009. Biopython: freely available python tools for computational molecular biology and bioinformatics. Bioinformatics. 25(11):1422–1423. DOI:10.1093/bioinformatics/btp163.

Tamura K, Stecher G, Kumar S. 2021. MEGA11: molecular evolutionary genetics analysis version 11. Molecular Biology and Evolution. 38(7):3022–3027. DOI:10.1093/molbev/msab120.

Xiong C, Huang Y, Li Z, Wu L, Liu Z, Zhu W, Li J, Xu R, Hong X. 2023. Comparative chloroplast genomics reveals the phylogeny and the adaptive evolution of Begonia in China. BMC Genomics. 24(1). DOI:10.1186/s12864-023-09563-3. https://pubmed.ncbi.nlm.nih.gov/37891463/ 

Yang Z. 2007. PAML 4: Phylogenetic analysis by maximum likelihood. Molecular Biology and Evolution. 24(8):1586–1591. DOI:10.1093/molbev/msm088. https://academic.oup.com/mbe/article/24/8/1586/1103731 

Yang Z, Wong WSW, Nielsen R. 2005. Bayes empirical Bayes inference of amino acid sites under positive selection. Molecular Biology and Evolution. 22(4):1107–1118. DOI:10.1093/molbev/msi097. https://academic.oup.com/mbe/article/22/4/1107/1083468 
