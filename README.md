# Definition of a genetic relatedness cutoff to exclude recent transmission of meticillin-resistant _Staphylococcus aureus_: a genomic epidemiology analysis

This GitHub project contains the data and code necessary to reproduce the findings of the study 'Definition of a genetic relatedness cutoff to exclude recent transmission of meticillin-resistant _Staphylococcus aureus_: a genomic epidemiology analysis', and includes:
* SupplementaryData1.xlsx: SNP and time distances between pairs of isolates from the same individual for cohort 1. Named as appendix 2 in the publication.
* SupplementaryData2.xlsx: SNP and time distances between pairs of isolates from the same individual for the external independent cohort. Named as appendix 3 in the publication.
* SupplementaryData3.xlsx: SNP distances among epidemiologically linked cases in cohort 2. Named as appendix 4 in the publication.
* snp-cutoff.within_host_diversity.R: The R code implementing approach A where data from cohort 1 (SupplementaryData1.xlsx) was used to calculate a theoretical SNP cutoff informed by the _S. aureus_ substitution rate and within-host diversity (approach A; figure 1A); and where a simulation model (approach B), parameterised with the same data, was used to estimate SNP distances between simulated transmission pairs (approach B; figure 1A).
* data_fit.R, data_fit_all_map.R and simu_transmission_fn.R: R code used in simulation model, in R script snp-cutoff.within_host_diversity.R.
* snp-cutoff.within_host_diversity.price2017.R: The R code implementing approach A where data from the external independent cohort (SupplementaryData2.xlsx) was used to calculate a theoretical SNP cutoff informed by the _S. aureus_ substitution rate and within-host diversity.
* snp-cutoff.between_hosts_links.R: SNP distances among epidemiologically linked cases in cohort 2 (SupplementaryData3.xlsx) used by Approach C.

# Required dependencies

See top lines in R scripts for required R packages.

# Citation

Coll F, Raven KE, Knight GM, Blane B, Harrison EM, Leek D, Enoch DA, Brown NM, Parkhill J, Peacock SJ. Definition of a genetic relatedness cutoff to exclude recent transmission of meticillin-resistant _Staphylococcus aureus_: a genomic epidemiology analysis. Lancet Microbe. 2020 Dec;1(8):e328-e335. doi: 10.1016/S2666-5247(20)30149-X. PMID: 33313577; PMCID:
PMC7721685. [URL](https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(20)30149-X/fulltext)
