# Shark-Chromosome-Synteny

Chromosome synteny analysis and comparative analysis of macro-chromosomes and micro-chromosomes.

This repository contains all scripts related to chromosome analysis from the paper. Updates will be made following the paper’s publication.

The flowchart below outlines the main steps in shark chromosome analysis.

![14291730896249_ pic](https://github.com/user-attachments/assets/7af9d366-cd4a-4456-b440-c2cd8b05a8bb)

## Chromosome synteny

We assume the following sequence of karyotype evolution:

​	1.	Separation into different chromosomes

​	2.	Fusion within the same chromosome without mixing

​	3.	Fusion within the same chromosome with mixing

We aim to determine if this sequence aligns with specific phylogenetic hypotheses; for more details, please see [Schultz et al., *Nature* 2023](https://www.nature.com/articles/s41586-023-05936-6).

We used the [odp pipeline](https://github.com/conchoecia/odp?tab=readme-ov-file) to reconstruct ancestral shark chromosomes for synteny analysis.
