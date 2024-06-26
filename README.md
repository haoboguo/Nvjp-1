# Nvjp-1 (PLoS One 2024)
This repository contains 1000 AF2 models of both Nvjp-1 (an intrinsically disordered protein) and T7RdhA (a globular protein).

R codes for residue-residue interaction network analysis, and secondary structure statistics, can be found in the folder /codes

Our models indicate that the Nvjp-1 models derivate from each other (mean RMSD ~6.5 Angs.) and all T7RdhA models show consistent structures (mean RMSD ~0.5 Angs).
However, the pLDDT scores of both proteins are highly consistent for all Nvjp-1 (381 amino acids, blue) models, or for all T7RdhA (406 amino acids, red) models, see the Figure. The shaded areas in the Figure are mean+/-sd (cyan for Nvjp-1 and pink for T7RdhA). Therefore, the pLDDT scores provided by AF2 may serve as a useful feature for protein structures and/or disorders.

![The pLDDT Profiles of T7RdhA (red) and Nvjp-1 (blue)](https://github.com/haoboguo/Nvjp-1/blob/main/t7rdha-nvjp1.plddt.all.png)

More info of T7RdhA could be found here: https://www.nature.com/articles/s41598-023-30310-x 

The residue-residue interaction networks (RINs) of Nvjp-1 (top) and T7RdhA (bottom) show distinct patterns, in which the residue-residue interactions of the IDP (Nvjp-1) are mostly transient, but those of the well-folded protein (T7RdhA) are persistent. In the RINs, each vertex represents one amino acid residue and each edge is an interaction. The red edges are persistent interactions that can be observed in more than 75% of all models, whereas the blue edges are transient interactions that can be observed in less than 25% of all models. Histograms (C and F) indicate interactions in the IDP's are mostly transient. A and D are circular representations of the RIN; B and E are optimized representations.

![Residue-residue interaction networks](https://github.com/haoboguo/Nvjp-1/blob/main/RIN.png)

Note: the RINs have been constructed based on the contact maps (physical adjacencies); the red bar in Fig C (persistent interactions in Nvjp-1) comes from the consecutive residues that are connected by the peptide bonds.

The secondary structure element (SSE) analysis shows that the Nvjp-1 models (top) possess little or no secondary structures; whereas T7RdhA models (bottom) have considerable alpha-helices and beta-strands. The number of residues in each SSE category is shown (mean+/-sd).

![Secondary Structure Elements](https://github.com/haoboguo/Nvjp-1/blob/main/SSE.png)

Molecular dynamics simulations were performed to the Nvjp-1 monomer and a hypothetical dimer. The RINs were constructed using configurations taken from the respective 500-ns MD trajectories, which also shows considerably more transient interactions than the persistent interactions. The hypothetical dimer has a doughnut shape with two monomeric chains tangled with each other, forming more stable interactions. A and D are circular representations of the RIN; the two circles in D represent the two chains of Nvjp-1 homodimer.

![MD RINs](https://github.com/haoboguo/Nvjp-1/blob/main/RIN-MD.png)
