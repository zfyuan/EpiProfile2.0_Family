# EpiProfile2.0_Family
EpiProfile 2.0: A Computational Platform for Processing Epi-Proteomics Mass Spectrometry Data

Summary

Epigenetic marks, mostly DNA methylation and histone post-translational modifications (PTMs), play important roles in chromatin structure and function. Accurate quantification of these marks is an ongoing challenge due to the variety of modifications and their wide dynamic range of abundance. We present EpiProfile 2.0, an extended version of our 2015 software (v1.0) for the enhanced quantification of histone peptides based on LC-MS/MS analysis. The unique features of EpiProfile 2.0 are discriminating the mixture of isobaric peptides and determining the retention time of modified peptides for all histones through data-independent acquisition (DIA). Various applications of EpiProfile 2.0 include different acquisition methods (i.e. DIA and DDA) and different sample preparations (e.g. label-free or isotopic labeling), source organisms, histone mutations, different derivatization strategies, and low-abundance PTMs (such as acyl-derived modifications). Being the first software of its kind we anticipate that EpiProfile 2.0 will play a fundamental role in epigenetic studies relevant to biology and translational medicine.

Methods

The details of each version of EpiProfile 2.0 are shown as below.

EpiProfile2.0_1Basic

	Support Human (norganism=1) and Mouse (norganism=2)
	Support label-free quantification (nsource=1), SILAC (Arg10 (nsource=2, nsubtype=0), Lys8Arg10 (nsource=2, nsubtype=1)), C13 on acetylation group (nsource=3), N15 labeling on amino acids (nsource=4), 13CD3 on methylation group and Methionine (nsource=5)
	Automatically determine the MS method (e.g. acquisition – DIA or DDA, fragmentation – CID, HCD, or ETD, resolution of MS1 – high or low, length of gradient)
	('norganism', 'nsource', and 'nsubtype' are parameters in the 'paras.txt' file under the folder of EpiProfile)

EpiProfile2.0_2Organisms

	Support 11 additional species (i.e. Bos taurus, Caenorhabditis elegans, Harpegnathos saltator, Heterocephalus glaber, Neurospora crassa, Oxytricha trifallax, Theileria annulata, Plasmodium falciparum, Saccharomyces cerevisiae, Saccharum officinarum, and Xenopus laevis)

EpiProfile2.0_3Mutations

	Support 30 known or predicted missense histone mutations (i.e. H33A29V_T32I, H33A15G, H33R17G, H33A29P, H33P121R, H33K27M, H33G34R, H33G34V, H33G34W, H33K36M, H31K27M, H33T45I, H33G90R, H33G33E, H33G34A, H33V35L, H33K36A, H33K36I, H33K36R, H33K36T, H33K36Q, H33K36E, H33K36Nle, H33K37E, H33K37Q, H33K37T, H33K37N, H33K37R, H31G34W, and H33K27R_G34R)

EpiProfile2.0_4Anhydrides

	Suppport 1 derivatized chemical anhydride (phenyl isocyanate – PIC)

EpiProfile2.0_5Low-abundancePTMs

	Support 4 low-abundance PTMs (i.e. H2AK13K15ub, H3K27acK36me, H3R17meR42me, and H3T3ph)
	Support 27 acyl-CoA species or combinations for eight H3 peptides (i.e. 3−8, 9−17, 18−26, 27−40, 54−63, 64−69, 73−83, and 117−128) and seven H4 peptides (i.e. 4−17, 20−23, 24−35, 41−45, 56−67, 68−78, and 79−92)

If more requests for EpiProfile 2.0, please feel free to contact us (emails in the 'User Manual of EpiProfile.pdf' file).
