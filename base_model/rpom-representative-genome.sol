=== Solution on 11-Dec-2024  17:28:04 for organism GCF_000011965 (version 28.5) using MetaFlux (Pathway Tools version 28.0) and solver SCIP (version 6.0.0)

==== --> This is a best solution (i.e., an optimal solution) <-- ====

==== --> This run was done in solving mode. <-- ====

==================== Model Statistics ====================

1) Reactions
------------

  1A) Specified reactions in model:   1604
  1B) Reactions from FBA input file after filtering:    807
      1C) Transport reactions:                            1
      1D) Spontaneous reactions:                         58
      1E) Enzymatic reactions:                          591
  1F) Reactions after instantiation of 1B:           915
  1G) After splitting reversible reactions of 1F:   1194
      1H) Carrying non-zero flux:                    185

2) Metabolites
--------------

  2A) Nutrients:               44
  2B) Secretions:              51
  2C) Biomass metabolites:     44
  2D) Metabolites that are reactants/products in reactions from 1B:    816
  2E) Metabolites that are reactants/products in reactions from 1F:    831

3) Genes/Proteins
-----------------

  3A) Enzymes for reactions from 1E:                   545
  3B) Genes coding for enzymes from 3A:                545
  3C) Transporters for transport reactions from 1C:      9
  3D) Genes coding for enzymes of 3C:                    9

Notes:
1G) It is these reactions that form the reaction set used in the FBA model.
See Section "Instantiation of Generic Reactions" in Chapter "MetaFlux: Flux Balance Analysis"
of the Pathway Tools' User Guide for a description of the reaction instantiation process.


==================== BIOMASS ====================


===== Produced Biomass Flux Metabolite(s) (biomass reaction flux of    30.000000000) (44 fixed biomass metabolites)

Flux:  30.000000000	ADP[CCO-CYTOSOL]	ADP
Flux:  30.000000000	ATP[CCO-CYTOSOL]	ATP
Flux:  30.000000000	CO+2[CCO-CYTOSOL]	Co2+
Flux:  30.000000000	FAD[CCO-CYTOSOL]	FAD
Flux:  30.000000000	FE+2[CCO-CYTOSOL]	Fe2+
Flux:  30.000000000	FE+3[CCO-CYTOSOL]	Fe3+
Flux:  30.000000000	GTP[CCO-CYTOSOL]	GTP
Flux:  30.000000000	PROTON[CCO-CYTOSOL]	H+
Flux:  30.000000000	WATER[CCO-CYTOSOL]	H2O
Flux:  30.000000000	L-ALPHA-ALANINE[CCO-CYTOSOL]	L-alanine
Flux:  30.000000000	ARG[CCO-CYTOSOL]	L-arginine
Flux:  30.000000000	ASN[CCO-CYTOSOL]	L-asparagine
Flux:  30.000000000	L-ASPARTATE[CCO-CYTOSOL]	L-aspartate
Flux:  30.000000000	CYS[CCO-CYTOSOL]	L-cysteine
Flux:  30.000000000	GLT[CCO-CYTOSOL]	L-glutamate
Flux:  30.000000000	GLN[CCO-CYTOSOL]	L-glutamine
Flux:  30.000000000	HIS[CCO-CYTOSOL]	L-histidine
Flux:  30.000000000	ILE[CCO-CYTOSOL]	L-isoleucine
Flux:  30.000000000	LEU[CCO-CYTOSOL]	L-leucine
Flux:  30.000000000	LYS[CCO-CYTOSOL]	L-lysine
Flux:  30.000000000	PHE[CCO-CYTOSOL]	L-phenylalanine
Flux:  30.000000000	PRO[CCO-CYTOSOL]	L-proline
Flux:  30.000000000	SER[CCO-CYTOSOL]	L-serine
Flux:  30.000000000	THR[CCO-CYTOSOL]	L-threonine
Flux:  30.000000000	TRP[CCO-CYTOSOL]	L-tryptophan
Flux:  30.000000000	TYR[CCO-CYTOSOL]	L-tyrosine
Flux:  30.000000000	VAL[CCO-CYTOSOL]	L-valine
Flux:  30.000000000	NAD[CCO-CYTOSOL]	NAD+
Flux:  30.000000000	NADH[CCO-CYTOSOL]	NADH
Flux:  30.000000000	NADP[CCO-CYTOSOL]	NADP+
Flux:  30.000000000	NADPH[CCO-CYTOSOL]	NADPH
Flux:  30.000000000	AMMONIUM[CCO-CYTOSOL]	ammonium
Flux:  30.000000000	DATP[CCO-CYTOSOL]	dATP
Flux:  30.000000000	DGTP[CCO-CYTOSOL]	dGTP
Flux:  30.000000000	TTP[CCO-CYTOSOL]	dTTP
Flux:  30.000000000	PPI[CCO-CYTOSOL]	diphosphate
Flux:  30.000000000	GLY[CCO-CYTOSOL]	glycine
Flux:  30.000000000	CPD-3[CCO-CYTOSOL]	molybdate
Flux:  30.000000000	Pi[CCO-CYTOSOL]	phosphate
Flux:  30.000000000	PROTOHEME[CCO-CYTOSOL]	protoheme
Flux:  30.000000000	PYRIDOXAL_PHOSPHATE[CCO-CYTOSOL]	pyridoxal 5'-phosphate
Flux:  30.000000000	RIBOFLAVIN[CCO-CYTOSOL]	riboflavin
Flux:  30.000000000	SULFATE[CCO-CYTOSOL]	sulfate
Flux:  30.000000000	THIAMINE-PYROPHOSPHATE[CCO-CYTOSOL]	thiamine diphosphate

==================== NUTRIENTS ====================

========== Consumed Nutrient(s) (16 such nutrients out of 44)

===== Consumed Fixed Nutrient(s) (16 such nutrients out of 44)

Flux:  30.000000000 	CO+2[CCO-CYTOSOL]	Co2+
Flux:  60.000000000 	FE+2[CCO-CYTOSOL]	Fe2+
Flux:  30.000000000 	FE+3[CCO-CYTOSOL]	Fe3+
Flux: 1046.400000000 	D-glucopyranose-6-phosphate[CCO-CYTOSOL]	D-glucopyranose 6-phosphate
Flux:  30.000000000 	CPD-3[CCO-CYTOSOL]	molybdate
Flux: 3000.000000000 	AMMONIUM[CCO-CYTOSOL]	ammonium
Flux:  90.000000000 	OXYGEN-MOLECULE[CCO-CYTOSOL]	dioxygen
Flux:   3.600000000 	Pi[CCO-CYTOSOL]	phosphate
Flux: 884.400000000 	PYRUVATE[CCO-CYTOSOL]	pyruvate
Flux:  30.000000000 	SULFATE[CCO-CYTOSOL]	sulfate
Flux: 120.000000000 	SUC[CCO-CYTOSOL]	succinate
Flux:  30.000000000 	SO3[CCO-CYTOSOL]	sulfite
Flux:  30.000000000 	THIAMINE[CCO-CYTOSOL]	thiamine
Flux: 120.000000000 	CIT[CCO-CYTOSOL]	citrate
Flux: 1539.600000000 	FUM[CCO-CYTOSOL]	fumarate
Flux:  30.000000000 	THYMIDINE[CCO-CYTOSOL]	thymidine

===== NOT Consumed Fixed Nutrient(s) (28 such metabolites out of 44)

Flux:   0.000000000 	ACET[CCO-CYTOSOL]	acetate
Flux:   0.000000000 	N-acetyl-D-glucosamine[CCO-CYTOSOL]	N-acetyl-D-glucosamine
Flux:   0.000000000 	CHITOBIOSE[CCO-CYTOSOL]	N,N'-diacetylchitobiose
Flux:   0.000000000 	CL-[CCO-CYTOSOL]	chloride
Flux:   0.000000000 	CARBON-DIOXIDE[CCO-CYTOSOL]	CO2
Flux:   0.000000000 	GALACTOSE[CCO-CYTOSOL]	beta-D-galactopyranose
Flux:   0.000000000 	GLYCEROL[CCO-CYTOSOL]	glycerol
Flux:   0.000000000 	GLYCOLLATE[CCO-CYTOSOL]	glycolate
Flux:  -0.000000000 	PROTON[CCO-CYTOSOL]	H+
Flux:   0.000000000 	WATER[CCO-CYTOSOL]	H2O
Flux:   0.000000000 	SELENATE[CCO-CYTOSOL]	selenate
Flux:   0.000000000 	SELENITE[CCO-CYTOSOL]	selenite
Flux:   0.000000000 	SPERMIDINE[CCO-CYTOSOL]	spermidine
Flux:   0.000000000 	BIOTIN[CCO-CYTOSOL]	biotin
Flux:   0.000000000 	HCO3[CCO-CYTOSOL]	hydrogencarbonate
Flux:   0.000000000 	CPD-335[CCO-CYTOSOL]	(R)-3-hydroxybutanoate
Flux:   0.000000000 	CPD-12693[CCO-CYTOSOL]	(2S)-2,3-dihydroxypropane-1-sulfonate
Flux:   0.000000000 	MAL[CCO-CYTOSOL]	(S)-malate
Flux:   0.000000000 	PUTRESCINE[CCO-CYTOSOL]	putrescine
Nutrient CA+2[CCO-CYTOSOL] is not in any reaction of the model
Nutrient FRU[CCO-CYTOSOL] is not in any reaction of the model
Nutrient Glucose[CCO-CYTOSOL] is not in any reaction of the model
Nutrient MANNOSE[CCO-CYTOSOL] is not in any reaction of the model
Nutrient MG+2[CCO-CYTOSOL] is not in any reaction of the model
Nutrient NA+[CCO-CYTOSOL] is not in any reaction of the model
Nutrient RIBOSE[CCO-CYTOSOL] is not in any reaction of the model
Nutrient D-L-Xylose[CCO-CYTOSOL] is not in any reaction of the model
Nutrient ZN+2[CCO-CYTOSOL] is not in any reaction of the model

==================== SECRETIONS ====================

===== Secreted Fixed Secretion(s) (4 such metabolites out of 51)

Flux: 2481.600000000 	ACET[CCO-CYTOSOL]	acetate
Flux: 1894.800000000 	PROTON[CCO-CYTOSOL]	H+
Flux: 3000.000000000 	WATER[CCO-CYTOSOL]	H2O
Flux: 1666.800000000 	HCO3[CCO-CYTOSOL]	hydrogencarbonate

===== NOT Secreted Fixed Secretion(s) (zero or near zero flux) (47 such metabolites)

Flux:   0.000000000 	N-acetyl-D-glucosamine[CCO-CYTOSOL]	N-acetyl-D-glucosamine
Flux:   0.000000000 	CHITOBIOSE[CCO-CYTOSOL]	N,N'-diacetylchitobiose
Flux:   0.000000000 	CL-[CCO-CYTOSOL]	chloride
Flux:   0.000000000 	CARBON-DIOXIDE[CCO-CYTOSOL]	CO2
Flux:   0.000000000 	CO+2[CCO-CYTOSOL]	Co2+
Flux:   0.000000000 	SS-DIMETHYL-BETA-PROPIOTHETIN[CCO-CYTOSOL]	dimethylsulfoniopropanoate
Flux:   0.000000000 	FE+2[CCO-CYTOSOL]	Fe2+
Flux:   0.000000000 	FE+3[CCO-CYTOSOL]	Fe3+
Flux:   0.000000000 	D-glucopyranose-6-phosphate[CCO-CYTOSOL]	D-glucopyranose 6-phosphate
Flux:   0.000000000 	GALACTOSE[CCO-CYTOSOL]	beta-D-galactopyranose
Flux:   0.000000000 	GLT[CCO-CYTOSOL]	L-glutamate
Flux:   0.000000000 	GLYCEROL[CCO-CYTOSOL]	glycerol
Flux:   0.000000000 	GLYCOLLATE[CCO-CYTOSOL]	glycolate
Flux:   0.000000000 	CPD-3[CCO-CYTOSOL]	molybdate
Flux:   0.000000000 	AMMONIUM[CCO-CYTOSOL]	ammonium
Flux:   0.000000000 	OXYGEN-MOLECULE[CCO-CYTOSOL]	dioxygen
Flux:   0.000000000 	Pi[CCO-CYTOSOL]	phosphate
Flux:   0.000000000 	PYRUVATE[CCO-CYTOSOL]	pyruvate
Flux:   0.000000000 	SELENATE[CCO-CYTOSOL]	selenate
Flux:   0.000000000 	SELENITE[CCO-CYTOSOL]	selenite
Flux:   0.000000000 	SULFATE[CCO-CYTOSOL]	sulfate
Flux:   0.000000000 	SPERMIDINE[CCO-CYTOSOL]	spermidine
Flux:   0.000000000 	SUC[CCO-CYTOSOL]	succinate
Flux:   0.000000000 	SO3[CCO-CYTOSOL]	sulfite
Flux:   0.000000000 	THIAMINE[CCO-CYTOSOL]	thiamine
Flux:   0.000000000 	BIOTIN[CCO-CYTOSOL]	biotin
Flux:   0.000000000 	CPD-335[CCO-CYTOSOL]	(R)-3-hydroxybutanoate
Flux:   0.000000000 	CPD-12693[CCO-CYTOSOL]	(2S)-2,3-dihydroxypropane-1-sulfonate
Flux:   0.000000000 	CARNITINE[CCO-CYTOSOL]	L-carnitine
Flux:   0.000000000 	CHOLINE[CCO-CYTOSOL]	choline
Flux:   0.000000000 	CIT[CCO-CYTOSOL]	citrate
Flux:   0.000000000 	L-CYSTEATE[CCO-CYTOSOL]	L-cysteate
Flux:   0.000000000 	FUM[CCO-CYTOSOL]	fumarate
Flux:   0.000000000 	MAL[CCO-CYTOSOL]	(S)-malate
Flux:   0.000000000 	PUTRESCINE[CCO-CYTOSOL]	putrescine
Flux:   0.000000000 	THYMIDINE[CCO-CYTOSOL]	thymidine
Flux:   0.000000000 	ADP[CCO-CYTOSOL]	ADP
Flux:  -0.000000000 	PPI[CCO-CYTOSOL]	diphosphate
Secretion CA+2[CCO-CYTOSOL] is not in any reaction of the model
Secretion FRU[CCO-CYTOSOL] is not in any reaction of the model
Secretion Glucose[CCO-CYTOSOL] is not in any reaction of the model
Secretion MANNOSE[CCO-CYTOSOL] is not in any reaction of the model
Secretion MG+2[CCO-CYTOSOL] is not in any reaction of the model
Secretion NA+[CCO-CYTOSOL] is not in any reaction of the model
Secretion RIBOSE[CCO-CYTOSOL] is not in any reaction of the model
Secretion D-L-Xylose[CCO-CYTOSOL] is not in any reaction of the model
Secretion ZN+2[CCO-CYTOSOL] is not in any reaction of the model

==================== REACTIONS ====================

Note: The reactions are written in the same direction as they carry flux.

====== Reactions from the PGDB (Fixed Reactions) with NON-zero Flux (185 such reactions)

Flux: 2964.000000000	(PEPSYNTH-RXN)	phosphoenolpyruvate + AMP + phosphate + 2 H+  ->  pyruvate + ATP + H2O

Flux: 2814.000000000	(ADENYL-KIN-RXN)	2 ADP  ->  ATP + AMP

Flux: 2451.600000000	(ACETATE--COA-LIGASE-RXN)	acetyl-CoA + AMP + diphosphate  ->  coenzyme A + acetate + ATP

Flux: 2176.800000000	(PYRUVFORMLY-RXN)	pyruvate + coenzyme A  ->  formate + acetyl-CoA

Flux: 2079.600000000	(ASPARTASE-RXN)	ammonium + fumarate  ->  L-aspartate

Flux: 1894.800000000	(FORMATETHFLIG-RXN *generic*)	a tetrahydrofolate + ATP + formate  ->  an N10-formyltetrahydrofolate + ADP + phosphate

Flux: 1696.800000000	(RXN0-5224)	CO2 + H2O  ->  hydrogencarbonate + H+

Flux: 1461.600000000	(PYRUVATEORTHOPHOSPHATE-DIKINASE-RXN)	pyruvate + ATP + phosphate  ->  phosphoenolpyruvate + AMP + diphosphate + H+

Flux: 1149.600000000	(PEPCARBOXYKIN-RXN)	oxaloacetate + ATP  ->  CO2 + phosphoenolpyruvate + ADP

Flux: 1144.800000000	(METHENYLTHFCYCLOHYDRO-RXN *generic*)	an N10-formyltetrahydrofolate + H+  ->  a 5,10-methenyltetrahydrofolate + H2O

Flux: 1144.800000000	(METHYLENETHFDEHYDROG-NADP-RXN *generic*)	a 5,10-methenyltetrahydrofolate + NADPH  ->  a 5,10-methylenetetrahydrofolate + NADP+

Flux: 1054.800000000	(GCVT-RXN *generic*)	a [glycine-cleavage complex H protein] N6-dihydrolipoyl-L-lysine + a 5,10-methylenetetrahydrofolate + ammonium  ->  a [glycine-cleavage complex H protein] N6-aminomethyldihydrolipoyl-L-lysine + a tetrahydrofolate

Flux: 1054.800000000	(GCVP-RXN *generic*)	a [glycine-cleavage complex H protein] N6-aminomethyldihydrolipoyl-L-lysine + CO2  ->  glycine + a [glycine-cleavage complex H protein] N6-[(R)-lipoyl]-L-lysine + H+

Flux: 1054.800000000	(RXN-8629 *generic*)	a [glycine-cleavage complex H protein] N6-[(R)-lipoyl]-L-lysine + NADH + H+  ->  a [glycine-cleavage complex H protein] N6-dihydrolipoyl-L-lysine + NAD+

Flux: 759.600000000	(Instantiation of 1.4.1.21-RXN)	L-aspartate + NADP+ + H2O  ->  oxaloacetate + ammonium + NADPH + H+

Flux: 750.000000000	(Instantiation of L-GLN-FRUCT-6-P-AMINOTRANS-RXN)	alpha-D-glucosamine 6-phosphate + L-glutamate  ->  beta-D-fructofuranose 6-phosphate + L-glutamine

Flux: 750.000000000	(GLUCOSAMINE-6-P-DEAMIN-RXN)	beta-D-fructofuranose 6-phosphate + ammonium  ->  alpha-D-glucosamine 6-phosphate + H2O

Flux: 540.000000000	(PRPPSYN-RXN *generic*)	D-ribose 5-phosphate + ATP  ->  5-phospho-alpha-D-ribose 1-diphosphate + AMP + H+

Flux: 540.000000000	(RXN-9952)	D-gluconate 6-phosphate + NADP+  ->  D-ribulose 5-phosphate + CO2 + NADPH

Flux: 540.000000000	(RIB5PISOM-RXN *generic*)	D-ribulose 5-phosphate  ->  D-ribose 5-phosphate

Flux: 540.000000000	(6PGLUCONOLACT-RXN)	6-phospho D-glucono-1,5-lactone + H2O  ->  D-gluconate 6-phosphate + H+

Flux: 540.000000000	(GLU6PDEHYDROG-RXN *generic*)	D-glucopyranose 6-phosphate + NADP+  ->  6-phospho D-glucono-1,5-lactone + NADPH + H+

Flux: 532.800000000	(PHOSGLYPHOS-RXN)	3-phospho-D-glyceroyl phosphate + ADP  ->  3-phospho-D-glycerate + ATP

Flux: 532.800000000	(3PGAREARR-RXN)	3-phospho-D-glycerate  ->  2-phospho-D-glycerate

Flux: 532.800000000	(2PGADEHYDRAT-RXN)	2-phospho-D-glycerate  ->  phosphoenolpyruvate + H2O

Flux: 532.800000000	(GAPOXNPHOSPHN-RXN)	D-glyceraldehyde 3-phosphate + NAD+ + phosphate  ->  3-phospho-D-glyceroyl phosphate + NADH + H+

Flux: 506.400000000	(RXN-16998 *generic*)	D-glucopyranose 6-phosphate + phosphorylated phosphoglucomutase  ->  alpha-D-glucose 1,6-bisphosphate + phosphoglucomutase

Flux: 506.400000000	(Instantiation of RXN-16998)	alpha-D-glucose 1,6-bisphosphate + phosphoglucomutase  ->  alpha-D-glucose 6-phosphate + phosphorylated phosphoglucomutase

Flux: 506.400000000	(PGLUCISOM-RXN)	alpha-D-glucose 6-phosphate  ->  beta-D-fructofuranose 6-phosphate

Flux: 462.000000000	(1.2.1.2-RXN)	formate + NAD+  ->  CO2 + NADH

Flux: 394.800000000	(THIOREDOXIN-REDUCT-NADPH-RXN *generic*)	oxidized thioredoxin + NADPH + H+  ->  reduced thioredoxin + NADP+

Flux: 390.000000000	(ASPAMINOTRANS-RXN)	L-aspartate + 2-oxoglutarate  ->  oxaloacetate + L-glutamate

Flux: 390.000000000	(IMPCYCLOHYDROLASE-RXN)	5-formamido-1-(5-phospho-D-ribosyl)-imidazole-4-carboxamide  ->  IMP + H2O

Flux: 390.000000000	(AICARTRANSFORM-RXN *generic*)	an N10-formyltetrahydrofolate + 5-amino-1-(5-phospho-D-ribosyl)imidazole-4-carboxamide  ->  a tetrahydrofolate + 5-formamido-1-(5-phospho-D-ribosyl)-imidazole-4-carboxamide

Flux: 386.400000000	(F16ALDOLASE-RXN)	beta-D-fructofuranose 1,6-bisphosphate  ->  glycerone phosphate + D-glyceraldehyde 3-phosphate

Flux: 386.400000000	(6PFRUCTPHOS-RXN)	ATP + beta-D-fructofuranose 6-phosphate  ->  ADP + beta-D-fructofuranose 1,6-bisphosphate + H+

Flux: 360.000000000	(AICARSYN-RXN)	5'-phosphoribosyl-4-(N-succinocarboxamide)-5-aminoimidazole  ->  fumarate + 5-amino-1-(5-phospho-D-ribosyl)imidazole-4-carboxamide

Flux: 360.000000000	(SAICARSYN-RXN)	ATP + 5-amino-1-(5-phospho-D-ribosyl)imidazole-4-carboxylate + L-aspartate  ->  ADP + 5'-phosphoribosyl-4-(N-succinocarboxamide)-5-aminoimidazole + phosphate + H+

Flux: 360.000000000	(AIRCARBOXY-RXN)	5-amino-1-(5-phospho-beta-D-ribosyl)imidazole + CO2  ->  5-amino-1-(5-phospho-D-ribosyl)imidazole-4-carboxylate + 2 H+

Flux: 360.000000000	(GDPKIN-RXN)	ATP + GDP  ->  ADP + GTP

Flux: 360.000000000	(PRPPAMIDOTRANS-RXN)	5-phospho-alpha-D-ribose 1-diphosphate + L-glutamine + H2O  ->  5-phospho-beta-D-ribosylamine + L-glutamate + diphosphate

Flux: 360.000000000	(AIRS-RXN)	ATP + 2-(formamido)-N1-(5-phospho-beta-D-ribosyl)acetamidine  ->  ADP + 5-amino-1-(5-phospho-beta-D-ribosyl)imidazole + phosphate + H+

Flux: 360.000000000	(GLYRIBONUCSYN-RXN)	5-phospho-beta-D-ribosylamine + ATP + glycine  ->  ADP + N1-(5-phospho-beta-D-ribosyl)glycinamide + phosphate + H+

Flux: 360.000000000	(FGAMSYN-RXN)	ATP + N2-formyl-N1-(5-phospho-beta-D-ribosyl)glycinamide + L-glutamine + H2O  ->  L-glutamate + ADP + 2-(formamido)-N1-(5-phospho-beta-D-ribosyl)acetamidine + phosphate + H+

Flux: 360.000000000	(GART-RXN *generic*)	an N10-formyltetrahydrofolate + N1-(5-phospho-beta-D-ribosyl)glycinamide  ->  a tetrahydrofolate + N2-formyl-N1-(5-phospho-beta-D-ribosyl)glycinamide + H+

Flux: 334.800000000	(RXN-7566 *generic*)	glycine + reduced thioredoxin + phosphate + H+  ->  acetyl phosphate + ammonium + oxidized thioredoxin + H2O

Flux: 334.800000000	(PHOSACETYLTRANS-RXN)	acetyl phosphate + coenzyme A  ->  acetyl-CoA + phosphate

Flux: 270.000000000	(SUCCCOASYN-RXN)	succinate + coenzyme A + ATP  ->  succinyl-CoA + ADP + phosphate

Flux: 270.000000000	(ADENYLOSUCCINATE-SYNTHASE-RXN)	L-aspartate + IMP + GTP  ->  N6-(1,2-dicarboxyethyl)AMP + GDP + phosphate + 2 H+

Flux: 270.000000000	(AMPSYN-RXN)	N6-(1,2-dicarboxyethyl)AMP  ->  AMP + fumarate

Flux: 266.400000000	(TRIOSEPISOMERIZATION-RXN)	glycerone phosphate  ->  D-glyceraldehyde 3-phosphate

Flux: 240.000000000	(5-AMINOLEVULINIC-ACID-SYNTHASE-RXN)	glycine + succinyl-CoA + H+  ->  CO2 + coenzyme A + 5-aminolevulinate

Flux: 210.000000000	(RXN-11811 *spontaneous*)	ammonium  ->  ammonia + H+

Flux: 120.000000000	(ACONITATEDEHYDR-RXN)	citrate  ->  cis-aconitate + H2O

Flux: 120.000000000	(RXN-9951)	D-threo-isocitrate + NADP+  ->  oxalosuccinate + NADPH + H+

Flux: 120.000000000	(2TRANSKETO-RXN)	beta-D-fructofuranose 6-phosphate + D-glyceraldehyde 3-phosphate  ->  D-erythrose 4-phosphate + D-xylulose 5-phosphate

Flux: 120.000000000	(PORPHOBILSYNTH-RXN)	2 5-aminolevulinate  ->  porphobilinogen + H+ + 2 H2O

Flux: 120.000000000	(IMP-DEHYDROG-RXN)	IMP + NAD+ + H2O  ->  XMP + NADH + H+

Flux: 120.000000000	(GUANYL-KIN-RXN)	ATP + GMP  ->  ADP + GDP

Flux: 120.000000000	(RIBULP3EPIM-RXN)	D-xylulose 5-phosphate  ->  D-ribulose 5-phosphate

Flux: 120.000000000	(RXN-8642)	oxalosuccinate + H+  ->  2-oxoglutarate + CO2

Flux: 120.000000000	(ACONITATEHYDR-RXN)	cis-aconitate + H2O  ->  D-threo-isocitrate

Flux: 120.000000000	(GMP-SYN-NH3-RXN)	XMP + ammonia + ATP  ->  GMP + AMP + diphosphate + H+

Flux: 120.000000000	(DIOHBUTANONEPSYN-RXN)	D-ribulose 5-phosphate  ->  formate + 1-deoxy-L-glycero-tetrulose 4-phosphate + H+

Flux: 120.000000000	(LUMAZINESYN-RXN)	5-amino-6-(D-ribitylamino)uracil + 1-deoxy-L-glycero-tetrulose 4-phosphate  ->  6,7-dimethyl-8-(1-D-ribityl)lumazine + phosphate + H+ + 2 H2O

Flux: 120.000000000	(RXN-9772)	L-aspartate + fumarate  ->  2-iminosuccinate + succinate + H+

Flux: 120.000000000	(NAD-SYNTH-NH3-RXN)	ammonium + ATP + nicotinate adenine dinucleotide  ->  AMP + NAD+ + diphosphate + H+

Flux: 120.000000000	(QUINOLINATE-SYNTHA-RXN)	2-iminosuccinate + glycerone phosphate  ->  quinolinate + phosphate + 2 H2O

Flux: 120.000000000	(QUINOPRIBOTRANS-RXN)	5-phospho-alpha-D-ribose 1-diphosphate + quinolinate + 2 H+  ->  beta-nicotinate D-ribonucleotide + CO2 + diphosphate

Flux: 120.000000000	(NICONUCADENYLYLTRAN-RXN)	beta-nicotinate D-ribonucleotide + ATP + H+  ->  nicotinate adenine dinucleotide + diphosphate

Flux:  90.000000000	(RXN-20673)	3-deoxy-D-arabino-heptulosonate 7-phosphate  ->  3-deoxy-D-arabino-heptulopyranuronate 7-phosphate

Flux:  90.000000000	(CHORISMATE-SYNTHASE-RXN)	5-enolpyruvoyl-shikimate 3-phosphate  ->  chorismate + phosphate

Flux:  90.000000000	(3-DEHYDROQUINATE-DEHYDRATASE-RXN)	3-dehydroquinate  ->  3-dehydroshikimate + H2O

Flux:  90.000000000	(GLYOHMETRANS-RXN *generic*)	glycine + a 5,10-methylenetetrahydrofolate + H2O  ->  L-serine + a tetrahydrofolate

Flux:  90.000000000	(ASPARTATE-SEMIALDEHYDE-DEHYDROGENASE-RXN)	L-aspartyl-4-phosphate + NADPH + H+  ->  L-aspartate 4-semialdehyde + NADP+ + phosphate

Flux:  90.000000000	(ASPARTATEKIN-RXN)	L-aspartate + ATP  ->  L-aspartyl-4-phosphate + ADP

Flux:  90.000000000	(SHIKIMATE-KINASE-RXN)	shikimate + ATP  ->  shikimate 3-phosphate + ADP + H+

Flux:  90.000000000	(DAHPSYN-RXN)	phosphoenolpyruvate + D-erythrose 4-phosphate + H2O  ->  3-deoxy-D-arabino-heptulosonate 7-phosphate + phosphate

Flux:  90.000000000	(2.5.1.19-RXN)	shikimate 3-phosphate + phosphoenolpyruvate  ->  5-enolpyruvoyl-shikimate 3-phosphate + phosphate

Flux:  90.000000000	(SHIKIMATE-5-DEHYDROGENASE-RXN)	3-dehydroshikimate + NADPH + H+  ->  shikimate + NADP+

Flux:  90.000000000	(RXN-20674)	3-deoxy-D-arabino-heptulopyranuronate 7-phosphate  ->  3-dehydroquinate + phosphate

Flux:  60.000000000	(Instantiation of HOMOSERDEHYDROG-RXN)	L-aspartate 4-semialdehyde + NADH + H+  ->  L-homoserine + NAD+

Flux:  60.000000000	(RIBOPHOSPHAT-RXN)	5-amino-6-(5-phospho-D-ribitylamino)uracil + H2O  ->  5-amino-6-(D-ribitylamino)uracil + phosphate

Flux:  60.000000000	(ACETOLACTSYN-RXN)	2 pyruvate + H+  ->  (S)-2-acetolactate + CO2

Flux:  60.000000000	(RIBOFLAVINSYNDEAM-RXN)	2,5-diamino-6-(5-phospho-D-ribosylamino)pyrimidin-4(3H)-one + H+ + H2O  ->  5-amino-6-(5-phospho-D-ribosylamino)uracil + ammonium

Flux:  60.000000000	(RIBOFLAVIN-SYN-RXN)	2 6,7-dimethyl-8-(1-D-ribityl)lumazine + H+  ->  5-amino-6-(D-ribitylamino)uracil + riboflavin

Flux:  60.000000000	(GLUTSEMIALDEHYDROG-RXN)	gamma-L-glutamyl 5-phosphate + NADPH + H+  ->  L-glutamate-5-semialdehyde + NADP+ + phosphate

Flux:  60.000000000	(GTP-CYCLOHYDRO-II-RXN)	GTP + 4 H2O  ->  2,5-diamino-6-(5-phospho-D-ribosylamino)pyrimidin-4(3H)-one + formate + 2 phosphate + 3 H+

Flux:  60.000000000	(CATAL-RXN)	2 hydrogen peroxide  ->  2 H2O + dioxygen

Flux:  60.000000000	(RIBOFLAVINSYNREDUC-RXN)	5-amino-6-(5-phospho-D-ribosylamino)uracil + NADPH + H+  ->  5-amino-6-(5-phospho-D-ribitylamino)uracil + NADP+

Flux:  60.000000000	(CHORISMATEMUT-RXN)	chorismate  ->  prephenate

Flux:  60.000000000	(DIHYDROXYISOVALDEHYDRAT-RXN)	(2R)-2,3-dihydroxy-3-methylbutanoate  ->  3-methyl-2-oxobutanoate + H2O

Flux:  60.000000000	(GLUTKIN-RXN)	L-glutamate + ATP  ->  gamma-L-glutamyl 5-phosphate + ADP

Flux:  60.000000000	(ACETOLACTREDUCTOISOM-RXN)	(S)-2-acetolactate + NADPH + H+  ->  (2R)-2,3-dihydroxy-3-methylbutanoate + NADP+

Flux:  60.000000000	(NADH-KINASE-RXN)	ATP + NADH  ->  ADP + NADPH + H+

Flux:  30.000000000	(SERINE-O-ACETTRAN-RXN)	L-serine + acetyl-CoA  ->  O-acetyl-L-serine + coenzyme A

Flux:  30.000000000	(HISTPRATPHYD-RXN)	1-(5-phospho-beta-D-ribosyl)-ATP + H2O  ->  1-(5-phospho-beta-D-ribosyl)-AMP + diphosphate + H+

Flux:  30.000000000	(PRAISOM-RXN)	N-(5-phosphoribosyl)-anthranilate  ->  1-(2-carboxyphenylamino)-1-deoxy-D-ribulose 5-phosphate

Flux:  30.000000000	(SUCCINYLDIAMINOPIMTRANS-RXN)	L-glutamate + N-succinyl-2-amino-6-ketopimelate  ->  2-oxoglutarate + N-succinyl-L,L-2,6-diaminopimelate

Flux:  30.000000000	(IGPSYN-RXN)	1-(2-carboxyphenylamino)-1-deoxy-D-ribulose 5-phosphate + H+  ->  (1S,2R)-1-C-(indol-3-yl)glycerol 3-phosphate + CO2 + H2O

Flux:  30.000000000	(DADPKIN-RXN)	ATP + dADP  ->  ADP + dATP

Flux:  30.000000000	(DIAMINOPIMEPIM-RXN)	L,L-diaminopimelate  ->  meso-diaminopimelate

Flux:  30.000000000	(RIBOFLAVINKIN-RXN)	ATP + riboflavin  ->  ADP + FMN + H+

Flux:  30.000000000	(Instantiation of RXN-16062)	(S)-2-aceto-2-hydroxybutanoate + NADH + H+  ->  (R)-2,3-dihydroxy-3-methylpentanoate + NAD+

Flux:  30.000000000	(DGDPKIN-RXN)	ATP + dGDP  ->  ADP + dGTP

Flux:  30.000000000	(ADPREDUCT-RXN *generic*)	ADP + reduced thioredoxin  ->  dADP + oxidized thioredoxin + H2O

Flux:  30.000000000	(DTDPKIN-RXN)	ATP + dTDP  ->  ADP + dTTP

Flux:  30.000000000	(HISTCYCLOHYD-RXN)	1-(5-phospho-beta-D-ribosyl)-AMP + H2O  ->  1-(5-phospho-beta-D-ribosyl)-5-[(5-phosphoribosylamino)methylideneamino]imidazole-4-carboxamide

Flux:  30.000000000	(ACETOOHBUTSYN-RXN)	pyruvate + 2-oxobutanoate + H+  ->  (S)-2-aceto-2-hydroxybutanoate + CO2

Flux:  30.000000000	(TETHYDPICSUCC-RXN)	(S)-2,3,4,5-tetrahydrodipicolinate + succinyl-CoA + H2O  ->  N-succinyl-2-amino-6-ketopimelate + coenzyme A

Flux:  30.000000000	(DIHYDROXYMETVALDEHYDRAT-RXN)	(R)-2,3-dihydroxy-3-methylpentanoate  ->  (3S)-3-methyl-2-oxopentanoate + H2O

Flux:  30.000000000	(RXN-13179)	4-phosphooxy-L-threonine + NAD+  ->  3-amino-1-hydroxyacetone 1-phosphate + CO2 + NADH

Flux:  30.000000000	(RXN-7800 *spontaneous*)	(2S)-2-isopropyl-3-oxosuccinate + H+  ->  4-methyl-2-oxopentanoate + CO2

Flux:  30.000000000	(BRANCHED-CHAINAMINOTRANSFERVAL-RXN)	L-glutamate + 3-methyl-2-oxobutanoate  ->  L-valine + 2-oxoglutarate

Flux:  30.000000000	(PSERTRANSAMPYR-RXN)	(3R)-3-hydroxy-2-oxo-4 phosphooxybutanoate + L-glutamate  ->  4-phosphooxy-L-threonine + 2-oxoglutarate

Flux:  30.000000000	(ACSERLY-RXN)	O-acetyl-L-serine + hydrogen sulfide  ->  L-cysteine + acetate + H+

Flux:  30.000000000	(UROGENIIISYN-RXN)	preuroporphyrinogen  ->  uroporphyrinogen-III + H2O

Flux:  30.000000000	(Instantiation of RXN-14014)	(2S,4S)-4-hydroxy-2,3,4,5-tetrahydrodipicolinate + NADH + H+  ->  (S)-2,3,4,5-tetrahydrodipicolinate + NAD+ + H2O

Flux:  30.000000000	(THYKI-RXN)	thymidine + ATP  ->  dTMP + ADP + H+

Flux:  30.000000000	(RXN-14196)	carbamate + ATP  ->  carbamoyl phosphate + ADP

Flux:  30.000000000	(RXN-15123 *spontaneous*)	2-iminobutanoate + H2O + H+  ->  2-oxobutanoate + ammonium

Flux:  30.000000000	(RXN-22729)	(2S)-2-amino-4-deoxy-chorismate  ->  anthranilate + pyruvate + H+

Flux:  30.000000000	(PDXJ-RXN)	1-deoxy-D-xylulose 5-phosphate + 3-amino-1-hydroxyacetone 1-phosphate  ->  pyridoxine 5'-phosphate + phosphate + H+ + 2 H2O

Flux:  30.000000000	(FADSYN-RXN)	ATP + FMN + H+  ->  FAD + diphosphate

Flux:  30.000000000	(THIAMIN-PYROPHOSPHOKINASE-RXN)	thiamine + ATP  ->  AMP + thiamine diphosphate + H+

Flux:  30.000000000	(RXN-20080 *generic*)	ATP + L-aspartyl-[tRNAAsn]  ->  ADP + 4-phosphooxy-L-aspartyl-[tRNAAsn]

Flux:  30.000000000	(PREPHENATEDEHYDROG-RXN)	prephenate + NAD+  ->  3-(4-hydroxyphenyl)pyruvate + CO2 + NADH

Flux:  30.000000000	(PRTRANS-RXN)	anthranilate + 5-phospho-alpha-D-ribose 1-diphosphate  ->  N-(5-phosphoribosyl)-anthranilate + diphosphate

Flux:  30.000000000	(RXN-15121 *spontaneous*)	(2Z)-2-aminobut-2-enoate  ->  2-iminobutanoate

Flux:  30.000000000	(SPONTPRO-RXN *spontaneous*)	L-glutamate-5-semialdehyde  ->  (S)-1-pyrroline-5-carboxylate + H+ + H2O

Flux:  30.000000000	(RXN-10814)	3-phenyl-2-oxopropanoate + L-glutamate  ->  L-phenylalanine + 2-oxoglutarate

Flux:  30.000000000	(PREPHENATEDEHYDRAT-RXN)	prephenate + H+  ->  3-phenyl-2-oxopropanoate + CO2 + H2O

Flux:  30.000000000	(RXN0-2382)	L-serine + indole  ->  L-tryptophan + H2O

Flux:  30.000000000	(RXN-21524 *generic*)	a [hydroxymethylbilane synthase] ES3 intermediate + porphobilinogen  ->  a [hydroxymethylbilane synthase] ES4 intermediate + ammonium

Flux:  30.000000000	(RXN-12460 *generic*)	L-asparaginyl-[tRNAAsn] + H2O  ->  L-asparagine + tRNAAsn + H+

Flux:  30.000000000	(RXN-16910)	carboxyphosphate + ammonium  ->  carbamate + phosphate + 2 H+

Flux:  30.000000000	(TYROSINE-AMINOTRANSFERASE-RXN)	3-(4-hydroxyphenyl)pyruvate + L-glutamate  ->  L-tyrosine + 2-oxoglutarate

Flux:  30.000000000	(RXN-21525 *generic*)	a [hydroxymethylbilane synthase] ES4 intermediate + H2O  ->  preuroporphyrinogen + a holo-[hydroxymethylbilane synthase]

Flux:  30.000000000	(BRANCHED-CHAINAMINOTRANSFERILEU-RXN)	L-glutamate + (3S)-3-methyl-2-oxopentanoate  ->  L-isoleucine + 2-oxoglutarate

Flux:  30.000000000	(THRESYN-RXN)	O-phospho-L-homoserine + H2O  ->  L-threonine + phosphate

Flux:  30.000000000	(SULFITE-REDUCT-RXN)	sulfite + 3 NADPH + 5 H+  ->  hydrogen sulfide + 3 NADP+ + 3 H2O

Flux:  30.000000000	(ARGSUCCINLYA-RXN)	L-arginino-succinate  ->  L-arginine + fumarate

Flux:  30.000000000	(RXN-16909)	hydrogencarbonate + ATP  ->  carboxyphosphate + ADP

Flux:  30.000000000	(RXN-21526 *generic*)	a holo-[hydroxymethylbilane synthase] + porphobilinogen  ->  a holo [hydroxymethylbilane synthase] ES intermediate + ammonium

Flux:  30.000000000	(RXN0-2381)	(1S,2R)-1-C-(indol-3-yl)glycerol 3-phosphate  ->  indole + D-glyceraldehyde 3-phosphate

Flux:  30.000000000	(HISTIDPHOS-RXN)	L-histidinol phosphate + H2O  ->  histidinol + phosphate

Flux:  30.000000000	(DIAMINOPIMDECARB-RXN)	meso-diaminopimelate + H+  ->  L-lysine + CO2

Flux:  30.000000000	(RXN-17900)	phosphoribulosylformimino-AICAR-phosphate + ammonia  ->  D-erythro-1-(imidazol-4-yl)-glycerol 3-phosphate + 5-amino-1-(5-phospho-D-ribosyl)imidazole-4-carboxamide + H2O

Flux:  30.000000000	(3-ISOPROPYLMALISOM-RXN)	(2S)-2-isopropylmalate  ->  2-isopropylmaleate + H2O

Flux:  30.000000000	(ORNCARBAMTRANSFER-RXN)	L-ornithine + carbamoyl phosphate  ->  L-citrulline + phosphate + H+

Flux:  30.000000000	(ERYTH4PDEHYDROG-RXN)	D-erythrose 4-phosphate + NAD+ + H2O  ->  4-phospho-D-erythronate + NADH + 2 H+

Flux:  30.000000000	(DTMPKI-RXN)	ATP + dTMP  ->  ADP + dTDP

Flux:  30.000000000	(PROTOPORGENOXI-RXN)	protoporphyrinogen IX + 3 dioxygen  ->  protoporphyrin IX + 3 hydrogen peroxide

Flux:  30.000000000	(HISTAMINOTRANS-RXN)	3-(imidazol-4-yl)-2-oxopropyl phosphate + L-glutamate  ->  L-histidinol phosphate + 2-oxoglutarate

Flux:  30.000000000	(HISTALDEHYD-RXN)	histidinal + NAD+ + H2O  ->  L-histidine + NADH + 2 H+

Flux:  30.000000000	(RXN-8991)	2-isopropylmaleate + H2O  ->  (2R,3S)-3-isopropylmalate

Flux:  30.000000000	(RXN-22711)	chorismate + ammonia + H+  ->  (2S)-2-amino-4-deoxy-chorismate + H2O

Flux:  30.000000000	(GDPREDUCT-RXN *generic*)	GDP + reduced thioredoxin  ->  dGDP + oxidized thioredoxin + H2O

Flux:  30.000000000	(RXN0-1461)	coproporphyrinogen III + dioxygen + 2 H+  ->  protoporphyrinogen IX + 2 CO2 + 2 H2O

Flux:  30.000000000	(UROGENDECARBOX-RXN)	uroporphyrinogen-III + 4 H+  ->  coproporphyrinogen III + 4 CO2

Flux:  30.000000000	(SUCCDIAMINOPIMDESUCC-RXN)	N-succinyl-L,L-2,6-diaminopimelate + H2O  ->  L,L-diaminopimelate + succinate

Flux:  30.000000000	(RXN-20081 *generic*)	4-phosphooxy-L-aspartyl-[tRNAAsn] + ammonia  ->  L-asparaginyl-[tRNAAsn] + phosphate

Flux:  30.000000000	(ARGSUCCINSYN-RXN)	L-aspartate + L-citrulline + ATP  ->  L-arginino-succinate + AMP + diphosphate + H+

Flux:  30.000000000	(PNPOXI-RXN)	pyridoxine 5'-phosphate + dioxygen  ->  hydrogen peroxide + pyridoxal 5'-phosphate

Flux:  30.000000000	(ATPPHOSPHORIBOSYLTRANS-RXN)	ATP + 5-phospho-alpha-D-ribose 1-diphosphate  ->  1-(5-phospho-beta-D-ribosyl)-ATP + diphosphate

Flux:  30.000000000	(Instantiation of PYRROLINECARBREDUCT-RXN)	(S)-1-pyrroline-5-carboxylate + NADH + 2 H+  ->  L-proline + NAD+

Flux:  30.000000000	(RXN-21527 *generic*)	a holo [hydroxymethylbilane synthase] ES intermediate + porphobilinogen  ->  a [hydroxymethylbilane synthase] ES2 intermediate + ammonium

Flux:  30.000000000	(ORNITHINE-GLU-AMINOTRANSFERASE-RXN)	L-glutamate + L-glutamate-5-semialdehyde  ->  L-ornithine + 2-oxoglutarate

Flux:  30.000000000	(ALANINE-DEHYDROGENASE-RXN)	pyruvate + ammonium + NADH + H+  ->  L-alanine + NAD+ + H2O

Flux:  30.000000000	(2-ISOPROPYLMALATESYN-RXN)	acetyl-CoA + 3-methyl-2-oxobutanoate + H2O  ->  (2S)-2-isopropylmalate + coenzyme A + H+

Flux:  30.000000000	(IMIDPHOSDEHYD-RXN)	D-erythro-1-(imidazol-4-yl)-glycerol 3-phosphate  ->  3-(imidazol-4-yl)-2-oxopropyl phosphate + H2O

Flux:  30.000000000	(DXS-RXN)	pyruvate + D-glyceraldehyde 3-phosphate + H+  ->  1-deoxy-D-xylulose 5-phosphate + CO2

Flux:  30.000000000	(RXN-21523 *generic*)	a [hydroxymethylbilane synthase] ES2 intermediate + porphobilinogen  ->  a [hydroxymethylbilane synthase] ES3 intermediate + ammonium

Flux:  30.000000000	(3-ISOPROPYLMALDEHYDROG-RXN)	(2R,3S)-3-isopropylmalate + NAD+  ->  (2S)-2-isopropyl-3-oxosuccinate + NADH + H+

Flux:  30.000000000	(RXN490-3616 *generic*)	tRNAAsn + L-aspartate + ATP  ->  L-aspartyl-[tRNAAsn] + AMP + diphosphate

Flux:  30.000000000	(BRANCHED-CHAINAMINOTRANSFERLEU-RXN)	L-glutamate + 4-methyl-2-oxopentanoate  ->  L-leucine + 2-oxoglutarate

Flux:  30.000000000	(RXN-14049)	L-homoserine  ->  (2Z)-2-aminobut-2-enoate + H+ + H2O

Flux:  30.000000000	(PRIBFAICARPISOM-RXN)	1-(5-phospho-beta-D-ribosyl)-5-[(5-phosphoribosylamino)methylideneamino]imidazole-4-carboxamide  ->  phosphoribulosylformimino-AICAR-phosphate

Flux:  30.000000000	(PROTOHEMEFERROCHELAT-RXN)	protoporphyrin IX + Fe2+  ->  protoheme + 2 H+

Flux:  30.000000000	(HISTOLDEHYD-RXN)	histidinol + NAD+  ->  histidinal + NADH + H+

Flux:  30.000000000	(ERYTHRON4PDEHYDROG-RXN)	4-phospho-D-erythronate + NAD+  ->  (3R)-3-hydroxy-2-oxo-4 phosphooxybutanoate + NADH + H+

Flux:  30.000000000	(HOMOSERKIN-RXN)	L-homoserine + ATP  ->  O-phospho-L-homoserine + ADP + H+

Flux:  30.000000000	(DIHYDRODIPICSYN-RXN)	pyruvate + L-aspartate 4-semialdehyde  ->  (2S,4S)-4-hydroxy-2,3,4,5-tetrahydrodipicolinate + H2O + H+

Flux:   0.000000000	(MALATE--COA-LIGASE-RXN)	(S)-malate + coenzyme A + ATP  ->  (S)-malyl-CoA + ADP + phosphate

Flux:   0.000000000	(FUMHYDR-RXN)	fumarate + H2O  ->  (S)-malate


====== Reactions from the PGDB (Fixed Reactions) with ZERO Flux (1009 such reactions)

(ACONITATEDEHYDR-RXN)	cis-aconitate + H2O  ->  citrate

(RXN0-366)	guanosine + H2O  ->  D-ribofuranose + guanine

(THREODEHYD-RXN)	L-threonine + NAD+  ->  L-2-amino-3-oxobutanoate + NADH + H+

(GLUCOSE-6-PHOSPHATE-1-EPIMERASE-RXN *spontaneous*)	alpha-D-glucose 6-phosphate  ->  beta-D-glucose 6-phosphate

(GLUCOSE-6-PHOSPHATE-1-EPIMERASE-RXN *spontaneous*)	beta-D-glucose 6-phosphate  ->  alpha-D-glucose 6-phosphate

(RXN-12865)	selenodiglutathione + NADPH + H+  ->  glutathioselenol + glutathione + NADP+

(NAD-SYNTH-GLN-RXN)	ATP + nicotinate adenine dinucleotide + L-glutamine + H2O  ->  L-glutamate + AMP + NAD+ + diphosphate + H+

(XYLULOKIN-RXN-CPD-24961/ATP//XYLULOSE-5-PHOSPHATE/ADP/PROTON.47. *instantiated*)	alpha-D-xylulofuranose + ATP  ->  D-xylulose 5-phosphate + ADP + H+

(RXN-17886 *spontaneous*)	S-(hydroxysulfenamide)glutathione + glutathione  ->  glutathione disulfide + hydroxylamine

(RXN-22316 *spontaneous*)	arsenite  ->  arsenite2- + H+

(RXN-22316 *spontaneous*)	arsenite2- + H+  ->  arsenite

(HOMOSERDEHYDROG-RXN-HOMO-SER/NADP//L-ASPARTATE-SEMIALDEHYDE/NADPH/PROTON.53. *instantiated*)	L-aspartate 4-semialdehyde + NADPH + H+  ->  L-homoserine + NADP+

(RXN-5721 *spontaneous*)	2-amino-3-carboxymuconate-6-semialdehyde  ->  quinolinate + H2O + H+

(DXPREDISOM-RXN)	1-deoxy-D-xylulose 5-phosphate + NADPH + H+  ->  2-C-methyl-D-erythritol 4-phosphate + NADP+

(RXN-13142)	(R)-NADPHX  ->  (S)-NADPHX

(RXN-17781)	(S)-3-hydroxy-(7Z)-hexadecenoyl-CoA + NAD+  ->  3-oxo-(7Z)-hexadecenoyl-CoA + NADH + H+

(CYSTEINE-AMINOTRANSFERASE-RXN)	L-cysteine + 2-oxoglutarate  ->  2-oxo-3-sulfanylpropanoate + L-glutamate

(CYSTEINE-AMINOTRANSFERASE-RXN)	2-oxo-3-sulfanylpropanoate + L-glutamate  ->  L-cysteine + 2-oxoglutarate

(RXN-14808 *spontaneous*)	aldehydo-L-arabinose  ->  L-arabinopyranose

(RXN-14808 *spontaneous*)	L-arabinopyranose  ->  aldehydo-L-arabinose

(RXN-10058)	2,5-diamino-6-(5-phospho-D-ribitylamino)pyrimidin-4(3H)-one + H+ + H2O  ->  5-amino-6-(5-phospho-D-ribitylamino)uracil + ammonium

(RXN-10058)	5-amino-6-(5-phospho-D-ribitylamino)uracil + ammonium  ->  2,5-diamino-6-(5-phospho-D-ribitylamino)pyrimidin-4(3H)-one + H+ + H2O

(RXN-22438)	a 5,10-methylenetetrahydrofolate + NADH + H+  ->  a 5-methyltetrahydrofolate + NAD+

(HEMN-RXN)	coproporphyrinogen III + 2 S-adenosyl-L-methionine  ->  protoporphyrinogen IX + 2 CO2 + 2 L-methionine + 2 5'-deoxyadenosine

(RXN-17782)	3-oxo-(7Z)-hexadecenoyl-CoA + coenzyme A  ->  (5Z)-tetradecenoyl-CoA + acetyl-CoA

(RXN-14048)	L-cystathionine  ->  L-cysteine + (2Z)-2-aminobut-2-enoate + H+

(GLUCISOM-RXN)	alpha-D-glucopyranose  ->  alpha-D-fructofuranose

(GLUCISOM-RXN)	alpha-D-fructofuranose  ->  alpha-D-glucopyranose

(RXN0-5391)	(2E,5Z)-tetradecenoyl-CoA  ->  (3E,5Z)-tetradeca-3,5-dienoyl-CoA

(RXN0-5391)	(3E,5Z)-tetradeca-3,5-dienoyl-CoA  ->  (2E,5Z)-tetradecenoyl-CoA

(RXN-9789)	[ThiS sulfur-carrier protein] + ATP + H+  ->  carboxy-adenylated-[ThiS sulfur-carrier protein] + diphosphate

(KYNURENINASE-RXN)	L-kynurenine + H2O  ->  L-alanine + anthranilate + H+

(KYNURENINASE-RXN)	L-alanine + anthranilate + H+  ->  L-kynurenine + H2O

(RXN-18208)	3-[(3'-methylsulfanyl)propyl]malate + NAD+  ->  3-carboxy-6-(methylsulfanyl)-2-oxohexanoate + NADH + H+

(RXN-18208)	3-carboxy-6-(methylsulfanyl)-2-oxohexanoate + NADH + H+  ->  3-[(3'-methylsulfanyl)propyl]malate + NAD+

(RXN-7594 *spontaneous*)	4-oxoglutaramate + H2O  ->  2-oxoglutarate + ammonium

(SUPEROX-DISMUT-RXN)	2 superoxide + 2 H+  ->  hydrogen peroxide + dioxygen

(ADCLY-RXN)	4-amino-4-deoxychorismate  ->  4-aminobenzoate + pyruvate + H+

(RXN-9951)	oxalosuccinate + NADPH + H+  ->  D-threo-isocitrate + NADP+

(THI-P-SYN-RXN)	4-methyl-5-(2-phosphooxyethyl)thiazole + 4-amino-2-methyl-5-(diphosphooxymethyl)pyrimidine + H+  ->  thiamine phosphate + diphosphate

(RXN-9550)	palmitoleoyl-[acp] + H2O  ->  palmitoleate + acyl-carrier protein + H+

(PHENYLSERINE-ALDOLASE-RXN)	L-threo-3-phenylserine  ->  glycine + benzaldehyde

(PHENYLSERINE-ALDOLASE-RXN)	glycine + benzaldehyde  ->  L-threo-3-phenylserine

(RXN-2881 *spontaneous*)	formaldehyde + a tetrahydrofolate  ->  a 5,10-methylenetetrahydrofolate + H2O

(RXN-22604 *spontaneous*)	(R)-lactate + H+  ->  (R)-lactic acid

(RXN-22604 *spontaneous*)	(R)-lactic acid  ->  (R)-lactate + H+

(RXN-11455)	7,8-dihydro-8-oxoguanine + H+ + H2O  ->  urate + ammonium

(RXN-11455)	urate + ammonium  ->  7,8-dihydro-8-oxoguanine + H+ + H2O

(RXN-22740)	L-ribulofuranose  ->  beta-L-ribofuranse

(RXN-22740)	beta-L-ribofuranse  ->  L-ribulofuranose

(RXN-14207)	phosphoenolpyruvate + dGDP + H+  ->  pyruvate + dGTP

(RXN-8141)	a CDP-diacylglycerol + a phosphatidylglycerol  ->  a cardiolipin + CMP + H+

(RXN-8141)	a cardiolipin + CMP + H+  ->  a CDP-diacylglycerol + a phosphatidylglycerol

(3.4.13.19-RXN-D-ALA-D-ALA/WATER//L-ALPHA-ALANINE.35. *instantiated*)	D-alanyl-D-alanine + H2O  ->  2 L-alanine

(ADDALT-RXN)	2'-deoxyadenosine + H2O + H+  ->  2'-deoxyinosine + ammonium

(2TRANSKETO-RXN)	D-erythrose 4-phosphate + D-xylulose 5-phosphate  ->  beta-D-fructofuranose 6-phosphate + D-glyceraldehyde 3-phosphate

(RXN0-1133)	coenzyme A + a [pyruvate dehydrogenase E2 protein] N6-S-acetyldihydrolipoyl-L-lysine  ->  acetyl-CoA + a [pyruvate dehydrogenase E2 protein] N6-dihydrolipoyl-L-lysine

(SEDOBISALDOL-RXN)	D-sedoheptulose 1,7-bisphosphate  ->  glycerone phosphate + D-erythrose 4-phosphate

(SEDOBISALDOL-RXN)	glycerone phosphate + D-erythrose 4-phosphate  ->  D-sedoheptulose 1,7-bisphosphate

(RXN-11375)	a ceramide + H2O  ->  a sphingoid base + a fatty acid

(RXN-11375)	a sphingoid base + a fatty acid  ->  a ceramide + H2O

(RXN-12134)	vernolate + H2O  ->  12,13-DiHOME

(RXN-12134)	12,13-DiHOME  ->  vernolate + H2O

(RXN-12197)	UDP + H2O  ->  UMP + phosphate + H+

(RXN-22739-CPD-12045//CPD-24959.21. *instantiated*)	alpha-L-arabinofuranose  ->  beta-L-ribulofuranose

(RXN-22739-CPD-12045//CPD-24959.21. *instantiated*)	beta-L-ribulofuranose  ->  alpha-L-arabinofuranose

(RXN-22739-CPD-12045//CPD-24938.21. *instantiated*)	alpha-L-arabinofuranose  ->  alpha-L-ribulofuranose

(RXN-22739-CPD-12045//CPD-24938.21. *instantiated*)	alpha-L-ribulofuranose  ->  alpha-L-arabinofuranose

(RXN-22739-CPD-12045//L-Ribulofuranose.28. *instantiated*)	alpha-L-arabinofuranose  ->  L-ribulofuranose

(RXN-22739-CPD-12045//L-Ribulofuranose.28. *instantiated*)	L-ribulofuranose  ->  alpha-L-arabinofuranose

(RXN-22739-CPD-12046//CPD-24959.21. *instantiated*)	beta-L-arabinofuranose  ->  beta-L-ribulofuranose

(RXN-22739-CPD-12046//CPD-24959.21. *instantiated*)	beta-L-ribulofuranose  ->  beta-L-arabinofuranose

(RXN-22739-CPD-12046//CPD-24938.21. *instantiated*)	beta-L-arabinofuranose  ->  alpha-L-ribulofuranose

(RXN-22739-CPD-12046//CPD-24938.21. *instantiated*)	alpha-L-ribulofuranose  ->  beta-L-arabinofuranose

(RXN-22739-CPD-12046//L-Ribulofuranose.28. *instantiated*)	beta-L-arabinofuranose  ->  L-ribulofuranose

(RXN-22739-CPD-12046//L-Ribulofuranose.28. *instantiated*)	L-ribulofuranose  ->  beta-L-arabinofuranose

(RXN-22739-L-arabinofuranose//CPD-24959.29. *instantiated*)	L-arabinofuranose  ->  beta-L-ribulofuranose

(RXN-22739-L-arabinofuranose//CPD-24959.29. *instantiated*)	beta-L-ribulofuranose  ->  L-arabinofuranose

(RXN-22739-L-arabinofuranose//CPD-24938.29. *instantiated*)	L-arabinofuranose  ->  alpha-L-ribulofuranose

(RXN-22739-L-arabinofuranose//CPD-24938.29. *instantiated*)	alpha-L-ribulofuranose  ->  L-arabinofuranose

(RXN-22739)	L-arabinofuranose  ->  L-ribulofuranose

(RXN-22739)	L-ribulofuranose  ->  L-arabinofuranose

(INOPHOSPHOR-RXN)	inosine + phosphate  ->  hypoxanthine + alpha-D-ribose-1-phosphate

(INOPHOSPHOR-RXN)	hypoxanthine + alpha-D-ribose-1-phosphate  ->  inosine + phosphate

(ARYLFORMAMIDASE-RXN)	N-Formyl-L-kynurenine + H2O  ->  L-kynurenine + formate + H+

(RXN-13141)	(S)-NADPHX + ADP  ->  AMP + NADPH + phosphate + H+

(RXN-12867)	glutathioselenol + NADPH + H+  ->  hydrogen selenide + glutathione + NADP+

(S-ADENMETSYN-RXN)	ATP + L-methionine + H2O  ->  S-adenosyl-L-methionine + phosphate + diphosphate

(RXN-22317 *spontaneous*)	arsenite2-  ->  arsenite3- + H+

(RXN-22317 *spontaneous*)	arsenite3- + H+  ->  arsenite2-

(RXN-11348)	2 {beta-D-GlcNAc-(1->4)-3-O-[L-Ala-gamma-D-iGln-(6-N-beta-D-Asn)-L-Lys-D-Ala]-Mur2Ac}-PP-Und  ->  a nascent peptidoglycan dimer (E. faecium, tetrapeptide) + di-trans,octa-cis-undecaprenyl diphosphate + H+

(RXN-11348)	a nascent peptidoglycan dimer (E. faecium, tetrapeptide) + di-trans,octa-cis-undecaprenyl diphosphate + H+  ->  2 {beta-D-GlcNAc-(1->4)-3-O-[L-Ala-gamma-D-iGln-(6-N-beta-D-Asn)-L-Lys-D-Ala]-Mur2Ac}-PP-Und

(2OXOGLUTDECARB-RXN)	2-oxoglutarate + a [2-oxoglutarate dehydrogenase E2 protein] N6-lipoyl-L-lysine + H+  ->  a [2-oxoglutarate dehydrogenase E2 protein] N6-S-succinyldihydrolipoyl-L-lysine + CO2

(RXN0-313)	beta-D-fructofuranose 6-phosphate  ->  dihydroxyacetone + D-glyceraldehyde 3-phosphate

(RXN0-313)	dihydroxyacetone + D-glyceraldehyde 3-phosphate  ->  beta-D-fructofuranose 6-phosphate

(SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31. *instantiated*)	succinate[in] + ubiquinone-8[generic membrane]  ->  fumarate[in] + ubiquinol-8[generic membrane]

(RXN-12610)	2-(2-carboxy-4-methylthiazol-5-yl)ethyl phosphate + 4-amino-2-methyl-5-(diphosphooxymethyl)pyrimidine + 2 H+  ->  thiamine phosphate + CO2 + diphosphate

(RXN-12609)	2-[(2R,5Z)-2-carboxy-4-methylthiazol-5(2H)-ylidene]ethyl phosphate  ->  2-(2-carboxy-4-methylthiazol-5-yl)ethyl phosphate

(RXN-22319 *spontaneous*)	acetic acid  ->  acetate + H+

(RXN-22319 *spontaneous*)	acetate + H+  ->  acetic acid

(RXN0-884)	(E)-4-hydroxy-3-methylbut-2-en-1-yl diphosphate + 2 reduced ferredoxin [iron-sulfur] cluster + 2 H+  ->  prenyl diphosphate + 2 oxidized ferredoxin [iron-sulfur] cluster + H2O

(GDPPYPHOSKIN-RXN)	ATP + GDP  ->  AMP + ppGpp + H+

(1.5.1.20-RXN)	a 5-methyltetrahydrofolate + NAD(P)+  ->  a 5,10-methylenetetrahydrofolate + NAD(P)H + H+

(1.5.1.20-RXN)	a 5,10-methylenetetrahydrofolate + NAD(P)H + H+  ->  a 5-methyltetrahydrofolate + NAD(P)+

(RXN0-2145)	(2E,5Z)-dodeca-3,5-dienoyl-[acp] + NADH + H+  ->  (5Z)-dodec-5-enoyl-[acp] + NAD+

(HYDROXYLAMINE-REDUCTASE-NADH-RXN)	hydroxylamine + NADH + 2 H+  ->  ammonium + NAD+ + H2O

(RXN-15889)	[EntF L-seryl-carrier protein] + coenzyme A  ->  holo-[EntF L-seryl-carrier protein] + adenosine 3',5'-bisphosphate + H+

(RXN-15889)	holo-[EntF L-seryl-carrier protein] + adenosine 3',5'-bisphosphate + H+  ->  [EntF L-seryl-carrier protein] + coenzyme A

(RXN-7609)	GMP + H2O  ->  guanosine + phosphate

(BIOTIN-CARBOXYL-RXN)	a [biotin carboxyl-carrier-protein dimer]-N6-biotinyl-L-lysine + hydrogencarbonate + ATP  ->  a [carboxyl-carrier protein dimer]-N6-carboxybiotinyl-L-lysine + ADP + phosphate + H+

(3-CH3-2-OXOBUTANOATE-OH-CH3-XFER-RXN)	a 5,10-methylenetetrahydrofolate + 3-methyl-2-oxobutanoate + H2O  ->  2-dehydropantoate + a tetrahydrofolate

(ACETATE--COA-LIGASE-RXN)	coenzyme A + acetate + ATP  ->  acetyl-CoA + AMP + diphosphate

(MALATE-DEH-RXN)	(S)-malate + NAD+  ->  oxaloacetate + NADH + H+

(MALATE-DEH-RXN)	oxaloacetate + NADH + H+  ->  (S)-malate + NAD+

(RXN-12864 *spontaneous*)	selenite + 4 glutathione + 2 H+  ->  selenodiglutathione + glutathione disulfide + 3 H2O

(NAG1P-URIDYLTRANS-RXN)	UTP + N-acetyl-alpha-D-glucosamine 1-phosphate + H+  ->  UDP-N-acetyl-alpha-D-glucosamine + diphosphate

(SULFATE-ADENYLYLTRANS-RXN)	sulfate + ATP + H+  ->  adenosine 5'-phosphosulfate + diphosphate

(PROTEIN-KINASE-RXN)	ATP + a [protein]-(L-serine/L-threonine)  ->  ADP + a [protein] (L-serine/L-threonine) phosphate + H+

(RXN-17824)	an O4-methylthymine in DNA + a [protein]-L-cysteine  ->  a thymine in DNA + a [protein]-S-methyl-L-cysteine

(RXN-17824)	a thymine in DNA + a [protein]-S-methyl-L-cysteine  ->  an O4-methylthymine in DNA + a [protein]-L-cysteine

(PANTOATE-BETA-ALANINE-LIG-RXN)	beta-alanine + (R)-pantoate + ATP  ->  (R)-pantothenate + AMP + diphosphate + H+

(1.4.1.21-RXN-L-ASPARTATE/NAD/WATER//OXALACETIC_ACID/AMMONIUM/NADH/PROTON.60. *instantiated*)	L-aspartate + NAD+ + H2O  ->  oxaloacetate + ammonium + NADH + H+

(1.2.7.4-RXN)	carbon monoxide + 2 oxidized ferredoxin [iron-sulfur] cluster + H2O  ->  CO2 + 2 reduced ferredoxin [iron-sulfur] cluster + 2 H+

(1.2.7.4-RXN)	CO2 + 2 reduced ferredoxin [iron-sulfur] cluster + 2 H+  ->  carbon monoxide + 2 oxidized ferredoxin [iron-sulfur] cluster + H2O

(KETOBUTFORMLY-RXN)	2-oxobutanoate + coenzyme A  ->  propanoyl-CoA + formate

(RXNQT-4191)	thiamine phosphate + H2O  ->  thiamine + phosphate

(RXN-14291)	all-trans-hexaprenyl diphosphate + 4-aminobenzoate  ->  3-all trans-hexaprenyl-4-aminobenzoate + diphosphate

(RXN-14291)	3-all trans-hexaprenyl-4-aminobenzoate + diphosphate  ->  all-trans-hexaprenyl diphosphate + 4-aminobenzoate

(N-FORMYLGLUTAMATE-DEFORMYLASE-RXN)	N-formyl-L-glutamate + H2O  ->  L-glutamate + formate

(RXN-10024)	butanoyl-[acp] + S-adenosyl-L-methionine  ->  acyl-carrier protein + S-methyl-5'-thioadenosine + PAI-1-2 + H+

(RXN-10024)	acyl-carrier protein + S-methyl-5'-thioadenosine + PAI-1-2 + H+  ->  butanoyl-[acp] + S-adenosyl-L-methionine

(RXN-14192)	phosphoenolpyruvate + dADP + H+  ->  pyruvate + dATP

(ACETYLGLUTKIN-RXN)	N-acetyl-L-glutamate + ATP  ->  N-acetylglutamyl-phosphate + ADP

(RXN0-5462)	GTP + H2O  ->  GDP + phosphate + H+

(RXN-18392 *spontaneous*)	2-dehydro-3-deoxy-D-octonate  ->  3-deoxy-alpha-D-manno-2-octulosonate

(RXN-18392 *spontaneous*)	3-deoxy-alpha-D-manno-2-octulosonate  ->  2-dehydro-3-deoxy-D-octonate

(RXN-9590)	acyl-[acyl-carrier protein] + phosphate  ->  an acyl phosphate + acyl-carrier protein

(RXN-9590)	an acyl phosphate + acyl-carrier protein  ->  acyl-[acyl-carrier protein] + phosphate

(1.6.99.5-RXN-NADH/UBIQUINONE-8/PROTON//NAD/CPD-9956.39. *instantiated*)	NADH[in] + ubiquinone-8[generic membrane] + H+[in]  ->  NAD+[in] + ubiquinol-8[generic membrane]

(ADENYL-KIN-RXN)	ATP + AMP  ->  2 ADP

(4-HYDROXYPHENYLPYRUVATE-DIOXYGENASE-RXN)	3-(4-hydroxyphenyl)pyruvate + dioxygen  ->  CO2 + homogentisate

(ARGININE-DEIMINASE-RXN)	L-arginine + H2O  ->  ammonium + L-citrulline

(RXN0-747)	dADP + oxidized NrdH glutaredoxin-like protein + H2O  ->  ADP + reduced NrdH glutaredoxin-like protein

(RXN0-747)	ADP + reduced NrdH glutaredoxin-like protein  ->  dADP + oxidized NrdH glutaredoxin-like protein + H2O

(IMIDAZOLONEPROPIONASE-RXN)	4-imidazolone-5-propanoate + H2O  ->  N-formimino-L-glutamate

(PSERTRANSAM-RXN)	O-phospho-L-serine + 2-oxoglutarate  ->  L-glutamate + 3-phosphooxypyruvate

(PSERTRANSAM-RXN)	L-glutamate + 3-phosphooxypyruvate  ->  O-phospho-L-serine + 2-oxoglutarate

(RXN0-5199)	guanosine + phosphate  ->  guanine + alpha-D-ribose-1-phosphate

(RXN0-5199)	guanine + alpha-D-ribose-1-phosphate  ->  guanosine + phosphate

(RXN-10655)	(7Z)-3-oxotetradec-7-enoyl-[acp] + NADPH + H+  ->  (3R,7Z)-3-hydroxytetradec-7-enoyl-[acp] + NADP+

(RXN-10025)	octanoyl-[acp] + S-adenosyl-L-methionine  ->  acyl-carrier protein + S-methyl-5'-thioadenosine + VAI-1-2 + H+

(RXN-10025)	acyl-carrier protein + S-methyl-5'-thioadenosine + VAI-1-2 + H+  ->  octanoyl-[acp] + S-adenosyl-L-methionine

(SPERMIDINESYN-RXN)	putrescine + S-adenosyl 3-(methylsulfanyl)propylamine  ->  spermidine + S-methyl-5'-thioadenosine + H+

(UDPGLUCEPIM-RXN)	UDP-alpha-D-glucose  ->  UDP-alpha-D-galactose

(UDPGLUCEPIM-RXN)	UDP-alpha-D-galactose  ->  UDP-alpha-D-glucose

(RXN-16391)	(2S)-ethylmalonyl-CoA  ->  (2R)-ethylmalonyl-CoA

(RXN-22)	canavaninosuccinate  ->  L-canavanine + fumarate

(RXN-22)	L-canavanine + fumarate  ->  canavaninosuccinate

(ACETYL-COA-CARBOXYLTRANSFER-RXN)	ATP + acetyl-CoA + hydrogencarbonate  ->  ADP + malonyl-CoA + phosphate + H+

(ACETYLORNDEACET-RXN)	N-acetyl-L-ornithine + H2O  ->  L-ornithine + acetate

(RXN-9788)	[ThiI sulfur-carrier protein]-S-sulfanyl-L-cysteine + carboxy-adenylated-[ThiS sulfur-carrier protein] + 2 reduced ferredoxin [iron-sulfur] cluster  ->  thiocarboxylated-[ThiS sulfur-carrier protein] + [ThiI sulfur-carrier protein]-L-cysteine + AMP + 2 oxidized ferredoxin [iron-sulfur] cluster

(F16ALDOLASE-RXN)	glycerone phosphate + D-glyceraldehyde 3-phosphate  ->  beta-D-fructofuranose 1,6-bisphosphate

(RXN-17889)	L-arginyl-[tRNAArg] + an N-terminal L-aspartyl-[protein]  ->  tRNAArg + an N-terminal L-arginiyl-L-aspartyl-[protein] + H+

(RXN-17889)	tRNAArg + an N-terminal L-arginiyl-L-aspartyl-[protein] + H+  ->  L-arginyl-[tRNAArg] + an N-terminal L-aspartyl-[protein]

(UNDECAPRENYL-DIPHOSPHATASE-RXN)	di-trans,octa-cis-undecaprenyl diphosphate + H2O  ->  di-trans,octa-cis-undecaprenyl phosphate + phosphate + H+

(RXN-14301)	L-methionine  ->  methanethiol + (2Z)-2-aminobut-2-enoate + H+

(HOMSUCTRAN-RXN)	L-homoserine + succinyl-CoA  ->  O-succinyl-L-homoserine + coenzyme A

(HOMSUCTRAN-RXN)	O-succinyl-L-homoserine + coenzyme A  ->  L-homoserine + succinyl-CoA

(SULFITE-DEHYDROGENASE-RXN-SO3/Oxidized-cytochromes-c551/WATER//SULFATE/Reduced-cytochromes-c551/PROTON.77. *instantiated*)	sulfite + 2 oxidized cytochrome c551 + H2O  ->  sulfate + 2 reduced cytochrome c551 + 2 H+

(SULFITE-DEHYDROGENASE-RXN-SO3/Oxidized-cytochromes-C4/WATER//SULFATE/Reduced-cytochromes-C4/PROTON.73. *instantiated*)	sulfite + 2 oxidized cytochrome c4 + H2O  ->  sulfate + 2 reduced cytochrome c4 + 2 H+

(SULFITE-DEHYDROGENASE-RXN-SO3/Oxidized-cytochromes-c553/WATER//SULFATE/Reduced-cytochromes-c553/PROTON.77. *instantiated*)	sulfite + 2 oxidized cytochrome c-553 + H2O  ->  sulfate + 2 reduced cytochrome c-553 + 2 H+

(SULFITE-DEHYDROGENASE-RXN-SO3/Oxidized-Cytochromes-C6/WATER//SULFATE/Reduced-Cytochromes-C6/PROTON.73. *instantiated*)	sulfite + 2 oxidized cytochrome c6 + H2O  ->  sulfate + 2 reduced cytochrome c6 + 2 H+

(SULFITE-DEHYDROGENASE-RXN-SO3/Oxidized-NapC-proteins/WATER//SULFATE/Reduced-NapC-proteins/PROTON.71. *instantiated*)	sulfite + 2 oxidized [NapC protein] + H2O  ->  sulfate + 2 reduced [NapC protein] + 2 H+

(SULFITE-DEHYDROGENASE-RXN-SO3/an-oxidized-NrfB-protein/WATER//SULFATE/a-reduced-NrfB-protein/PROTON.74. *instantiated*)	sulfite + 2 oxidized [NrfB protein] + H2O  ->  sulfate + 2 reduced [NrfB protein] + 2 H+

(SULFITE-DEHYDROGENASE-RXN-SO3/Oxidized-CycA1-cytochromes/WATER//SULFATE/Reduced-CycA1-cytochromes/PROTON.79. *instantiated*)	sulfite + 2 oxidized CycA1 cytochrome + H2O  ->  sulfate + 2 reduced CycA1 cytochrome + 2 H+

(SULFITE-DEHYDROGENASE-RXN)	sulfite + 2 oxidized c-type cytochrome + H2O  ->  sulfate + 2 reduced c-type cytochrome + 2 H+

(OHACYL-COA-DEHYDROG-RXN-POLYMER-INST-L-3-HYDROXYACYL-COA-C16-H32/NAD//POLYMER-INST-3-KETOACYL-COA-C16-H32/NADH/PROTON.94. *instantiated*)	POLYMER-INST-L-3-HYDROXYACYL-COA-C16-H32 + NAD+  ->  POLYMER-INST-3-KETOACYL-COA-C16-H32 + NADH + H+

(OHACYL-COA-DEHYDROG-RXN-POLYMER-INST-L-3-HYDROXYACYL-COA-C14-H28/NAD//CPD-10260/NADH/PROTON.68. *instantiated*)	POLYMER-INST-L-3-HYDROXYACYL-COA-C14-H28 + NAD+  ->  3-oxooctadecanoyl-CoA + NADH + H+

(OHACYL-COA-DEHYDROG-RXN-POLYMER-INST-L-3-HYDROXYACYL-COA-C12-H24/NAD//POLYMER-INST-3-KETOACYL-COA-C12-H24/NADH/PROTON.94. *instantiated*)	POLYMER-INST-L-3-HYDROXYACYL-COA-C12-H24 + NAD+  ->  POLYMER-INST-3-KETOACYL-COA-C12-H24 + NADH + H+

(OHACYL-COA-DEHYDROG-RXN-POLYMER-INST-L-3-HYDROXYACYL-COA-C10-H20/NAD//POLYMER-INST-3-KETOACYL-COA-C10-H20/NADH/PROTON.94. *instantiated*)	POLYMER-INST-L-3-HYDROXYACYL-COA-C10-H20 + NAD+  ->  POLYMER-INST-3-KETOACYL-COA-C10-H20 + NADH + H+

(RXN-20676)	(S)-3-hydroxydodecanoyl-CoA + NAD+  ->  3-oxododecanoyl-CoA + NADH + H+

(RXN-11662)	(S)-3-hydroxybutanoyl-CoA + NAD+  ->  acetoacetyl-CoA + NADH + H+

(RXN-11662)	acetoacetyl-CoA + NADH + H+  ->  (S)-3-hydroxybutanoyl-CoA + NAD+

(KDO-8PSYNTH-RXN)	D-arabinofuranose 5-phosphate + phosphoenolpyruvate + H2O  ->  3-deoxy-D-manno-octulosonate 8-phosphate + phosphate

(PHOSGLYPHOS-RXN)	3-phospho-D-glycerate + ATP  ->  3-phospho-D-glyceroyl phosphate + ADP

(RXN-10851 *spontaneous*)	2 S-sulfanylglutathione + reduced two electron carrier  ->  bisorganyltrisulfane + hydrogen sulfide + oxidized electron carrier

(RXN-10851 *spontaneous*)	bisorganyltrisulfane + hydrogen sulfide + oxidized electron carrier  ->  2 S-sulfanylglutathione + reduced two electron carrier

(RXN-18031 *spontaneous*)	carbonic acid  ->  hydrogencarbonate + H+

(RXN-18031 *spontaneous*)	hydrogencarbonate + H+  ->  carbonic acid

(HYDROXYMETHYLGLUTARYL-COA-REDUCTASE-RXN)	(R)-mevalonate + coenzyme A + 2 NAD+  ->  (S)-3-hydroxy-3-methylglutaryl-CoA + 2 NADH + 2 H+

(HYDROXYMETHYLGLUTARYL-COA-REDUCTASE-RXN)	(S)-3-hydroxy-3-methylglutaryl-CoA + 2 NADH + 2 H+  ->  (R)-mevalonate + coenzyme A + 2 NAD+

(THI-P-KIN-RXN)	thiamine phosphate + ATP  ->  ADP + thiamine diphosphate

(RXN-17779)	(7Z)-hexadecenoyl-CoA + oxidized electron-transfer flavoprotein + H+  ->  (2E,7Z)-hexadecenoyl-CoA + reduced electron-transfer flavoprotein

(RXN0-2141)	(3Z)-dec-3-enoyl-[acp] + malonyl-[acp] + H+  ->  (5Z)-3-oxododec-5-enoyl-[acp] + acyl-carrier protein + CO2

(BRANCHED-CHAINAMINOTRANSFERVAL-RXN)	L-valine + 2-oxoglutarate  ->  L-glutamate + 3-methyl-2-oxobutanoate

(RXN0-5330-NADH/UBIQUINONE-8/PROTON//NAD/CPD-9956.39. *instantiated*)	NADH[in] + ubiquinone-8[generic membrane] + H+[in]  ->  NAD+[in] + ubiquinol-8[generic membrane]

(RXN-12140)	(9Z,12Z)-15,16-epoxyoctadeca-9,12-dienoate + H2O  ->  15,16-DiHODE

(RXN-12140)	15,16-DiHODE  ->  (9Z,12Z)-15,16-epoxyoctadeca-9,12-dienoate + H2O

(PSERTRANSAMPYR-RXN)	4-phosphooxy-L-threonine + 2-oxoglutarate  ->  (3R)-3-hydroxy-2-oxo-4 phosphooxybutanoate + L-glutamate

(RXN-14857)	pyrithiamine + oxidized electron carrier + H2O  ->  4-amino-2-methyl-5-pyrimidinemethanol + 2-(2-methylpyridin-3-yl)ethanol + reduced two electron carrier

(RXN-14857)	4-amino-2-methyl-5-pyrimidinemethanol + 2-(2-methylpyridin-3-yl)ethanol + reduced two electron carrier  ->  pyrithiamine + oxidized electron carrier + H2O

(RXN-11347)	UDP-N-acetyl-alpha-D-muramoyl-L-alanyl-gamma-D-glutamyl-L-lysyl-D-alanine + di-trans,octa-cis-undecaprenyl phosphate  ->  Und-PP-Mur2Ac-L-Ala-gamma-D-Glu-L-Lys-D-Ala + UMP

(RXN-11347)	Und-PP-Mur2Ac-L-Ala-gamma-D-Glu-L-Lys-D-Ala + UMP  ->  UDP-N-acetyl-alpha-D-muramoyl-L-alanyl-gamma-D-glutamyl-L-lysyl-D-alanine + di-trans,octa-cis-undecaprenyl phosphate

(RXN-15733)	6-[1-hydroxyethyl]-7-methyl-7,8-dihydropterin + ATP  ->  [1-(2-amino-7-methyl-4-oxo-7,8-dihydro-3H-pteridin-6-yl)]ethyl diphosphate + AMP + H+

(RXN-15733)	[1-(2-amino-7-methyl-4-oxo-7,8-dihydro-3H-pteridin-6-yl)]ethyl diphosphate + AMP + H+  ->  6-[1-hydroxyethyl]-7-methyl-7,8-dihydropterin + ATP

(MALATE-DEHYDROGENASE-ACCEPTOR-RXN-MAL/UBIQUINONE-8//OXALACETIC_ACID/CPD-9956.43. *instantiated*)	(S)-malate[in] + ubiquinone-8[generic membrane]  ->  oxaloacetate[in] + ubiquinol-8[generic membrane]

(MALATE-DEHYDROGENASE-ACCEPTOR-RXN-MAL/UBIQUINONE-8//OXALACETIC_ACID/CPD-9956.43. *instantiated*)	oxaloacetate[in] + ubiquinol-8[generic membrane]  ->  (S)-malate[in] + ubiquinone-8[generic membrane]

(DCTP-DEAM-RXN)	dCTP + H+ + H2O  ->  ammonium + dUTP

(IMP-DEHYDROG-RXN)	XMP + NADH + H+  ->  IMP + NAD+ + H2O

(CITRAMALATE-LYASE-RXN)	(S)-citramalate  ->  pyruvate + acetate

(DCDPKIN-RXN)	ATP + dCDP  ->  ADP + dCTP

(ACETATEKIN-RXN)	ATP + acetate  ->  ADP + acetyl phosphate

(ACETATEKIN-RXN)	ADP + acetyl phosphate  ->  ATP + acetate

(ALDXANAU-RXN)	chloroacetaldehyde + NAD+ + H2O  ->  2-chloroacetate + NADH + 2 H+

(ALDXANAU-RXN)	2-chloroacetate + NADH + 2 H+  ->  chloroacetaldehyde + NAD+ + H2O

(RXN-19921 *spontaneous*)	2''-O-acetyl-ADP-ribose  ->  3''-O-acetyl-ADP-ribose

(RXN-19921 *spontaneous*)	3''-O-acetyl-ADP-ribose  ->  2''-O-acetyl-ADP-ribose

(RXN-15261)	7,8-dihydropterin + H2O + H+  ->  7,8-dihydrolumazine + ammonium

(RXN-15261)	7,8-dihydrolumazine + ammonium  ->  7,8-dihydropterin + H2O + H+

(RXN-17783)	(5Z)-tetradecenoyl-CoA + oxidized electron-transfer flavoprotein + H+  ->  (2E,5Z)-tetradecenoyl-CoA + reduced electron-transfer flavoprotein

(PPPGPPHYDRO-RXN)	pppGpp + H2O  ->  ppGpp + phosphate + H+

(3.1.2.21-RXN)	dodecanoyl-[acp] + H2O  ->  acyl-carrier protein + laurate + H+

(RXN-20084)	ammonia + 2-oxoglutarate + NADPH + 2 H+  ->  L-glutamate + NADP+ + H2O

(RXN0-4641)	N-acetyl-D-muramate 6-phosphate + H2O  ->  N-acetyl-D-glucosamine 6-phosphate + (R)-lactate

(RXN0-4641)	N-acetyl-D-glucosamine 6-phosphate + (R)-lactate  ->  N-acetyl-D-muramate 6-phosphate + H2O

(MANDELATE-RACEMASE-RXN)	(S)-mandelate  ->  (R)-mandelate

(MANDELATE-RACEMASE-RXN)	(R)-mandelate  ->  (S)-mandelate

(LTAA-RXN)	L-allo-threonine  ->  glycine + acetaldehyde

(LTAA-RXN)	glycine + acetaldehyde  ->  L-allo-threonine

(3PGAREARR-RXN)	2-phospho-D-glycerate  ->  3-phospho-D-glycerate

(METHYLASPARTATE-AMMONIA-LYASE-RXN)	(2S, 3S)-3-methylaspartate  ->  ammonium + mesaconate

(HYDROXY-MANDELATE-RACEMASE-RXN)	(S)-4-hydroxymandelate  ->  (R)-4-hydroxymandelate

(HYDROXY-MANDELATE-RACEMASE-RXN)	(R)-4-hydroxymandelate  ->  (S)-4-hydroxymandelate

(PTAALT-RXN)	propanoyl-CoA + phosphate  ->  propanoyl phosphate + coenzyme A

(RXN-17887 *spontaneous*)	a [protein] 3-nitrosothio-L-alanine + glutathione  ->  a [protein]-L-cysteine + S-nitrosoglutathione

(RXN-17887 *spontaneous*)	a [protein]-L-cysteine + S-nitrosoglutathione  ->  a [protein] 3-nitrosothio-L-alanine + glutathione

(RXN0-901)	xanthine + NAD+ + H2O  ->  urate + NADH + H+

(RXN0-901)	urate + NADH + H+  ->  xanthine + NAD+ + H2O

(3-HYDROXYISOBUTYRATE-DEHYDROGENASE-RXN)	(S)-3-hydroxy-isobutanoate + NAD+  ->  (S)-methylmalonate-semialdehyde + NADH + H+

(3-HYDROXYISOBUTYRATE-DEHYDROGENASE-RXN)	(S)-methylmalonate-semialdehyde + NADH + H+  ->  (S)-3-hydroxy-isobutanoate + NAD+

(1.10.2.2-RXN-Oxidized-CycA1-cytochromes/CPD-9956//Reduced-CycA1-cytochromes/UBIQUINONE-8/PROTON.83. *instantiated*)	2 oxidized CycA1 cytochrome[out] + ubiquinol-8[generic membrane]  ->  2 reduced CycA1 cytochrome[out] + ubiquinone-8[generic membrane] + 2 H+[out]

(RXN-15816-Oxidized-CycA1-cytochromes/CPD-9956//Reduced-CycA1-cytochromes/UBIQUINONE-8/PROTON.83. *instantiated*)	2 oxidized CycA1 cytochrome[out] + ubiquinol-8[T]  ->  2 reduced CycA1 cytochrome[out] + ubiquinone-8[T] + 2 H+[in]

(PYRUVATEORTHOPHOSPHATE-DIKINASE-RXN)	phosphoenolpyruvate + AMP + diphosphate + H+  ->  pyruvate + ATP + phosphate

(1.5.1.9-RXN)	L-saccharopine + NAD+ + H2O  ->  L-glutamate + (S)-2-amino-6-oxohexanoate + NADH + H+

(1.5.1.9-RXN)	L-glutamate + (S)-2-amino-6-oxohexanoate + NADH + H+  ->  L-saccharopine + NAD+ + H2O

(RXN-12570)	(S)-3-hydroxyhexanoyl-CoA + NAD+  ->  3-oxohexanoyl-CoA + NADH + H+

(ARGDECARBOX-RXN)	L-arginine + H+  ->  agmatine + CO2

(RXN-14014-DELTA1-PIPERIDEINE-2-6-DICARBOXYLATE/NADP/WATER//CPD-14443/NADPH/PROTON.72. *instantiated*)	(2S,4S)-4-hydroxy-2,3,4,5-tetrahydrodipicolinate + NADPH + H+  ->  (S)-2,3,4,5-tetrahydrodipicolinate + NADP+ + H2O

(RXN-20459 *spontaneous*)	S-adenosyl-L-methionine  ->  (R)-S-adenosyl-L-methionine

(RXN-20459 *spontaneous*)	(R)-S-adenosyl-L-methionine  ->  S-adenosyl-L-methionine

(RXN-7931)	(3Z)-dodec-3-enoyl-CoA  ->  (2E)-dodec-2-enoyl-CoA

(RXN-7931)	(2E)-dodec-2-enoyl-CoA  ->  (3Z)-dodec-3-enoyl-CoA

(RXN0-7139)	glucosamine 1,6-diphosphate + phosphoglucosamine mutase  ->  alpha-D-glucosamine 1-phosphate + phosphorylated phosphoglucosamine mutase

(RXN0-7139)	alpha-D-glucosamine 1-phosphate + phosphorylated phosphoglucosamine mutase  ->  glucosamine 1,6-diphosphate + phosphoglucosamine mutase

(F16BDEPHOS-RXN)	beta-D-fructofuranose 1,6-bisphosphate + H2O  ->  beta-D-fructofuranose 6-phosphate + phosphate

(RXN-14196)	carbamoyl phosphate + ADP  ->  carbamate + ATP

(OXALODECARB-RXN)	oxaloacetate + H+  ->  pyruvate + CO2

(RXN-11832)	ATP + CMP  ->  ADP + CDP

(RXN-11832)	ADP + CDP  ->  ATP + CMP

(CYSTEINE--TRNA-LIGASE-RXN)	tRNACys + L-cysteine + ATP  ->  L-cysteinyl-[tRNACys] + AMP + diphosphate

(CITSYN-RXN)	acetyl-CoA + oxaloacetate + H2O  ->  citrate + coenzyme A + H+

(RXN-13037)	octanoyl-[acp] + [glycine cleavage system lipoyl-carrier protein]-L-lysine  ->  a [glycine-cleavage complex H protein] N6-octanoyl-L-lysine + acyl-carrier protein + H+

(RXN-13037)	a [glycine-cleavage complex H protein] N6-octanoyl-L-lysine + acyl-carrier protein + H+  ->  octanoyl-[acp] + [glycine cleavage system lipoyl-carrier protein]-L-lysine

(GUANYL-KIN-RXN)	ADP + GDP  ->  ATP + GMP

(RXN-15583)	L-cysteate  ->  sulfite + 2-aminoprop-2-enoate + H+

(RXN-15583)	sulfite + 2-aminoprop-2-enoate + H+  ->  L-cysteate

(RXN-13733)	a dihydroceramide + H2O  ->  sphinganine + a carboxylate

(RXN-13733)	sphinganine + a carboxylate  ->  a dihydroceramide + H2O

(RXN-17778)	3-oxo-(9Z)-octadecenoyl-CoA + coenzyme A  ->  (7Z)-hexadecenoyl-CoA + acetyl-CoA

(DIHYDROURACIL-DEHYDROGENASE-NAD+-RXN)	5,6-dihydrouracil + NAD+  ->  uracil + NADH + H+

(DIHYDROURACIL-DEHYDROGENASE-NAD+-RXN)	uracil + NADH + H+  ->  5,6-dihydrouracil + NAD+

(1.13.11.6-RXN)	3-hydroxyanthranilate + dioxygen  ->  2-amino-3-carboxymuconate-6-semialdehyde

(THIAMINASE-RXN)	thiamine + H2O  ->  5-(2-hydroxyethyl)-4-methylthiazole + 4-amino-2-methyl-5-pyrimidinemethanol + H+

(DEOXYINOPHOSPHOR-RXN)	2'-deoxyinosine + phosphate  ->  hypoxanthine + 2-deoxy-alpha-D-ribose 1-phosphate

(DEOXYINOPHOSPHOR-RXN)	hypoxanthine + 2-deoxy-alpha-D-ribose 1-phosphate  ->  2'-deoxyinosine + phosphate

(PHOSPHASERSYN-RXN)	a CDP-diacylglycerol + L-serine  ->  CMP + a phosphatidylserine + H+

(FADSYN-RXN)	FAD + diphosphate  ->  ATP + FMN + H+

(RXN-7743)	acetyl-CoA + pyruvate + H2O  ->  (R)-citramalate + coenzyme A + H+

(RXN-12754 *spontaneous*)	NADH + H2O  ->  (R)-NADHX

(3.1.3.16-RXN)	a [protein] (L-serine/L-threonine) phosphate + H2O  ->  a [protein]-(L-serine/L-threonine) + phosphate

(RIBULP3EPIM-RXN)	D-ribulose 5-phosphate  ->  D-xylulose 5-phosphate

(RXN-21315)	5-deoxy-alpha-D-ribose 1-phosphate  ->  5-deoxy-D-ribulose 1-phosphate

(RXN-21315)	5-deoxy-D-ribulose 1-phosphate  ->  5-deoxy-alpha-D-ribose 1-phosphate

(RXN-10020)	3-oxooctanoyl-[acp] + S-adenosyl-L-methionine  ->  acyl-carrier protein + S-methyl-5'-thioadenosine + AAI-1 + H+

(RXN-10020)	acyl-carrier protein + S-methyl-5'-thioadenosine + AAI-1 + H+  ->  3-oxooctanoyl-[acp] + S-adenosyl-L-methionine

(CHD-RXN)	choline + oxidized electron carrier  ->  betaine aldehyde + reduced two electron carrier

(CHD-RXN)	betaine aldehyde + reduced two electron carrier  ->  choline + oxidized electron carrier

(2PGADEHYDRAT-RXN)	phosphoenolpyruvate + H2O  ->  2-phospho-D-glycerate

(RXN-18384)	acryloyl-CoA + CO2 + NADPH  ->  (R)-methylmalonyl-CoA + NADP+

(3-HYDROXYBUTYRATE-DEHYDROGENASE-RXN)	(R)-3-hydroxybutanoate + NAD+  ->  acetoacetate + NADH + H+

(3-HYDROXYBUTYRATE-DEHYDROGENASE-RXN)	acetoacetate + NADH + H+  ->  (R)-3-hydroxybutanoate + NAD+

(3-HYDROXBUTYRYL-COA-DEHYDRATASE-RXN)	(3R)-3-hydroxybutanoyl-CoA  ->  crotonyl-CoA + H2O

(3-HYDROXBUTYRYL-COA-DEHYDRATASE-RXN)	crotonyl-CoA + H2O  ->  (3R)-3-hydroxybutanoyl-CoA

(CYSTATHIONINE-BETA-SYNTHASE-RXN)	L-serine + L-homocysteine  ->  L-cystathionine + H2O

(CYSTATHIONINE-BETA-SYNTHASE-RXN)	L-cystathionine + H2O  ->  L-serine + L-homocysteine

(PEPCARBOXYKIN-RXN)	CO2 + phosphoenolpyruvate + ADP  ->  oxaloacetate + ATP

(RXN-22388 *spontaneous*)	dihydrogen phosphate  ->  phosphate + H+

(RXN-22388 *spontaneous*)	phosphate + H+  ->  dihydrogen phosphate

(RXN-8850)	dUMP + a 5,10-methylenetetrahydrofolate + NADPH + H+  ->  dTMP + a tetrahydrofolate + NADP+

(RXN-16949)	3-hydroxy-4-methyl-benzoate + NADH + dioxygen + H+  ->  4-methylgentisate + NAD+ + H2O

(RXN-16949)	4-methylgentisate + NAD+ + H2O  ->  3-hydroxy-4-methyl-benzoate + NADH + dioxygen + H+

(OHMETPYRKIN-RXN)	ATP + 4-amino-2-methyl-5-pyrimidinemethanol  ->  ADP + 4-amino-2-methyl-5-(phosphooxymethyl)pyrimidine + H+

(RXN-8961)	(2R,3S)-beta-methylmalyl-CoA  ->  glyoxylate + propanoyl-CoA

(RXN0-2144)	(3R,5Z)-3-hydroxydodec-5-enoyl-[acp]  ->  (2E,5Z)-dodeca-3,5-dienoyl-[acp] + H2O

(ACETYLHOMOSER-CYS-RXN)	O-acetyl-L-homoserine + hydrogen sulfide  ->  L-homocysteine + acetate + H+

(3.5.1.80-RXN)	N-acetyl-D-galactosamine 6-phosphate + H2O  ->  D-galactosamine 6-phosphate + acetate

(3.5.1.80-RXN)	D-galactosamine 6-phosphate + acetate  ->  N-acetyl-D-galactosamine 6-phosphate + H2O

(ASPAMINOTRANS-RXN)	oxaloacetate + L-glutamate  ->  L-aspartate + 2-oxoglutarate

(GALACTURIDYLYLTRANS-RXN)	UDP-alpha-D-glucose + alpha-D-galactose 1-phosphate  ->  UDP-alpha-D-galactose + alpha-D-glucopyranose 1-phosphate

(RXN-16112)	(3R,7Z,10Z,13Z,16Z)-3-hydroxydocosatetraenoyl-CoA + NADP+  ->  (7Z,10Z,13Z,16Z)-3-oxodocosatetraenoyl-CoA + NADPH + H+

(RXN-16112)	(7Z,10Z,13Z,16Z)-3-oxodocosatetraenoyl-CoA + NADPH + H+  ->  (3R,7Z,10Z,13Z,16Z)-3-hydroxydocosatetraenoyl-CoA + NADP+

(RXN-9540)	3-oxohexadecanoyl-[acp] + NADPH + H+  ->  (3R)-3-hydroxyhexadecanoyl-[acp] + NADP+

(RXN-8958)	(2R)-ethylmalonyl-CoA  ->  (2S)-methylsuccinyl-CoA

(BADH-RXN)	betaine aldehyde + NAD+ + H2O  ->  glycine betaine + NADH + 2 H+

(RXN-14950)	a [glycine-cleavage complex H protein] N6-octanoyl-L-lysine + 2 sulfurated [sulfur carrier] + 2 reduced [2Fe-2S] ferredoxin + 2 S-adenosyl-L-methionine  ->  a [glycine-cleavage complex H protein] N6-[(R)-lipoyl]-L-lysine + 2 L-methionine + 2 5'-deoxyadenosine + 2 unsulfurated [sulfur carrier] + 2 oxidized [2Fe-2S] ferredoxin

(RXN-7913)	ATP + dCMP  ->  ADP + dCDP

(RXN-14118)	CTP + H+ + H2O  ->  UTP + ammonium

(KDO-8PPHOSPHAT-RXN)	3-deoxy-D-manno-octulosonate 8-phosphate + H2O  ->  3-deoxy-alpha-D-manno-2-octulosonate + phosphate

(NQOR-RXN-UBIQUINONE-8/NADH/PROTON//CPD-9956/NAD.39. *instantiated*)	ubiquinone-8[generic membrane] + NADH[in] + H+[in]  ->  ubiquinol-8[generic membrane] + NAD+[in]

(SARCOX-RXN)	sarcosine + dioxygen + H2O  ->  glycine + formaldehyde + hydrogen peroxide

(KDPGALDOL-RXN)	2-dehydro-3-deoxy-D-gluconate 6-phosphate  ->  D-glyceraldehyde 3-phosphate + pyruvate

(325-BISPHOSPHATE-NUCLEOTIDASE-RXN)	adenosine 3',5'-bisphosphate + H2O  ->  AMP + phosphate

(RXN-14393)	(S)-3-hydroxy-(5Z)-tetradecenoyl-CoA + NAD+  ->  3-oxo-(5Z)-tetradecenoyl-CoA + NADH + H+

(O-ACETYLHOMOSERINE-THIOL-LYASE-RXN)	O-acetyl-L-homoserine + methanethiol  ->  L-methionine + acetate + H+

(RXN-11217)	hydantoin + H2O  ->  N-carbamoylglycine + H+

(RXN-11217)	N-carbamoylglycine + H+  ->  hydantoin + H2O

(1.2.4.4-RXN)	3-methyl-2-oxobutanoate + a [BCAA dehydrogenase E2 protein] N6-lipoyl-L-lysine + H+  ->  an [apo BCAA dehydrogenase E2 protein] N6-S-[2-methylpropanoyl]dihydrolipoyl-L-lysine + CO2

(1.2.4.4-RXN)	an [apo BCAA dehydrogenase E2 protein] N6-S-[2-methylpropanoyl]dihydrolipoyl-L-lysine + CO2  ->  3-methyl-2-oxobutanoate + a [BCAA dehydrogenase E2 protein] N6-lipoyl-L-lysine + H+

(RXN-9380)	carboxynorspermidine + H+  ->  norspermidine + CO2

(RXN-9380)	norspermidine + CO2  ->  carboxynorspermidine + H+

(ACETOOHBUTREDUCTOISOM-RXN)	(S)-2-aceto-2-hydroxybutanoate + NADPH + H+  ->  (R)-2,3-dihydroxy-3-methylpentanoate + NADP+

(XANPRIBOSYLTRAN-RXN)	xanthine + 5-phospho-alpha-D-ribose 1-diphosphate  ->  XMP + diphosphate

(RXN-18703)	a [3-mercaptopyruvate sulfurtransferase]-S-sulfanyl-L-cysteine + reduced thioredoxin  ->  a [3-mercaptopyruvate sulfurtransferase]-L-cysteine + hydrogen sulfide + oxidized thioredoxin

(RXN-18702)	2-oxo-3-sulfanylpropanoate + a [3-mercaptopyruvate sulfurtransferase]-L-cysteine  ->  a [3-mercaptopyruvate sulfurtransferase]-S-sulfanyl-L-cysteine + pyruvate

(S-FORMYLGLUTATHIONE-HYDROLASE-RXN)	S-formylglutathione + H2O  ->  formate + glutathione + H+

(RIB5PISOM-RXN)	D-ribose 5-phosphate  ->  D-ribulose 5-phosphate

(RXN-18206)	3-[(4'-methylsulfanyl)butyl]malate + NAD+  ->  3-carboxy-7-(methylsulfanyl)-2-oxoheptanoate + NADH + H+

(RXN-18206)	3-carboxy-7-(methylsulfanyl)-2-oxoheptanoate + NADH + H+  ->  3-[(4'-methylsulfanyl)butyl]malate + NAD+

(PPGPPSYN-RXN)	ppGpp + H2O  ->  GDP + diphosphate

(RXN0-302)	2-phospho-4-(cytidine 5'-diphospho)-2-C-methyl-D-erythritol  ->  2-C-methyl-D-erythritol-2,4-cyclodiphosphate + CMP

(FORMIMINOGLUTAMATE-DEIMINASE-RXN)	N-formimino-L-glutamate + H2O  ->  ammonium + N-formyl-L-glutamate

(N-ACETYLTRANSFER-RXN)	L-glutamate + acetyl-CoA  ->  N-acetyl-L-glutamate + coenzyme A + H+

(HISTIDINE-AMMONIA-LYASE-RXN)	L-histidine  ->  ammonium + urocanate

(RXN-20155)	a [protein] S-acetyl-L-cysteine + acetyl-CoA  ->  acetoacetyl-CoA + a [protein]-L-cysteine

(RXN-20155)	acetoacetyl-CoA + a [protein]-L-cysteine  ->  a [protein] S-acetyl-L-cysteine + acetyl-CoA

(ACETALD-DEHYDROG-RXN)	acetaldehyde + coenzyme A + NAD+  ->  acetyl-CoA + NADH + H+

(NADH-DEHYDROG-A-RXN-NADH/UBIQUINONE-8/PROTON//NAD/CPD-9956/PROTON.46. *instantiated*)	NADH[in] + ubiquinone-8[T] + 5 H+[in]  ->  NAD+[in] + ubiquinol-8[T] + 4 H+[out]

(NADH-DEHYDROG-A-RXN-NADH/UBIQUINONE-8/PROTON//NAD/CPD-9956/PROTON.46. *instantiated*)	NAD+[in] + ubiquinol-8[T] + 4 H+[out]  ->  NADH[in] + ubiquinone-8[T] + 5 H+[in]

(RXN-18210)	3-ethylmalate + NAD+  ->  3-ethyl-2-oxosuccinate + NADH + H+

(RXN-18210)	3-ethyl-2-oxosuccinate + NADH + H+  ->  3-ethylmalate + NAD+

(RXN-16955)	3-hydroxy-4-methylbenzaldehyde + NADP+ + H2O  ->  3-hydroxy-4-methyl-benzoate + NADPH + 2 H+

(RXN-16955)	3-hydroxy-4-methyl-benzoate + NADPH + 2 H+  ->  3-hydroxy-4-methylbenzaldehyde + NADP+ + H2O

(RXN-10659)	(9Z)-3-oxohexadec-9-enoyl-[acp] + NADPH + H+  ->  (3R,9Z)-3-hydroxyhexadec-9-enoyl-[acp] + NADP+

(MALYL-COA-LYASE-RXN)	glyoxylate + acetyl-CoA  ->  (S)-malyl-CoA

(RXN-8957)	crotonyl-CoA + CO2 + NADPH  ->  (2S)-ethylmalonyl-CoA + NADP+

(RXN-20457 *spontaneous*)	NADPH + H2O  ->  (S)-NADPHX

(D-2-HYDROXY-ACID-DEHYDROGENASE-RXN)	(R)-lactate + oxidized electron carrier  ->  pyruvate + reduced two electron carrier

(DEHYDDEOXGALACTKIN-RXN)	2-dehydro-3-deoxy-D-galactonate + ATP  ->  2-dehydro-3-deoxy-D-galactonate 6-phosphate + ADP + H+

(DEHYDDEOXGALACTKIN-RXN)	2-dehydro-3-deoxy-D-galactonate 6-phosphate + ADP + H+  ->  2-dehydro-3-deoxy-D-galactonate + ATP

(RXN-14382)	an [L-cysteine desulfurase]-S-sulfanyl-L-cysteine + [ThiI sulfur-carrier protein]-L-cysteine  ->  an [L-cysteine desulfurase]-L-cysteine + [ThiI sulfur-carrier protein]-S-sulfanyl-L-cysteine

(RXN-14382)	an [L-cysteine desulfurase]-L-cysteine + [ThiI sulfur-carrier protein]-S-sulfanyl-L-cysteine  ->  an [L-cysteine desulfurase]-S-sulfanyl-L-cysteine + [ThiI sulfur-carrier protein]-L-cysteine

(RXN0-4401)	NADH + H2O  ->  AMP + reduced beta-nicotinamide D-ribonucleotide + 2 H+

(RXN0-4401)	AMP + reduced beta-nicotinamide D-ribonucleotide + 2 H+  ->  NADH + H2O

(ALCOHOL-DEHYDROG-RXN)	ethanol + NAD+  ->  acetaldehyde + NADH + H+

(ALCOHOL-DEHYDROG-RXN)	acetaldehyde + NADH + H+  ->  ethanol + NAD+

(3-OXOADIPATE-COA-TRANSFERASE-RXN)	succinyl-CoA + 3-oxoadipate  ->  succinate + 3-oxoadipyl-CoA

(GCVT-RXN)	a [glycine-cleavage complex H protein] N6-aminomethyldihydrolipoyl-L-lysine + a tetrahydrofolate  ->  a [glycine-cleavage complex H protein] N6-dihydrolipoyl-L-lysine + a 5,10-methylenetetrahydrofolate + ammonium

(RXN-17919)	adenosine-5'-diphospho-5'-[DNA] + [DNA]-3'-hydroxyl  ->  [DNA] + AMP + H+

(RXN-14278)	hexanoyl-CoA + oxidized electron-transfer flavoprotein + H+  ->  (2E)-hexenoyl-CoA + reduced electron-transfer flavoprotein

(RXN-14229)	octanoyl-CoA + oxidized electron-transfer flavoprotein + H+  ->  (2E)-oct-2-enoyl-CoA + reduced electron-transfer flavoprotein

(RXN-13615)	decanoyl-CoA + oxidized electron-transfer flavoprotein + H+  ->  (2E)-dec-2-enoyl-CoA + reduced electron-transfer flavoprotein

(ACYLCOADEHYDROG-RXN-LAUROYLCOA-CPD/ETF-Oxidized/PROTON//CPD-7222/ETF-Reduced.57. *instantiated*)	lauroyl-CoA + oxidized electron-transfer flavoprotein + H+  ->  (2E)-dodec-2-enoyl-CoA + reduced electron-transfer flavoprotein

(ACYLCOADEHYDROG-RXN-POLYMER-INST-Saturated-Fatty-Acyl-CoA-C10-H20/ETF-Oxidized/PROTON//POLYMER-INST-TRANS-D2-ENOYL-COA-C10-H20/ETF-Reduced.119. *instantiated*)	POLYMER-INST-Saturated-Fatty-Acyl-CoA-C10-H20 + oxidized electron-transfer flavoprotein + H+  ->  POLYMER-INST-TRANS-D2-ENOYL-COA-C10-H20 + reduced electron-transfer flavoprotein

(ACYLCOADEHYDROG-RXN-PALMITYL-COA/ETF-Oxidized/PROTON//POLYMER-INST-TRANS-D2-ENOYL-COA-C12-H24/ETF-Reduced.86. *instantiated*)	palmitoyl-CoA + oxidized electron-transfer flavoprotein + H+  ->  POLYMER-INST-TRANS-D2-ENOYL-COA-C12-H24 + reduced electron-transfer flavoprotein

(ACYLCOADEHYDROG-RXN-STEAROYL-COA/ETF-Oxidized/PROTON//POLYMER-INST-TRANS-D2-ENOYL-COA-C14-H28/ETF-Reduced.86. *instantiated*)	stearoyl-CoA + oxidized electron-transfer flavoprotein + H+  ->  POLYMER-INST-TRANS-D2-ENOYL-COA-C14-H28 + reduced electron-transfer flavoprotein

(METHYLGLUTACONYL-COA-HYDRATASE-RXN)	3-methylglutaconyl-CoA + H2O  ->  (S)-3-hydroxy-3-methylglutaryl-CoA

(RXN-12614)	glycine + dioxygen  ->  2-iminoacetate + hydrogen peroxide + H+

(MALATE-DEHYDROGENASE-NADP+-RXN)	(S)-malate + NADP+  ->  oxaloacetate + NADPH + H+

(MALATE-DEHYDROGENASE-NADP+-RXN)	oxaloacetate + NADPH + H+  ->  (S)-malate + NADP+

(RXN0-882)	2-C-methyl-D-erythritol-2,4-cyclodiphosphate + 2 reduced ferredoxin [iron-sulfur] cluster + H+  ->  (E)-4-hydroxy-3-methylbut-2-en-1-yl diphosphate + 2 oxidized ferredoxin [iron-sulfur] cluster + H2O

(GUANPRIBOSYLTRAN-RXN)	guanine + 5-phospho-alpha-D-ribose 1-diphosphate  ->  GMP + diphosphate

(GLURS-RXN)	tRNAglu + L-glutamate + ATP  ->  L-glutamyl-[tRNAGlu] + AMP + diphosphate

(RXN-21316)	5-deoxy-D-ribulose 1-phosphate  ->  glycerone phosphate + acetaldehyde

(RXN-21316)	glycerone phosphate + acetaldehyde  ->  5-deoxy-D-ribulose 1-phosphate

(RXN-12721)	adenosine 5'-phosphoselenate + 2 glutathione  ->  selenite + AMP + glutathione disulfide + 2 H+

(RXN-12611)	2-[(2R,5Z)-2-carboxy-4-methylthiazol-5(2H)-ylidene]ethyl phosphate + 4-amino-2-methyl-5-(diphosphooxymethyl)pyrimidine + 2 H+  ->  thiamine phosphate + CO2 + diphosphate

(XANTHOSINEPHOSPHORY-RXN)	xanthosine + phosphate  ->  xanthine + alpha-D-ribose-1-phosphate

(XANTHOSINEPHOSPHORY-RXN)	xanthine + alpha-D-ribose-1-phosphate  ->  xanthosine + phosphate

(RXN-15635)	3-methyl-2-oxobutanoate + 5,10-methylene-tetrahydromethanopterin + H2O  ->  2-dehydropantoate + tetrahydromethanopterin

(RXN-15635)	2-dehydropantoate + tetrahydromethanopterin  ->  3-methyl-2-oxobutanoate + 5,10-methylene-tetrahydromethanopterin + H2O

(RXN-9633)	3-oxooctadecanoyl-[acp] + NADPH + H+  ->  (3R)-3-hydroxyoctadecanoyl-[acp] + NADP+

(3-OXOACYL-ACP-REDUCT-RXN-POLYMER-INST-OH-ACYL-ACP-C12-H24/NADP//3-oxo-palmitoyl-ACPs/NADPH/PROTON.73. *instantiated*)	3-oxohexadecanoyl-[acp] + NADPH + H+  ->  POLYMER-INST-OH-ACYL-ACP-C12-H24 + NADP+

(RXN-9532)	3-oxododecanoyl-[acp] + NADPH + H+  ->  (3R)-3-hydroxydodecanoyl-[acp] + NADP+

(RXN-9518)	3-oxohexanoyl-[acp] + NADPH + H+  ->  (3R)-3-hydroxyhexanoyl-[acp] + NADP+

(RXN-9514)	acetoacetyl-[acp] + NADPH + H+  ->  (3R)-3-hydroxybutanoyl-[acp] + NADP+

(SPONTPRO-RXN *spontaneous*)	(S)-1-pyrroline-5-carboxylate + H+ + H2O  ->  L-glutamate-5-semialdehyde

(RXN-8315 *spontaneous*)	sulfite + H+  ->  hydrogensulfite

(RXN-8315 *spontaneous*)	hydrogensulfite  ->  sulfite + H+

(PHENYLPYRUVATE-TAUTOMERASE-RXN *spontaneous*)	3-phenyl-2-oxopropanoate  ->  enol-phenylpyruvate

(PHENYLPYRUVATE-TAUTOMERASE-RXN *spontaneous*)	enol-phenylpyruvate  ->  3-phenyl-2-oxopropanoate

(RXN-22915)	(R)-4'-phosphopantothenoyl-cytidylate + L-cysteine  ->  (R)-4'-phosphopantothenoyl-L-cysteine + CMP + 2 H+

(GAPOXNPHOSPHN-RXN)	3-phospho-D-glyceroyl phosphate + NADH + H+  ->  D-glyceraldehyde 3-phosphate + NAD+ + phosphate

(3-HYDROXYBUTYRYL-COA-DEHYDROGENASE-RXN)	acetoacetyl-CoA + NADPH + H+  ->  (S)-3-hydroxybutanoyl-CoA + NADP+

(RXN-12625)	N,N'-diacetylchitobiose + H2O  ->  2 N-acetyl-D-glucosamine

(4.2.1.61-RXN)	(3R)-3-hydroxyhexadecanoyl-[acp]  ->  (2E)-hexadec-2-enoyl-[acp] + H2O

(RXN-22706)	N,N-dimethylglycine + a tetrahydrofolate + dioxygen  ->  sarcosine + a 5,10-methylenetetrahydrofolate + hydrogen peroxide

(3-OH-BENZALDEHYDE-DEHYDROG-NADP+-RXN)	3-hydroxybenzaldehyde + NADP+ + H2O  ->  3-hydroxybenzoate + NADPH + 2 H+

(3-OH-BENZALDEHYDE-DEHYDROG-NADP+-RXN)	3-hydroxybenzoate + NADPH + 2 H+  ->  3-hydroxybenzaldehyde + NADP+ + H2O

(RXN-9527)	octanoyl-[acp] + malonyl-[acp] + H+  ->  3-oxodecanoyl-[acp] + CO2 + acyl-carrier protein

(RXN-9535)	dodecanoyl-[acp] + malonyl-[acp] + H+  ->  3-oxotetradecanoyl-[acp] + CO2 + acyl-carrier protein

(RXN-9632)	palmitoyl-[acp] + malonyl-[acp] + H+  ->  3-oxooctadecanoyl-[acp] + CO2 + acyl-carrier protein

(1-ACYLGLYCEROL-3-P-ACYLTRANSFER-RXN)	acyl-[acyl-carrier protein] + a 2-lysophosphatidate  ->  a phosphatidate + acyl-carrier protein

(SAMDECARB-RXN)	S-adenosyl-L-methionine + H+  ->  CO2 + S-adenosyl 3-(methylsulfanyl)propylamine

(RXN-12196)	UTP + H2O  ->  UDP + phosphate + H+

(RXN-10814)	L-phenylalanine + 2-oxoglutarate  ->  3-phenyl-2-oxopropanoate + L-glutamate

(PHOSMANMUT-RXN)	alpha-D-mannose 1-phosphate  ->  D-mannopyranose 6-phosphate

(PHOSMANMUT-RXN)	D-mannopyranose 6-phosphate  ->  alpha-D-mannose 1-phosphate

(DIHYDROPYRIMIDINASE-RXN)	5,6-dihydrouracil + H2O  ->  3-ureidopropanoate + H+

(NAD-KIN-RXN)	ATP + NAD+  ->  ADP + NADP+ + H+

(4.1.1.81-RXN)	L-threonine 3-O-phosphate + H+  ->  (R)-1-amino-2-propanol O-2-phosphate + CO2

(4.1.1.81-RXN)	(R)-1-amino-2-propanol O-2-phosphate + CO2  ->  L-threonine 3-O-phosphate + H+

(AMP-DEPHOSPHORYLATION-RXN)	AMP + H2O  ->  adenosine + phosphate

(RXN-12195)	CTP + H2O  ->  CDP + phosphate + H+

(RXN-15714)	5-aminolevulinate + coenzyme A + ATP  ->  5-aminolevulinyl-CoA + AMP + diphosphate

(RXN-15714)	5-aminolevulinyl-CoA + AMP + diphosphate  ->  5-aminolevulinate + coenzyme A + ATP

(GLYOHMETRANS-RXN)	L-serine + a tetrahydrofolate  ->  glycine + a 5,10-methylenetetrahydrofolate + H2O

(P-PANTOCYSDECARB-RXN)	(R)-4'-phosphopantothenoyl-L-cysteine + H+  ->  CO2 + 4'-phosphopantetheine

(RXN-14971-SUC/UBIQUINONE-8//FUM/CPD-9956.31. *instantiated*)	succinate[in] + ubiquinone-8[generic membrane]  ->  fumarate[in] + ubiquinol-8[generic membrane]

(RXN-14971-SUC/UBIQUINONE-8//FUM/CPD-9956.31. *instantiated*)	fumarate[in] + ubiquinol-8[generic membrane]  ->  succinate[in] + ubiquinone-8[generic membrane]

(GLUTATHIONE-REDUCT-NADPH-RXN)	glutathione disulfide + NADPH + H+  ->  2 glutathione + NADP+

(RXN-14025)	UMP + H2O  ->  uridine + phosphate

(RXN-14858)	oxythiamine + H2O  ->  5-(2-hydroxyethyl)-4-methylthiazole + 5-(hydroxymethyl)-2-methyl-4(1H)-pyrimidinone + H+

(RXN-14858)	5-(2-hydroxyethyl)-4-methylthiazole + 5-(hydroxymethyl)-2-methyl-4(1H)-pyrimidinone + H+  ->  oxythiamine + H2O

(ALDOSE-1-EPIMERASE-RXN *spontaneous*)	alpha-D-glucopyranose  ->  beta-D-glucopyranose

(ALDOSE-1-EPIMERASE-RXN *spontaneous*)	beta-D-glucopyranose  ->  alpha-D-glucopyranose

(RXN-16804 *spontaneous*)	aldehydo-D-arabinose 5-phosphate  ->  D-arabinofuranose 5-phosphate

(RXN-16804 *spontaneous*)	D-arabinofuranose 5-phosphate  ->  aldehydo-D-arabinose 5-phosphate

(CDPKIN-RXN)	ATP + CDP  ->  ADP + CTP

(RXN-22914)	(R)-4'-phosphopantothenate + CTP + H+  ->  (R)-4'-phosphopantothenoyl-cytidylate + diphosphate

(RXN-8665)	L-tryptophan + dioxygen  ->  N-Formyl-L-kynurenine

(ORNITHINE-CYCLODEAMINASE-RXN)	L-ornithine  ->  L-proline + ammonium

(RXN-9516)	butanoyl-[acp] + malonyl-[acp] + H+  ->  3-oxohexanoyl-[acp] + CO2 + acyl-carrier protein

(RXN-7716)	a [2-oxoglutarate dehydrogenase E2 protein] N6-dihydrolipoyl-L-lysine + NAD+  ->  a [2-oxoglutarate dehydrogenase E2 protein] N6-lipoyl-L-lysine + NADH + H+

(RXN-22713)	chorismate + ammonia + H+  ->  4-amino-4-deoxychorismate + H2O

(RIBONUCLEOSIDE-DIP-REDUCTII-RXN)	dCDP + oxidized NrdH glutaredoxin-like protein + H2O  ->  CDP + reduced NrdH glutaredoxin-like protein

(RIBONUCLEOSIDE-DIP-REDUCTII-RXN)	CDP + reduced NrdH glutaredoxin-like protein  ->  dCDP + oxidized NrdH glutaredoxin-like protein + H2O

(ASPARAGINE--TRNA-LIGASE-RXN)	tRNAAsn + L-asparagine + ATP  ->  L-asparaginyl-[tRNAAsn] + AMP + diphosphate

(RXN0-6565)	5,6-dihydrothymine + NAD+  ->  thymine + NADH + H+

(RXN0-6565)	thymine + NADH + H+  ->  5,6-dihydrothymine + NAD+

(ADENPRIBOSYLTRAN-RXN)	adenine + 5-phospho-alpha-D-ribose 1-diphosphate  ->  AMP + diphosphate

(AIRCARBOXY-RXN)	5-amino-1-(5-phospho-D-ribosyl)imidazole-4-carboxylate + 2 H+  ->  5-amino-1-(5-phospho-beta-D-ribosyl)imidazole + CO2

(PEPDEPHOS-RXN)	phosphoenolpyruvate + ADP + H+  ->  pyruvate + ATP

(RXN-18200)	3-[(7'-methylsulfanyl)heptyl]malate + NAD+  ->  3-carboxy-10-(methylsulfanyl)-2-oxodecanoate + NADH + H+

(RXN-18200)	3-carboxy-10-(methylsulfanyl)-2-oxodecanoate + NADH + H+  ->  3-[(7'-methylsulfanyl)heptyl]malate + NAD+

(ADENPHOSPHOR-RXN)	adenosine + phosphate  ->  adenine + alpha-D-ribose-1-phosphate

(ADENPHOSPHOR-RXN)	adenine + alpha-D-ribose-1-phosphate  ->  adenosine + phosphate

(CARDIOLIPSYN-RXN)	2 a phosphatidylglycerol  ->  a cardiolipin + glycerol

(RXN-10811)	coenzyme A + H2O  ->  adenosine 3',5'-bisphosphate + 4'-phosphopantetheine + 2 H+

(RXN-10086)	(2S)-methylsuccinyl-CoA + acetyl-CoA  ->  4-methyl-3-oxoadipyl-CoA + coenzyme A

(RXN-10086)	4-methyl-3-oxoadipyl-CoA + coenzyme A  ->  (2S)-methylsuccinyl-CoA + acetyl-CoA

(GLUTATHIONE-SYN-RXN)	glycine + gamma-L-glutamyl-L-cysteine + ATP  ->  ADP + glutathione + phosphate + H+

(RXN-13697)	L-aspartate + 2-oxoglutarate  ->  oxaloacetate + L-glutamate

(RXN-13697)	oxaloacetate + L-glutamate  ->  L-aspartate + 2-oxoglutarate

(GLYC3PDEHYDROGBIOSYN-RXN-GLYCEROL-3P/NAD//DIHYDROXY-ACETONE-PHOSPHATE/NADH/PROTON.57. *instantiated*)	glycerone phosphate + NADH + H+  ->  sn-glycerol 3-phosphate + NAD+

(GLYC3PDEHYDROGBIOSYN-RXN-GLYCEROL-3P/NADP//DIHYDROXY-ACETONE-PHOSPHATE/NADPH/PROTON.59. *instantiated*)	glycerone phosphate + NADPH + H+  ->  sn-glycerol 3-phosphate + NADP+

(RXN-13908)	4a-hydroxy-N10-formyltetrahydrofolate  ->  an N10-formyl-7,8-dihydrofolate + H2O

(RXN-13908)	an N10-formyl-7,8-dihydrofolate + H2O  ->  4a-hydroxy-N10-formyltetrahydrofolate

(TRIOSEPISOMERIZATION-RXN)	D-glyceraldehyde 3-phosphate  ->  glycerone phosphate

(TYROSINE-AMINOTRANSFERASE-RXN)	L-tyrosine + 2-oxoglutarate  ->  3-(4-hydroxyphenyl)pyruvate + L-glutamate

(RXN-8667)	hydrogen peroxide + reduced two electron carrier  ->  oxidized electron carrier + 2 H2O

(RXN-8667)	oxidized electron carrier + 2 H2O  ->  hydrogen peroxide + reduced two electron carrier

(XMPXAN-RXN)	XMP + H2O  ->  xanthosine + phosphate

(RXN-8642)	2-oxoglutarate + CO2  ->  oxalosuccinate + H+

(RXN-22623)	benzyl alcohol + NADP+  ->  benzaldehyde + NADPH + H+

(RXN-22623)	benzaldehyde + NADPH + H+  ->  benzyl alcohol + NADP+

(GLUCOKIN-RXN-ALPHA-GLUCOSE/ATP//ALPHA-GLC-6-P/ADP/PROTON.44. *instantiated*)	alpha-D-glucopyranose + ATP  ->  alpha-D-glucose 6-phosphate + ADP + H+

(GLUCOKIN-RXN-ALPHA-GLUCOSE/ATP//GLC-6-P/ADP/PROTON.38. *instantiated*)	alpha-D-glucopyranose + ATP  ->  beta-D-glucose 6-phosphate + ADP + H+

(GLUCOKIN-RXN-ALPHA-GLUCOSE/ATP//D-glucopyranose-6-phosphate/ADP/PROTON.58. *instantiated*)	alpha-D-glucopyranose + ATP  ->  D-glucopyranose 6-phosphate + ADP + H+

(GLUCOKIN-RXN-GLC/ATP//ALPHA-GLC-6-P/ADP/PROTON.34. *instantiated*)	beta-D-glucopyranose + ATP  ->  alpha-D-glucose 6-phosphate + ADP + H+

(GLUCOKIN-RXN-GLC/ATP//GLC-6-P/ADP/PROTON.28. *instantiated*)	beta-D-glucopyranose + ATP  ->  beta-D-glucose 6-phosphate + ADP + H+

(GLUCOKIN-RXN-GLC/ATP//D-glucopyranose-6-phosphate/ADP/PROTON.48. *instantiated*)	beta-D-glucopyranose + ATP  ->  D-glucopyranose 6-phosphate + ADP + H+

(RXN-15829-Oxidized-CycA1-cytochromes/CPD-9956//Reduced-CycA1-cytochromes/UBIQUINONE-8/PROTON.83. *instantiated*)	2 oxidized CycA1 cytochrome[out] + ubiquinol-8[T]  ->  2 reduced CycA1 cytochrome[out] + ubiquinone-8[T] + 2 H+[in]

(RXN-15829-Oxidized-CycA1-cytochromes/CPD-9956//Reduced-CycA1-cytochromes/UBIQUINONE-8/PROTON.83. *instantiated*)	2 reduced CycA1 cytochrome[out] + ubiquinone-8[T] + 2 H+[in]  ->  2 oxidized CycA1 cytochrome[out] + ubiquinol-8[T]

(RXN-20895)	L-leucine  ->  D-leucine

(RXN-20895)	D-leucine  ->  L-leucine

(ASPARAGHYD-RXN)	L-asparagine + H2O  ->  L-aspartate + ammonium

(RXN-10656)	(3R,7Z)-3-hydroxytetradec-7-enoyl-[acp]  ->  (2E,7Z)-tetradeca-2,7-dienoyl-[acp] + H2O

(RXN-14056 *spontaneous*)	(R)-methylmalonate-semialdehyde  ->  (S)-methylmalonate-semialdehyde

(RXN-14056 *spontaneous*)	(S)-methylmalonate-semialdehyde  ->  (R)-methylmalonate-semialdehyde

(DEOXYGUANPHOSPHOR-RXN)	2'-deoxyguanosine + phosphate  ->  guanine + 2-deoxy-alpha-D-ribose 1-phosphate

(DEOXYGUANPHOSPHOR-RXN)	guanine + 2-deoxy-alpha-D-ribose 1-phosphate  ->  2'-deoxyguanosine + phosphate

(RXN-15682)	3-benzyl-3,6 -bis(cysteinylglycine)- 6-(hydroxymethyl)-diketopiperazine + 2 H2O  ->  3-benzyl-3,6 -bis(cysteinyl)- 6-(hydroxymethyl)-diketopiperazine + 2 glycine

(RXN-15682)	3-benzyl-3,6 -bis(cysteinyl)- 6-(hydroxymethyl)-diketopiperazine + 2 glycine  ->  3-benzyl-3,6 -bis(cysteinylglycine)- 6-(hydroxymethyl)-diketopiperazine + 2 H2O

(3.5.1.88-RXN)	a [protein] N-terminal-formyl-L-methionine + H2O  ->  an N-terminal L-methionyl-[protein] + formate

(3.5.1.88-RXN)	an N-terminal L-methionyl-[protein] + formate  ->  a [protein] N-terminal-formyl-L-methionine + H2O

(RXN0-5055)	acetyl-CoA + a [carboxyl-carrier protein dimer]-N6-carboxybiotinyl-L-lysine  ->  malonyl-CoA + a [biotin carboxyl-carrier-protein dimer]-N6-biotinyl-L-lysine

(RXN-9922 *spontaneous*)	3-hydroxy-cis,cis-muconate  ->  2-maleylacetate

(RXN-9922 *spontaneous*)	2-maleylacetate  ->  3-hydroxy-cis,cis-muconate

(RXN-7744)	citraconate + H2O  ->  (2R,3S)-3-methylmalate

(RXN0-6427)	pppGpp + H2O  ->  GTP + diphosphate

(PGLYCDEHYDROG-RXN)	3-phospho-D-glycerate + NAD+  ->  3-phosphooxypyruvate + NADH + H+

(PGLYCDEHYDROG-RXN)	3-phosphooxypyruvate + NADH + H+  ->  3-phospho-D-glycerate + NAD+

(METHANETHIOL-OXIDASE-RXN)	methanethiol + dioxygen + H2O  ->  formaldehyde + hydrogen sulfide + hydrogen peroxide

(RXN-11209)	thymine + NADPH + H+  ->  5,6-dihydrothymine + NADP+

(RXN-19381)	D-cysteine  ->  hydrogen sulfide + 2-aminoprop-2-enoate

(RXN-12541 *spontaneous*)	2 Fe2+ + 2 dioxygen  ->  2 superoxide + 2 Fe3+

(ISOVALERYLCOA-DHLIPOAMIDE-RXN)	3-methylbutanoyl-CoA + an [apo BCAA dehydrogenase E2 protein] N6-dihydrolipoyl-L-lysine  ->  a [apo BCAA dehydrogenase E2 protein] N6-S-[3-methylbutanoyl]dihydrolipoyl-L-lysine + coenzyme A

(ISOVALERYLCOA-DHLIPOAMIDE-RXN)	a [apo BCAA dehydrogenase E2 protein] N6-S-[3-methylbutanoyl]dihydrolipoyl-L-lysine + coenzyme A  ->  3-methylbutanoyl-CoA + an [apo BCAA dehydrogenase E2 protein] N6-dihydrolipoyl-L-lysine

(1.8.4.8-RXN)	3'-phosphoadenylyl-sulfate + reduced thioredoxin  ->  adenosine 3',5'-bisphosphate + sulfite + oxidized thioredoxin + 2 H+

(RXN-14812-FRUCTOSE-6P//FRUCTOSE-6P.25. *instantiated* *spontaneous*)	beta-D-fructofuranose 6-phosphate  ->  beta-D-fructofuranose 6-phosphate

(RXN-14812-FRUCTOSE-6P//FRUCTOSE-6P.25. *instantiated* *spontaneous*)	beta-D-fructofuranose 6-phosphate  ->  beta-D-fructofuranose 6-phosphate

(RXN-14812 *spontaneous*)	keto-D-fructose 6-phosphate  ->  beta-D-fructofuranose 6-phosphate

(RXN-14812 *spontaneous*)	beta-D-fructofuranose 6-phosphate  ->  keto-D-fructose 6-phosphate

(THIOESTER-RXN-STEAROYL-COA/WATER//STEARIC_ACID/CO-A/PROTON.45. *instantiated*)	stearoyl-CoA + H2O  ->  stearate + coenzyme A + H+

(THIOESTER-RXN-BUTYRYL-COA/WATER//BUTYRIC_ACID/CO-A/PROTON.44. *instantiated*)	butanoyl-CoA + H2O  ->  butanoate + coenzyme A + H+

(THIOESTER-RXN-PALMITYL-COA/WATER//PALMITATE/CO-A/PROTON.42. *instantiated*)	palmitoyl-CoA + H2O  ->  palmitate + coenzyme A + H+

(THIOESTER-RXN-PROPIONYL-COA/WATER//PROPIONATE/CO-A/PROTON.44. *instantiated*)	propanoyl-CoA + H2O  ->  propanoate + coenzyme A + H+

(THIOESTER-RXN-LAUROYLCOA-CPD/WATER//DODECANOATE/CO-A/PROTON.46. *instantiated*)	lauroyl-CoA + H2O  ->  laurate + coenzyme A + H+

(THIOESTER-RXN-OLEOYL-COA/WATER//OLEATE-CPD/CO-A/PROTON.41. *instantiated*)	oleoyl-CoA + H2O  ->  oleate + coenzyme A + H+

(THIOESTER-RXN-CPD-19144/WATER//CPD-9245/CO-A/PROTON.38. *instantiated*)	(7Z)-hexadecenoyl-CoA + H2O  ->  palmitoleate + coenzyme A + H+

(THIOESTER-RXN-CPD-196/WATER//CPD-195/CO-A/PROTON.35. *instantiated*)	octanoyl-CoA + H2O  ->  octanoate + coenzyme A + H+

(RXN-14904 *spontaneous*)	alpha-D-ribofuranose  ->  beta-D-ribofuranose

(RXN-14904 *spontaneous*)	beta-D-ribofuranose  ->  alpha-D-ribofuranose

(RXN-17884)	S-nitrosoglutathione + NADH + H+  ->  S-(hydroxysulfenamide)glutathione + NAD+

(PANTOTHENATE-KIN-RXN)	(R)-pantothenate + ATP  ->  (R)-4'-phosphopantothenate + ADP + H+

(CDPDIGLYSYN-RXN)	CTP + a phosphatidate + H+  ->  a CDP-diacylglycerol + diphosphate

(PYRUVATE-CARBOXYLASE-RXN)	pyruvate + hydrogencarbonate + ATP  ->  oxaloacetate + ADP + phosphate + H+

(PEPCARBOX-RXN)	phosphoenolpyruvate + hydrogencarbonate  ->  oxaloacetate + phosphate

(RXN-16998)	alpha-D-glucose 1,6-bisphosphate + phosphoglucomutase  ->  D-glucopyranose 6-phosphate + phosphorylated phosphoglucomutase

(RXN-16998-ALPHA-GLUCOSE-16-BISPHOSPHATE/Phosphoglucomutase//GLC-6-P/Phosphorylated-phosphoglucomutase.92. *instantiated*)	alpha-D-glucose 1,6-bisphosphate + phosphoglucomutase  ->  beta-D-glucose 6-phosphate + phosphorylated phosphoglucomutase

(RXN-16998-ALPHA-GLUCOSE-16-BISPHOSPHATE/Phosphoglucomutase//GLC-6-P/Phosphorylated-phosphoglucomutase.92. *instantiated*)	beta-D-glucose 6-phosphate + phosphorylated phosphoglucomutase  ->  alpha-D-glucose 1,6-bisphosphate + phosphoglucomutase

(RXN-16998-ALPHA-GLUCOSE-16-BISPHOSPHATE/Phosphoglucomutase//ALPHA-GLC-6-P/Phosphorylated-phosphoglucomutase.98. *instantiated*)	alpha-D-glucose 6-phosphate + phosphorylated phosphoglucomutase  ->  alpha-D-glucose 1,6-bisphosphate + phosphoglucomutase

(RXN-15122)	L-threonine  ->  (2Z)-2-aminobut-2-enoate + H2O + H+

(D-ALANINE-AMINOTRANSFERASE-RXN)	2-oxoglutarate + D-alanine  ->  D-glutamate + pyruvate

(D-ALANINE-AMINOTRANSFERASE-RXN)	D-glutamate + pyruvate  ->  2-oxoglutarate + D-alanine

(2-DEHYDROPANTOATE-REDUCT-RXN)	2-dehydropantoate + NADPH + H+  ->  (R)-pantoate + NADP+

(RXN-12456)	a 5,6-dihydrouracil20 in tRNA + NAD(P)+  ->  a uracil20 in tRNA + NAD(P)H + H+

(RXN-12456)	a uracil20 in tRNA + NAD(P)H + H+  ->  a 5,6-dihydrouracil20 in tRNA + NAD(P)+

(RXN-10660)	(3R,9Z)-3-hydroxyhexadec-9-enoyl-[acp]  ->  (2E,9Z)-hexadeca-2,9-dienoyl-[acp] + H2O

(2.3.1.168-RXN)	isobutanoyl-CoA + an [apo BCAA dehydrogenase E2 protein] N6-dihydrolipoyl-L-lysine  ->  an [apo BCAA dehydrogenase E2 protein] N6-S-[2-methylpropanoyl]dihydrolipoyl-L-lysine + coenzyme A

(2.3.1.168-RXN)	an [apo BCAA dehydrogenase E2 protein] N6-S-[2-methylpropanoyl]dihydrolipoyl-L-lysine + coenzyme A  ->  isobutanoyl-CoA + an [apo BCAA dehydrogenase E2 protein] N6-dihydrolipoyl-L-lysine

(RXN0-7012)	a phosphatidylethanolamine + a phosphatidylglycerol  ->  a cardiolipin + ethanolamine

(RXN-8668)	a [protein]-L-methionine + oxidized thioredoxin + H2O  ->  a protein-L-methionine-(S)-S-oxide + reduced thioredoxin

(RXN-8668)	a protein-L-methionine-(S)-S-oxide + reduced thioredoxin  ->  a [protein]-L-methionine + oxidized thioredoxin + H2O

(PROPIONYL-COA-CARBOXY-RXN)	ATP + propanoyl-CoA + hydrogencarbonate  ->  (S)-methylmalonyl-CoA + ADP + phosphate + H+

(PROPIONYL-COA-CARBOXY-RXN)	(S)-methylmalonyl-CoA + ADP + phosphate + H+  ->  ATP + propanoyl-CoA + hydrogencarbonate

(1.3.1.2-RXN)	uracil + NADPH + H+  ->  5,6-dihydrouracil + NADP+

(RXN-9549)	palmitoyl-[acp] + H2O  ->  palmitate + acyl-carrier protein + H+

(R125-RXN)	(S)-5-amino-3-oxohexanoate + acetyl-CoA  ->  (S)-3-aminobutanoyl-CoA + acetoacetate

(R125-RXN)	(S)-3-aminobutanoyl-CoA + acetoacetate  ->  (S)-5-amino-3-oxohexanoate + acetyl-CoA

(SUCCCOASYN-RXN)	succinyl-CoA + ADP + phosphate  ->  succinate + coenzyme A + ATP

(RXN-22715)	D-erythrose 4-phosphate  ->  D-threose 4-phosphate

(RXN-22715)	D-threose 4-phosphate  ->  D-erythrose 4-phosphate

(METHYLCROTONYL-COA-CARBOXYLASE-RXN)	3-methylcrotonyl-CoA + hydrogencarbonate + ATP  ->  3-methylglutaconyl-CoA + ADP + phosphate + H+

(RXN-15200)	L-tyrosine + 3-phenyl-2-oxopropanoate  ->  3-(4-hydroxyphenyl)pyruvate + L-phenylalanine

(RXN-15200)	3-(4-hydroxyphenyl)pyruvate + L-phenylalanine  ->  L-tyrosine + 3-phenyl-2-oxopropanoate

(GSADENYLATION-RXN)	a [glutamine-synthetase]-L-tyrosine + ATP  ->  a [glutamine synthetase]-O4-(5'-adenylyl)-L-tyrosine + diphosphate

(GSADENYLATION-RXN)	a [glutamine synthetase]-O4-(5'-adenylyl)-L-tyrosine + diphosphate  ->  a [glutamine-synthetase]-L-tyrosine + ATP

(R101-RXN)	L-2,4-diaminobutanoate + 2-oxoglutarate  ->  L-aspartate 4-semialdehyde + L-glutamate

(R101-RXN)	L-aspartate 4-semialdehyde + L-glutamate  ->  L-2,4-diaminobutanoate + 2-oxoglutarate

(RXN-15125)	L-serine  ->  2-aminoprop-2-enoate + H2O

(RXN-8960)	2-methylfumaryl-CoA + H2O  ->  (2R,3S)-beta-methylmalyl-CoA

(DUDPKIN-RXN)	ATP + dUDP  ->  ADP + dUTP

(BUTYRYL-COA-DEHYDROGENASE-RXN)	butanoyl-CoA + oxidized electron-transfer flavoprotein + H+  ->  crotonyl-CoA + reduced electron-transfer flavoprotein

(BUTYRYL-COA-DEHYDROGENASE-RXN)	crotonyl-CoA + reduced electron-transfer flavoprotein  ->  butanoyl-CoA + oxidized electron-transfer flavoprotein + H+

(RXN-10661)	(2E,9Z)-hexadeca-2,9-dienoyl-[acp] + NADH + H+  ->  palmitoleoyl-[acp] + NAD+

(RXN-18427)	L-aspartate 4-semialdehyde + NAD+ + H2O  ->  L-aspartate + NADH + 2 H+

(RXN-10658)	(7Z)-tetradec-7-enoyl-[acp] + malonyl-[acp] + H+  ->  (9Z)-3-oxohexadec-9-enoyl-[acp] + acyl-carrier protein + CO2

(GLUTRNAREDUCT-RXN)	L-glutamyl-[tRNAGlu] + NADPH + H+  ->  (S)-4-amino-5-oxopentanoate + tRNAglu + NADP+

(BENZALDEHYDE-DEHYDROGENASE-NADP+-RXN)	benzaldehyde + NADP+ + H2O  ->  benzoate + NADPH + 2 H+

(BENZALDEHYDE-DEHYDROGENASE-NADP+-RXN)	benzoate + NADPH + 2 H+  ->  benzaldehyde + NADP+ + H2O

(RXN-16956)	3-hydroxy-5-methylbenzaldehyde + NADP+ + H2O  ->  3-hydroxy-5-methyl-benzoate + NADPH + 2 H+

(RXN-16956)	3-hydroxy-5-methyl-benzoate + NADPH + 2 H+  ->  3-hydroxy-5-methylbenzaldehyde + NADP+ + H2O

(RXN-4821 *spontaneous*)	(S)-2,3,4,5-tetrahydrodipicolinate + H+ + H2O  ->  L-alpha-amino-epsilon-keto-pimelate

(RXN-4821 *spontaneous*)	L-alpha-amino-epsilon-keto-pimelate  ->  (S)-2,3,4,5-tetrahydrodipicolinate + H+ + H2O

(RXN1G01-46)	L-glutamine  ->  D-glutamine

(RXN1G01-46)	D-glutamine  ->  L-glutamine

(METHENYLTHFCYCLOHYDRO-RXN)	a 5,10-methenyltetrahydrofolate + H2O  ->  an N10-formyltetrahydrofolate + H+

(XYLISOM-RXN)	alpha-D-xylopyranose  ->  alpha-D-xylulofuranose

(XYLISOM-RXN)	alpha-D-xylulofuranose  ->  alpha-D-xylopyranose

(PROPKIN-RXN)	propanoate + ATP  ->  propanoyl phosphate + ADP

(PROPKIN-RXN)	propanoyl phosphate + ADP  ->  propanoate + ATP

(GLUTAMATE-N-ACETYLTRANSFERASE-RXN)	L-glutamate + N-acetyl-L-ornithine  ->  N-acetyl-L-glutamate + L-ornithine

(CHMS-DEHYDROGENASE-RXN)	2-hydroxy-5-carboxymethylmuconate semialdehyde + NAD+ + H2O  ->  5-carboxymethyl-2-hydroxymuconate + NADH + 2 H+

(CHMS-DEHYDROGENASE-RXN)	5-carboxymethyl-2-hydroxymuconate + NADH + 2 H+  ->  2-hydroxy-5-carboxymethylmuconate semialdehyde + NAD+ + H2O

(DCTP-PYROPHOSPHATASE-RXN)	dCTP + H2O  ->  dCMP + diphosphate + H+

(3-ISOPROPYLMALISOM-RXN)	2-isopropylmaleate + H2O  ->  (2S)-2-isopropylmalate

(CYTOCHROME-C-PEROXIDASE-RXN-HYDROGEN-PEROXIDE/Reduced-cytochromes-c553/PROTON//Oxidized-cytochromes-c553/WATER.83. *instantiated*)	hydrogen peroxide + 2 reduced cytochrome c-553 + 2 H+  ->  2 oxidized cytochrome c-553 + 2 H2O

(CYTOCHROME-C-PEROXIDASE-RXN-HYDROGEN-PEROXIDE/Reduced-CycA1-cytochromes/PROTON//Oxidized-CycA1-cytochromes/WATER.85. *instantiated*)	hydrogen peroxide + 2 reduced CycA1 cytochrome + 2 H+  ->  2 oxidized CycA1 cytochrome + 2 H2O

(CYTOCHROME-C-PEROXIDASE-RXN-HYDROGEN-PEROXIDE/a-reduced-NrfB-protein/PROTON//an-oxidized-NrfB-protein/WATER.80. *instantiated*)	hydrogen peroxide + 2 reduced [NrfB protein] + 2 H+  ->  2 oxidized [NrfB protein] + 2 H2O

(CYTOCHROME-C-PEROXIDASE-RXN-HYDROGEN-PEROXIDE/Reduced-Cytochromes-C6/PROTON//Oxidized-Cytochromes-C6/WATER.79. *instantiated*)	hydrogen peroxide + 2 reduced cytochrome c6 + 2 H+  ->  2 oxidized cytochrome c6 + 2 H2O

(CYTOCHROME-C-PEROXIDASE-RXN-HYDROGEN-PEROXIDE/Reduced-NapC-proteins/PROTON//Oxidized-NapC-proteins/WATER.77. *instantiated*)	hydrogen peroxide + 2 reduced [NapC protein] + 2 H+  ->  2 oxidized [NapC protein] + 2 H2O

(CYTOCHROME-C-PEROXIDASE-RXN-HYDROGEN-PEROXIDE/Reduced-cytochromes-c551/PROTON//Oxidized-cytochromes-c551/WATER.83. *instantiated*)	hydrogen peroxide + 2 reduced cytochrome c551 + 2 H+  ->  2 oxidized cytochrome c551 + 2 H2O

(CYTOCHROME-C-PEROXIDASE-RXN-HYDROGEN-PEROXIDE/Reduced-cytochromes-C4/PROTON//Oxidized-cytochromes-C4/WATER.79. *instantiated*)	hydrogen peroxide + 2 reduced cytochrome c4 + 2 H+  ->  2 oxidized cytochrome c4 + 2 H2O

(CYTOCHROME-C-PEROXIDASE-RXN)	hydrogen peroxide + 2 reduced c-type cytochrome + 2 H+  ->  2 oxidized c-type cytochrome + 2 H2O

(4-HYDROXYPROLINE-EPIMERASE-RXN)	trans-4-hydroxy-L-proline  ->  (4R)-4-hydroxy-D-proline

(4-HYDROXYPROLINE-EPIMERASE-RXN)	(4R)-4-hydroxy-D-proline  ->  trans-4-hydroxy-L-proline

(RXN-9644)	oleate + coenzyme A + ATP  ->  oleoyl-CoA + AMP + diphosphate

(RXN-7904-CPD-9245/CO-A/ATP//CPD-19144/AMP/PPI.37. *instantiated*)	palmitoleate + coenzyme A + ATP  ->  (7Z)-hexadecenoyl-CoA + AMP + diphosphate

(1TRANSKETO-RXN)	D-sedoheptulose 7-phosphate + D-glyceraldehyde 3-phosphate  ->  D-ribose 5-phosphate + D-xylulose 5-phosphate

(1TRANSKETO-RXN)	D-ribose 5-phosphate + D-xylulose 5-phosphate  ->  D-sedoheptulose 7-phosphate + D-glyceraldehyde 3-phosphate

(MALATE--COA-LIGASE-RXN)	(S)-malyl-CoA + ADP + phosphate  ->  (S)-malate + coenzyme A + ATP

(RXN-11667)	crotonyl-CoA + H2O  ->  (S)-3-hydroxybutanoyl-CoA

(RXN-13616)	(2E)-dec-2-enoyl-CoA + H2O  ->  (S)-3-hydroxydecanoyl-CoA

(ENOYL-COA-HYDRAT-RXN-POLYMER-INST-L-3-HYDROXYACYL-COA-C10-H20//POLYMER-INST-TRANS-D2-ENOYL-COA-C10-H20/WATER.88. *instantiated*)	POLYMER-INST-TRANS-D2-ENOYL-COA-C10-H20 + H2O  ->  POLYMER-INST-L-3-HYDROXYACYL-COA-C10-H20

(ENOYL-COA-HYDRAT-RXN-POLYMER-INST-L-3-HYDROXYACYL-COA-C12-H24//POLYMER-INST-TRANS-D2-ENOYL-COA-C12-H24/WATER.88. *instantiated*)	POLYMER-INST-TRANS-D2-ENOYL-COA-C12-H24 + H2O  ->  POLYMER-INST-L-3-HYDROXYACYL-COA-C12-H24

(ENOYL-COA-HYDRAT-RXN-POLYMER-INST-L-3-HYDROXYACYL-COA-C14-H28//POLYMER-INST-TRANS-D2-ENOYL-COA-C14-H28/WATER.88. *instantiated*)	POLYMER-INST-TRANS-D2-ENOYL-COA-C14-H28 + H2O  ->  POLYMER-INST-L-3-HYDROXYACYL-COA-C14-H28

(RXN-15127 *spontaneous*)	2-iminopropanoate + H2O  ->  pyruvate + ammonium

(RXN0-5222 *spontaneous*)	carbamate + 2 H+  ->  CO2 + ammonium

(ADENOSINE-KINASE-RXN)	adenosine + ATP  ->  AMP + ADP + H+

(ORNCARBAMTRANSFER-RXN)	L-citrulline + phosphate + H+  ->  L-ornithine + carbamoyl phosphate

(RXN-22736)	carboxyphosphate + ammonia  ->  carbamate + phosphate + H+

(5.3.3.14-RXN)	(2E)-dec-2-enoyl-[acp]  ->  (3Z)-dec-3-enoyl-[acp]

(5.3.3.14-RXN)	(3Z)-dec-3-enoyl-[acp]  ->  (2E)-dec-2-enoyl-[acp]

(DEOXYRIBOSE-P-ALD-RXN)	2-deoxy-D-ribose 5-phosphate  ->  D-glyceraldehyde 3-phosphate + acetaldehyde

(RXN-22285 *spontaneous*)	hydrogen selenide  ->  selenide + 2 H+

(RXN-22285 *spontaneous*)	selenide + 2 H+  ->  hydrogen selenide

(GALACTOKIN-RXN)	alpha-D-galactopyranose + ATP  ->  alpha-D-galactose 1-phosphate + ADP + H+

(RXN-10023)	3-oxododecanoyl-[acp] + S-adenosyl-L-methionine  ->  acyl-carrier protein + S-methyl-5'-thioadenosine + PAI-1 + H+

(RXN-10023)	acyl-carrier protein + S-methyl-5'-thioadenosine + PAI-1 + H+  ->  3-oxododecanoyl-[acp] + S-adenosyl-L-methionine

(RXN-15511)	a [protein]-L-histidine + 2,3-diphospho-D-glycerate  ->  a [protein]-Npi-phospho-L-histidine + 3-phospho-D-glycerate

(RXN-15511)	a [protein]-Npi-phospho-L-histidine + 3-phospho-D-glycerate  ->  a [protein]-L-histidine + 2,3-diphospho-D-glycerate

(ACONITATEHYDR-RXN)	D-threo-isocitrate  ->  cis-aconitate + H2O

(RXN-16997)	phosphorylated phosphoglucomutase + alpha-D-glucopyranose 1-phosphate  ->  alpha-D-glucose 1,6-bisphosphate + phosphoglucomutase

(RXN-16997)	alpha-D-glucose 1,6-bisphosphate + phosphoglucomutase  ->  phosphorylated phosphoglucomutase + alpha-D-glucopyranose 1-phosphate

(RXN-8631)	beta-D-fructofuranose 1-phosphate  ->  glycerone phosphate + D-glyceraldehyde

(RXN-8631)	glycerone phosphate + D-glyceraldehyde  ->  beta-D-fructofuranose 1-phosphate

(PYRIMSYN3-RXN)	4-amino-2-methyl-5-(phosphooxymethyl)pyrimidine + ATP  ->  4-amino-2-methyl-5-(diphosphooxymethyl)pyrimidine + ADP

(URACIL-PRIBOSYLTRANS-RXN)	5-phospho-alpha-D-ribose 1-diphosphate + uracil  ->  UMP + diphosphate

(RXN-18032 *spontaneous*)	hydrogencarbonate  ->  carbonate + H+

(RXN-18032 *spontaneous*)	carbonate + H+  ->  hydrogencarbonate

(GMP-SYN-NH3-RXN)	GMP + AMP + diphosphate + H+  ->  XMP + ammonia + ATP

(RXN0-5114)	O-phospho-L-serine + H2O  ->  L-serine + phosphate

(RXN-10022)	3-oxohexanoyl-[acp] + S-adenosyl-L-methionine  ->  acyl-carrier protein + S-methyl-5'-thioadenosine + VAI-1 + H+

(RXN-10022)	acyl-carrier protein + S-methyl-5'-thioadenosine + VAI-1 + H+  ->  3-oxohexanoyl-[acp] + S-adenosyl-L-methionine

(RXN0-1147)	coenzyme A + a [2-oxoglutarate dehydrogenase E2 protein] N6-S-succinyldihydrolipoyl-L-lysine  ->  succinyl-CoA + a [2-oxoglutarate dehydrogenase E2 protein] N6-dihydrolipoyl-L-lysine

(RXN0-6563)	4-hydroxy-L-threonine  ->  glycine + glycolaldehyde

(RXN0-6563)	glycine + glycolaldehyde  ->  4-hydroxy-L-threonine

(RXN-17777)	(S)-3-hydroxy-(9Z)-octadecenoyl-CoA + NAD+  ->  3-oxo-(9Z)-octadecenoyl-CoA + NADH + H+

(RXN-7933)	N-acetyl-L-citrulline + H2O  ->  L-citrulline + acetate

(RXN-7933)	L-citrulline + acetate  ->  N-acetyl-L-citrulline + H2O

(RXN-20679)	(3S)-3-hydroxyoctanoyl-CoA + NAD+  ->  3-oxooctanoyl-CoA + NADH + H+

(PEPSYNTH-RXN)	pyruvate + ATP + H2O  ->  phosphoenolpyruvate + AMP + phosphate + 2 H+

(2.3.1.78-RXN)	alpha-D-glucosaminide-[heparan sulfate] + acetyl-CoA  ->  N-acetyl-alpha-D-glucosaminide-[heparan sulfate] + coenzyme A

(2.3.1.78-RXN)	N-acetyl-alpha-D-glucosaminide-[heparan sulfate] + coenzyme A  ->  alpha-D-glucosaminide-[heparan sulfate] + acetyl-CoA

(HISTAMINOTRANS-RXN)	L-histidinol phosphate + 2-oxoglutarate  ->  3-(imidazol-4-yl)-2-oxopropyl phosphate + L-glutamate

(RXN-14883 *spontaneous*)	aldehydo-D-ribose  ->  D-ribopyranose

(RXN-14883 *spontaneous*)	D-ribopyranose  ->  aldehydo-D-ribose

(RXN-20678)	(2E)-oct-2-enoyl-CoA + H2O  ->  (3S)-3-hydroxyoctanoyl-CoA

(RXN-16424)	phosphorylated phosphoglucosamine mutase + D-glucosamine 6-phosphate  ->  glucosamine 1,6-diphosphate + phosphoglucosamine mutase

(RXN-16424)	glucosamine 1,6-diphosphate + phosphoglucosamine mutase  ->  phosphorylated phosphoglucosamine mutase + D-glucosamine 6-phosphate

(RXN-16424-Phosphorylated-Phosphoglucosamine-Mutase/CPD-13469//CPD0-1096/Phosphoglucosamine-Mutase.88. *instantiated*)	phosphorylated phosphoglucosamine mutase + alpha-D-glucosamine 6-phosphate  ->  glucosamine 1,6-diphosphate + phosphoglucosamine mutase

(RXN-16424-Phosphorylated-Phosphoglucosamine-Mutase/CPD-13469//CPD0-1096/Phosphoglucosamine-Mutase.88. *instantiated*)	glucosamine 1,6-diphosphate + phosphoglucosamine mutase  ->  phosphorylated phosphoglucosamine mutase + alpha-D-glucosamine 6-phosphate

(RXN-12490)	(S)-3-hydroxydecanoyl-CoA + NAD+  ->  3-oxodecanoyl-CoA + NADH + H+

(METHYLENETHFDEHYDROG-NADP-RXN)	a 5,10-methylenetetrahydrofolate + NADP+  ->  a 5,10-methenyltetrahydrofolate + NADPH

(CARNDETRU-RXN)	(R)-carnitinyl-CoA  ->  crotonobetainyl-CoA + H2O

(CARNDETRU-RXN)	crotonobetainyl-CoA + H2O  ->  (R)-carnitinyl-CoA

(ARGININE--TRNA-LIGASE-RXN)	tRNAArg + L-arginine + ATP  ->  L-arginyl-[tRNAArg] + AMP + diphosphate

(RXN0-2142)	(5Z)-3-oxododec-5-enoyl-[acp] + NADPH + H+  ->  (3R,5Z)-3-hydroxydodec-5-enoyl-[acp] + NADP+

(2.3.1.180-RXN)	acetyl-CoA + malonyl-[acp] + H+  ->  acetoacetyl-[acp] + coenzyme A + CO2

(3-HYDROXY-KYNURENINASE-RXN)	3-hydroxy-L-kynurenine + H2O  ->  3-hydroxyanthranilate + L-alanine + H+

(RXN-9531)	decanoyl-[acp] + malonyl-[acp] + H+  ->  3-oxododecanoyl-[acp] + CO2 + acyl-carrier protein

(ASPARTASE-RXN)	L-aspartate  ->  ammonium + fumarate

(PANTEPADENYLYLTRAN-RXN)	ATP + 4'-phosphopantetheine + H+  ->  3'-dephospho-CoA + diphosphate

(PEPTIDYLPROLYL-ISOMERASE-RXN)	a [protein]-L-proline (omega = 180)  ->  a [protein]-L-proline (omega = 0)

(PEPTIDYLPROLYL-ISOMERASE-RXN)	a [protein]-L-proline (omega = 0)  ->  a [protein]-L-proline (omega = 180)

(RXN0-5393)	(2E,5Z)-tetradecenoyl-CoA + H2O  ->  (S)-3-hydroxy-(5Z)-tetradecenoyl-CoA

(RXN-10994)	apo-[VibB aryl-carrier protein] + coenzyme A  ->  holo-[VibB aryl-carrier protein] + adenosine 3',5'-bisphosphate + H+

(RXN-10994)	holo-[VibB aryl-carrier protein] + adenosine 3',5'-bisphosphate + H+  ->  apo-[VibB aryl-carrier protein] + coenzyme A

(RXN-9384)	O-succinyl-L-homoserine + hydrogen sulfide  ->  L-homocysteine + succinate + H+

(RXN-9384)	L-homocysteine + succinate + H+  ->  O-succinyl-L-homoserine + hydrogen sulfide

(RXN-9539)	myristoyl-[acp] + malonyl-[acp] + H+  ->  3-oxohexadecanoyl-[acp] + CO2 + acyl-carrier protein

(PYRNUTRANSHYDROGEN-RXN)	NAD+ + NADPH  ->  NADH + NADP+

(PYRNUTRANSHYDROGEN-RXN)	NADH + NADP+  ->  NAD+ + NADPH

(RXN-8991)	(2R,3S)-3-isopropylmalate  ->  2-isopropylmaleate + H2O

(RXN0-2301)	3-methylbutanoyl-CoA + oxidized electron-transfer flavoprotein + H+  ->  3-methylcrotonyl-CoA + reduced electron-transfer flavoprotein

(N-ACETYLGLUTPREDUCT-RXN)	N-acetylglutamyl-phosphate + NADPH + H+  ->  N-acetyl-L-glutamate 5-semialdehyde + NADP+ + phosphate

(GCVP-RXN)	glycine + a [glycine-cleavage complex H protein] N6-[(R)-lipoyl]-L-lysine + H+  ->  a [glycine-cleavage complex H protein] N6-aminomethyldihydrolipoyl-L-lysine + CO2

(GLUCONOKIN-RXN)	ATP + D-gluconate  ->  D-gluconate 6-phosphate + ADP + H+

(RXN-9528)	3-oxodecanoyl-[acp] + NADPH + H+  ->  (3R)-3-hydroxydecanoyl-[acp] + NADP+

(GLYCOLATEDEHYDRO-RXN)	glycolate + oxidized electron carrier  ->  glyoxylate + reduced two electron carrier

(FUMHYDR-RXN)	(S)-malate  ->  fumarate + H2O

(RXN-8899 *spontaneous*)	2-aminoprop-2-enoate + H2O  ->  pyruvate + ammonium

(GLUTCYSLIG-RXN)	L-glutamate + L-cysteine + ATP  ->  gamma-L-glutamyl-L-cysteine + ADP + phosphate + H+

(THREONINE-ALDOLASE-RXN)	L-threonine  ->  glycine + acetaldehyde

(THREONINE-ALDOLASE-RXN)	glycine + acetaldehyde  ->  L-threonine

(RXN-14957)	a [pyruvate dehydrogenase E2 protein] N6-octanoyl-L-lysine + 2 sulfurated [sulfur carrier] + 2 reduced [2Fe-2S] ferredoxin + 2 S-adenosyl-L-methionine  ->  a [pyruvate dehydrogenase E2 protein] N6-lipoyl-L-lysine + 2 L-methionine + 2 5'-deoxyadenosine + 2 unsulfurated [sulfur carrier] + 2 oxidized [2Fe-2S] ferredoxin

(RXN-14957)	a [pyruvate dehydrogenase E2 protein] N6-lipoyl-L-lysine + 2 L-methionine + 2 5'-deoxyadenosine + 2 unsulfurated [sulfur carrier] + 2 oxidized [2Fe-2S] ferredoxin  ->  a [pyruvate dehydrogenase E2 protein] N6-octanoyl-L-lysine + 2 sulfurated [sulfur carrier] + 2 reduced [2Fe-2S] ferredoxin + 2 S-adenosyl-L-methionine

(RXN-22603 *spontaneous*)	(R)-lactic acid  ->  (S)-lactic acid

(RXN-22603 *spontaneous*)	(S)-lactic acid  ->  (R)-lactic acid

(RXNI-2)	succinyl-CoA + acetoacetate  ->  succinate + acetoacetyl-CoA

(RXNI-2)	succinate + acetoacetyl-CoA  ->  succinyl-CoA + acetoacetate

(ACETOACETATE--COA-LIGASE-RXN)	coenzyme A + acetoacetate + ATP  ->  acetoacetyl-CoA + AMP + diphosphate

(RXN-12198)	CDP + H2O  ->  CMP + phosphate + H+

(4.2.1.58-RXN)	(3R)-3-hydroxybutanoyl-[acp]  ->  (2E)-but-2-enoyl-[acp] + H2O

(4.2.1.59-RXN)	(3R)-3-hydroxyoctanoyl-[acp]  ->  (2E)-oct-2-enoyl-[acp] + H2O

(RXN-9655)	(3R)-3-hydroxydecanoyl-[acp]  ->  (2E)-dec-2-enoyl-[acp] + H2O

(RXN-9533)	(3R)-3-hydroxydodecanoyl-[acp]  ->  (2E)-dodec-2-enoyl-[acp] + H2O

(RXN-9537)	(3R)-3-hydroxytetradecanoyl-[acp]  ->  (2E)-tetradec-2-enoyl-[acp] + H2O

(3-HYDROXYDECANOYL-ACP-DEHYDR-RXN-POLYMER-INST-OH-ACYL-ACP-C12-H24//2-Hexadecenoyl-ACPs/WATER.60. *instantiated*)	POLYMER-INST-OH-ACYL-ACP-C12-H24  ->  (2E)-hexadec-2-enoyl-[acp] + H2O

(RXN-9634)	(3R)-3-hydroxyoctadecanoyl-[acp]  ->  (2E)-octadec-2-enoyl-[acp] + H2O

(ADENODEAMIN-RXN)	adenosine + H2O + H+  ->  inosine + ammonium

(ASPARTATE-SEMIALDEHYDE-DEHYDROGENASE-RXN)	L-aspartate 4-semialdehyde + NADP+ + phosphate  ->  L-aspartyl-4-phosphate + NADPH + H+

(RXN-9520)	(3R)-3-hydroxyhexanoyl-[acp]  ->  (2E)-hex-2-enoyl-[acp] + H2O

(RXN0-743)	N5-carboxyaminoimidazole ribonucleotide  ->  5-amino-1-(5-phospho-D-ribosyl)imidazole-4-carboxylate

(L-ASPARTATE-OXID-RXN)	L-aspartate + dioxygen  ->  2-iminosuccinate + hydrogen peroxide + H+

(RXN-9524)	3-oxooctanoyl-[acp] + NADPH + H+  ->  (3R)-3-hydroxyoctanoyl-[acp] + NADP+

(PPENTOMUT-RXN)	alpha-D-ribose-1-phosphate  ->  D-ribose 5-phosphate

(RXN-12753 *spontaneous*)	NADH + H2O  ->  (S)-NADHX

(RXN0-3962)	acetaldehyde + NADP+ + H2O  ->  acetate + NADPH + 2 H+

(RXN-10462)	acyl-[acyl-carrier protein] + sn-glycerol 3-phosphate  ->  a 2-lysophosphatidate + acyl-carrier protein

(METHYLMALONYL-COA-EPIM-RXN)	(S)-methylmalonyl-CoA  ->  (R)-methylmalonyl-CoA

(RXN-10026)	(3R)-3-hydroxybutanoyl-[acp] + S-adenosyl-L-methionine  ->  acyl-carrier protein + S-methyl-5'-thioadenosine + HAI-1 + H+

(RXN-10026)	acyl-carrier protein + S-methyl-5'-thioadenosine + HAI-1 + H+  ->  (3R)-3-hydroxybutanoyl-[acp] + S-adenosyl-L-methionine

(RXN-9544)	(3R)-3-hydroxyoctadecanoyl-CoA + NADP+  ->  3-oxooctadecanoyl-CoA + NADPH + H+

(RXN-9544)	3-oxooctadecanoyl-CoA + NADPH + H+  ->  (3R)-3-hydroxyoctadecanoyl-CoA + NADP+

(RXN-11350)	2 {beta-D-GlcNAc-(1->4)-3-O-[L-Ala-gamma-D-iGln-(6-N-beta-D-Asn)-L-Lys-D-Ala-D-Ala]-Mur2Ac}-PP-Und  ->  peptidoglycan dimer (E. faeciums) + di-trans,octa-cis-undecaprenyl diphosphate + H+

(RXN-11350)	peptidoglycan dimer (E. faeciums) + di-trans,octa-cis-undecaprenyl diphosphate + H+  ->  2 {beta-D-GlcNAc-(1->4)-3-O-[L-Ala-gamma-D-iGln-(6-N-beta-D-Asn)-L-Lys-D-Ala-D-Ala]-Mur2Ac}-PP-Und

(RXN-15348 *spontaneous*)	S-sulfanylglutathione + glutathione  ->  hydrogen sulfide + glutathione disulfide

(RXN-15348 *spontaneous*)	hydrogen sulfide + glutathione disulfide  ->  S-sulfanylglutathione + glutathione

(RXN-9548)	stearoyl-[acp] + H2O  ->  stearate + acyl-carrier protein + H+

(INORGPYROPHOSPHAT-RXN)	diphosphate + H2O  ->  2 phosphate + H+

(RXN-9657)	(2E)-but-2-enoyl-[acp] + NADH + H+  ->  butanoyl-[acp] + NAD+

(RXN-9658)	(2E)-hex-2-enoyl-[acp] + NADH + H+  ->  hexanoyl-[acp] + NAD+

(RXN-9659)	(2E)-oct-2-enoyl-[acp] + NADH + H+  ->  octanoyl-[acp] + NAD+

(RXN-9660)	(2E)-dec-2-enoyl-[acp] + NADH + H+  ->  decanoyl-[acp] + NAD+

(ENOYL-ACP-REDUCT-NADH-RXN-Stearoyl-ACPs/NAD//Octadec-2-enoyl-ACPs/NADH/PROTON.52. *instantiated*)	(2E)-octadec-2-enoyl-[acp] + NADH + H+  ->  stearoyl-[acp] + NAD+

(RXN-15510)	a [protein]-Npi-phospho-L-histidine + 2-phospho-D-glycerate  ->  a [protein]-L-histidine + 2,3-diphospho-D-glycerate

(RXN-15510)	a [protein]-L-histidine + 2,3-diphospho-D-glycerate  ->  a [protein]-Npi-phospho-L-histidine + 2-phospho-D-glycerate

(RXN-2962-S-HYDROXYMETHYLGLUTATHIONE/NAD//CPD-548/NADH/PROTON.52. *instantiated*)	S-(hydroxymethyl)glutathione + NAD+  ->  S-formylglutathione + NADH + H+

(RXN-2962-S-HYDROXYMETHYLGLUTATHIONE/NADP//CPD-548/NADPH/PROTON.54. *instantiated*)	S-(hydroxymethyl)glutathione + NADP+  ->  S-formylglutathione + NADPH + H+

(RXN-18204)	3-[(5'-methylsulfanyl)pentyl]malate + NAD+  ->  3-carboxy-8-(methylsulfanyl)-2-oxooctanoate + NADH + H+

(RXN-18204)	3-carboxy-8-(methylsulfanyl)-2-oxooctanoate + NADH + H+  ->  3-[(5'-methylsulfanyl)pentyl]malate + NAD+

(GLUTATHIONE-PEROXIDASE-RXN)	hydrogen peroxide + 2 glutathione  ->  glutathione disulfide + 2 H2O

(RXN-17823)	DNA methyl phosphotriester + a [protein]-L-cysteine  ->  [DNA] + a [protein]-S-methyl-L-cysteine + H+

(RXN-17823)	[DNA] + a [protein]-S-methyl-L-cysteine + H+  ->  DNA methyl phosphotriester + a [protein]-L-cysteine

(RXN-12565)	butanoyl-CoA + acetyl-CoA  ->  3-oxohexanoyl-CoA + coenzyme A

(RXN-12565)	3-oxohexanoyl-CoA + coenzyme A  ->  butanoyl-CoA + acetyl-CoA

(RXN-14277)	3-oxooctanoyl-CoA + coenzyme A  ->  hexanoyl-CoA + acetyl-CoA

(RXN-13617)	3-oxodecanoyl-CoA + coenzyme A  ->  octanoyl-CoA + acetyl-CoA

(KETOACYLCOATHIOL-RXN-LAUROYLCOA-CPD/ACETYL-COA//POLYMER-INST-3-KETOACYL-COA-C10-H20/CO-A.68. *instantiated*)	POLYMER-INST-3-KETOACYL-COA-C10-H20 + coenzyme A  ->  lauroyl-CoA + acetyl-CoA

(KETOACYLCOATHIOL-RXN-POLYMER-INST-Saturated-Fatty-Acyl-CoA-C10-H20/ACETYL-COA//POLYMER-INST-3-KETOACYL-COA-C12-H24/CO-A.99. *instantiated*)	POLYMER-INST-3-KETOACYL-COA-C12-H24 + coenzyme A  ->  POLYMER-INST-Saturated-Fatty-Acyl-CoA-C10-H20 + acetyl-CoA

(KETOACYLCOATHIOL-RXN-PALMITYL-COA/ACETYL-COA//CPD-10260/CO-A.40. *instantiated*)	3-oxooctadecanoyl-CoA + coenzyme A  ->  palmitoyl-CoA + acetyl-CoA

(KETOACYLCOATHIOL-RXN-STEAROYL-COA/ACETYL-COA//POLYMER-INST-3-KETOACYL-COA-C16-H32/CO-A.66. *instantiated*)	POLYMER-INST-3-KETOACYL-COA-C16-H32 + coenzyme A  ->  stearoyl-CoA + acetyl-CoA

(GLUCONATE-5-DEHYDROGENASE-RXN-GLUCONATE/NAD//5-DEHYDROGLUCONATE/NADH/PROTON.46. *instantiated*)	D-gluconate + NAD+  ->  5-dehydro-D-gluconate + NADH + H+

(GLUCONATE-5-DEHYDROGENASE-RXN-GLUCONATE/NAD//5-DEHYDROGLUCONATE/NADH/PROTON.46. *instantiated*)	5-dehydro-D-gluconate + NADH + H+  ->  D-gluconate + NAD+

(GLUCONATE-5-DEHYDROGENASE-RXN-GLUCONATE/NADP//5-DEHYDROGLUCONATE/NADPH/PROTON.48. *instantiated*)	D-gluconate + NADP+  ->  5-dehydro-D-gluconate + NADPH + H+

(GLUCONATE-5-DEHYDROGENASE-RXN-GLUCONATE/NADP//5-DEHYDROGLUCONATE/NADPH/PROTON.48. *instantiated*)	5-dehydro-D-gluconate + NADPH + H+  ->  D-gluconate + NADP+

(ACYLCOASYN-RXN-BUTYRIC_ACID/CO-A/ATP//BUTYRYL-COA/PPI/AMP.43. *instantiated*)	butanoate + coenzyme A + ATP  ->  butanoyl-CoA + diphosphate + AMP

(R223-RXN)	octanoate + ATP + coenzyme A  ->  octanoyl-CoA + AMP + diphosphate

(RXN-16393)	laurate + coenzyme A + ATP  ->  lauroyl-CoA + AMP + diphosphate

(RXN-9623)	palmitate + coenzyme A + ATP  ->  palmitoyl-CoA + AMP + diphosphate

(RXN-16380)	stearate + coenzyme A + ATP  ->  stearoyl-CoA + AMP + diphosphate

(RXN-8959)	(2S)-methylsuccinyl-CoA + oxidized electron-transfer flavoprotein + H+  ->  2-methylfumaryl-CoA + reduced electron-transfer flavoprotein

(RXN0-748)	dGDP + oxidized NrdH glutaredoxin-like protein + H2O  ->  GDP + reduced NrdH glutaredoxin-like protein

(RXN0-748)	GDP + reduced NrdH glutaredoxin-like protein  ->  dGDP + oxidized NrdH glutaredoxin-like protein + H2O

(METHYLMALONYL-COA-MUT-RXN)	(R)-methylmalonyl-CoA  ->  succinyl-CoA

(RXN0-742)	5-amino-1-(5-phospho-beta-D-ribosyl)imidazole + ATP + hydrogencarbonate  ->  N5-carboxyaminoimidazole ribonucleotide + ADP + phosphate + 2 H+

(RXN-7745)	(2R,3S)-3-methylmalate + NAD+  ->  2-oxobutanoate + CO2 + NADH

(RXN-8173 *spontaneous*)	(S)-2-amino-6-oxohexanoate  ->  1-piperideine 6-carboxylate + H2O + H+

(RXN-8173 *spontaneous*)	1-piperideine 6-carboxylate + H2O + H+  ->  (S)-2-amino-6-oxohexanoate

(RXN-14959)	[2-oxoglutarate-dehydrogenase E2 protein] N6-octanoyl-L-lysine + 2 sulfurated [sulfur carrier] + 2 reduced [2Fe-2S] ferredoxin + 2 S-adenosyl-L-methionine  ->  a [2-oxoglutarate dehydrogenase E2 protein] N6-lipoyl-L-lysine + 2 L-methionine + 2 5'-deoxyadenosine + 2 unsulfurated [sulfur carrier] + 2 oxidized [2Fe-2S] ferredoxin

(RXN-14959)	a [2-oxoglutarate dehydrogenase E2 protein] N6-lipoyl-L-lysine + 2 L-methionine + 2 5'-deoxyadenosine + 2 unsulfurated [sulfur carrier] + 2 oxidized [2Fe-2S] ferredoxin  ->  [2-oxoglutarate-dehydrogenase E2 protein] N6-octanoyl-L-lysine + 2 sulfurated [sulfur carrier] + 2 reduced [2Fe-2S] ferredoxin + 2 S-adenosyl-L-methionine

(RXN-20456 *spontaneous*)	NADPH + H2O  ->  (R)-NADPHX

(2.7.7.60-RXN)	2-C-methyl-D-erythritol 4-phosphate + CTP + H+  ->  4-(cytidine 5'-diphospho)-2-C-methyl-D-erythritol + diphosphate

(CYTIDEAM2-RXN)	cytidine + H+ + H2O  ->  uridine + ammonium

(CYTIDEAM2-RXN)	uridine + ammonium  ->  cytidine + H+ + H2O

(CERAMIDASE-RXN)	an N-acylsphingosine + H2O  ->  sphingosine + a fatty acid

(CERAMIDASE-RXN)	sphingosine + a fatty acid  ->  an N-acylsphingosine + H2O

(RXN-14304)	5'-deoxyadenosine + phosphate  ->  5-deoxy-alpha-D-ribose 1-phosphate + adenine

(HOMOSERINE-O-ACETYLTRANSFERASE-RXN)	L-homoserine + acetyl-CoA  ->  O-acetyl-L-homoserine + coenzyme A

(RXN-5901)	acetoacetyl-CoA + NADPH + H+  ->  (3R)-3-hydroxybutanoyl-CoA + NADP+

(DALADALALIG-RXN)	2 D-alanine + ATP  ->  D-alanyl-D-alanine + ADP + phosphate + H+

(R-2-METHYLMALATE-DEHYDRATASE-RXN)	(R)-citramalate  ->  citraconate + H2O

(R-2-METHYLMALATE-DEHYDRATASE-RXN)	citraconate + H2O  ->  (R)-citramalate

(RXN-9536)	3-oxotetradecanoyl-[acp] + NADPH + H+  ->  (3R)-3-hydroxytetradecanoyl-[acp] + NADP+

(RXN-12752)	(R)-NADHX  ->  (S)-NADHX

(HYDROXYMETHYLGLUTARYL-COA-LYASE-RXN)	(S)-3-hydroxy-3-methylglutaryl-CoA  ->  acetoacetate + acetyl-CoA

(DUTP-PYROP-RXN)	dUTP + H2O  ->  dUMP + diphosphate + H+

(MALEYLACETOACETATE-ISOMERASE-RXN)	4-maleyl-acetoacetate  ->  4-fumaryl-acetoacetate

(KETOISOCAPROATE-RXN)	4-methyl-2-oxopentanoate + a [BCAA dehydrogenase E2 protein] N6-lipoyl-L-lysine + H+  ->  a [apo BCAA dehydrogenase E2 protein] N6-S-[3-methylbutanoyl]dihydrolipoyl-L-lysine + CO2

(ATPPHOSPHORIBOSYLTRANS-RXN)	1-(5-phospho-beta-D-ribosyl)-ATP + diphosphate  ->  ATP + 5-phospho-alpha-D-ribose 1-diphosphate

(PYRROLINECARBREDUCT-RXN-PRO/NADP//L-DELTA1-PYRROLINE_5-CARBOXYLATE/NADPH/PROTON.56. *instantiated*)	(S)-1-pyrroline-5-carboxylate + NADPH + 2 H+  ->  L-proline + NADP+

(RXN-9662)	(2E)-tetradec-2-enoyl-[acp] + NADH + H+  ->  myristoyl-[acp] + NAD+

(RXN-8751)	hypoxanthine + NAD+ + H2O  ->  6,8-dihydroxypurine + NADH + H+

(RXN-8751)	6,8-dihydroxypurine + NADH + H+  ->  hypoxanthine + NAD+ + H2O

(RXN-8629)	a [glycine-cleavage complex H protein] N6-dihydrolipoyl-L-lysine + NAD+  ->  a [glycine-cleavage complex H protein] N6-[(R)-lipoyl]-L-lysine + NADH + H+

(RXN-14809 *spontaneous*)	L-arabinofuranose  ->  L-arabinopyranose

(RXN-14809 *spontaneous*)	L-arabinopyranose  ->  L-arabinofuranose

(RXN-14809-CPD-12046//L-arabinopyranose.29. *instantiated* *spontaneous*)	beta-L-arabinofuranose  ->  L-arabinopyranose

(RXN-14809-CPD-12046//L-arabinopyranose.29. *instantiated* *spontaneous*)	L-arabinopyranose  ->  beta-L-arabinofuranose

(RXN-14809-CPD-12045//L-arabinopyranose.29. *instantiated* *spontaneous*)	alpha-L-arabinofuranose  ->  L-arabinopyranose

(RXN-14809-CPD-12045//L-arabinopyranose.29. *instantiated* *spontaneous*)	L-arabinopyranose  ->  alpha-L-arabinofuranose

(PGLUCONDEHYDRAT-RXN)	D-gluconate 6-phosphate  ->  2-dehydro-3-deoxy-D-gluconate 6-phosphate + H2O

(METHIONYL-TRNA-FORMYLTRANSFERASE-RXN)	an N10-formyltetrahydrofolate + L-methionyl-[initiator tRNAMet]  ->  a tetrahydrofolate + N-formyl-L-methionyl-[initiator tRNAmet] + H+

(METHIONYL-TRNA-FORMYLTRANSFERASE-RXN)	a tetrahydrofolate + N-formyl-L-methionyl-[initiator tRNAmet] + H+  ->  an N10-formyltetrahydrofolate + L-methionyl-[initiator tRNAMet]

(DHLBXANAU-RXN)	2-chloroacetate + H2O  ->  glycolate + chloride + H+

(DHLBXANAU-RXN)	glycolate + chloride + H+  ->  2-chloroacetate + H2O

(ADENINE-DEAMINASE-RXN)	adenine + H+ + H2O  ->  ammonium + hypoxanthine

(GSAAMINOTRANS-RXN)	(S)-4-amino-5-oxopentanoate  ->  5-aminolevulinate

(GLYCEROL-KIN-RXN)	glycerol + ATP  ->  sn-glycerol 3-phosphate + ADP + H+

(RXN-8752)	6,8-dihydroxypurine + NAD+ + H2O  ->  urate + NADH + H+

(GLUTAMIN-RXN)	L-glutamine + H2O  ->  L-glutamate + ammonia + H+

(RXN-12508)	2-(alpha-hydroxyethyl)thiamine diphosphate + a [pyruvate dehydrogenase E2 protein] N6-lipoyl-L-lysine  ->  a [pyruvate dehydrogenase E2 protein] N6-S-acetyldihydrolipoyl-L-lysine + thiamine diphosphate

(PRODISULFREDUCT-A-RXN-Oxidized-NrdH-Proteins/GLUTATHIONE//OXIDIZED-GLUTATHIONE/Reduced-NrdH-Proteins.79. *instantiated* *spontaneous*)	oxidized NrdH glutaredoxin-like protein + 2 glutathione  ->  glutathione disulfide + reduced NrdH glutaredoxin-like protein

(AMINO-ACID-RACEMASE-RXN-L-ASPARTATE//CPD-302.21. *instantiated*)	L-aspartate  ->  D-aspartate

(AMINO-ACID-RACEMASE-RXN-L-ASPARTATE//CPD-302.21. *instantiated*)	D-aspartate  ->  L-aspartate

(AMINO-ACID-RACEMASE-RXN-TYR//D-TYROSINE.16. *instantiated*)	L-tyrosine  ->  D-tyrosine

(AMINO-ACID-RACEMASE-RXN-TYR//D-TYROSINE.16. *instantiated*)	D-tyrosine  ->  L-tyrosine

(AMINO-ACID-RACEMASE-RXN-CYS//D-CYSTEINE.16. *instantiated*)	L-cysteine  ->  D-cysteine

(AMINO-ACID-RACEMASE-RXN-CYS//D-CYSTEINE.16. *instantiated*)	D-cysteine  ->  L-cysteine

(AMINO-ACID-RACEMASE-RXN-TRP//D-TRYPTOPHAN.18. *instantiated*)	L-tryptophan  ->  D-tryptophan

(AMINO-ACID-RACEMASE-RXN-TRP//D-TRYPTOPHAN.18. *instantiated*)	D-tryptophan  ->  L-tryptophan

(RXN-20896)	L-asparagine  ->  D-asparagine

(RXN-20896)	D-asparagine  ->  L-asparagine

(GLUTRACE-RXN)	L-glutamate  ->  D-glutamate

(RXN-18202)	3-[(6'-methylsulfanyl)hexyl]malate + NAD+  ->  3-carboxy-9-(methylsulfanyl)-2-oxononanoate + NADH + H+

(RXN-18202)	3-carboxy-9-(methylsulfanyl)-2-oxononanoate + NADH + H+  ->  3-[(6'-methylsulfanyl)hexyl]malate + NAD+

(L-GLN-FRUCT-6-P-AMINOTRANS-RXN-FRUCTOSE-6P/GLN//CPD-13469/GLT.31. *instantiated*)	beta-D-fructofuranose 6-phosphate + L-glutamine  ->  alpha-D-glucosamine 6-phosphate + L-glutamate

(L-GLN-FRUCT-6-P-AMINOTRANS-RXN)	beta-D-fructofuranose 6-phosphate + L-glutamine  ->  D-glucosamine 6-phosphate + L-glutamate

(L-GLN-FRUCT-6-P-AMINOTRANS-RXN)	D-glucosamine 6-phosphate + L-glutamate  ->  beta-D-fructofuranose 6-phosphate + L-glutamine

(TRYPTOPHAN-AMINOTRANSFERASE-RXN)	L-tryptophan + 2-oxoglutarate  ->  (indol-3-yl)pyruvate + L-glutamate

(TRYPTOPHAN-AMINOTRANSFERASE-RXN)	(indol-3-yl)pyruvate + L-glutamate  ->  L-tryptophan + 2-oxoglutarate

(PGLUCISOM-RXN)	beta-D-fructofuranose 6-phosphate  ->  alpha-D-glucose 6-phosphate

(HOMOGENTISATE-12-DIOXYGENASE-RXN)	homogentisate + dioxygen  ->  4-maleyl-acetoacetate + H+

(1.2.2.3-RXN)	2 reduced cytochrome c-553 + CO2 + H+  ->  2 oxidized cytochrome c-553 + formate

(RXN-12583)	pyruvate + thiamine diphosphate + H+  ->  2-(alpha-hydroxyethyl)thiamine diphosphate + CO2

(D-PPENTOMUT-RXN)	2-deoxy-alpha-D-ribose 1-phosphate  ->  2-deoxy-D-ribose 5-phosphate

(ISPH2-RXN)	(E)-4-hydroxy-3-methylbut-2-en-1-yl diphosphate + 2 reduced ferredoxin [iron-sulfur] cluster + 2 H+  ->  3-methylbut-3-en-1-yl diphosphate + 2 oxidized ferredoxin [iron-sulfur] cluster + H2O

(RXN-20154)	acetyl-CoA + a [protein]-L-cysteine  ->  a [protein] S-acetyl-L-cysteine + coenzyme A

(RXN-20154)	a [protein] S-acetyl-L-cysteine + coenzyme A  ->  acetyl-CoA + a [protein]-L-cysteine

(THIAZOLSYN3-RXN)	ATP + 5-(2-hydroxyethyl)-4-methylthiazole  ->  ADP + 4-methyl-5-(2-phosphooxyethyl)thiazole + H+

(RXN-7719)	an [apo BCAA dehydrogenase E2 protein] N6-dihydrolipoyl-L-lysine + NAD+  ->  a [BCAA dehydrogenase E2 protein] N6-lipoyl-L-lysine + NADH + H+

(RXN-7719)	a [BCAA dehydrogenase E2 protein] N6-lipoyl-L-lysine + NADH + H+  ->  an [apo BCAA dehydrogenase E2 protein] N6-dihydrolipoyl-L-lysine + NAD+

(GLUTAMINESYN-RXN)	ammonium + L-glutamate + ATP  ->  L-glutamine + ADP + phosphate + H+

(UDPKIN-RXN)	ATP + UDP  ->  ADP + UTP

(RXN-2961 *spontaneous*)	formaldehyde + glutathione  ->  S-(hydroxymethyl)glutathione

(RXN-14116)	L-glutamate-5-semialdehyde + NAD+ + H2O  ->  L-glutamate + NADH + 2 H+

(RXN-16950)	3-hydroxy-5-methyl-benzoate + NADH + dioxygen + H+  ->  3-methylgentisate + NAD+ + H2O

(RXN-16950)	3-methylgentisate + NAD+ + H2O  ->  3-hydroxy-5-methyl-benzoate + NADH + dioxygen + H+

(AKBLIG-RXN)	glycine + acetyl-CoA  ->  L-2-amino-3-oxobutanoate + coenzyme A

(AKBLIG-RXN)	L-2-amino-3-oxobutanoate + coenzyme A  ->  glycine + acetyl-CoA

(RXN-6268 *spontaneous*)	betaine aldehyde hydrate  ->  betaine aldehyde + H2O

(RXN-6268 *spontaneous*)	betaine aldehyde + H2O  ->  betaine aldehyde hydrate

(GLUCOSAMINE-6-P-DEAMIN-RXN)	alpha-D-glucosamine 6-phosphate + H2O  ->  beta-D-fructofuranose 6-phosphate + ammonium

(RXN-12002)	ATP + UMP  ->  ADP + UDP

(RXN-9515)	(2E)-but-2-enoyl-[acp] + NADPH + H+  ->  butanoyl-[acp] + NADP+

(RXN-9521)	(2E)-hex-2-enoyl-[acp] + NADPH + H+  ->  hexanoyl-[acp] + NADP+

(RXN-9526)	(2E)-oct-2-enoyl-[acp] + NADPH + H+  ->  octanoyl-[acp] + NADP+

(RXN-9530)	(2E)-dec-2-enoyl-[acp] + NADPH + H+  ->  decanoyl-[acp] + NADP+

(RXN-9534)	(2E)-dodec-2-enoyl-[acp] + NADPH + H+  ->  dodecanoyl-[acp] + NADP+

(RXN-9538)	(2E)-tetradec-2-enoyl-[acp] + NADPH + H+  ->  myristoyl-[acp] + NADP+

(RXN-9542)	(2E)-hexadec-2-enoyl-[acp] + NADPH + H+  ->  palmitoyl-[acp] + NADP+

(RXN-14882 *spontaneous*)	aldehydo-D-ribose  ->  D-ribofuranose

(RXN-14882 *spontaneous*)	D-ribofuranose  ->  aldehydo-D-ribose

(RXN-8975)	UDP-N-acetyl-alpha-D-muramoyl-L-alanyl-gamma-D-glutamyl-L-lysyl-D-alanyl-D-alanine + di-trans,octa-cis-undecaprenyl phosphate  ->  Und-PP-Mur2Ac-L-Ala-gamma-D-Glu-L-Lys-D-Ala-D-Ala + UMP

(RXN-8975)	Und-PP-Mur2Ac-L-Ala-gamma-D-Glu-L-Lys-D-Ala-D-Ala + UMP  ->  UDP-N-acetyl-alpha-D-muramoyl-L-alanyl-gamma-D-glutamyl-L-lysyl-D-alanyl-D-alanine + di-trans,octa-cis-undecaprenyl phosphate

(PYRUFLAVREDUCT-RXN)	pyruvate + coenzyme A + 2 oxidized ferredoxin [iron-sulfur] cluster  ->  acetyl-CoA + CO2 + 2 reduced ferredoxin [iron-sulfur] cluster + H+

(PYRUFLAVREDUCT-RXN)	acetyl-CoA + CO2 + 2 reduced ferredoxin [iron-sulfur] cluster + H+  ->  pyruvate + coenzyme A + 2 oxidized ferredoxin [iron-sulfur] cluster

(ETHAMLY-RXN)	ethanolamine  ->  ammonium + acetaldehyde

(RXN-34)	L-canavanine + H2O  ->  L-canaline + urea + H+

(RXN-34)	L-canaline + urea + H+  ->  L-canavanine + H2O

(HYDROXYHEPTA-DIENEDIOATE-HYDROXY-RXN)	(4Z)-2-oxohept-4-enedioate + H2O  ->  (4S)-4-hydroxy-2-oxoheptanedioate

(HYDROXYHEPTA-DIENEDIOATE-HYDROXY-RXN)	(4S)-4-hydroxy-2-oxoheptanedioate  ->  (4Z)-2-oxohept-4-enedioate + H2O

(DEOXYADENPHOSPHOR-RXN)	2'-deoxyadenosine + phosphate  ->  adenine + 2-deoxy-alpha-D-ribose 1-phosphate

(DEOXYADENPHOSPHOR-RXN)	adenine + 2-deoxy-alpha-D-ribose 1-phosphate  ->  2'-deoxyadenosine + phosphate

(RXN-7682)	hypoxanthine + NAD+ + H2O  ->  xanthine + NADH + H+

(RXN0-6727)	(S)-NADHX + ADP  ->  AMP + NADH + phosphate + H+

(2.3.1.157-RXN)	acetyl-CoA + alpha-D-glucosamine 1-phosphate  ->  N-acetyl-alpha-D-glucosamine 1-phosphate + coenzyme A + H+

(PHOSACETYLTRANS-RXN)	acetyl-CoA + phosphate  ->  acetyl phosphate + coenzyme A

(RXN-11291)	2 {beta-D-GlcNAc-(1->4)-3-O-[L-Ala-gamma-D-iGln-(6-N-Gly5)-L-Lys-D-Ala-D-Ala]-Mur2Ac}-PP-Und  ->  di-trans,octa-cis-undecaprenyl diphosphate + peptidoglycan dimer (S. aureus) + H+

(RXN-11291)	di-trans,octa-cis-undecaprenyl diphosphate + peptidoglycan dimer (S. aureus) + H+  ->  2 {beta-D-GlcNAc-(1->4)-3-O-[L-Ala-gamma-D-iGln-(6-N-Gly5)-L-Lys-D-Ala-D-Ala]-Mur2Ac}-PP-Und

(RXN-6383)	3-hydroxypropanoyl-CoA  ->  acryloyl-CoA + H2O

(RXN-6383)	acryloyl-CoA + H2O  ->  3-hydroxypropanoyl-CoA

(RXN0-5305 *spontaneous*)	beta-D-ribofuranose  ->  aldehydo-D-ribose

(RXN0-5305 *spontaneous*)	aldehydo-D-ribose  ->  beta-D-ribofuranose

(MALSYN-RXN)	acetyl-CoA + glyoxylate + H2O  ->  (S)-malate + coenzyme A + H+

(RXN0-5224)	hydrogencarbonate + H+  ->  CO2 + H2O

(UROCANATE-HYDRATASE-RXN)	urocanate + H2O  ->  4-imidazolone-5-propanoate

(ORNITHINE-GLU-AMINOTRANSFERASE-RXN)	L-ornithine + 2-oxoglutarate  ->  L-glutamate + L-glutamate-5-semialdehyde

(ALANINE-DEHYDROGENASE-RXN)	L-alanine + NAD+ + H2O  ->  pyruvate + ammonium + NADH + H+

(RXN-9523)	hexanoyl-[acp] + malonyl-[acp] + H+  ->  3-oxooctanoyl-[acp] + CO2 + acyl-carrier protein

(RXN-22455)	(2E)-octadec-2-enoyl-[acp] + NADPH + H+  ->  stearoyl-[acp] + NADP+

(RXN-6002)	L-malic semialdehyde + NADP+ + H2O  ->  (S)-malate + NADPH + 2 H+

(RXN-6002)	(S)-malate + NADPH + 2 H+  ->  L-malic semialdehyde + NADP+ + H2O

(RXNN-386)	2-hydroxychromene-2-carboxylate  ->  trans-O-hydroxybenzylidenepyruvate

(RXNN-386)	trans-O-hydroxybenzylidenepyruvate  ->  2-hydroxychromene-2-carboxylate

(KYNURENINE-3-MONOOXYGENASE-RXN)	L-kynurenine + NADPH + dioxygen + H+  ->  3-hydroxy-L-kynurenine + NADP+ + H2O

(PHOSPHASERDECARB-RXN)	a phosphatidylserine + H+  ->  a phosphatidylethanolamine + CO2

(BETA-UREIDOPROPIONASE-RXN)	3-ureidopropanoate + H2O + 2 H+  ->  beta-alanine + ammonium + CO2

(RXN-14394)	3-oxo-(5Z)-tetradecenoyl-CoA + coenzyme A  ->  (3Z)-dodec-3-enoyl-CoA + acetyl-CoA

(RXN-17780)	(2E,7Z)-hexadecenoyl-CoA + H2O  ->  (S)-3-hydroxy-(7Z)-hexadecenoyl-CoA

(2.1.1.5-RXN)	L-homocysteine + glycine betaine  ->  L-methionine + N,N-dimethylglycine

(GUANINE-DEAMINASE-RXN)	guanine + H+ + H2O  ->  xanthine + ammonium

(RXN-12019)	adenosine 5'-phosphosulfate + reduced thioredoxin  ->  sulfite + AMP + oxidized thioredoxin + 2 H+

(THIAZOLSYN2-RXN)	1-deoxy-D-xylulose 5-phosphate + 2-iminoacetate + thiocarboxylated-[ThiS sulfur-carrier protein]  ->  2-[(2R,5Z)-2-carboxy-4-methylthiazol-5(2H)-ylidene]ethyl phosphate + [ThiS sulfur-carrier protein] + 2 H2O

(RXN-10657)	(2E,7Z)-tetradeca-2,7-dienoyl-[acp] + NADH + H+  ->  (7Z)-tetradec-7-enoyl-[acp] + NAD+

(MALEYLACETATE-REDUCTASE-RXN-3-KETO-ADIPATE/NAD//CPD-294/NADH/PROTON.40. *instantiated*)	2-maleylacetate + NADH + H+  ->  3-oxoadipate + NAD+

(MALEYLACETATE-REDUCTASE-RXN-3-KETO-ADIPATE/NADP//CPD-294/NADPH/PROTON.42. *instantiated*)	2-maleylacetate + NADPH + H+  ->  3-oxoadipate + NADP+

(GTPPYPHOSKIN-RXN)	ATP + GTP  ->  pppGpp + AMP + H+

(RXN-3641)	3-oxoadipyl-CoA + coenzyme A  ->  succinyl-CoA + acetyl-CoA

(RXN-11811 *spontaneous*)	ammonia + H+  ->  ammonium

(GLYCERATE-DEHYDROGENASE-RXN)	D-glycerate + NAD+  ->  hydroxypyruvate + NADH + H+

(GLYCERATE-DEHYDROGENASE-RXN)	hydroxypyruvate + NADH + H+  ->  D-glycerate + NAD+

(RXN-16637)	L-cysteinyl-[tRNACys] + H2O  ->  L-cysteine + tRNACys + H+

(RXN-16637)	L-cysteine + tRNACys + H+  ->  L-cysteinyl-[tRNACys] + H2O

(RXN-12621)	an [L-cysteine desulfurase]-S-sulfanyl-L-cysteine + carboxy-adenylated-[ThiS sulfur-carrier protein] + reduced two electron carrier  ->  thiocarboxylated-[ThiS sulfur-carrier protein] + an [L-cysteine desulfurase]-L-cysteine + AMP + oxidized electron carrier + 2 H+

(RXN-11852)	acetyl-CoA + (3Z)-hex-3-en-1-ol  ->  (3Z)-hex-3-en-1-yl acetate + coenzyme A

(RXN-11852)	(3Z)-hex-3-en-1-yl acetate + coenzyme A  ->  acetyl-CoA + (3Z)-hex-3-en-1-ol

(RXN66-1)	ethanol + hydrogen peroxide  ->  acetaldehyde + 2 H2O

(RXN0-300)	D-glycerate + NADP+  ->  hydroxypyruvate + NADPH + H+

(RXN0-300)	hydroxypyruvate + NADPH + H+  ->  D-glycerate + NADP+

(GLU6PDEHYDROG-RXN-ALPHA-GLC-6-P/NADP//D-6-P-GLUCONO-DELTA-LACTONE/NADPH/PROTON.61. *instantiated*)	alpha-D-glucose 6-phosphate + NADP+  ->  6-phospho D-glucono-1,5-lactone + NADPH + H+

(GLU6PDEHYDROG-RXN-GLC-6-P/NADP//D-6-P-GLUCONO-DELTA-LACTONE/NADPH/PROTON.55. *instantiated*)	beta-D-glucose 6-phosphate + NADP+  ->  6-phospho D-glucono-1,5-lactone + NADPH + H+

(RXN-12587)	an [L-cysteine desulfurase]-S-sulfanyl-L-cysteine + unsulfurated [sulfur carrier]  ->  an [L-cysteine desulfurase]-L-cysteine + sulfurated [sulfur carrier]

(RXN-12587)	an [L-cysteine desulfurase]-L-cysteine + sulfurated [sulfur carrier]  ->  an [L-cysteine desulfurase]-S-sulfanyl-L-cysteine + unsulfurated [sulfur carrier]

(ARGINASE-RXN)	L-arginine + H2O  ->  urea + L-ornithine

(DARAB5PISOM-RXN)	D-ribulose 5-phosphate  ->  aldehydo-D-arabinose 5-phosphate

(ADENYLYLSULFKIN-RXN)	adenosine 5'-phosphosulfate + ATP  ->  3'-phosphoadenylyl-sulfate + ADP + H+

(S-2-METHYLMALATE-DEHYDRATASE-RXN)	mesaconate + H2O  ->  (S)-citramalate

(RXN-13329 *spontaneous*)	2-iminoacetate + H+ + H2O  ->  glyoxylate + ammonium

(DEPHOSPHOCOAKIN-RXN)	3'-dephospho-CoA + ATP  ->  coenzyme A + ADP + H+

(RXN-14274)	3-oxododecanoyl-CoA + coenzyme A  ->  decanoyl-CoA + acetyl-CoA

(RXN-14117)	phosphoenolpyruvate + GDP + H+  ->  pyruvate + GTP

(AGMATIN-RXN)	agmatine + H2O  ->  urea + putrescine

(RXN-12720)	selenate + ATP + H+  ->  adenosine 5'-phosphoselenate + diphosphate

(RXN-14325)	ATP + UTP + ammonia  ->  ADP + CTP + phosphate + H+

(RXN0-1132)	a [pyruvate dehydrogenase E2 protein] N6-dihydrolipoyl-L-lysine + NAD+  ->  a [pyruvate dehydrogenase E2 protein] N6-lipoyl-L-lysine + NADH + H+

(RXN-8976)	Und-PP-Mur2Ac-L-Ala-gamma-D-Glu-L-Lys-D-Ala-D-Ala + UDP-N-acetyl-alpha-D-glucosamine  ->  Und-PP-beta-D-GlcNAc-(1->4)-MurNAc-L-Ala-gamma-D-Glu-L-Lys-D-Ala-D-Ala + UDP + H+

(RXN-8976)	Und-PP-beta-D-GlcNAc-(1->4)-MurNAc-L-Ala-gamma-D-Glu-L-Lys-D-Ala-D-Ala + UDP + H+  ->  Und-PP-Mur2Ac-L-Ala-gamma-D-Glu-L-Lys-D-Ala-D-Ala + UDP-N-acetyl-alpha-D-glucosamine

(BRANCHED-CHAINAMINOTRANSFERLEU-RXN)	L-leucine + 2-oxoglutarate  ->  L-glutamate + 4-methyl-2-oxopentanoate

(RXN-22286 *spontaneous*)	hydrosulfide + H+  ->  hydrogen sulfide

(RXN-22286 *spontaneous*)	hydrogen sulfide  ->  hydrosulfide + H+

(RXN-20675)	(2E)-dodec-2-enoyl-CoA + H2O  ->  (S)-3-hydroxydodecanoyl-CoA

(RXN0-722)	UDP + reduced NrdH glutaredoxin-like protein  ->  dUDP + oxidized NrdH glutaredoxin-like protein + H2O

(RXN-15129)	L-cysteine  ->  2-aminoprop-2-enoate + hydrogen sulfide

(RXN-11564)	carboxyspermidine + H+  ->  spermidine + CO2

(RXN-11564)	spermidine + CO2  ->  carboxyspermidine + H+

(ACETYLORNTRANSAM-RXN)	N-acetyl-L-ornithine + 2-oxoglutarate  ->  L-glutamate + N-acetyl-L-glutamate 5-semialdehyde

(ACETYLORNTRANSAM-RXN)	L-glutamate + N-acetyl-L-glutamate 5-semialdehyde  ->  N-acetyl-L-ornithine + 2-oxoglutarate

(RXN-9087)	acryloyl-CoA + NADPH + H+  ->  propanoyl-CoA + NADP+

(NAG6PDEACET-RXN-N-ACETYL-D-GLUCOSAMINE-6-P/WATER//CPD-13469/ACET.49. *instantiated*)	N-acetyl-D-glucosamine 6-phosphate + H2O  ->  alpha-D-glucosamine 6-phosphate + acetate

(NAG6PDEACET-RXN)	N-acetyl-D-glucosamine 6-phosphate + H2O  ->  D-glucosamine 6-phosphate + acetate

(IPPISOM-RXN)	3-methylbut-3-en-1-yl diphosphate  ->  prenyl diphosphate

(IPPISOM-RXN)	prenyl diphosphate  ->  3-methylbut-3-en-1-yl diphosphate

(RXN-17776)	(2E,9Z)-octadecenoyl-CoA + H2O  ->  (S)-3-hydroxy-(9Z)-octadecenoyl-CoA

(RXN-12567)	(2E)-hexenoyl-CoA + H2O  ->  (S)-3-hydroxyhexanoyl-CoA

(RXN-9663)	(2E)-hexadec-2-enoyl-[acp] + NADH + H+  ->  palmitoyl-[acp] + NAD+

(HISTOLDEHYD-RXN)	histidinal + NADH + H+  ->  histidinol + NAD+

(ORNITHINE-RACEMASE-RXN)	L-ornithine  ->  D-ornithine

(ORNITHINE-RACEMASE-RXN)	D-ornithine  ->  L-ornithine

(HYPOXANPRIBOSYLTRAN-RXN)	hypoxanthine + 5-phospho-alpha-D-ribose 1-diphosphate  ->  IMP + diphosphate

(ALDOSE1EPIM-RXN)	beta-D-galactopyranose  ->  alpha-D-galactopyranose

(RXN-19774)	an (S-diacyl-sn-glyceryl)-L-cysteinyl-[apolipoprotein] + a phosphatidylethanolamine  ->  an N-acyl-(S-diacyl-sn-glyceryl)-L-cysteinyl-[lipoprotein] + a 1-lysophosphatidylethanolamine + H+

(RXN-19774)	an N-acyl-(S-diacyl-sn-glyceryl)-L-cysteinyl-[lipoprotein] + a 1-lysophosphatidylethanolamine + H+  ->  an (S-diacyl-sn-glyceryl)-L-cysteinyl-[apolipoprotein] + a phosphatidylethanolamine

(RXN-10654)	(5Z)-dodec-5-enoyl-[acp] + malonyl-[acp] + H+  ->  (7Z)-3-oxotetradec-7-enoyl-[acp] + acyl-carrier protein + CO2

(FUMARYLACETOACETASE-RXN)	4-fumaryl-acetoacetate + H2O  ->  fumarate + acetoacetate + H+

(THYM-PHOSPH-RXN)	thymidine + phosphate  ->  2-deoxy-alpha-D-ribose 1-phosphate + thymine

(THYM-PHOSPH-RXN)	2-deoxy-alpha-D-ribose 1-phosphate + thymine  ->  thymidine + phosphate

(OHBUTYRYL-COA-EPIM-RXN-CPD-650//S-3-HYDROXYBUTANOYL-COA.33. *instantiated*)	(3R)-3-hydroxybutanoyl-CoA  ->  (S)-3-hydroxybutanoyl-CoA

(OHBUTYRYL-COA-EPIM-RXN-CPD-650//S-3-HYDROXYBUTANOYL-COA.33. *instantiated*)	(S)-3-hydroxybutanoyl-CoA  ->  (3R)-3-hydroxybutanoyl-CoA

(OHBUTYRYL-COA-EPIM-RXN-POLYMER-INST-D-3-HYDROXYACYL-COA-C2-H4//OH-HEXANOYL-COA.56. *instantiated*)	POLYMER-INST-D-3-HYDROXYACYL-COA-C2-H4  ->  (S)-3-hydroxyhexanoyl-CoA

(OHBUTYRYL-COA-EPIM-RXN-POLYMER-INST-D-3-HYDROXYACYL-COA-C2-H4//OH-HEXANOYL-COA.56. *instantiated*)	(S)-3-hydroxyhexanoyl-CoA  ->  POLYMER-INST-D-3-HYDROXYACYL-COA-C2-H4

(OHBUTYRYL-COA-EPIM-RXN-POLYMER-INST-D-3-HYDROXYACYL-COA-C4-H8//CPD-22313.50. *instantiated*)	POLYMER-INST-D-3-HYDROXYACYL-COA-C4-H8  ->  (3S)-3-hydroxyoctanoyl-CoA

(OHBUTYRYL-COA-EPIM-RXN-POLYMER-INST-D-3-HYDROXYACYL-COA-C4-H8//CPD-22313.50. *instantiated*)	(3S)-3-hydroxyoctanoyl-CoA  ->  POLYMER-INST-D-3-HYDROXYACYL-COA-C4-H8

(OHBUTYRYL-COA-EPIM-RXN-POLYMER-INST-D-3-HYDROXYACYL-COA-C6-H12//CPD0-2244.51. *instantiated*)	POLYMER-INST-D-3-HYDROXYACYL-COA-C6-H12  ->  (S)-3-hydroxydecanoyl-CoA

(OHBUTYRYL-COA-EPIM-RXN-POLYMER-INST-D-3-HYDROXYACYL-COA-C6-H12//CPD0-2244.51. *instantiated*)	(S)-3-hydroxydecanoyl-CoA  ->  POLYMER-INST-D-3-HYDROXYACYL-COA-C6-H12

(OHBUTYRYL-COA-EPIM-RXN-POLYMER-INST-D-3-HYDROXYACYL-COA-C8-H16//CPD0-2107.51. *instantiated*)	POLYMER-INST-D-3-HYDROXYACYL-COA-C8-H16  ->  (S)-3-hydroxydodecanoyl-CoA

(OHBUTYRYL-COA-EPIM-RXN-POLYMER-INST-D-3-HYDROXYACYL-COA-C8-H16//CPD0-2107.51. *instantiated*)	(S)-3-hydroxydodecanoyl-CoA  ->  POLYMER-INST-D-3-HYDROXYACYL-COA-C8-H16

(OHBUTYRYL-COA-EPIM-RXN-POLYMER-INST-D-3-HYDROXYACYL-COA-C10-H20//POLYMER-INST-L-3-HYDROXYACYL-COA-C10-H20.83. *instantiated*)	POLYMER-INST-D-3-HYDROXYACYL-COA-C10-H20  ->  POLYMER-INST-L-3-HYDROXYACYL-COA-C10-H20

(OHBUTYRYL-COA-EPIM-RXN-POLYMER-INST-D-3-HYDROXYACYL-COA-C10-H20//POLYMER-INST-L-3-HYDROXYACYL-COA-C10-H20.83. *instantiated*)	POLYMER-INST-L-3-HYDROXYACYL-COA-C10-H20  ->  POLYMER-INST-D-3-HYDROXYACYL-COA-C10-H20

(OHBUTYRYL-COA-EPIM-RXN-POLYMER-INST-D-3-HYDROXYACYL-COA-C12-H24//POLYMER-INST-L-3-HYDROXYACYL-COA-C12-H24.83. *instantiated*)	POLYMER-INST-D-3-HYDROXYACYL-COA-C12-H24  ->  POLYMER-INST-L-3-HYDROXYACYL-COA-C12-H24

(OHBUTYRYL-COA-EPIM-RXN-POLYMER-INST-D-3-HYDROXYACYL-COA-C12-H24//POLYMER-INST-L-3-HYDROXYACYL-COA-C12-H24.83. *instantiated*)	POLYMER-INST-L-3-HYDROXYACYL-COA-C12-H24  ->  POLYMER-INST-D-3-HYDROXYACYL-COA-C12-H24

(OHBUTYRYL-COA-EPIM-RXN-CPD-10261//POLYMER-INST-L-3-HYDROXYACYL-COA-C14-H28.52. *instantiated*)	(3R)-3-hydroxyoctadecanoyl-CoA  ->  POLYMER-INST-L-3-HYDROXYACYL-COA-C14-H28

(OHBUTYRYL-COA-EPIM-RXN-CPD-10261//POLYMER-INST-L-3-HYDROXYACYL-COA-C14-H28.52. *instantiated*)	POLYMER-INST-L-3-HYDROXYACYL-COA-C14-H28  ->  (3R)-3-hydroxyoctadecanoyl-CoA

(OHBUTYRYL-COA-EPIM-RXN-POLYMER-INST-D-3-HYDROXYACYL-COA-C16-H32//POLYMER-INST-L-3-HYDROXYACYL-COA-C16-H32.83. *instantiated*)	POLYMER-INST-D-3-HYDROXYACYL-COA-C16-H32  ->  POLYMER-INST-L-3-HYDROXYACYL-COA-C16-H32

(OHBUTYRYL-COA-EPIM-RXN-POLYMER-INST-D-3-HYDROXYACYL-COA-C16-H32//POLYMER-INST-L-3-HYDROXYACYL-COA-C16-H32.83. *instantiated*)	POLYMER-INST-L-3-HYDROXYACYL-COA-C16-H32  ->  POLYMER-INST-D-3-HYDROXYACYL-COA-C16-H32

(RXN-22287 *spontaneous*)	hydrosulfide  ->  S2- + H+

(RXN-22287 *spontaneous*)	S2- + H+  ->  hydrosulfide

(THIOREDOXIN-RXN)	reduced thioredoxin + oxidized electron carrier  ->  oxidized thioredoxin + reduced two electron carrier

(MALONYL-COA-ACP-TRANSACYL-RXN)	malonyl-CoA + acyl-carrier protein  ->  coenzyme A + malonyl-[acp]

(RXN-9661)	(2E)-dodec-2-enoyl-[acp] + NADH + H+  ->  dodecanoyl-[acp] + NAD+

(UREASE-RXN)	urea + 2 H+ + H2O  ->  2 ammonium + CO2

(RXN-10674)	AsbD acyl-carrier protein + coenzyme A  ->  adenosine 3',5'-bisphosphate + holo-[AsbD acyl-carrier protein] + H+

(RXN-10674)	adenosine 3',5'-bisphosphate + holo-[AsbD acyl-carrier protein] + H+  ->  AsbD acyl-carrier protein + coenzyme A

(GLYOXIII-RXN)	methylglyoxal + H2O  ->  (R)-lactate + H+

(ALARACECAT-RXN)	L-alanine  ->  D-alanine

(ALARACECAT-RXN)	D-alanine  ->  L-alanine

(RXN-14026)	CMP + H2O  ->  cytidine + phosphate

(RXN-19916)	D-sedoheptulose 1-phosphate  ->  glycerone phosphate + D-erythrose

(RXN-19916)	glycerone phosphate + D-erythrose  ->  D-sedoheptulose 1-phosphate

(UDPREDUCT-RXN)	UDP + reduced thioredoxin  ->  dUDP + oxidized thioredoxin + H2O

(RXN0-308)	L-cysteine + an [L-cysteine desulfurase]-L-cysteine  ->  L-alanine + an [L-cysteine desulfurase]-S-sulfanyl-L-cysteine

(RXN-17775)	oleoyl-CoA + oxidized electron-transfer flavoprotein + H+  ->  (2E,9Z)-octadecenoyl-CoA + reduced electron-transfer flavoprotein

(METHYLASPARTATE-MUTASE-RXN)	L-glutamate  ->  (2S, 3S)-3-methylaspartate

(RXN-22315 *spontaneous*)	arsenous acid  ->  arsenite + H+

(RXN-22315 *spontaneous*)	arsenite + H+  ->  arsenous acid

(METHGLYSYN-RXN)	glycerone phosphate  ->  methylglyoxal + phosphate

(RXN-7593 *spontaneous*)	4-imidazolone-5-propanoate + reduced two electron carrier + dioxygen + H2O  ->  4-oxoglutaramate + ammonium + formate + oxidized electron carrier

(RXN-22742)	a tetrahydrofolate + sarcosine + dioxygen  ->  glycine + a 5,10-methylenetetrahydrofolate + hydrogen peroxide

(RXN-15124 *spontaneous*)	2-aminoprop-2-enoate  ->  2-iminopropanoate

(CDPREDUCT-RXN)	CDP + reduced thioredoxin  ->  dCDP + oxidized thioredoxin + H2O

(4OH2OXOGLUTARALDOL-RXN)	(4R)-4-hydroxy-2-oxoglutarate  ->  glyoxylate + pyruvate

(4OH2OXOGLUTARALDOL-RXN)	glyoxylate + pyruvate  ->  (4R)-4-hydroxy-2-oxoglutarate

(2.7.1.148-RXN)	4-(cytidine 5'-diphospho)-2-C-methyl-D-erythritol + ATP  ->  2-phospho-4-(cytidine 5'-diphospho)-2-C-methyl-D-erythritol + ADP + H+

(4.1.1.75-RXN)	5-guanidino-2-oxopentanoate + H+  ->  CO2 + 4-guanidinobutyraldehyde

(4.1.1.75-RXN)	CO2 + 4-guanidinobutyraldehyde  ->  5-guanidino-2-oxopentanoate + H+

(ATPSYN-RXN)	ATP + H2O + 3 H+[in]  ->  ADP + phosphate + 4 H+[out]

(ATPSYN-RXN)	ADP + phosphate + 4 H+[out]  ->  ATP + H2O + 3 H+[in]


====== End of the solution file.
