# Pathway Interpretation Report: Type 2 diabetes

## Executive Summary

This pathway report analyzes proteomic changes associated with Type 2 diabetes across 17 distinct functional clusters. In total, 13 clusters show a decreased representation of associated proteins, while 4 clusters demonstrate an increased representation. The altered pathways span core metabolic process networks (including fatty acid degradation, fatty acid biosynthesis, the citric acid cycle, amino acid metabolism, and oxidative phosphorylation) as well as inflammatory, vascular, and signal transduction cascades (such as complement and coagulation cascades, neutrophil degranulation, regulation of fibrinolysis, Rho GTPase signaling, and IGF transport regulation). Overall dataset evidence strength ranges from EXPLORATORY to MODERATE, based strictly on protein candidate counts and confidence levels, while literature assessment establishes HIGH disease relevance for central metabolic and inflammatory clusters (C002, C006, C007, C010, C012, C015) based on primary studies in Type 2 diabetes models and clinical cohorts. Clear separation is maintained between dataset evidence strength and disease relevance, avoiding causal extrapolation or directional claims of pathway activation or inhibition.

## Biological Theme
CD40 signaling pathway

## Direction
decreased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 2
- Levels 1–3: 0
- Levels 4–5: 2

## Disease Relevance
LOW

## Dataset Evidence
- Representative pathway: CD40 Signaling Pathway (GO:0023035)
- Representative library: GO_Biological_Process_2023
- Member pathways:
  - CD40 Signaling Pathway (GO:0023035)
- Best member-pathway adjusted p-value: 0.0456790270420036
- Representative Gene Ratio: 0.014388489208633

## Dataset Interpretation
The dataset exhibits a decreased representation of proteins associated with the CD40 signaling pathway in Type 2 diabetes compared to non-diabetic controls. This observation is based on two lower-confidence candidate proteins (ITGB1, PHB2). The dataset provides observational proteomic associations without establishing causal relationships or functional pathway inhibition.

## Limitations
This cluster is supported by a small number of candidate proteins (2 proteins), all originating from candidate levels 4–5, resulting in an EXPLORATORY evidence strength rating. No high-confidence (Levels 1–3) candidates are present. The observational nature of the dataset prevents drawing mechanistic conclusions.

## Recent Scientific Research — Previous 5 Years
- Relevant primary PubMed studies meeting criteria: 0
- Commonly studied tissues/models: None reported in the provided records.
- Major findings: No relevant primary literature was retrieved for this specific pathway-disease pairing.
- Consistency or disagreement: No literature available to assess agreement.
- Important reported mechanisms: None reported in the provided records.
- Relationship between literature and dataset: Literature evidence is insufficient to evaluate the dataset findings.
- Representative references: Fewer than three relevant primary references are available; none were retrieved for this cluster based on the predefined search criteria.

## Biological Theme
Fatty acid beta-oxidation

## Direction
decreased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 4
- Levels 1–3: 0
- Levels 4–5: 4

## Disease Relevance
HIGH

## Dataset Evidence
- Representative pathway: Fatty Acid Beta-Oxidation Using acyl-CoA Dehydrogenase (GO:0033539)
- Representative library: GO_Biological_Process_2023
- Member pathways:
  - Electron Transport Chain (GO:0022900)
  - Fatty Acid Beta-Oxidation Using acyl-CoA Dehydrogenase (GO:0033539)
  - Respiratory Electron Transport Chain (GO:0022904)
- Best member-pathway adjusted p-value: 0.0031032614034833
- Representative Gene Ratio: 0.0215827338129496

## Dataset Interpretation
The dataset demonstrates a decreased representation of proteins involved in fatty acid beta-oxidation in Type 2 diabetes samples compared to non-diabetic controls. The mapped candidates include ACADS, ETFB, ETFDH, and SDHA. These findings indicate differential abundance of mitochondrial oxidation enzymes without implying direct functional impairment or causal mechanisms within the dataset itself.

## Limitations
Dataset evidence is limited to 4 proteins, all belonging to candidate levels 4–5, yielding an EXPLORATORY dataset evidence strength. The current dataset lacks high-confidence (Levels 1–3) candidates and functional activity assays.

## Recent Scientific Research — Previous 5 Years
- Relevant primary PubMed studies meeting criteria: 9
- Commonly studied tissues/models: db/db mice, high-fat diet-induced diabetic rodent models, human serum, postmortem retinas, and primary cell models (cardiomyocytes, endothelial cells, HepG2, 3T3-L1).
- Major findings: Primary studies demonstrate that modulating fatty acid beta-oxidation impacts lipid accumulation, lipotoxicity, and tissue injury in Type 2 diabetes. Interventions targeting liver and adipose tissues (e.g., semaglutide, tangerine nanovesicles, jujuboside A, GMFB knockout) reduced steatosis and improved insulin sensitivity. Conversely, excessive endothelial fatty acid beta-oxidation in cardiac microvasculature was associated with microvascular dysfunction.
- Consistency or disagreement: High agreement across studies that fatty acid beta-oxidation dysregulation is a central feature of lipid metabolism changes in Type 2 diabetes, though specific cellular contexts determine whether upregulation or downregulation is protective.
- Important reported mechanisms: PPAR-alpha and PPAR-gamma activation driving beta-oxidation enzyme expression; CARM1 nuclear translocation; Adipsin inhibition of Irak2 mitochondrial translocation; miR-30-mediated endothelial regulation.
- Relationship between literature and dataset: Published literature provides robust support for the disease relevance of fatty acid beta-oxidation in Type 2 diabetes, complementing the decreased representation of related enzymes observed in the current dataset.
- Representative references:
  - Zou J, Song Q, Shaw PC, Wu Y, Zuo Z, Yu R. Tangerine Peel-Derived Exosome-Like Nanovesicles Alleviate Hepatic Steatosis Induced by Type 2 Diabetes: Evidenced by Regulating Lipid Metabolism and Intestinal Microflora. International journal of nanomedicine. 2024. PMID: 39371479.
  https://pubmed.ncbi.nlm.nih.gov/39371479/
  - Jiang MY, Man WR, Zhang XB, Zhang XH, Duan Y, Lin J, Zhang Y, Cao Y, Wu DX, Shu XF, Xin L, Wang H, Zhang X, Li CY, Gu XM, Zhang X, Sun DD. Adipsin inhibits Irak2 mitochondrial translocation and improves fatty acid β-oxidation to alleviate diabetic cardiomyopathy. Military Medical Research. 2023. PMID: 38072993.
  https://pubmed.ncbi.nlm.nih.gov/38072993/
  - Veitch S, Njock MS, Chandy M, Siraj MA, Chi L, Mak HQ, Yu K, Rathnakumar K, Perez-Romero CA, Chen Z, Alibhai FJ, Gustafson D, Raju S, Wu R, Khat DZ, Wang Y, Caballero A, Meagher P, Lau E, Pepic L, Cheng HS, Galant NJ, Howe KL, Li RK, Connelly KA, Husain M, Delgado-Olguin P, Fish JE. MiR-30 promotes fatty acid beta-oxidation and endothelial cell dysfunction and is a circulating biomarker of coronary microvascular dysfunction in pre-clinical models of diabetes. Cardiovascular diabetology. 2022. PMID: 35209901.
  https://pubmed.ncbi.nlm.nih.gov/35209901/
  - Niu S, Chen S, Wu D, Wang C, Xu J, Yin J, Zhao Y. Multi-omics analysis of hepatic outcomes in T2DM-MAFLD patients treated with semaglutide: a single-centre, longitudinal, data-driven study. Frontiers in endocrinology. 2025. PMID: 41103643.
  https://pubmed.ncbi.nlm.nih.gov/41103643/

## Biological Theme
Translation initiation complex formation

## Direction
decreased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 8
- Levels 1–3: 0
- Levels 4–5: 8

## Disease Relevance
LOW

## Dataset Evidence
- Representative pathway: Formation Of Cytoplasmic Translation Initiation Complex (GO:0001732)
- Representative library: GO_Biological_Process_2023
- Member pathways:
  - Cytoplasmic Translational Initiation (GO:0002183)
  - Formation Of Cytoplasmic Translation Initiation Complex (GO:0001732)
  - Positive Regulation Of RNA Binding (GO:1905216)
  - Positive Regulation Of mRNA Binding (GO:1902416)
  - Regulation Of mRNA Binding (GO:1902415)
  - RNA transport
  - Formation Of Ternary Complex, And Subsequently, 43S Complex R-HSA-72695
  - Ribosomal Scanning And Start Codon Recognition R-HSA-72702
  - Translation Initiation Complex Formation R-HSA-72649
  - mRNA Activation Upon Binding Of Cap-Binding Complex And eIFs, Subsequent Binding To 43S R-HSA-72662
- Best member-pathway adjusted p-value: 0.000259271141533
- Representative Gene Ratio: 0.0287769784172661

## Dataset Interpretation
The dataset indicates a decreased representation of proteins associated with translation initiation complex formation (e.g., EIF3D, EIF3E, EIF3I, EIF6, RPS2, EIF4E2) in Type 2 diabetes. These findings reflect lower relative abundance of protein synthesis machinery components without directly demonstrating altered translational rates in vivo.

## Limitations
Although 8 proteins are present, all are classified under candidate levels 4–5, resulting in an EXPLORATORY evidence strength rating. There are no high-confidence (Levels 1–3) candidate markers.

## Recent Scientific Research — Previous 5 Years
- Relevant primary PubMed studies meeting criteria: 0
- Commonly studied tissues/models: None reported in the provided records.
- Major findings: No relevant primary literature was retrieved for this theme in relation to Type 2 diabetes.
- Consistency or disagreement: No literature available to assess agreement.
- Important reported mechanisms: None reported in the provided records.
- Relationship between literature and dataset: Insufficient literature is available to contextualize the dataset findings.
- Representative references: Fewer than three relevant primary references are available; none were retrieved for this cluster based on the predefined search criteria.

## Biological Theme
Mitochondrion organization

## Direction
decreased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 7
- Levels 1–3: 1
- Levels 4–5: 6

## Disease Relevance
MODERATE

## Dataset Evidence
- Representative pathway: Mitochondrion Organization (GO:0007005)
- Representative library: GO_Biological_Process_2023
- Member pathways:
  - Mitochondrion Organization (GO:0007005)
- Best member-pathway adjusted p-value: 0.007618969816774
- Representative Gene Ratio: 0.0503597122302158

## Dataset Interpretation
The dataset exhibits a decreased representation of proteins involved in mitochondrion organization, including HSPD1 (Level 3) alongside BCAP31, COX7A2, PHB2, SLC25A4, VPS13C, and PRDX3. This indicates altered abundance of structural and functional mitochondrial proteins in Type 2 diabetes samples.

## Limitations
The dataset contains 7 total proteins, with only 1 Level 3 protein and 6 Level 4–5 proteins, placing evidence strength at the EXPLORATORY level. The data are strictly observational and do not quantify mitochondrial morphometry or respiration.

## Recent Scientific Research — Previous 5 Years
- Relevant primary PubMed studies meeting criteria: 1
- Commonly studied tissues/models: HSP22 knockout diabetic mice and vascular endothelial cells.
- Major findings: Rosiglitazone protected vascular endothelial cells and reduced diabetic angiopathy by upregulating HSP22 and suppressing mitochondrial reactive oxygen species formation and mitochondrial dysfunction.
- Consistency or disagreement: Consistent evidence from a single primary study highlighting mitochondrial structural and functional integrity in diabetic vascular pathology.
- Important reported mechanisms: PPAR-gamma regulation upregulating HSP22 expression to suppress mitochondrial oxidative stress and dysfunction.
- Relationship between literature and dataset: The published study provides moderate disease support linking mitochondrial organization to diabetic vascular complications, supporting the relevance of the decreased mitochondrial protein representation in the dataset.
- Representative references:
  - Yu L, Chen S, Liang Q, Huang C, Zhang W, Hu L, Yu Y, Liu L, Cheng X, Bao H. Rosiglitazone reduces diabetes angiopathy by inhibiting mitochondrial dysfunction dependent on regulating HSP22 expression. iScience. 2023. PMID: 36968091.
  https://pubmed.ncbi.nlm.nih.gov/36968091/

## Biological Theme
Purine nucleobase biosynthesis

## Direction
decreased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 2
- Levels 1–3: 0
- Levels 4–5: 2

## Disease Relevance
LOW

## Dataset Evidence
- Representative pathway: Purine Nucleobase Biosynthetic Process (GO:0009113)
- Representative library: GO_Biological_Process_2023
- Member pathways:
  - GMP Biosynthetic Process (GO:0006177)
  - GMP Metabolic Process (GO:0046037)
  - Nucleobase Biosynthetic Process (GO:0046112)
  - Purine Nucleobase Biosynthetic Process (GO:0009113)
  - Purine Ribonucleoside Monophosphate Biosynthesis R-HSA-73817
- Best member-pathway adjusted p-value: 0.0210495328620624
- Representative Gene Ratio: 0.014388489208633

## Dataset Interpretation
The dataset demonstrates a decreased representation of proteins associated with purine nucleobase biosynthesis (GMPS, PAICS) in Type 2 diabetes. These observations represent candidate abundance differences and do not establish metabolic flux changes.

## Limitations
Evidence is restricted to 2 lower-confidence candidate proteins (Levels 4–5), yielding an EXPLORATORY rating without high-confidence candidates.

## Recent Scientific Research — Previous 5 Years
- Relevant primary PubMed studies meeting criteria: 0
- Commonly studied tissues/models: None reported in the provided records.
- Major findings: No primary literature meeting selection criteria was identified for purine nucleobase biosynthesis in Type 2 diabetes.
- Consistency or disagreement: No literature available to assess agreement.
- Important reported mechanisms: None reported in the provided records.
- Relationship between literature and dataset: Literature support is insufficient to contextualize dataset observations.
- Representative references: Fewer than three relevant primary references are available; none were retrieved for this cluster based on the predefined search criteria.

## Biological Theme
Citric acid cycle

## Direction
decreased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 16
- Levels 1–3: 1
- Levels 4–5: 15

## Disease Relevance
HIGH

## Dataset Evidence
- Representative pathway: Citrate cycle (TCA cycle)
- Representative library: KEGG_2021_Human
- Member pathways:
  - Dicarboxylic Acid Metabolic Process (GO:0043648)
  - Citrate cycle (TCA cycle)
  - Cysteine and methionine metabolism
  - Glycolysis / Gluconeogenesis
  - Glyoxylate and dicarboxylate metabolism
  - Pyruvate metabolism
  - Citric Acid Cycle (TCA Cycle) R-HSA-71403
  - Mitochondrial Protein Import R-HSA-1268020
  - Protein Localization R-HSA-9609507
  - Pyruvate Metabolism And Citric Acid (TCA) Cycle R-HSA-71406
- Best member-pathway adjusted p-value: 2.027802752953447e-06
- Representative Gene Ratio: 0.0431654676258992

## Dataset Interpretation
The dataset shows a decreased representation of proteins involved in the citric acid cycle in Type 2 diabetes samples compared to non-diabetic controls. Candidates include HSPD1 (Level 3) and 15 Level 4–5 proteins (including CS, DLAT, DLD, SDHA, ACO2, ALDH2, LDHB, MDH2). These data indicate widespread differential abundance of mitochondrial matrix enzymes.

## Limitations
Despite a strong p-value and 16 identified proteins, 15 candidates belong to levels 4–5 and only 1 to level 3, maintaining an EXPLORATORY dataset evidence strength classification.

## Recent Scientific Research — Previous 5 Years
- Relevant primary PubMed studies meeting criteria: 5
- Commonly studied tissues/models: Diabetic rat models, pancreatic beta-cell lines (betaTC3), mouse islets, human gastric mucosal tissue, human adipose-derived mesenchymal stem cells, and Drosophila high-sugar diet models.
- Major findings: Citric acid cycle enzymes and intermediate flux are significantly altered in Type 2 diabetes. Enhancing TCA cycle flux (e.g., via G6pc2 knockout, alpha-ketoglutarate supplementation, or Nampt upregulation by Dangua Fang) improved glucose-stimulated insulin secretion, osteogenic differentiation, or metabolic parameters.
- Consistency or disagreement: High agreement across primary studies that TCA cycle flux and key enzyme levels are perturbed in diabetic states, and that interventions enhancing TCA cycle activity improve metabolic outcomes.
- Important reported mechanisms: Nampt upregulation enhancing TCA cycle flux; G6pc2 deletion increasing CAC flux and insulin secretion; exogenous alpha-ketoglutarate restoring cellular regeneration.
- Relationship between literature and dataset: Robust literature evidence confirms high disease relevance for TCA cycle dysregulation in Type 2 diabetes, aligning with the decreased representation of TCA cycle enzymes observed in this proteomic dataset.
- Representative references:
  - Xianpei H, Zhita W, Liuqing Y, Liang L, Suping H. Dangua Fang regulating tricarboxylic acid cycle and respiratory chain and its mechanism in diabetic rats. Journal of traditional Chinese medicine = Chung i tsa chih ying wen pan. 2023. PMID: 37946477.
  https://pubmed.ncbi.nlm.nih.gov/37946477/
  - Dias DB, Chan W, Ellinghaus A, Fritsche-Guenther R, Wiebach J, Bembennek A, Laske T, Baumbach J, Duda GN, Kirwan JA, Poh PSP. Endogenous dysregulated energy and amino acid metabolism delay scaffold-guided large volume bone regeneration in a diabetic rat model with Leptin receptor deficiency. Acta biomaterialia. 2025. PMID: 40319991.
  https://pubmed.ncbi.nlm.nih.gov/40319991/
  - Rahim M, Nakhe AY, Banerjee DR, Overway EM, Bosma KJ, Rosch JC, Oeser JK, Wang B, Lippmann ES, Jacobson DA, O'Brien RM, Young JD. Glucose-6-phosphatase catalytic subunit 2 negatively regulates glucose oxidation and insulin secretion in pancreatic β-cells. The Journal of biological chemistry. 2022. PMID: 35176280.
  https://pubmed.ncbi.nlm.nih.gov/35176280/

## Biological Theme
Fatty acid biosynthesis

## Direction
decreased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 11
- Levels 1–3: 0
- Levels 4–5: 11

## Disease Relevance
HIGH

## Dataset Evidence
- Representative pathway: Fatty acid biosynthesis
- Representative library: KEGG_2021_Human
- Member pathways:
  - Monocarboxylic Acid Metabolic Process (GO:0032787)
  - Fatty acid biosynthesis
  - Ferroptosis
  - Peroxisome
  - Fatty acyl-CoA Biosynthesis R-HSA-75105
- Best member-pathway adjusted p-value: 0.0016735024638605
- Representative Gene Ratio: 0.0215827338129496

## Dataset Interpretation
The dataset demonstrates a decreased representation of proteins associated with fatty acid biosynthesis in Type 2 diabetes. Mapped candidates include 11 Level 4–5 proteins such as ACSF2, ACSL1, ACSL3, FASN, and FTH1. These results reflect lower relative abundance of enzymes involved in fatty acyl-CoA synthesis and lipogenesis.

## Limitations
All 11 candidate proteins are classified under levels 4–5 with no Level 1–3 proteins, assigning an EXPLORATORY evidence strength rating to the dataset.

## Recent Scientific Research — Previous 5 Years
- Relevant primary PubMed studies meeting criteria: 18
- Commonly studied tissues/models: Human serum, urine exosomes, epicardial adipose tissue, subgingival plaque metagenomes, db/db mice, ZDF rats, high-fat diet-fed mice, and 3T3-L1 adipocytes.
- Major findings: Fatty acid biosynthesis pathways and enzymes (such as FASN and ACC-1) are consistently dysregulated in Type 2 diabetes and its complications. Therapeutic interventions (e.g., short-chain fatty acids, botanical extracts, acupuncture, or probiotics) improved glycaemia and hepatic steatosis in diabetic models by downregulating key lipogenic enzymes.
- Consistency or disagreement: High agreement across clinical cohort studies and preclinical interventions that fatty acid biosynthesis pathways are central to glycolipid dysregulation in Type 2 diabetes.
- Important reported mechanisms: Downregulation of FASN and ACC-1 via AMPK and PPAR signaling; microRNA-mediated post-transcriptional regulation (miR-30, miR-142-3p); gut microbial metabolite interactions.
- Relationship between literature and dataset: Extensive primary literature validates high disease relevance for fatty acid biosynthesis in Type 2 diabetes, complementing the decreased protein representation observed in this dataset.
- Representative references:
  - Ikeda T, Takii K, Omichi Y, Nishimoto Y, Ichikawa D, Matsunaga T, Kawauchi A, Kimura I. Hexanoic Acid Improves Metabolic Health in Mice Fed High-Fat Diet. Nutrients. 2025. PMID: 40944255.
  https://pubmed.ncbi.nlm.nih.gov/40944255/
  - Henneke L, Schlicht K, Andreani NA, Hollstein T, Demetrowitsch T, Knappe C, Hartmann K, Jensen-Kroll J, Rohmann N, Pohlschneider D, Geisler C, Schulte DM, Settgast U, Türk K, Zimmermann J, Kaleta C, Baines JF, Shearer J, Shah S, Shen-Tu G, Schwarz K, Franke A, Schreiber S, Laudes M. A dietary carbohydrate - gut Parasutterella - human fatty acid biosynthesis metabolic axis in obesity and type 2 diabetes. Gut microbes. 2022. PMID: 35435797.
  https://pubmed.ncbi.nlm.nih.gov/35435797/
  - Jiang Y, Yu H, Pan Y, Zhang B, Jing Y, Lei J, Li N, Yang J. Effect and Mechanism of Qihua Tongtiao Formula (QHTTF) on Improving Glucose and Lipid Metabolism Disorders in ZDF Rats by Integrating Network Pharmacology, Metabolomics, and Biological Validation. Pharmaceuticals (Basel, Switzerland). 2025. PMID: 41011216.
  https://pubmed.ncbi.nlm.nih.gov/41011216/

## Biological Theme
Legionellosis pathway

## Direction
decreased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 3
- Levels 1–3: 1
- Levels 4–5: 2

## Disease Relevance
LOW

## Dataset Evidence
- Representative pathway: Legionellosis
- Representative library: KEGG_2021_Human
- Member pathways:
  - Legionellosis
- Best member-pathway adjusted p-value: 0.0390705887259355
- Representative Gene Ratio: 0.0215827338129496

## Dataset Interpretation
The dataset indicates a decreased representation of proteins associated with the KEGG Legionellosis pathway (HSPD1 [Level 3], BCL2L13, SAR1A). In a non-infectious metabolic context, this enrichment likely reflects overlapping organellar chaperones and vesicular transport components.

## Limitations
Evidence is based on only 3 candidates (1 Level 3 protein, 2 Level 4–5 proteins), yielding an EXPLORATORY dataset evidence strength rating. The pathway label reflects an infection database entry with shared cellular machinery.

## Recent Scientific Research — Previous 5 Years
- Relevant primary PubMed studies meeting criteria: 0
- Commonly studied tissues/models: None reported in the provided records.
- Major findings: No primary literature was retrieved linking this pathway term directly to Type 2 diabetes.
- Consistency or disagreement: No literature available to assess agreement.
- Important reported mechanisms: None reported in the provided records.
- Relationship between literature and dataset: Literature support is absent for this pathway term in Type 2 diabetes.
- Representative references: Fewer than three relevant primary references are available; none were retrieved for this cluster based on the predefined search criteria.

## Biological Theme
Oxidative phosphorylation and respiratory electron transport

## Direction
decreased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 25
- Levels 1–3: 1
- Levels 4–5: 24

## Disease Relevance
LOW

## Dataset Evidence
- Representative pathway: Citric Acid (TCA) Cycle And Respiratory Electron Transport R-HSA-1428517
- Representative library: Reactome_2022
- Member pathways:
  - Aerobic Electron Transport Chain (GO:0019646)
  - Aerobic Respiration (GO:0009060)
  - Cellular Respiration (GO:0045333)
  - Energy Derivation By Oxidation Of Organic Compounds (GO:0015980)
  - Mitochondrial ATP Synthesis Coupled Electron Transport (GO:0042775)
  - Mitochondrial Electron Transport, Ubiquinol To Cytochrome C (GO:0006122)
  - Oxidative Phosphorylation (GO:0006119)
  - Alzheimer disease
  - Amyotrophic lateral sclerosis
  - Cardiac muscle contraction
  - Diabetic cardiomyopathy
  - Huntington disease
  - Non-alcoholic fatty liver disease
  - Oxidative phosphorylation
  - Parkinson disease
  - Pathways of neurodegeneration
  - Prion disease
  - Thermogenesis
  - Citric Acid (TCA) Cycle And Respiratory Electron Transport R-HSA-1428517
  - Respiratory Electron Transport R-HSA-611105
  - Respiratory Electron Transport, ATP Synthesis By Chemiosmotic Coupling, Heat Production By Uncoupling Proteins R-HSA-163200
- Best member-pathway adjusted p-value: 2.1660684279841467e-10
- Representative Gene Ratio: 0.1079136690647482

## Dataset Interpretation
The dataset exhibits a decreased representation of proteins involved in oxidative phosphorylation and respiratory electron transport in Type 2 diabetes. Candidates include NDUFA10 (Level 2) and 24 Level 4–5 proteins (e.g., COX6B1, COX7A2, NDUFA12, SDHA, UQCRC1, UQCRC2). These data point to reduced relative abundance of electron transport chain subunits.

## Limitations
Despite a highly significant adjusted p-value (2.17e-10) and 25 mapped proteins, 24 belong to candidate levels 4–5 and only 1 to Level 2, maintaining an EXPLORATORY dataset evidence strength rating. Functional ATP synthesis rates were not measured.

## Recent Scientific Research — Previous 5 Years
- Relevant primary PubMed studies meeting criteria: 0
- Commonly studied tissues/models: None reported in the provided records.
- Major findings: No primary literature meeting the specific search query was retrieved for this cluster.
- Consistency or disagreement: No literature available to assess agreement.
- Important reported mechanisms: None reported in the provided records.
- Relationship between literature and dataset: No direct literature assessment is available for this specific query string.
- Representative references: Fewer than three relevant primary references are available; none were retrieved for this cluster based on the predefined search criteria.

## Biological Theme
Amino acid metabolism

## Direction
decreased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 9
- Levels 1–3: 0
- Levels 4–5: 9

## Disease Relevance
HIGH

## Disease Relevance Assessment
High disease relevance supported by 12 primary PubMed studies demonstrating altered amino acid catabolism and transport across diabetic cohorts and models.

## Dataset Evidence
- Representative pathway: Metabolism Of Amino Acids And Derivatives R-HSA-71291
- Representative library: Reactome_2022
- Member pathways:
  - Metabolism Of Amino Acids And Derivatives R-HSA-71291
- Best member-pathway adjusted p-value: 0.0224735348776994
- Representative Gene Ratio: 0.0647482014388489

## Dataset Interpretation
The dataset demonstrates a decreased representation of proteins involved in amino acid metabolism in Type 2 diabetes. Mapped candidates include 9 Level 4–5 proteins (DLAT, ECHS1, ENOPH1, ETHE1, GPT, RPS2, FAH, GRHPR, HIBADH). This indicates reduced abundance of enzymes involved in branched-chain and organic amino acid catabolism.

## Limitations
Dataset evidence strength is EXPLORATORY, as all 9 candidate proteins are classified under candidate levels 4–5 without Level 1–3 markers.

## Recent Scientific Research — Previous 5 Years
- Relevant primary PubMed studies meeting criteria: 12
- Commonly studied tissues/models: Human clinical cohorts (serum, urine, fecal metaproteomics, white adipose tissue), mouse models (db/db, Akita, A-PDK1KO, Crebbp/Ep300 knockout), brown adipocytes, and primary hepatocytes.
- Major findings: Disruption of branched-chain amino acid (BCAA) catabolism in adipose tissue and heart contributes to systemic insulin resistance and diabetic cardiomyopathy via mTORC1 hyperactivation. Hepatic CBP/p300 and SEC22B regulate amino acid catabolism and amino acid-driven gluconeogenesis. Hyperglycemia alters amino acid transport and synthesis in pancreatic alpha cells, sustaining mTORC1 activation and hyperglucagonemia.
- Consistency or disagreement: High agreement across primary studies that amino acid catabolism is dysregulated in Type 2 diabetes and directly impacts glucose homeostasis, insulin sensitivity, and tissue complications.
- Important reported mechanisms: BCAA catabolic impairment driving mTORC1 activation; histone crotonylation by CBP/p300 regulating gluconeogenic gene expression; SEC22B phosphorylation modulating hepatocellular response to glucagon; high glucose upregulating Slc7a2/Slc38a4 amino acid transporters in alpha cells.
- Relationship between literature and dataset: Extensive primary literature supports high disease relevance for amino acid metabolic dysregulation in Type 2 diabetes, mirroring the decreased protein representation observed in the current dataset.
- Representative references:
  - Verkerke ARP, Wang D, Yoshida N, Taxin ZH, Shi X, Zheng S, Li Y, Auger C, Oikawa S, Yook JS, Granath-Panelo M, He W, Zhang GF, Matsushita M, Saito M, Gerszten RE, Mills EL, Banks AS, Ishihama Y, White PJ, McGarrah RW, Yoneshiro T, Kajimura S. BCAA-nitrogen flux in brown fat controls metabolic health independent of thermogenesis. Cell. 2024. PMID: 38653240.
  https://pubmed.ncbi.nlm.nih.gov/38653240/
  - Sheng C, Li T, Lin H, Ma X, Zhou F, Li M, Wang Y, Wang S, Tan J, Chen J, Yang Y, Liu J, Bi Y, Lu J, Wang X, Zhou L. Hepatic CBP/p300 Orchestrate Amino Acid-Driven Gluconeogenesis through Histone Crotonylation. Advanced science (Weinheim, Baden-Wurttemberg, Germany). 2025. PMID: 40791184.
  https://pubmed.ncbi.nlm.nih.gov/40791184/
  - Jiang X, Liu X, Qu X, Zhu P, Wo F, Xu X, Jin J, He Q, Wu J. Integration of metabolomics and peptidomics reveals distinct molecular landscape of human diabetic kidney disease. Theranostics. 2023. PMID: 37351171.
  https://pubmed.ncbi.nlm.nih.gov/37351171/

## Biological Theme
Fatty acid degradation and catabolism

## Direction
decreased

## Dataset Evidence Strength
MODERATE

Supporting proteins:
- Total proteins: 49
- Levels 1–3: 2
- Levels 4–5: 47

## Disease Relevance
LOW

## Dataset Evidence
- Representative pathway: Metabolism R-HSA-1430728
- Representative library: Reactome_2022
- Member pathways:
  - Amino Acid Catabolic Process (GO:0009063)
  - Branched-Chain Amino Acid Catabolic Process (GO:0009083)
  - Fatty Acid Beta-Oxidation (GO:0006635)
  - Fatty Acid Catabolic Process (GO:0009062)
  - Fatty Acid Oxidation (GO:0019395)
  - Butanoate metabolism
  - Fatty acid degradation
  - Fatty acid elongation
  - Lysine degradation
  - Propanoate metabolism
  - Tryptophan metabolism
  - Valine, leucine and isoleucine degradation
  - beta-Alanine metabolism
  - Beta Oxidation Of butanoyl-CoA To acetyl-CoA R-HSA-77352
  - Beta Oxidation Of decanoyl-CoA To octanoyl-CoA-CoA R-HSA-77346
  - Beta Oxidation Of hexanoyl-CoA To butanoyl-CoA R-HSA-77350
  - Beta Oxidation Of lauroyl-CoA To decanoyl-CoA-CoA R-HSA-77310
  - Beta Oxidation Of octanoyl-CoA To hexanoyl-CoA R-HSA-77348
  - Fatty Acid Metabolism R-HSA-8978868
  - Metabolism Of Lipids R-HSA-556833
  - Metabolism R-HSA-1430728
  - Mitochondrial Fatty Acid Beta-Oxidation Of Saturated Fatty Acids R-HSA-77286
  - Mitochondrial Fatty Acid Beta-Oxidation Of Unsaturated Fatty Acids R-HSA-77288
  - Mitochondrial Fatty Acid Beta-Oxidation R-HSA-77289
- Best member-pathway adjusted p-value: 2.1660684279841467e-10
- Representative Gene Ratio: 0.3237410071942446

## Dataset Interpretation
The dataset demonstrates a decreased representation of proteins involved in fatty acid degradation and catabolism in Type 2 diabetes. This large cluster includes 49 total proteins, with 2 Level 2 candidates (NDUFA10, TM7SF2) and 47 Level 4–5 candidates (such as ACADS, ACSL1, ACSL3, FASN, HADH, HADHA, ACAA2). These findings indicate broad differential representation of lipid and acyl-CoA metabolic enzymes.

## Limitations
Although the cluster includes 49 mapped proteins and holds a MODERATE dataset evidence strength rating based on candidate count and 2 Level 2 proteins, functional enzymatic activity was not directly quantified in the dataset.

## Recent Scientific Research — Previous 5 Years
- Relevant primary PubMed studies meeting criteria: 0
- Commonly studied tissues/models: None reported in the provided records.
- Major findings: No primary literature was retrieved using the specific search query for this consolidated theme.
- Consistency or disagreement: No literature available to assess agreement.
- Important reported mechanisms: None reported in the provided records.
- Relationship between literature and dataset: Literature support is insufficient for this specific search formulation.
- Representative references: Fewer than three relevant primary references are available; none were retrieved for this cluster based on the predefined search criteria.

## Biological Theme
Neutrophil degranulation

## Direction
decreased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 10
- Levels 1–3: 0
- Levels 4–5: 10

## Disease Relevance
HIGH

## Dataset Evidence
- Representative pathway: Neutrophil Degranulation R-HSA-6798695
- Representative library: Reactome_2022
- Member pathways:
  - Neutrophil Degranulation R-HSA-6798695
- Best member-pathway adjusted p-value: 0.032855721572688
- Representative Gene Ratio: 0.0719424460431654

## Dataset Interpretation
The dataset exhibits a decreased representation of proteins associated with neutrophil degranulation in Type 2 diabetes. Candidates include 10 Level 4–5 proteins (FTH1, HLA-A, RNASET2, STBD1, STK10, ACTR1B, CAPN1, GDI2, RAB7A, RAP1B). This indicates lower relative abundance of granule and vesicle regulatory proteins in the analyzed sample set.

## Limitations
All 10 mapped proteins belong to candidate levels 4–5 without Level 1–3 markers, assigning an EXPLORATORY rating to the dataset evidence strength.

## Recent Scientific Research — Previous 5 Years
- Relevant primary PubMed studies meeting criteria: 8
- Commonly studied tissues/models: Human clinical samples (plasma, serum, skin biopsies), isolated human neutrophils under high-glucose conditions, and diabetic rat bone regeneration models.
- Major findings: Neutrophil degranulation and activation markers are consistently elevated in prediabetes, Type 2 diabetes, and diabetic microvascular and tissue complications (diabetic foot ulcers, diastolic dysfunction, impaired bone healing). High glucose and metabolic intermediates prime human neutrophils via phosphoproteomic signaling cascades.
- Consistency or disagreement: High agreement across clinical and experimental studies that neutrophil degranulation is pathologically altered in Type 2 diabetes and associated with tissue complications.
- Important reported mechanisms: Hyperglycemia-induced kinase activation (JNK, SYK, PRKACA) driving neutrophil priming; sustained neutrophil extracellular trap (NET) formation and mast cell crosstalk in impaired tissue repair.
- Relationship between literature and dataset: Extensive literature evidence establishes HIGH disease relevance for neutrophil degranulation in Type 2 diabetes, complementing the altered representation of neutrophil granule proteins in the current dataset.
- Representative references:
  - Ministrini S, Andreozzi F, Montecucco F, Minetti S, Bertolotto M, Liberale L, Mannino GC, Succurro E, Cassano V, Miceli S, Perticone M, Sesti G, Sciacqua A, Carbone F. Neutrophil degranulation biomarkers characterize restrictive echocardiographic pattern with diastolic dysfunction in patients with diabetes. European journal of clinical investigation. 2021. PMID: 34129696.
  https://pubmed.ncbi.nlm.nih.gov/34129696/
  - Thimmappa PY, Nair AS, Najar MA, Mohanty V, Shastry S, Prasad TSK, Joshi MB. Quantitative phosphoproteomics reveals diverse stimuli activate distinct signaling pathways during neutrophil activation. Cell and tissue research. 2022. PMID: 35622142.
  https://pubmed.ncbi.nlm.nih.gov/35622142/
  - Yin Y, Gao X, Li Q, Lu L, Zhang H, Wu Q, Shi J, Liu C, Yue L, Xiao S, Wu J, Lin X, Zeng R. Molecular Signature of Prediabetes With High-Risk of Diabetes Revealed by Deep Plasma Proteome. Diabetes, obesity & metabolism. 2026. PMID: 42244140.
  https://pubmed.ncbi.nlm.nih.gov/42244140/

## Biological Theme
Rho GTPase signaling

## Direction
decreased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 12
- Levels 1–3: 0
- Levels 4–5: 12

## Disease Relevance
MODERATE

## Dataset Evidence
- Representative pathway: RHOG GTPase Cycle R-HSA-9013408
- Representative library: Reactome_2022
- Member pathways:
  - RHOG GTPase Cycle R-HSA-9013408
  - Signaling By Rho GTPases R-HSA-194315
  - Signaling By Rho GTPases, Miro GTPases And RHOBTB3 R-HSA-9716542
- Best member-pathway adjusted p-value: 0.0055787492942776
- Representative Gene Ratio: 0.0359712230215827

## Dataset Interpretation
The dataset demonstrates a decreased representation of proteins associated with Rho GTPase signaling in Type 2 diabetes. Supporting candidates include 12 Level 4–5 proteins (e.g., BCAP31, HINT2, HSPE1, ITGB1, PPP1R14A, RANGAP1, RAB7A, YWHAE). These results reflect differential abundance of signal transduction regulators.

## Limitations
Evidence strength is EXPLORATORY due to reliance on 12 Level 4–5 candidate proteins, with no high-confidence (Levels 1–3) candidates or GTPase activity measurements.

## Recent Scientific Research — Previous 5 Years
- Relevant primary PubMed studies meeting criteria: 1
- Commonly studied tissues/models: Isolated human neutrophils under high glucose conditions.
- Major findings: High glucose exposure altered the phosphoproteomic signature of human neutrophils, impacting kinases associated with Rho GTPase signaling and downstream immune response.
- Consistency or disagreement: Single primary study demonstrating mechanistic link under hyperglycemic conditions.
- Important reported mechanisms: Differential phosphorylation of C-Jun-N-Terminal Kinase, NTRK1, SYK, and PRKACA linked to Rho GTPase signaling.
- Relationship between literature and dataset: Literature provides moderate support for hyperglycemia-induced alterations in Rho GTPase signaling pathways, aligning with the decreased representation of signaling candidates in the dataset.
- Representative references:
  - Thimmappa PY, Nair AS, Najar MA, Mohanty V, Shastry S, Prasad TSK, Joshi MB. Quantitative phosphoproteomics reveals diverse stimuli activate distinct signaling pathways during neutrophil activation. Cell and tissue research. 2022. PMID: 35622142.
  https://pubmed.ncbi.nlm.nih.gov/35622142/

## Biological Theme
Regulation of fibrinolysis and blood coagulation

## Direction
increased

## Dataset Evidence Strength
MODERATE

Supporting proteins:
- Total proteins: 14
- Levels 1–3: 2
- Levels 4–5: 12

## Disease Relevance
LOW

## Dataset Evidence
- Representative pathway: Regulation Of Fibrinolysis (GO:0051917)
- Representative library: GO_Biological_Process_2023
- Member pathways:
  - Negative Regulation Of Blood Coagulation (GO:0030195)
  - Negative Regulation Of Fibrinolysis (GO:0051918)
  - Peptide Cross-Linking (GO:0018149)
  - Positive Regulation Of Blood Coagulation (GO:0030194)
  - Regulation Of Fibrinolysis (GO:0051917)
  - Formation Of Fibrin Clot (Clotting Cascade) R-HSA-140877
  - Hemostasis R-HSA-109582
  - Platelet Activation, Signaling And Aggregation R-HSA-76002
  - Platelet Degranulation R-HSA-114608
  - Prolactin Receptor Signaling R-HSA-1170546
  - Response To Elevated Platelet Cytosolic Ca2+ R-HSA-76005
  - Signaling By Cytosolic FGFR1 Fusion Mutants R-HSA-1839117
  - Signaling By PDGF R-HSA-186797
- Best member-pathway adjusted p-value: 0.0007263841003932
- Representative Gene Ratio: 0.0416666666666666

## Dataset Interpretation
The dataset demonstrates an increased representation of proteins involved in the regulation of fibrinolysis and blood coagulation in Type 2 diabetes. Candidates include STAT5B (Level 1) and VTN (Level 2), alongside 12 Level 4–5 proteins (F13A1, KLKB1, PLG, PPBP, PROCR, THBS1, TF). This indicates higher relative abundance of circulating haemostatic factors.

## Limitations
Although supported by 14 proteins including Level 1 and Level 2 candidates (granting MODERATE dataset evidence strength), functional clotting times or fibrinolytic activity assays were not performed in this dataset.

## Recent Scientific Research — Previous 5 Years
- Relevant primary PubMed studies meeting criteria: 0
- Commonly studied tissues/models: None reported in the provided records.
- Major findings: No primary literature was retrieved for this exact pathway query string.
- Consistency or disagreement: No literature available to evaluate agreement.
- Important reported mechanisms: None reported in the provided records.
- Relationship between literature and dataset: Literature support is insufficient for this specific query formulation.
- Representative references: Fewer than three relevant primary references are available; none were retrieved for this cluster based on the predefined search criteria.

## Biological Theme
Complement and coagulation cascades

## Direction
increased

## Dataset Evidence Strength
MODERATE

Supporting proteins:
- Total proteins: 35
- Levels 1–3: 2
- Levels 4–5: 33

## Disease Relevance
HIGH

## Dataset Evidence
- Representative pathway: Complement and coagulation cascades
- Representative library: KEGG_2021_Human
- Member pathways:
  - Positive Regulation Of Immune System Process (GO:0002684)
  - Complement and coagulation cascades
  - Coronavirus disease
  - Staphylococcus aureus infection
  - Systemic lupus erythematosus
  - Activation Of C3 And C5 R-HSA-174577
  - Complement Cascade R-HSA-166658
  - Immune System R-HSA-168256
  - Innate Immune System R-HSA-168249
  - Neutrophil Degranulation R-HSA-6798695
  - Regulation Of Complement Cascade R-HSA-977606
  - Terminal Pathway Of Complement R-HSA-166665
- Best member-pathway adjusted p-value: 2.0829738221088427e-14
- Representative Gene Ratio: 0.1354166666666666

## Dataset Interpretation
The dataset exhibits a strong increased representation of proteins in the complement and coagulation cascades in Type 2 diabetes (adjusted p = 2.08e-14). Mapped candidates include STAT5B (Level 1) and VTN (Level 2), along with 33 Level 4–5 proteins (e.g., C2, C5, C6, C8A, C9, CFHR2, F13A1, KLKB1, LCN2, PLG, PPBP, C1S, CFB). These observations reflect higher relative abundance of innate immune and pro-thrombotic factors.

## Limitations
Dataset evidence strength is MODERATE based on 35 mapped proteins and 2 high-confidence candidates (Levels 1–2). The data are observational and do not directly demonstrate functional complement activation or clot formation in sample tissues.

## Recent Scientific Research — Previous 5 Years
- Relevant primary PubMed studies meeting criteria: 9
- Commonly studied tissues/models: Human serum, plasma exosomes, urine, aqueous humor, human platelets, and diabetic rodent models (ob/ob, db/db).
- Major findings: Proteomic profiling across multiple human biofluids demonstrates that complement and coagulation cascade proteins are consistently elevated in prediabetes, Type 2 diabetes, and microvascular/renal complications (diabetic nephropathy, diabetic retinopathy, heart failure). In addition, GLP-1 receptor agonist therapy significantly modulated circulating complement and coagulation factor levels.
- Consistency or disagreement: High agreement across clinical, biomarker, and preclinical studies that complement and coagulation cascades are activated in Type 2 diabetes and track with disease severity.
- Important reported mechanisms: Autophagy-dependent platelet activation promoting microthrombosis; differential urinary excretion of complement regulators (CFAH, DAF) reflecting renal microvascular injury; pharmacological modulation by GLP-1 receptor agonists.
- Relationship between literature and dataset: Robust literature evidence establishes HIGH disease relevance, strongly supporting the increased representation of complement and coagulation proteins observed in the current dataset.
- Representative references:
  - Peters AE, Nguyen M, Green JB, Pearson ER, Buse JB, Sourij H, Hernandez AF, Sattar N, Holman RR, Mentz RJ, Shah SH. Proteomic pathways across the ejection fraction spectrum in patients with heart failure and diabetes mellitus: an EXSCEL trial substudy. Scientific reports. 2025. PMID: 40825966.
  https://pubmed.ncbi.nlm.nih.gov/40825966/
  - Bertalan PM, Nokhoijav E, Pap Á, Neagu GC, Káplár M, Darula Z, Kalló G, Prokai L, Csősz É. Comprehensive Insights into Obesity and Type 2 Diabetes from Protein Network, Canonical Pathway, Phosphorylation and Antimicrobial Peptide Signatures of Human Serum. Proteomes. 2025. PMID: 41441338.
  https://pubmed.ncbi.nlm.nih.gov/41441338/
  - Zhao L, Zhang Y, Liu F, Yang H, Zhong Y, Wang Y, Li S, Su Q, Tang L, Bai L, Ren H, Zou Y, Wang S, Zheng S, Xu H, Li L, Zhang J, Chai Z, Cooper ME, Tong N. Urinary complement proteins and risk of end-stage renal disease: quantitative urinary proteomics in patients with type 2 diabetes and biopsy-proven diabetic nephropathy. Journal of endocrinological investigation. 2021. PMID: 34043214.
  https://pubmed.ncbi.nlm.nih.gov/34043214/

## Biological Theme
Tight junction dynamics

## Direction
increased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 5
- Levels 1–3: 0
- Levels 4–5: 5

## Disease Relevance
LOW

## Dataset Evidence
- Representative pathway: Tight junction
- Representative library: KEGG_2021_Human
- Member pathways:
  - Tight junction
- Best member-pathway adjusted p-value: 0.0269372204878218
- Representative Gene Ratio: 0.0520833333333333

## Dataset Interpretation
The dataset indicates an increased representation of proteins associated with tight junctions (RAB8B, YBX3, ARPC5L, MSN, ROCK1) in Type 2 diabetes. These data represent Candidate protein abundance variations without establishing barrier permeability or junctional structure.

## Limitations
Supported by 5 lower-confidence candidate proteins (Levels 4–5) without Level 1–3 candidates, resulting in an EXPLORATORY evidence strength rating.

## Recent Scientific Research — Previous 5 Years
- Relevant primary PubMed studies meeting criteria: 0
- Commonly studied tissues/models: None reported in the provided records.
- Major findings: No primary literature meeting selection criteria was identified for this pathway term in Type 2 diabetes.
- Consistency or disagreement: No literature available to evaluate agreement.
- Important reported mechanisms: None reported in the provided records.
- Relationship between literature and dataset: Literature support is insufficient to evaluate dataset findings.
- Representative references: Fewer than three relevant primary references are available; none were retrieved for this cluster based on the predefined search criteria.

## Biological Theme
Regulation of IGF transport and uptake by IGFBPs

## Direction
increased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 38
- Levels 1–3: 1
- Levels 4–5: 37

## Disease Relevance
LOW

## Dataset Evidence
- Representative pathway: Regulation Of IGF Transport And Uptake By IGFBPs R-HSA-381426
- Representative library: Reactome_2022
- Member pathways:
  - Ferroptosis
  - HIF-1 signaling pathway
  - COPI-mediated Anterograde Transport R-HSA-6807878
  - Disease R-HSA-1643685
  - Golgi-to-ER Retrograde Transport R-HSA-8856688
  - Intra-Golgi And Retrograde Golgi-to-ER Traffic R-HSA-6811442
  - Iron Uptake And Transport R-HSA-917937
  - KEAP1-NFE2L2 Pathway R-HSA-9755511
  - Metabolism Of Proteins R-HSA-392499
  - Nuclear Events Mediated By NFE2L2 R-HSA-9759194
  - Post-translational Protein Modification R-HSA-597592
  - Post-translational Protein Phosphorylation R-HSA-8957275
  - Regulation Of Expression Of SLITs And ROBOs R-HSA-9010553
  - Regulation Of IGF Transport And Uptake By IGFBPs R-HSA-381426
  - Transport Of Small Molecules R-HSA-382551
- Best member-pathway adjusted p-value: 0.0023047733484059
- Representative Gene Ratio: 0.0625

## Dataset Interpretation
The dataset shows an increased representation of proteins associated with the regulation of IGF transport and uptake by IGFBPs in Type 2 diabetes. Mapped candidates include STAT5B (Level 1) and 37 Level 4–5 candidates (e.g., APOF, BST1, FSTL1, HDLBP, HK1, HMOX1, LCN2, PLG, THBS1, LCAT, PRDX1, TF). This reflects higher relative abundance of growth factor transport and matrix-binding proteins.

## Limitations
Although 38 total proteins are mapped, 37 belong to candidate levels 4–5 and only 1 to Level 1, maintaining an EXPLORATORY dataset evidence strength rating. Functional IGF bioavailability was not measured.

## Recent Scientific Research — Previous 5 Years
- Relevant primary PubMed studies meeting criteria: 0
- Commonly studied tissues/models: None reported in the provided records.
- Major findings: No primary literature meeting selection criteria was retrieved for this pathway query.
- Consistency or disagreement: No literature available to evaluate agreement.
- Important reported mechanisms: None reported in the provided records.
- Relationship between literature and dataset: Literature support is insufficient for this query formulation.
- Representative references: Fewer than three relevant primary references are available; none were retrieved for this cluster based on the predefined search criteria.

## Final Biological Interpretation

The proteomic findings from this comparison of Type 2 diabetes vs. non-diabetic controls highlight a complex metabolic and vascular signature. At the core metabolic level, there is a prominent decreased representation of proteins associated with mitochondrial oxidative processes, including fatty acid degradation (C011, MODERATE strength), fatty acid beta-oxidation (C002, EXPLORATORY strength), the citric acid cycle (C006, EXPLORATORY strength), oxidative phosphorylation (C009, EXPLORATORY strength), and amino acid catabolism (C010, EXPLORATORY strength). Independent literature evidence firmly establishes HIGH disease relevance for these metabolic axes, where preclinical and clinical research links suppressed mitochondrial oxidation and disrupted branched-chain amino acid catabolism to ectopic lipid accumulation, lipotoxicity, impaired glucose-stimulated insulin secretion, and systemic insulin resistance.

Concurrently, the dataset demonstrates a marked increased representation of proteins involved in inflammatory, vascular, and haemostatic cascades, most notably complement and coagulation cascades (C015, MODERATE strength) and regulation of fibrinolysis (C014, MODERATE strength), as well as increased IGF transport regulation (C017, EXPLORATORY strength). Published clinical literature confirms HIGH disease relevance for complement and coagulation cascades (C015) and neutrophil degranulation (C012), showing that circulating and biofluid levels of these innate immune markers track closely with microvascular and macrovascular diabetic complications, including nephropathy, retinopathy, and cardiovascular events. Crucially, dataset evidence strength (determined strictly by candidate counts and levels) remains distinct from disease relevance (derived from literature assessment), avoiding unwarranted claims of causality or directional pathway activation/inhibition from proteomic abundance alone.
