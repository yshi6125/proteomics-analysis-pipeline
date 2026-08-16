# Pathway Interpretation Report: Type 2 diabetes

## Executive Summary

This report provides a systematic biological assessment of differential pathway enrichment from a proteomics investigation comparing individuals with Type 2 Diabetes (DM) to non-diabetic controls (NDM). Across 16 evaluated biological clusters, the dataset demonstrates widespread shifts in metabolic, structural, and immune pathways. The most prominent statistical signature is characterized by increased representation of proteins involved in complement and coagulation cascades (Cluster C013), which achieves strong dataset evidence strength supported by robust representation of candidate proteins. In addition, multiple energy metabolism pathways exhibit consistent decreased representation, including cellular respiration (C001), fatty acid beta-oxidation (C002), the citric acid cycle (C007), fatty acid degradation (C008), and peroxisomal metabolism (C009). These metabolic clusters show high-to-moderate disease relevance based on extensive clinical and experimental literature in Type 2 diabetes, although their dataset evidence strength remains exploratory due to reliance on candidate-level protein measurements. Other exploratory clusters with increased representation include regulation of fibrinolysis (C012), tight junction organization (C014), Golgi-to-ER retrograde transport (C015), and regulation of IGF transport and uptake by IGFBPs (C016). Independent evaluation of dataset evidence strength versus published disease relevance ensures that strong literature support does not artificially inflate dataset confidence, nor does high statistical dataset strength imply clinical causation.

## Biological Theme
Cellular respiration

## Direction
decreased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 25
- Levels 1–3: 1
- Levels 4–5: 24

## Disease Relevance
HIGH

## Dataset Evidence
- Representative pathway: Cellular Respiration (GO:0045333)
- Representative library: GO_Biological_Process_2023
- Member pathways:
  - GO_Biological_Process_2023: Aerobic Electron Transport Chain (GO:0019646)
  - GO_Biological_Process_2023: Aerobic Respiration (GO:0009060)
  - GO_Biological_Process_2023: Cellular Respiration (GO:0045333)
  - GO_Biological_Process_2023: Energy Derivation By Oxidation Of Organic Compounds (GO:0015980)
  - GO_Biological_Process_2023: Mitochondrial ATP Synthesis Coupled Electron Transport (GO:0042775)
  - GO_Biological_Process_2023: Mitochondrial Electron Transport, Ubiquinol To Cytochrome C (GO:0006122)
  - GO_Biological_Process_2023: Oxidative Phosphorylation (GO:0006119)
  - KEGG_2021_Human: Alzheimer disease
  - KEGG_2021_Human: Amyotrophic lateral sclerosis
  - KEGG_2021_Human: Cardiac muscle contraction
  - KEGG_2021_Human: Diabetic cardiomyopathy
  - KEGG_2021_Human: Huntington disease
  - KEGG_2021_Human: Non-alcoholic fatty liver disease
  - KEGG_2021_Human: Oxidative phosphorylation
  - KEGG_2021_Human: Parkinson disease
  - KEGG_2021_Human: Pathways of neurodegeneration
  - KEGG_2021_Human: Prion disease
  - KEGG_2021_Human: Thermogenesis
  - Reactome_2022: Citric Acid (TCA) Cycle And Respiratory Electron Transport R-HSA-1428517
  - Reactome_2022: Respiratory Electron Transport R-HSA-611105
  - Reactome_2022: Respiratory Electron Transport, ATP Synthesis By Chemiosmotic Coupling, Heat Production By Uncoupling Proteins R-HSA-163200
- Best member-pathway adjusted p-value: 1.3971876455356703e-09
- Representative Gene Ratio: 0.0816326530612244

## Dataset Interpretation
The present dataset shows a decreased representation of proteins involved in cellular respiration in the DM versus NDM comparison. Twenty-five proteins contribute to this cluster, including electron transport chain components and mitochondrial respiratory complexes (e.g., NDUFA10, NDUFA12, SDHA, UQCRC1, UQCRC2, COX6B1, COX7A2). These proteomic variations indicate lower detected levels of key enzymes involved in oxidative phosphorylation and respiratory electron transport in the diabetic group. This dataset observation is descriptive and does not directly prove functional suppression or causality in mitochondrial flux.

## Limitations
Dataset evidence strength is classified as EXPLORATORY because 24 out of 25 supporting proteins belong to candidate levels 4–5, with only a single protein meeting level 1–3 criteria. The observational design prevents inferring whether decreased protein representation is a driver or consequence of metabolic dysfunction in Type 2 diabetes.

## Recent Scientific Research — Previous 5 Years
3 relevant primary PubMed studies meeting predefined criteria were identified.
- Commonly studied tissues/models: 3D engineered beige adipose tissue microvascular fragments, rodent skeletal (C2C12) and cardiac (H9C2) myotubes, and human Huh7 hepatocytes/liver cohorts.
- Major findings: Primary research demonstrates that cellular respiration is modulated in metabolic and T2D conditions. Skeletal myotubes under diabetic stress display reduced ATP-linked respiration, whereas cardiac myotubes show compensatory respiratory increases. Modulation of metabolic pathways (such as HIBCH/3-HIB) directly alters cellular respiration and lipid storage in hepatocytes.
- Consistency/disagreement: High agreement across studies that cellular respiration is altered under metabolic stress, displaying tissue-specific respiratory capacity and adaptive responses.
- Important reported mechanisms: Valine/3-HIB pathway feedback regulation modulating respiration, and suppression of ATP-linked respiration under diabetic conditions.
- Relationship to current dataset: Literature confirms that cellular respiration components are heavily impacted in diabetic tissues, reinforcing the clinical relevance of the decreased representation observed in this dataset.

Representative references:
Acosta et al. Engineering Functional Vascularized Beige Adipose Tissue from Microvascular Fragments of Models of Healthy and Type II Diabetes Conditions. Journal of tissue engineering. 2022. PMID: 35782994.
https://pubmed.ncbi.nlm.nih.gov/35782994/

Kopp et al. Modeling and Phenotyping Acute and Chronic Type 2 Diabetes Mellitus In Vitro in Rodent Heart and Skeletal Muscle Cells. Cells. 2023. PMID: 38132105.
https://pubmed.ncbi.nlm.nih.gov/38132105/

Bjune et al. Metabolic role of the hepatic valine/3-hydroxyisobutyrate (3-HIB) pathway in fatty liver disease. EBioMedicine. 2023. PMID: 37084480.
https://pubmed.ncbi.nlm.nih.gov/37084480/

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
  - GO_Biological_Process_2023: Electron Transport Chain (GO:0022900)
  - GO_Biological_Process_2023: Fatty Acid Beta-Oxidation Using acyl-CoA Dehydrogenase (GO:0033539)
  - GO_Biological_Process_2023: Respiratory Electron Transport Chain (GO:0022904)
- Best member-pathway adjusted p-value: 0.0039304901500761
- Representative Gene Ratio: 0.0204081632653061

## Dataset Interpretation
The dataset demonstrates a decreased representation of proteins involved in fatty acid beta-oxidation via acyl-CoA dehydrogenase in DM compared to NDM. Four candidate proteins (ACADS, ETFB, ETFDH, SDHA) are mapped to this cluster. These findings reflect decreased abundance of specific enzymes facilitating mitochondrial fatty acid catabolism in the diabetic state. No causal claims regarding metabolic rate or flux can be made from these protein abundance differences alone.

## Limitations
Evidence strength is EXPLORATORY because all 4 supporting proteins are categorized under Level 4 candidate criteria, lacking Level 1–3 confirmed protein support. The limited candidate number and exploratory nature require caution when generalizing findings.

## Recent Scientific Research — Previous 5 Years
8 relevant primary PubMed studies meeting predefined criteria were identified.
- Commonly studied tissues/models: db/db and HFD-induced diabetic mice, primary cardiomyocytes, human serum, human postmortem retinas, and Gmfb knockout rodents.
- Major findings: Studies demonstrate that promoting mitochondrial fatty acid beta-oxidation in metabolic organs facilitates lipid clearance and improves insulin sensitivity. Clinical treatment with semaglutide in T2DM-MAFLD patients induces metabolic shifts supporting enhanced beta-oxidation. Impaired beta-oxidation in retinas correlates with diabetic retinopathy progression, whereas excessive endothelial beta-oxidation in diabetic hearts can promote microvascular dysfunction.
- Consistency/disagreement: Strong consensus across studies that fatty acid beta-oxidation is central to T2D pathophysiology, though tissue-specific nuances exist between metabolic tissues and endothelium.
- Important reported mechanisms: Inhibition of Irak2 translocation by Adipsin preserving mitochondrial integrity, GMFB depletion driving PPARγ/CARM1-mediated beta-oxidation, and YY1/CYP2E1/PPARα signaling activation.
- Relationship to current dataset: Literature supports the high disease relevance of mitochondrial beta-oxidation in T2D, aligning with the decreased protein representation observed in the current proteomics analysis.

Representative references:
Niu et al. Multi-omics analysis of hepatic outcomes in T2DM-MAFLD patients treated with semaglutide: a single-centre, longitudinal, data-driven study. Frontiers in endocrinology. 2025. PMID: 41103643.
https://pubmed.ncbi.nlm.nih.gov/41103643/

Jiang et al. Adipsin inhibits Irak2 mitochondrial translocation and improves fatty acid β-oxidation to alleviate diabetic cardiomyopathy. Military Medical Research. 2023. PMID: 38072993.
https://pubmed.ncbi.nlm.nih.gov/38072993/

Veitch et al. MiR-30 promotes fatty acid beta-oxidation and endothelial cell dysfunction and is a circulating biomarker of coronary microvascular dysfunction in pre-clinical models of diabetes. Cardiovascular diabetology. 2022. PMID: 35209901.
https://pubmed.ncbi.nlm.nih.gov/35209901/

Fort et al. Diminished retinal complex lipid synthesis and impaired fatty acid β-oxidation associated with human diabetic retinopathy. JCI insight. 2021. PMID: 34437304.
https://pubmed.ncbi.nlm.nih.gov/34437304/

Zou et al. Tangerine Peel-Derived Exosome-Like Nanovesicles Alleviate Hepatic Steatosis Induced by Type 2 Diabetes: Evidenced by Regulating Lipid Metabolism and Intestinal Microflora. International journal of nanomedicine. 2024. PMID: 39371479.
https://pubmed.ncbi.nlm.nih.gov/39371479/

## Biological Theme
Cytoplasmic translation initiation

## Direction
decreased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 4
- Levels 1–3: 0
- Levels 4–5: 4

## Disease Relevance
LOW

## Dataset Evidence
- Representative pathway: Formation Of Cytoplasmic Translation Initiation Complex (GO:0001732)
- Representative library: GO_Biological_Process_2023
- Member pathways:
  - GO_Biological_Process_2023: Cytoplasmic Translational Initiation (GO:0002183)
  - GO_Biological_Process_2023: Formation Of Cytoplasmic Translation Initiation Complex (GO:0001732)
  - GO_Biological_Process_2023: Positive Regulation Of RNA Binding (GO:1905216)
  - GO_Biological_Process_2023: Positive Regulation Of mRNA Binding (GO:1902416)
  - GO_Biological_Process_2023: Regulation Of mRNA Binding (GO:1902415)
- Best member-pathway adjusted p-value: 0.0003511970191567
- Representative Gene Ratio: 0.0272108843537414

## Dataset Interpretation
The dataset indicates a decreased representation of proteins involved in the formation of the cytoplasmic translation initiation complex in DM relative to NDM. Four proteins (EIF3D, EIF3E, EIF3I, EIF6) map to this initiation assembly process. The findings reflect altered levels of eukaryotic initiation factors in the analyzed dataset. Functional impairment of global protein synthesis cannot be inferred from these abundance data.

## Limitations
Dataset evidence strength is EXPLORATORY, supported by 4 proteins all categorized as Level 4 candidates. Disease relevance is LOW due to the absence of direct literature linking this specific pathway to Type 2 diabetes in the assessed query.

## Recent Scientific Research — Previous 5 Years
0 relevant primary PubMed studies meeting predefined criteria were identified.
- Commonly studied tissues/models: None reported in the retrieved literature.
- Major findings: No primary research studies linking cytoplasmic translation initiation directly to Type 2 diabetes were retrieved for the evaluated period.
- Consistency/disagreement: Not applicable due to lack of primary studies.
- Important reported mechanisms: None reported.
- Relationship to current dataset: The current dataset finding remains unsupported by recent primary literature in the context of T2D.

No primary references are available.

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
  - GO_Biological_Process_2023: Mitochondrion Organization (GO:0007005)
- Best member-pathway adjusted p-value: 0.0118742251437276
- Representative Gene Ratio: 0.0476190476190476

## Dataset Interpretation
A decreased representation of proteins associated with mitochondrion organization is present in DM compared to NDM. Seven proteins (HSPD1, BCAP31, COX7A2, PHB2, SLC25A4, VPS13C, PRDX3) contribute to this cluster. These differences indicate reduced detected levels of structural and regulatory mitochondrial proteins in diabetic samples. These changes do not establish causal structural degradation or organelle dysfunction.

## Limitations
Dataset evidence strength is EXPLORATORY, with 1 Level 3 protein and 6 Level 4–5 proteins. The sample size and candidate distribution limit firm dataset-level conclusions.

## Recent Scientific Research — Previous 5 Years
1 relevant primary PubMed study meeting predefined criteria was identified.
- Commonly studied tissues/models: HSP22 knockout T2DM mice and vascular endothelial cells.
- Major findings: Primary research indicates that targeting mitochondrial organization and structural stability reduces diabetic vascular complications. Rosiglitazone preserves endothelial function by upregulating HSP22 and PPAR-γ signaling, mitigating mitochondrial oxidative stress and structural disruption.
- Consistency/disagreement: Consistent evidence from a single study demonstrating the importance of mitochondrial structural preservation in diabetic angiopathy.
- Important reported mechanisms: HSP22 upregulation and PPAR-γ activation suppressing mitochondrial reactive oxygen species formation.
- Relationship to current dataset: The literature supports moderate disease relevance by highlighting mitochondrial structural maintenance in diabetic vascular pathology, complementing the lower protein representation observed in this dataset.

Fewer than three relevant primary references are available:
Yu et al. Rosiglitazone reduces diabetes angiopathy by inhibiting mitochondrial dysfunction dependent on regulating HSP22 expression. iScience. 2023. PMID: 36968091.
https://pubmed.ncbi.nlm.nih.gov/36968091/

## Biological Theme
Regulation of ubiquitin-protein transferase activity

## Direction
decreased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 3
- Levels 1–3: 0
- Levels 4–5: 3

## Disease Relevance
LOW

## Dataset Evidence
- Representative pathway: Positive Regulation Of Ubiquitin-Protein Transferase Activity (GO:0051443)
- Representative library: GO_Biological_Process_2023
- Member pathways:
  - GO_Biological_Process_2023: Positive Regulation Of Ubiquitin-Protein Transferase Activity (GO:0051443)
- Best member-pathway adjusted p-value: 0.0388122493171869
- Representative Gene Ratio: 0.0204081632653061

## Dataset Interpretation
The dataset reveals a decreased representation of proteins associated with the positive regulation of ubiquitin-protein transferase activity in DM versus NDM. Three candidate proteins (PIN1, RPS2, SKP1) map to this pathway. This pattern reflects lower abundance of specific ubiquitination regulatory components in DM. Causality regarding proteasomal degradation rates cannot be inferred.

## Limitations
Dataset evidence strength is EXPLORATORY, supported by 3 candidate proteins (Levels 4–5). Disease relevance is classified as LOW due to the lack of direct primary literature retrieved for this theme in T2D.

## Recent Scientific Research — Previous 5 Years
0 relevant primary PubMed studies meeting predefined criteria were identified.
- Commonly studied tissues/models: None reported in the retrieved literature.
- Major findings: No primary studies investigating regulation of ubiquitin-protein transferase activity in Type 2 diabetes were retrieved.
- Consistency/disagreement: Not applicable.
- Important reported mechanisms: None reported.
- Relationship to current dataset: Unverified by published primary literature in the specific T2D context.

No primary references are available.

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
  - GO_Biological_Process_2023: GMP Biosynthetic Process (GO:0006177)
  - GO_Biological_Process_2023: Nucleobase Biosynthetic Process (GO:0046112)
  - GO_Biological_Process_2023: Purine Nucleobase Biosynthetic Process (GO:0009113)
- Best member-pathway adjusted p-value: 0.0248889763070487
- Representative Gene Ratio: 0.0136054421768707

## Dataset Interpretation
Decreased representation of proteins involved in purine nucleobase biosynthesis is observed in DM compared to NDM. Two candidate proteins (GMPS, PAICS) are associated with this cluster. These data show decreased detected levels of enzymes in de novo purine synthesis in diabetic subjects. Alterations in nucleotide pool sizes or metabolic fluxes cannot be directly determined from protein abundances alone.

## Limitations
Dataset evidence strength is EXPLORATORY (2 Level 4 candidate proteins). Disease relevance is LOW based on absence of direct primary literature for this specific term.

## Recent Scientific Research — Previous 5 Years
0 relevant primary PubMed studies meeting predefined criteria were identified.
- Commonly studied tissues/models: None reported in the retrieved literature.
- Major findings: No primary literature was retrieved for purine nucleobase biosynthesis in T2D.
- Consistency/disagreement: Not applicable.
- Important reported mechanisms: None reported.
- Relationship to current dataset: Unsubstantiated by primary literature in this specific disease comparison.

No primary references are available.

## Biological Theme
Citric acid cycle

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
- Representative pathway: Citrate cycle (TCA cycle)
- Representative library: KEGG_2021_Human
- Member pathways:
  - GO_Biological_Process_2023: Dicarboxylic Acid Metabolic Process (GO:0043648)
  - KEGG_2021_Human: Citrate cycle (TCA cycle)
  - KEGG_2021_Human: Glyoxylate and dicarboxylate metabolism
  - KEGG_2021_Human: Pyruvate metabolism
  - Reactome_2022: Citric Acid Cycle (TCA Cycle) R-HSA-71403
  - Reactome_2022: Pyruvate Metabolism And Citric Acid (TCA) Cycle R-HSA-71406
- Best member-pathway adjusted p-value: 4.401184107912802e-05
- Representative Gene Ratio: 0.0408163265306122

## Dataset Interpretation
The current dataset exhibits a decreased representation of proteins associated with the citric acid cycle in DM versus NDM. Eleven candidate proteins (ALDH5A1, CS, DLAT, DLD, MRPS36, SDHA, ACO2, ALDH2, GRHPR, LDHB, MDH2) are included. These data indicate lower detected levels of key enzymes involved in tricarboxylic acid oxidation and pyruvate handling in DM. Functional suppression of mitochondrial ATP generation cannot be asserted without enzymatic assays.

## Limitations
Dataset evidence strength is EXPLORATORY because all 11 supporting proteins reside in candidate levels 4–5. Observational data preclude causal interpretation.

## Recent Scientific Research — Previous 5 Years
4 relevant primary PubMed studies meeting predefined criteria were identified.
- Commonly studied tissues/models: High-glucose/high-fat diet diabetic rats, Leptin receptor-deficient T2D rats, human T2D stem cells, and human gastric mucosal biopsies post-bariatric surgery.
- Major findings: Primary research confirms that TCA cycle intermediates, flux, and enzyme expression are altered across clinical samples and diabetic models. Exogenous supplementation of the TCA cycle intermediate alpha-ketoglutarate improves matrix formation in T2D stem cells, while traditional formulations like Dangua Fang restore TCA flux via Nampt upregulation.
- Consistency/disagreement: High agreement across primary studies that TCA cycle activity is compromised in T2D contexts and responsive to metabolic intervention.
- Important reported mechanisms: Nampt upregulation enhancing TCA cycle metabolic flux, and anaplerotic supplementation of alpha-ketoglutarate.
- Relationship to current dataset: Published literature confirms high disease relevance, validating the centrality of TCA cycle alterations in T2D, while dataset evidence strength remains independent and exploratory.

Representative references:
Heng et al. Dangua Fang regulating tricarboxylic acid cycle and respiratory chain and its mechanism in diabetic rats. Journal of traditional Chinese medicine = Chung i tsa chih ying wen pan. 2023. PMID: 37946477.
https://pubmed.ncbi.nlm.nih.gov/37946477/

Dias et al. Endogenous dysregulated energy and amino acid metabolism delay scaffold-guided large volume bone regeneration in a diabetic rat model with Leptin receptor deficiency. Acta biomaterialia. 2025. PMID: 40319991.
https://pubmed.ncbi.nlm.nih.gov/40319991/

Larson et al. Novel insights on remnant stomach following Roux-en-Y gastric bypass surgery based on histological evaluation and quantitative proteomics analysis. Scientific reports. 2025. PMID: 40652086.
https://pubmed.ncbi.nlm.nih.gov/40652086/

Yang et al. Identification of Shared Pathways and Molecules Between Type 2 Diabetes and Lung Adenocarcinoma and the Impact of High Glucose Environment on Lung Adenocarcinoma. International journal of endocrinology. 2025. PMID: 40212965.
https://pubmed.ncbi.nlm.nih.gov/40212965/

## Biological Theme
Fatty acid degradation

## Direction
decreased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 21
- Levels 1–3: 0
- Levels 4–5: 21

## Disease Relevance
HIGH

## Dataset Evidence
- Representative pathway: Fatty acid degradation
- Representative library: KEGG_2021_Human
- Member pathways:
  - GO_Biological_Process_2023: Amino Acid Catabolic Process (GO:0009063)
  - GO_Biological_Process_2023: Branched-Chain Amino Acid Catabolic Process (GO:0009083)
  - GO_Biological_Process_2023: Fatty Acid Beta-Oxidation (GO:0006635)
  - GO_Biological_Process_2023: Fatty Acid Catabolic Process (GO:0009062)
  - GO_Biological_Process_2023: Fatty Acid Oxidation (GO:0019395)
  - KEGG_2021_Human: Butanoate metabolism
  - KEGG_2021_Human: Fatty acid degradation
  - KEGG_2021_Human: Fatty acid elongation
  - KEGG_2021_Human: Lysine degradation
  - KEGG_2021_Human: Propanoate metabolism
  - KEGG_2021_Human: Tryptophan metabolism
  - KEGG_2021_Human: Valine, leucine and isoleucine degradation
  - KEGG_2021_Human: beta-Alanine metabolism
  - Reactome_2022: Beta Oxidation Of butanoyl-CoA To acetyl-CoA R-HSA-77352
  - Reactome_2022: Beta Oxidation Of decanoyl-CoA To octanoyl-CoA-CoA R-HSA-77346
  - Reactome_2022: Beta Oxidation Of hexanoyl-CoA To butanoyl-CoA R-HSA-77350
  - Reactome_2022: Beta Oxidation Of lauroyl-CoA To decanoyl-CoA-CoA R-HSA-77310
  - Reactome_2022: Beta Oxidation Of octanoyl-CoA To hexanoyl-CoA R-HSA-77348
  - Reactome_2022: Fatty Acid Metabolism R-HSA-8978868
  - Reactome_2022: Mitochondrial Fatty Acid Beta-Oxidation Of Saturated Fatty Acids R-HSA-77286
  - Reactome_2022: Mitochondrial Fatty Acid Beta-Oxidation R-HSA-77289
- Best member-pathway adjusted p-value: 1.9604197393568753e-07
- Representative Gene Ratio: 0.0612244897959183

## Dataset Interpretation
A decreased representation of proteins associated with fatty acid degradation is observed in DM relative to NDM. Twenty-one proteins map to this broad metabolic cluster (e.g., ACADS, ACSL1, ACSL3, ECHS1, ECI1, ETFB, ETFDH, FASN, HADH, OXCT1, ACAA2, HADHA). These changes indicate lower measured levels of enzymes involved in fatty acid breakdown and branched-chain amino acid catabolism in DM. The findings describe differential protein abundance without demonstrating dynamic metabolic inhibition.

## Limitations
Dataset evidence strength is EXPLORATORY because all 21 supporting proteins belong to Level 4–5 candidate tiers. Causal direction cannot be inferred.

## Recent Scientific Research — Previous 5 Years
10 relevant primary PubMed studies meeting predefined criteria were identified.
- Commonly studied tissues/models: Human patient cohorts (blood/plasma/serum), db/db and HFD/STZ-induced diabetic mice, and HepG2 cells.
- Major findings: Primary clinical and experimental studies consistently show that fatty acid degradation pathways are dysregulated in T2D and its vascular/cognitive complications. Pharmacological interventions (e.g., dapagliflozin, epimedin C, Rosa roxburghii Tratt ellagitannin, and JTTZR) enhance fatty acid degradation to reduce lipid accumulation and lipotoxicity.
- Consistency/disagreement: High agreement across studies that fatty acid degradation is compromised in diabetic states and responsive to lipid-lowering therapies.
- Important reported mechanisms: AMPK-SIRT1 activation upregulating CPT1alpha, PPAR alpha/gamma activation, and gut microbiota-mediated metabolic shifts.
- Relationship to current dataset: Literature confirms the strong disease relevance of fatty acid degradation in T2D, providing context for the decreased protein representation observed in this dataset.

Representative references:
Lee et al. Sodium-Glucose Cotransporter-2 Inhibitor Enhances Hepatic Gluconeogenesis and Reduces Lipid Accumulation via AMPK-SIRT1 Activation and Autophagy Induction. Endocrinology and metabolism (Seoul, Korea). 2025. PMID: 40350800.
https://pubmed.ncbi.nlm.nih.gov/40350800/

Zhou et al. The Mechanism Underlying the Hypoglycemic Effect of Epimedin C on Mice with Type 2 Diabetes Mellitus Based on Proteomic Analysis. Nutrients. 2023. PMID: 38201855.
https://pubmed.ncbi.nlm.nih.gov/38201855/

Tan et al. Transcriptomics Reveals the Mechanism of Rosa roxburghii Tratt Ellagitannin in Improving Hepatic Lipid Metabolism Disorder in db/db Mice. Nutrients. 2023. PMID: 37836471.
https://pubmed.ncbi.nlm.nih.gov/37836471/

Zhang et al. Decoding cardiovascular risks: analyzing type 2 diabetes mellitus and ASCVD gene expression. Frontiers in endocrinology. 2024. PMID: 38715799.
https://pubmed.ncbi.nlm.nih.gov/38715799/

Bao et al. Exploring the Regulation of Jiangtang Tiaozhi Formula on the Biological Network of Obese T2DM Complicated With Dyslipidemia Based on Clinical Transcriptomics. Frontiers in endocrinology. 2022. PMID: 35957821.
https://pubmed.ncbi.nlm.nih.gov/35957821/

## Biological Theme
Peroxisomal metabolism

## Direction
decreased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 11
- Levels 1–3: 0
- Levels 4–5: 11

## Disease Relevance
MODERATE

## Dataset Evidence
- Representative pathway: Peroxisome
- Representative library: KEGG_2021_Human
- Member pathways:
  - GO_Biological_Process_2023: Monocarboxylic Acid Metabolic Process (GO:0032787)
  - KEGG_2021_Human: Fatty acid biosynthesis
  - KEGG_2021_Human: Peroxisome
- Best member-pathway adjusted p-value: 0.0033797115361891
- Representative Gene Ratio: 0.0408163265306122

## Dataset Interpretation
The dataset indicates a decreased representation of proteins associated with peroxisomal metabolism in DM compared to NDM. Eleven candidate proteins (ACSF2, ACSL1, ACSL3, ECHDC3, FASN, PEX19, ECH1, GRHPR, LDHB, PRDX5, SOD1) map to this cluster. These results show lower relative abundance of peroxisomal biogenesis and lipid-processing proteins in DM. Peroxisomal functional loss cannot be asserted directly from these data.

## Limitations
Dataset evidence strength is EXPLORATORY (11 Level 4–5 candidate proteins). Observational data prevent drawing mechanistic inferences.

## Recent Scientific Research — Previous 5 Years
1 relevant primary PubMed study meeting predefined criteria was identified.
- Commonly studied tissues/models: Human epicardial adipose tissue from T2D patients.
- Major findings: Clinical research shows that GPR146 expression in epicardial adipose tissue from T2D patients correlates with peroxisomal metabolism genes, serum lipid levels, and response to SGLT2 inhibitor therapy.
- Consistency/disagreement: Consistent clinical observation linking peroxisomal transcriptomic signatures in adipose tissue to diabetic dyslipidemia.
- Important reported mechanisms: GPR146 co-expression with peroxisomal lipid metabolism genes.
- Relationship to current dataset: Literature supports moderate disease relevance, aligning with the observed lower representation of peroxisomal proteins in DM.

Fewer than three relevant primary references are available:
Ryk et al. Cholesin receptor signalling is active in cardiovascular system-associated adipose tissue and correlates with SGLT2i treatment in patients with diabetes. Cardiovascular diabetology. 2024. PMID: 38902687.
https://pubmed.ncbi.nlm.nih.gov/38902687/

## Biological Theme
RHOG GTPase cycle

## Direction
decreased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 5
- Levels 1–3: 0
- Levels 4–5: 5

## Disease Relevance
LOW

## Dataset Evidence
- Representative pathway: RHOG GTPase Cycle R-HSA-9013408
- Representative library: Reactome_2022
- Member pathways:
  - Reactome_2022: RHOG GTPase Cycle R-HSA-9013408
- Best member-pathway adjusted p-value: 0.0360539855655348
- Representative Gene Ratio: 0.0340136054421768

## Dataset Interpretation
A decreased representation of proteins associated with the RHOG GTPase cycle is observed in DM compared to NDM. Five proteins (HSPE1, ITGB1, LETM1, STBD1, RAB7A) contribute to this cluster. These data reflect lower relative abundance of Small GTPase signaling elements in diabetic samples. Changes in GTPase signaling activity cannot be inferred.

## Limitations
Evidence strength is EXPLORATORY, based on 5 Level 4–5 candidate proteins. Disease relevance is LOW due to lack of retrieved primary literature.

## Recent Scientific Research — Previous 5 Years
0 relevant primary PubMed studies meeting predefined criteria were identified.
- Commonly studied tissues/models: None reported in the retrieved literature.
- Major findings: No primary research was retrieved connecting the RHOG GTPase cycle to T2D.
- Consistency/disagreement: Not applicable.
- Important reported mechanisms: None reported.
- Relationship to current dataset: Unvalidated by primary literature in this specific context.

No primary references are available.

## Biological Theme
Regulation of pinocytosis

## Direction
increased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 2
- Levels 1–3: 0
- Levels 4–5: 2

## Disease Relevance
LOW

## Dataset Evidence
- Representative pathway: Positive Regulation Of Pinocytosis (GO:0048549)
- Representative library: GO_Biological_Process_2023
- Member pathways:
  - GO_Biological_Process_2023: Positive Regulation Of Pinocytosis (GO:0048549)
- Best member-pathway adjusted p-value: 0.0445514661793356
- Representative Gene Ratio: 0.019047619047619

## Dataset Interpretation
The dataset demonstrates an increased representation of proteins involved in the regulation of pinocytosis in DM relative to NDM. Two candidate proteins (APPL1, PPT1) map to this pathway. These findings indicate higher measured levels of fluid-phase endocytosis regulators in DM. Endocytic flux changes cannot be inferred.

## Limitations
Dataset evidence strength is EXPLORATORY (2 Level 4–5 proteins). Disease relevance is LOW due to lack of supporting literature.

## Recent Scientific Research — Previous 5 Years
0 relevant primary PubMed studies meeting predefined criteria were identified.
- Commonly studied tissues/models: None reported in the retrieved literature.
- Major findings: No primary studies were retrieved for pinocytosis regulation in T2D.
- Consistency/disagreement: Not applicable.
- Important reported mechanisms: None reported.
- Relationship to current dataset: Currently unsupported by literature.

No primary references are available.

## Biological Theme
Regulation of fibrinolysis

## Direction
increased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 9
- Levels 1–3: 1
- Levels 4–5: 8

## Disease Relevance
LOW

## Dataset Evidence
- Representative pathway: Regulation Of Fibrinolysis (GO:0051917)
- Representative library: GO_Biological_Process_2023
- Member pathways:
  - GO_Biological_Process_2023: Negative Regulation Of Blood Coagulation (GO:0030195)
  - GO_Biological_Process_2023: Negative Regulation Of Fibrinolysis (GO:0051918)
  - GO_Biological_Process_2023: Peptide Cross-Linking (GO:0018149)
  - GO_Biological_Process_2023: Positive Regulation Of Blood Coagulation (GO:0030194)
  - GO_Biological_Process_2023: Regulation Of Fibrinolysis (GO:0051917)
  - Reactome_2022: Platelet Degranulation R-HSA-114608
  - Reactome_2022: Response To Elevated Platelet Cytosolic Ca2+ R-HSA-76005
- Best member-pathway adjusted p-value: 0.0010742123216896
- Representative Gene Ratio: 0.0380952380952381

## Dataset Interpretation
An increased representation of proteins involved in the regulation of fibrinolysis is present in DM compared to NDM. Nine supporting proteins (VTN, F13A1, ITIH3, KLKB1, KRT1, PLG, PPBP, THBS1, TF) are enriched in this cluster. These data show higher detected levels of coagulation and fibrinolytic regulators in DM. Alterations in clot clearance rates or systemic thrombosis cannot be inferred directly from candidate abundance.

## Limitations
Dataset evidence strength is EXPLORATORY (1 Level 2 protein, 8 Level 4–5 proteins). Disease relevance is LOW as no primary literature was retrieved for this precise search query.

## Recent Scientific Research — Previous 5 Years
0 relevant primary PubMed studies meeting predefined criteria were identified.
- Commonly studied tissues/models: None reported in the retrieved literature.
- Major findings: No primary studies were retrieved for regulation of fibrinolysis in T2D.
- Consistency/disagreement: Not applicable.
- Important reported mechanisms: None reported.
- Relationship to current dataset: Lacks published primary literature validation under the current search parameters.

No primary references are available.

## Biological Theme
Complement and coagulation cascades

## Direction
increased

## Dataset Evidence Strength
STRONG

Supporting proteins:
- Total proteins: 37
- Levels 1–3: 3
- Levels 4–5: 34

## Disease Relevance
HIGH

## Dataset Evidence
- Representative pathway: Complement and coagulation cascades
- Representative library: KEGG_2021_Human
- Member pathways:
  - GO_Biological_Process_2023: Positive Regulation Of Immune System Process (GO:0002684)
  - KEGG_2021_Human: Complement and coagulation cascades
  - KEGG_2021_Human: Coronavirus disease
  - KEGG_2021_Human: Staphylococcus aureus infection
  - KEGG_2021_Human: Systemic lupus erythematosus
  - Reactome_2022: Activation Of C3 And C5 R-HSA-174577
  - Reactome_2022: Complement Cascade R-HSA-166658
  - Reactome_2022: Immune System R-HSA-168256
  - Reactome_2022: Innate Immune System R-HSA-168249
  - Reactome_2022: Regulation Of Complement Cascade R-HSA-977606
  - Reactome_2022: Terminal Pathway Of Complement R-HSA-166665
- Best member-pathway adjusted p-value: 1.239292998867119e-12
- Representative Gene Ratio: 0.1238095238095238

## Dataset Interpretation
The dataset demonstrates an increased representation of proteins involved in complement and coagulation cascades in DM versus NDM. Thirty-seven proteins support this cluster, including 3 Level 1–2 proteins (STAT5B, DCTN6, VTN) and 34 Level 4–5 proteins (e.g., C2, C5, C6, C8A, C9, CFHR2, F13A1, KLKB1, PLG, PROCR, C1S, CFB). This robust signal reflects increased systemic abundance of innate immune and pro-coagulant proteins in diabetic subjects. Pathway activation or causal role in vascular injury cannot be assumed without functional testing.

## Limitations
Although dataset evidence strength is STRONG due to protein counts and inclusion of Level 1–2 targets, observational design restricts conclusions to biomarker representation rather than causal etiology.

## Recent Scientific Research — Previous 5 Years
9 relevant primary PubMed studies meeting predefined criteria were identified.
- Commonly studied tissues/models: Human serum, plasma, urine, aqueous humor, circulating exosomes, platelets, and ob/ob / db/db mice.
- Major findings: Extensive primary literature demonstrates that complement and coagulation cascade proteins are elevated across serum, plasma, exosomes, urine, and aqueous humor in T2D and its vascular/renal complications. Urinary complement proteins (such as CFAH and DAF) independently predict end-stage renal disease in diabetic nephropathy. In platelets, autophagy-driven processes promote complement activation and microthrombosis.
- Consistency/disagreement: High consistency across independent human cohort and animal model studies.
- Important reported mechanisms: Autophagy-dependent platelet activation, urinary complement dysregulation driving renal injury, and hepatic secretion of complement factors into circulation.
- Relationship to current dataset: Literature confirms the high disease relevance of complement and coagulation activation in T2D, providing strong external validation for the statistical signature observed in this dataset.

Representative references:
Bertalan et al. Comprehensive Insights into Obesity and Type 2 Diabetes from Protein Network, Canonical Pathway, Phosphorylation and Antimicrobial Peptide Signatures of Human Serum. Proteomes. 2025. PMID: 41441338.
https://pubmed.ncbi.nlm.nih.gov/41441338/

Wu et al. Autophagy-enabled protein degradation: Key to platelet activation and ANGII production in patients with type 2 diabetes mellitus. Heliyon. 2024. PMID: 39253219.
https://pubmed.ncbi.nlm.nih.gov/39253219/

Jin et al. Stage-dependent proteomic alterations in aqueous humor of diabetic retinopathy patients based on data-independent acquisition and parallel reaction monitoring. Journal of translational medicine. 2025. PMID: 40281624.
https://pubmed.ncbi.nlm.nih.gov/40281624/

Zhao et al. Urinary complement proteins and risk of end-stage renal disease: quantitative urinary proteomics in patients with type 2 diabetes and biopsy-proven diabetic nephropathy. Journal of endocrinological investigation. 2021. PMID: 34043214.
https://pubmed.ncbi.nlm.nih.gov/34043214/

Liu et al. Plasma exosome proteomics in different glucose statuses: a cross-sectional study on type 2 diabetes pathogenesis. Clinical proteomics. 2026. PMID: 41864872.
https://pubmed.ncbi.nlm.nih.gov/41864872/

## Biological Theme
Tight junction organization

## Direction
increased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 7
- Levels 1–3: 0
- Levels 4–5: 7

## Disease Relevance
LOW

## Dataset Evidence
- Representative pathway: Tight junction
- Representative library: KEGG_2021_Human
- Member pathways:
  - GO_Biological_Process_2023: Actin Nucleation (GO:0045010)
  - KEGG_2021_Human: Tight junction
- Best member-pathway adjusted p-value: 0.0264055584513487
- Representative Gene Ratio: 0.0571428571428571

## Dataset Interpretation
An increased representation of proteins associated with tight junction organization is observed in DM compared to NDM. Seven candidate proteins (LMOD1, RAB8B, YBX3, ARPC2, ARPC5L, MSN, ROCK1) contribute to this cluster. These data indicate higher levels of structural and cytoskeletal regulators linked to tight junctions in diabetic samples. Changes in barrier integrity or cellular permeability cannot be determined from these static levels.

## Limitations
Evidence strength is EXPLORATORY (7 Level 4–5 proteins). Disease relevance is LOW due to lack of retrieved literature.

## Recent Scientific Research — Previous 5 Years
0 relevant primary PubMed studies meeting predefined criteria were identified.
- Commonly studied tissues/models: None reported in the retrieved literature.
- Major findings: No primary studies were retrieved for tight junction organization under the search criteria.
- Consistency/disagreement: Not applicable.
- Important reported mechanisms: None reported.
- Relationship to current dataset: Unverified by published primary literature in this specific context.

No primary references are available.

## Biological Theme
Golgi-to-ER retrograde transport

## Direction
increased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 5
- Levels 1–3: 1
- Levels 4–5: 4

## Disease Relevance
LOW

## Dataset Evidence
- Representative pathway: Golgi-to-ER Retrograde Transport R-HSA-8856688
- Representative library: Reactome_2022
- Member pathways:
  - Reactome_2022: Golgi-to-ER Retrograde Transport R-HSA-8856688
- Best member-pathway adjusted p-value: 0.0392235166053615
- Representative Gene Ratio: 0.0476190476190476

## Dataset Interpretation
The dataset shows an increased representation of proteins involved in Golgi-to-ER retrograde transport in DM relative to NDM. Five proteins (DCTN6, CENPE, COPB2, COPZ1, DCTN3) support this cluster. These results indicate increased abundance of coatomer and dynactin transport components in DM. Intracellular trafficking rates cannot be inferred.

## Limitations
Dataset evidence strength is EXPLORATORY (1 Level 2 protein, 4 Level 4 proteins). Disease relevance is LOW due to absent literature.

## Recent Scientific Research — Previous 5 Years
0 relevant primary PubMed studies meeting predefined criteria were identified.
- Commonly studied tissues/models: None reported in the retrieved literature.
- Major findings: No primary literature was retrieved for Golgi-to-ER retrograde transport in T2D.
- Consistency/disagreement: Not applicable.
- Important reported mechanisms: None reported.
- Relationship to current dataset: Unsubstantiated by primary literature.

No primary references are available.

## Biological Theme
Regulation of IGF transport and uptake by IGFBPs

## Direction
increased

## Dataset Evidence Strength
EXPLORATORY

Supporting proteins:
- Total proteins: 26
- Levels 1–3: 1
- Levels 4–5: 25

## Disease Relevance
LOW

## Dataset Evidence
- Representative pathway: Regulation Of IGF Transport And Uptake By IGFBPs R-HSA-381426
- Representative library: Reactome_2022
- Member pathways:
  - Reactome_2022: Iron Uptake And Transport R-HSA-917937
  - Reactome_2022: Metabolism Of Proteins R-HSA-392499
  - Reactome_2022: Post-translational Protein Phosphorylation R-HSA-8957275
  - Reactome_2022: Regulation Of IGF Transport And Uptake By IGFBPs R-HSA-381426
- Best member-pathway adjusted p-value: 0.0031419709015453
- Representative Gene Ratio: 0.0666666666666666

## Dataset Interpretation
An increased representation of proteins associated with the regulation of IGF transport and uptake by IGFBPs is observed in DM compared to NDM. Twenty-six proteins support this cluster (e.g., DCTN6, BST1, COPB2, COPZ1, FSTL1, HMOX1, KLKB1, LCN2, PLG, THBS1, CP, TF). These data reflect higher detected levels of extracellular matrix, carrier, and regulatory proteins in diabetic samples. IGF-1 bioavailability or signaling activity cannot be established from protein abundance alone.

## Limitations
Dataset evidence strength is EXPLORATORY (1 Level 2 protein, 25 Level 4–5 candidate proteins). Disease relevance is LOW as no primary literature was retrieved for this search string.

## Recent Scientific Research — Previous 5 Years
0 relevant primary PubMed studies meeting predefined criteria were identified.
- Commonly studied tissues/models: None reported in the retrieved literature.
- Major findings: No primary studies were retrieved for this specific pathway term in T2D.
- Consistency/disagreement: Not applicable.
- Important reported mechanisms: None reported.
- Relationship to current dataset: Unvalidated by primary literature under the evaluated search terms.

No primary references are available.

## Final Biological Interpretation

This comprehensive analysis synthesizes proteomic pathway enrichment with published scientific literature in Type 2 Diabetes (DM vs. NDM comparison). A major finding of the study is the strong dataset evidence strength observed for complement and coagulation cascades (C013), characterized by a large cluster size (37 proteins), high gene ratio (0.1238), and low adjusted p-value (1.24e-12), including Level 1–2 confirmed proteins. This statistical signature is heavily reinforced by literature demonstrating systemic activation of complement factors across diabetic serum, plasma, exosome, and urinary cohorts, where specific complement proteins independently predict diabetic microvascular complications. In parallel, the dataset demonstrates a widespread decreased representation of key mitochondrial energy metabolic pathways, including cellular respiration (C001), fatty acid beta-oxidation (C002), the citric acid cycle (C007), fatty acid degradation (C008), and peroxisomal metabolism (C009). Although these metabolic clusters are categorized as having exploratory dataset evidence strength due to reliance on candidate-tier (Levels 4–5) proteins, their disease relevance is consistently high to moderate across extensive independent literature in diabetic models and human cohorts. Strict separation between independent dataset evidence strength and disease relevance ensures that strong statistical findings (such as C013) are evaluated based on internal data quality, while highly relevant metabolic themes (C001, C002, C007, C008) are recognized for their established disease context without overstating dataset-level certainty. Overall, these findings highlight a dual proteomic profile in Type 2 diabetes marked by heightened representation of pro-inflammatory/coagulation cascade proteins and lowered representation of central mitochondrial and fatty acid oxidative machinery.
