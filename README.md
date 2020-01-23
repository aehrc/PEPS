## For PEPS2 see below

# PEPS: Polygenic Epistatic Phenotype Simulation

**Based on real genotype data.**
Polygenic Risk Score (PRS) models are proven to find better risk indicators for polygenic diseases. However, these models are limited to the additive effect of individual SNPs. What if a polygenic trait is formed in a more complex way where a superset of individual SNPs along with sets of interactive SNPs (Higher-Order Complex Epistasis Interactions) contribute to the phenotype. PEPS is developed to generate such complex phenotypes. The phenotype is simulated for real genotype data (i.e. 1000 Genomes Project) where all genomic pattern exists in the data (only the phenotype is simulated).

## Why Simulation?

Whether such a complex phenotype exists in the real world is a research question. Before answering this question, we should find out the software that can utilize all individual and interactive SNPs all together. Once such software is identified, it can be used to process real-world datasets to form more accurate risk indicators for complex diseases.

## What are Challenges?

There are two mathematical challenges to form such complex phenotypes:
Producing a phenotype that depends on the additive effect of many variables. In the PEPS world, a variable could be either an SNP or a set of SNPs interact in a complex way (i.e. 2-way, 3-way ... n-way epistasis interaction).
Producing each of the n-way epistasis interactions is a complex problem. Existing epistasis simulation software (EpiSim and GAMETES) simulate both phenotype and genotype. Yet GAMETES could not guarantee to form the phenotype.

## How PEPS Works?

**For more details see "PEPS Workflow" below.**
PEPS bypasses the above challenges by using a feedback loop simulation. In the first step, PEPS forms the variable that the user asked (i.e. 20 individual SNPs, 15 2-way interactive SNPs and 35 3-way interactive SNPs that is 20+15+35 variables in total). Then a randomly generated phenotype is assigned to the samples such that half of the sample are case and the other half are controls.

Next PEPS identifies the probability of being a case for each genotype of each variable. For SNP variables there are 3 genotypes 0/0 ( R), 0/1 (H) and 1/1 (A). For epistasis variables, the genotype is the concatenation of individual SNPs in the interaction. HRAH, HHAR, and RRAR are examples of genotype for a 4-way interaction variable. Thus an n-way interaction can have 3^n different genotype.

Finally, PEPS computes the probability of being a case for each sample by adding the probability of being a case for the genotype of that sample in all variables. Samples with the probability of being a case higher than average are assigned to the case group. Other samples are assigned to the control group.

## What are the Drawbacks and How to Address them?

The feedback loop PEPS used for simulation does not guarantee to associate all variables with the phenotype. Thus the resulting phenotype would be a function of a set of variables PEPS starts with. To address this situation PEPS computes the chi2 p-value of each variable with the resulting phenotype. Only variables that exceed the p-value threshold (given by the user) are considered as truth variables. Subsequently, SNPs that form truth variables are called truth SNPs.

When using PEPS, the input genotype should not include SNPs that are associated with any population structure. For example, in the case of 1000 Genomes project we observe that if we include SNPs that can be a predictor of ethnicity, the resulting synthetic phenotype mimics the ethnicity of the population. To address this problem, you should exclude all SNPs associated with any population structure from the genotype file given to PEPS. Such SNP should not be a truth SNP for a synthetic phenotype.

## Config file

This file is in JSON format and includes all the parameter for the simulation

- inputType `str`: "csv" or "vcf"
- dumpCSV `boolean`: In case the input is vcf save it in csv format (for faster access in later simulation)
- shuffleSnps `boolean`: Shuffle SNPs before assigning them to variables.
- outputPrefix `str`: Output file prefix
- inputPrefix `str`: Input file prefix. Depending on inputType ".csv" or ".vcf" is added to the prefix.
- pvalueThr `float`: Threshold for chi2 p-value of truth variable.
- numTree `int`: Number of tree in RandomForest model to compute AUC
- numLoop `int`: Number of times the simulation feedback loop is repeated
- variables: `array` of
  - numSnpsInVar `int`: Type of variable. 1 for SNP variable, 2 for 2-way (2-SNP) interaction Variable and n for n-way (n-SNP) interaction Variable
  - numVar `int`: Number of variable of this type

## Sample Data

The [SampleData](SampleData) directory includes examples of input, output and config file as well as processed PEPS notebook in HTML format.

**Input:**

There are two input file both of them are a subset of 1000Genomes data available in vcf and csv format. The `small` input with ~162 SNPs and the `large` input with 4969 SNPs.

**config:**

There are 2 config files. [config-small.json](SampleData/config-small.json) and [config-large.json](SampleData/config-large.json) are made to process small vcf file and large csv file respectively.

**Output:**

Output files related to both config files above.

**PEPS processed notebook**

[PEPS-small.html](SampleData/PEPS-small.html) and [PEPS-large.html](SampleData/PEPS-large.html) are examples that show how PEPS notebook look likes after processing small and large input.

## VCF to CSV conversion

To parse VCF files, PEPS uses a library called `pdbio`. This library is slow and only works for a tiny vcf file. **We strongly recommend to prepare input data in CSV format.** To do so you can use the [VCF_2_CSV.sh](VCF_2_CSV.sh) script.

## How to Run PEPS

You may pass the config file as an argument to [PEPS.py](PEPS.py) like below or use the PEPS Jupyter Notebook [PEPS.ipynb](PEPS.ipynb). Since the program is in the early stage of development, we strongly recommend to use the Jupyter Notebook and keep eye on all the charts and intermediate data plotted in the notebook. The path to the config file is set in the first cell of the notebook.

```sh
$ python3 PEPS.py SampleData/config-large.json
```

# PEPS2

[PEPS2.ipynb](PEPS2.ipynb) uses a different way to simulate phenotype. Instead of assigning random phenotype and computing Risk genotypes based on random phenotype, PEPS2 computes the Risk genotype using mathematics.

PEPS2 requier "seed" parameter in the config file.

An example config file: [config-peps2.json](SampleData/config-peps.json)

An example processed notebook: [PEPS2.html](SampleData/PEPS2.html)

PEPS2 is available in notebook format only.
