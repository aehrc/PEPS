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

The [SampleData](SampleData) directory includes examples of input, output and config files. The [HTML_Notebook](HTML_Notebook) contains processed PEPS notebooks in HTML format.

**Input:**

There are two input file both of them are a subset of 1000Genomes data available in vcf and csv format. The `small` input with ~162 SNPs and the `large` input with 4969 SNPs.

**config:**

There are 2 config files. [config-small.json](SampleData/config-small.json) and [config-large.json](SampleData/config-large.json) are made to process small vcf file and large csv file respectively.

**Output:**

Output files related to both config files above.

**PEPS processed notebook**

[PEPS-small.html](HTML_Notebook/PEPS-small.html) and [PEPS-large.html](HTML_Notebook/PEPS-large.html) are examples that show how PEPS notebook look likes after processing small and large input.

## VCF to CSV conversion

To parse VCF files, PEPS uses a library called `pdbio`. This library is slow and only works for a tiny vcf file. **We strongly recommend to prepare input data in CSV format.** To do so you can use the [VCF_2_CSV.sh](VCF_2_CSV.sh) script.

## How to Run PEPS3

Use the PEPS3 Jupyter Notebook [PEPS3.ipynb](PEPS.ipynb). We strongly recommend to keep eye on all the charts and intermediate data plotted in the notebook. The path to the config file is set in the first cell of the notebook.


# PEPS3

[PEPS3.ipynb](PEPS3.ipynb) uses a probabilistic model to calculate a genotype
for a set of variables within a population.

Let the desired frequency of cases be _Q_ and, similarly, the
frequency of controls is _P_ = 1-_Q_. If we let the
_k_<sup>th</sup> variable on average decrease the probability of a sample being a
control by the fraction _q_<sub>_k_</sub> then for _g_ variables we arrive at:

_Q_ = (1-_q_<sub>1</sub>)(1-_q_<sub>2</sub>) ... (1-_q_<sub>_g_</sub>)

In the interests of making each variable similar in importance, we set each
_q_<sub>_k_</sub> equal to the same value _q_, which leaves us with
_Q_ = (1-_q_)<sup>_g_</sup>, and rearrange to find _q_ = 1 - _Q_<sup>1/_g_</sup>.
As _q_ is an average over all the values for a variable, it can be expressed as

_q_ = _p_<sub>1</sub>_f_<sub>1</sub> + _p_<sub>2</sub>_f_<sub>2</sub> + ... +
_p_<sub>_n+1_</sub>_f_<sub>_n+1_</sub>

for _n+1_ possible values of the variable, where _f<sub>k</sub>_ is the
frequency of the _k_<sup>th</sup> variable. What remains is to choose the
effect on the phenotype _p_<sub>_k_</sub>
for each value of the variable. We choose to set the value that is most
frequent in the population to 0, and set the others such that they all
contribute an equal amount, i.e. _p_<sub>1</sub>_f_<sub>1</sub> = _p_<sub>2</sub>_f_<sub>2</sub> = ... =
_p_<sub>_n_</sub>_f_<sub>_n_</sub>. So we have for _p_<sub>_k_</sub>

_q_ = _np_<sub>_k_</sub>_f_<sub>_k_</sub>

_p_<sub>_k_</sub> = _q_/(_nf<sub>k</sub>_)

_p_<sub>_k_</sub> = (1-_Q_<sup>1/_g_</sup>)/_nf<sub>k</sub>_

This _p_ is the calculated for each value of each variable
(with a maximum value of 1) and represents the fraction by which the
probability of a sample with that value being in the control group is reduced.

The probability that a sample _s_ is in the control group is then

_Q<sub>s</sub>_ = (1-*p<sub>s1</sub>)(1-*p<sub>s2</sub>) ... (1-_p_<sub>_sg_</sub>)

where _p_<sub>_sk_</sub> is the _p_ associated with the value that sample _s_ has for variable _k_.

The cases and controls could then be calculated stochastically, but in the
interests of minimising noise and maximising the importance of each variable,
the samples are sorted by increasing _Q<sub>s</sub>_ and the first _P_ fraction
are selected as cases.

PEPS3 uses the "seed" parameter in the config file to ensure that a phenotype can be
reproduced. If the "shuffleSnps" parameter is true the seed is ignored and the
phenotype will be random.

PEPS2 does not require "numLoop" and "pvalueThr" parameter in the config file.

An example config file: [config-peps3.json](SampleData/config-peps3.json)

An example processed notebook: [PEPS3.html](HTML_Notebook/PEPS3.html)

It also considers the SNPs included in Epistasis interactions also appear indivudually as well.
For exampel if O3V5 made of 3 SNP (rs123, rs456, rs789) then each of these SNPs will form a Variable (O3V5S1, O3V5S2, O3V5S3) and independantly affect the phenotype
