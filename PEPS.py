#!/usr/bin/env python
# coding: utf-8

# # =========================================================
# # PolyEpi: Polygenic Phenotype with Higher-Order Epistasis Interactios
# # =========================================================

# # Path to the config file (Modify this before running the notebook.

# In[ ]:


# # Initialisation

# In[ ]:


import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import json

from pprint import pprint
from pdbio.vcfdataframe import VcfDataFrame

import itertools

from sklearn.feature_selection import SelectKBest, chi2
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from scipy import stats
import sklearn

TA = [1, .5, 1e-1, 1e-2, 1e-3, 1e-8, 1e-20]
configFilePath = sys.argv[1]


# # Read the config file

# In[ ]:


with open(configFilePath, 'r') as f:
    config = json.load(f)


# In[ ]:


shuffleSnps = config['shuffleSnps']
inputType = config['inputType']
dumpCSV = config['dumpCSV']

vcfInputPath = config['inputPrefix']+'.vcf'
csvInputPath = config['inputPrefix']+'.csv'

outputPrefix = config['outputPrefix']

pvalueThr = config['pvalueThr']
numTree = config['numTree']
numLoop = config['numLoop']


# # Compute total numebr of variables and number of requierd SNPs to form the variables

# In[ ]:


numVariables = 0
numSnpsNeeded = 0
maxOrder = 0

for v in config['variables']:
    v['numSNPs'] = v['numVar'] * v['numSnpsInVar']
    numVariables += v['numVar']
    numSnpsNeeded += v['numSNPs']
    if v['numSnpsInVar'] > maxOrder:
        maxOrder = v['numSnpsInVar']

config['numVariables'] = numVariables
config['numSnpsNeeded'] = numSnpsNeeded
config['maxOrder'] = maxOrder


# # Print config data and write it in "outputPrefix.config.json"

# In[ ]:


print("========== Configuration >>>")
pprint(config)
print("============================")
with open(outputPrefix+'.config.json', 'w') as outfile:
    json.dump(config, outfile)


# # Parse input genotype data from a VCF or CSV file
# ## If read from VCF file, the SNP id is set to CHROM:POS:REF:ALT

# In[ ]:


if inputType == 'vcf':
    vcfdf = VcfDataFrame(path=vcfInputPath)
    df = vcfdf.df
    df['SNP'] = df['#CHROM'].astype(str) + ':' + df['POS'].astype(
        str) + ':' + df['REF'].astype(str) + ':' + df['ALT'].astype(str)
    df = df.set_index('SNP')
    snpData = df.iloc[:, 9:].replace(['0/0', '0/1', '1/1'], [0, 1, 2])
    if dumpCSV:
        snpData.to_csv(csvInputPath)
elif inputType == 'csv':
    snpData = pd.read_csv(csvInputPath)
    snpData = snpData.set_index('SNP')
else:
    print("Incorrect inputType (should be 'vcf' or 'csv')")
    exit()


# In[ ]:


snpData.iloc[:5, :5]


# # There should be enough SNPs in the input file to create all variables

# In[ ]:


print("Number of SNPs in the input file: ", snpData.shape[0])
print("Number of SNPs needed: ", numSnpsNeeded)

if snpData.shape[0] < numSnpsNeeded:
    print("There are not enough SNPs in the input file")
    exit(1)
else:
    print("There are enough SNPs in the input file")


# # Suffle SNPs if asked in the config file.
# ## When SNPs are shuffled different set of SNPs used to form each variables each time

# In[ ]:


if shuffleSnps:
    snpData = snpData.sample(frac=1)


# # Transpose the genotype data and print number of snps and samples
# ## Also rename 0/0, 0/1 and 1/1 to R, H and A

# In[ ]:


snpData = snpData.T


# In[ ]:


df = snpData.replace([0, 1, 2], ['R', 'H', 'A'])
numSampels = df.shape[0]
numSNPs = df.shape[1]
print("number of sample", numSampels)
print("number of snp", numSNPs)


# In[ ]:


df.iloc[:5, :5]


# # Form variables from SNPs
# ## A variable could be a SNPs or a set of Interactive SNPs
# ## First identify whcih SNPs belong to each variable and then form the variables
# ## Naming of variables: O3V4 is the 4th variable with 3-interactive SNPs
# ## Write Variables SNPs infor in "outputPrefix.varData.csv"

# In[ ]:


colNames = list()  # to store variable names
for o, v in enumerate(config['variables']):
    for i in range(0, v['numVar']):
        colNames.append('O'+str(o+1)+'V'+str(i+1))


# In[ ]:


rowNames = ['order']
for o in range(maxOrder):
    rowNames.append('snp_'+str(o+1))


# In[ ]:


varData = pd.DataFrame(index=rowNames, columns=colNames)


# In[ ]:


idx = 0
for o, v in enumerate(config['variables']):
    for i in range(0, v['numVar']):
        name = 'O'+str(o+1)+'V'+str(i+1)
        varData.at['order', name] = str(o+1)
        for k in range(0, v['numSnpsInVar']):
            snp = 'snp_'+str(k+1)
            varData.at[snp, name] = df.columns[idx]
            idx += 1
varData = varData.fillna('---')


# In[ ]:


varData.to_csv(outputPrefix+'.varData.csv')


# In[ ]:


varData.iloc[:6, -5:]


# # Form Variable Genotype and write it to "outputPrefix.varGT.csv"
# ## For variables with more than one SNPs the genotype is the concatination of all SNPs involved
# ## For example RHA, ARH and AAR could be genotype value of a variable with 3 snps

# In[ ]:


varGT = df.iloc[:, -1:0].copy()


# In[ ]:


for o, v in enumerate(config['variables']):
    for i in range(0, v['numVar']):
        name = 'O'+str(o+1)+'V'+str(i+1)
        varGT[name] = ''
        for k in range(0, v['numSnpsInVar']):
            snp = 'snp_'+str(k+1)
            varGT[name] = varGT[name] + df[varData.loc[snp, name]]


# In[ ]:


varGT.to_csv(outputPrefix+'.varGT.csv')


# In[ ]:


varGT.iloc[:5, -5:]


# # This function compute Risks based on given phenotype

# In[ ]:


def FindRisks(df):
    varNames = df.columns[0:-1]

    case, ctrl = df[df['lbl'] == 1], df[df['lbl'] == 0]

    GS = dict()

    for i, v in enumerate(varNames):

        caseC = case[v].value_counts()
        ctrlC = ctrl[v].value_counts()

        count = caseC.to_frame().join(ctrlC.to_frame(), lsuffix='_case',
                                      rsuffix='_ctrl', how='outer').fillna(0)
        #count = count.iloc[:,:-1]
        count.columns = ['case', 'ctrl']
        count['p'] = (count['case']+1)/(count['ctrl']+1)
        #count['p'] /= count['p'].sum()

        GS[v] = count
    return GS


# # This function update phenotype based on given Risks

# In[ ]:


def Predict(df, GS):
    df['prob'] = 0.0
    for sample in df.index:
        prob = 0.0
        for var in df.drop(['lbl', 'prob'], axis=1).columns:
            G = GS[var]
            g = df.loc[sample, var]
            if g in G.index:
                prob += G.loc[g, 'p']
        df.at[sample, 'prob'] = prob

    mid = df['prob'].mean()
    df['lbl'] = df['prob'].apply(lambda x: 1 if(x > mid) else 0)


# # Assign random phenotype to samples
# # 0 for control and 1 for case

# In[ ]:


varGT['lbl'] = np.random.randint(0, 2, varGT.shape[0])


# In[ ]:


varGT[['lbl']].head(5)


# # Compute and plot chi2-pvalue (log10) of the variables for the random phenotype

# In[ ]:


features = varGT.columns[:-1]
corrDict = dict()
for v in features:
    corrDict[v] = stats.chi2_contingency(
        pd.crosstab(varGT['lbl'], varGT[v]).values)[1]
a = np.asarray(list(corrDict.values()))
b = - np.log10(a)
plt.plot(np.sort(b))
nsat = list()
for t in TA:
    nsat.append([t, np.where(a < t)[0].shape[0]])
x = pd.DataFrame(nsat)
x.columns = ['p-value', 'number of SNPs exceed the p-value']
x.set_index('p-value')


# # >>>> Phenotype Simulation Feedback loop <<<<

# In[ ]:


for i in range(numLoop):
    print('loop', i+1)
    Risks = FindRisks(varGT)
    Predict(varGT, Risks)
    varGT = varGT.drop(['prob'], axis=1)
print("Done")


# # Write Phenotype into a file outputPrefix.pheno.csv

# In[ ]:


phen = varGT[['lbl']].copy()
phen.index.name = 'sample'
phen.to_csv(outputPrefix+'.pheno.csv')


# In[ ]:


phen.head()


# # Compute and plot chi2-pvalue (log10) of the variables for the simulated phenotype

# In[ ]:


features = varGT.columns[:-1]
corrDict = dict()
for v in features:
    corrDict[v] = stats.chi2_contingency(
        pd.crosstab(varGT['lbl'], varGT[v]).values)[1]
a = np.asarray(list(corrDict.values()))
b = - np.log10(a)
plt.plot(np.sort(b))
nsat = list()
for t in TA:
    nsat.append([t, np.where(a < t)[0].shape[0]])
x = pd.DataFrame(nsat)
x.columns = ['p-value', 'number of SNPs exceed the p-value']
x.set_index('p-value')


# # Filter Variable that satisfy p-value treshold as Truth Variable and write it to "outputPrefix.varDataTruth.csv"

# In[ ]:


variables = list()
for v in corrDict:
    if corrDict[v] < pvalueThr:
        variables.append(v)

varDataTruth = varData[variables]
varDataTruth.to_csv(outputPrefix+'.varDataTruth.csv')


# In[ ]:


print("Number of Variables: ", varData.shape[1])
print("Number of Truth Variables: ", varDataTruth.shape[1])


# In[ ]:


varDataTruth.iloc[:6, -5:]


# # Filter SNPs included in All and Truth Variables
# ## write the Truth SNP names in  "outputPrefix.TruthSNP.csv"

# In[ ]:


snps = np.unique(varDataTruth.replace(
    np.nan, '', regex=True).drop('order').values.ravel())[1:]
snpDataTruth = snpData.loc[:, snpData.columns.isin(snps)]

snps2 = np.unique(varData.replace(
    np.nan, '', regex=True).drop('order').values.ravel())[1:]
snpDataVar = snpData.loc[:, snpData.columns.isin(snps2)]

pd.DataFrame(snps).rename(columns={0: 'v'}).to_csv(
    outputPrefix+'.TruthSNP.csv', index=False)


# In[ ]:


print("Number of SNPs used to form Variables: ", snpDataVar.shape[1])
print("Number of Truth SNPs (SNPs in Truth Variables): ",
      snpDataTruth.shape[1])


# In[ ]:


snpDataTruth.iloc[:5, -5:]


# # Add Phenotype to SNP Data

# In[ ]:


snpDataTruth['lbl'] = varGT['lbl']
snpDataVar['lbl'] = varGT['lbl']
snpData['lbl'] = varGT['lbl']


# # This function used to predict lable using RandomForest
# ## 75% training and 25% test

# In[ ]:


def RF_AUC(dfx, nTree):
    df = dfx.copy()
    features = df.columns[:-1]

    df['is_train'] = np.random.uniform(0, 1, len(df)) <= .75
    train, test = df[df['is_train'] == True], df[df['is_train'] == False]
    clf = RandomForestClassifier(n_jobs=2, n_estimators=nTree, random_state=0)
    clf.fit(train[features], train['lbl'])
    # clf.predict(test[features])
    prob = clf.predict_proba(test[features])
    y_true = test['lbl']
    y_scores = prob[:, 1]
    return clf, roc_auc_score(y_true, y_scores)


# # Train and test RandomForest for TruthSNP as well as for all SNP in the input file
# ## Print AUC
# ## Plot Importance Score

# In[ ]:


clfTruth, AucTruth = RF_AUC(snpDataTruth, nTree=numTree)
clfVar,   AucVar = RF_AUC(snpDataVar, nTree=numTree)
clf,      Auc = RF_AUC(snpData, nTree=numTree)


# In[ ]:


print("AUC Truth SNPs            : ", AucTruth)
print("AUC All SNPs in Variables : ", AucVar)
print("AUC All SNPs in input file: ", Auc)


# In[ ]:


pd.DataFrame(np.sort(clfTruth.feature_importances_)).plot(title="Truth SNP")
pd.DataFrame(np.sort(clfVar.feature_importances_)
             ).plot(title="All SNP in Variables")
pd.DataFrame(np.sort(clf.feature_importances_)).plot(
    title="All SNP in input file")
