#!/usr/bin/env python
# coding: utf-8

# In[4]:


import pandas as pd
import numpy as np
from pyplink import PyPlink
import statsmodels.api as sm
from IPython.display import clear_output


# In[5]:


expression = pd.read_table("GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz", delim_whitespace=True)
expression.set_index('TargetID', inplace=True)
expression.index = expression.index.str.replace("\.\d+", "", regex=True)
annotation = pd.read_csv("gene_annot.txt", sep="\t")
annotation.set_index('SYM', inplace=True)


# In[10]:


bim = PyPlink("LDREF/1000G.EUR.1").get_bim()
#bim


# In[13]:


pyplinkdict = {i: PyPlink("LDREF/1000G.EUR.{}".format(i)) for i in range(1,23)}


# In[21]:


#subset expression data for those w
expression = expression[pyplinkdict[1].get_fam()['fid'].values]


# In[20]:


for curr_chromosome in range(2,23):

    #expression = expression[pyplinkdict[curr_chromosome].get_fam()['fid'].values]
    progress = 0
    for symbol in expression.index:

        total_genes = len(expression.index)
        progress += 1

        if symbol not in annotation.index:
            continue
        gene_name, chromosome, gene_start = annotation.loc[symbol][["ID", "CHR", "START"]].values
        cis_start, cis_end = gene_start - 500000, gene_start + 500000
        if cis_start < 0:
            cis_start = 0
        max_stop = max(annotation["STOP"])
        if cis_end > max_stop:
            cis_end = max_stop

        curr_bim = pyplinkdict[chromosome].get_bim()
        curr_expression = expression.loc[symbol].values

        #select snps within gene body
        snp_subset = curr_bim[(curr_bim['pos'] > cis_start) & (curr_bim['pos'] < cis_end)].reset_index()

        #display(snp_subset)
        summary_matrix = []

        N = pyplinkdict[chromosome].get_nb_samples()
        snp_subset['N'] = [N] * snp_subset.shape[0]
        snp_subset['type'] = ['cc'] * snp_subset.shape[0]

        genotypes = []
        mafs = []
        for _, [snp, a1] in snp_subset[['snp', 'a1']].iterrows():
            #print(snp)
            alleles = pyplinkdict[chromosome].get_geno_marker(snp)
            unique_cts = np.unique(alleles, return_counts=True)[1]
            maf = min(unique_cts/sum(unique_cts))
            genotypes += [alleles]
            mafs += [maf]
        snp_subset['genotypes'] = genotypes
        snp_subset['MAF'] = mafs

        #print(len(curr_expression), len(genotypes[0]))
        snp_subset['results'] = snp_subset['genotypes'].apply(lambda x: sm.OLS(curr_expression, x).fit())

        snp_subset['pval'] = snp_subset['results'].apply(lambda x: x.pvalues[0])
        snp_subset['bse'] = snp_subset['results'].apply(lambda x: x.bse[0])
        snp_subset['varbeta'] = snp_subset['bse']**2
        snp_subset['effect_size'] = np.e**snp_subset['results'].apply(lambda x: x.params[0])
        snp_subset['beta'] = snp_subset['results'].apply(lambda x: x.params[0])

        snp_subset[['beta', 'varbeta', 'snp', 'pos', 'type', 'N', 'MAF', 'pval']].to_csv("eQTL_subsets/{}.csv".format(symbol), index=False)
        #print("{}".format(symbol))
        #clear_output(wait=True)
        print("Wrote {}/{} summary stats".format(progress, total_genes, curr_chromosome))

# In[ ]:




