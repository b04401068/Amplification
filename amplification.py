import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import numpy as np

#import sys
#sys.path.insert(0, './utility/')

def save_obj( data, name ):
    with open(name,'wb') as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL);


def load_obj( name ):
    with open(name,'rb') as f:
        return pickle.load(f);


###############################################################
#1-extracting cbioportal data
###############################################################
'''from cbio_py import cbio_mod as cb
study = cb.getAllStudies()
study = pd.DataFrame(study)
study.to_pickle('data/study.pkl')
#'''

'''from cbio_py import cbio_mod as cb
study = pd.read_pickle('data/study.pkl')
#print(study.iloc[0,:])
study = study.loc[ study.loc[:,'cancerTypeId'] == 'brca',:]
study = study.loc[ study.loc[:,'referenceGenome'] == 'hg19',:]
#print(study)
#print(study.iloc[0,:])
study.to_csv('data/study_breast.tsv',sep = '\t')
for i in study.index:
    print(study.loc[i,'name'])
    print(study.loc[i,'description'])
#'''

'''from cbio_py import cbio_mod as cb
study = pd.read_csv('data/study_breast.tsv',sep = '\t', index_col = 0)
study = study.loc[ study.loc[:,'select'] == 1,:]
print(study)
clinical_info = []
for i in study.index:
    print(study.loc[i,'name'])
    print(study.loc[i,'studyId'])
    print(study.loc[i,'description'])
    clinical = cb.getAllClinicalDataInStudy(study.loc[i,'studyId'], return_type = 'dict')
    clinical = pd.DataFrame(clinical)
    print(clinical.iloc[0,:])
    clinical_info.append(clinical)
clinical_info = pd.concat(clinical_info, axis = 0)
print(clinical_info)
study.to_pickle('data/study_info_0614.pkl')
clinical_info.to_pickle('data/clinical_info_0614.pkl')
#'''

'''from cbio_py import cbio_mod as cb
study = pd.read_csv('data/study_breast.tsv',sep = '\t', index_col = 0)
study = study.loc[ study.loc[:,'select'] == 1,:]
print(study)
clinical_info = []
for i in study.index:
    print(study.loc[i,'name'])
    print(study.loc[i,'studyId'])
    print(study.loc[i,'description'])
    clinical = cb.getClinicalAttributesByStudyId(study.loc[i,'studyId'], return_type = 'dict')
    clinical = pd.DataFrame(clinical)
    print(clinical)
    if False:
        for j in clinical.index:
            clinicalattr = cb.getClinicalAttributeInStudy(study.loc[i,'studyId'], clinical.loc[j,'clinicalAttributeId'], return_type = 'dict')
            print(clinicalattr)
            clinicalattr = pd.DataFrame(clinicalattr)
            clinical_info.append(clinicalattr)
#clinical_info = pd.concat(clinical_info, axis = 0)
#print(clinical_info)
#study.to_pickle('data/study_info_0614.pkl')
#clinical_info.to_pickle('data/clinical_att_0614.pkl')
#'''

'''clinical_info = pd.read_pickle('data/clinical_info_0614.pkl')
clinical_info = clinical_info.pivot(columns = 'clinicalAttributeId', index = 'sampleId', values = 'value')
print(clinical_info.iloc[0,:])
print(clinical_info)
clinical_info.to_pickle('data/clinical_info_0614_wide.pkl')
clinical_info.to_csv('data/clinical_info_0614_wide.tsv', sep = '\t')
#'''

'''
from cbio_py import cbio_mod as cb
study = pd.read_csv('data/study_breast.tsv',sep = '\t', index_col = 0)
study = study.loc[ study.loc[:,'select'] == 1,:]
print(study)
sample = pd.read_csv('data/clinical_data_from_cbioportal.tsv', sep = '\t')
print(sample)
for i in study.index:
    print(study.loc[i,'studyId'])
    tmp = sample.loc[ sample.loc[:,'Study ID'] == study.loc[i,'studyId'],:]
    tmp = tmp.dropna( how = 'all', axis = 1)
    print(tmp.iloc[0,:])
    print(tmp.columns.tolist())
    if 'Subtype' in tmp.columns:
        print(tmp.iloc[0,:].tolist())
        print(tmp.loc[:,'Subtype'].value_counts())
#'''

###############################################################
#CNA and mRNA
###############################################################
'''study = pd.read_csv('data/study_breast.tsv',sep = '\t', index_col = 0)
study = study.loc[ study.loc[:,'select'] == 1,:]
print(study)
for i in study.loc[:,'studyId']:
    print(i)
    cna = pd.read_csv('R_data/'+i+'_cna.tsv', sep = '\t', index_col = 0)
    print(cna)
    cna.to_pickle('data/'+i+'_cna.pkl')
    try:
        mRNA = pd.read_csv('R_data/'+i+'_mRNA.tsv', sep = '\t', index_col = 0)
        print(mRNA)
        mRNA.to_pickle('data/'+i+'_mRNA.pkl')
    except:
        print('no mRNA')
#'''

'''study = pd.read_csv('data/study_breast.tsv',sep = '\t', index_col = 0)
study = study.loc[ study.loc[:,'select'] == 1,:]
print(study)
cna_all = []
mRNA_all = []
for i in study.loc[:,'studyId']:
    print(i)
    cna = pd.read_pickle('data/'+i+'_cna.pkl')
    cna = cna.drop(columns = ['group','group_name'])
    cna_all.append(cna)
    try:
        mRNA = pd.read_pickle('data/'+i+'_mRNA.pkl')
        mRNA = mRNA.drop(columns = ['group','group_name'])
        if i == 'brca_mbcproject_2022':
            mRNA = mRNA.reindex(index = cna.index)
        mRNA_all.append(mRNA)
    except:
        print('no mRNA')
cna_all = pd.concat(cna_all, axis = 1)
print(cna_all)
cna_all.to_pickle('data/cna_all.pkl')
mRNA_all = pd.concat(mRNA_all, axis = 1)
print(mRNA_all)
mRNA_all.to_pickle('data/mRNA_all.pkl')
#'''

#remove NA for furthur processing
'''cna = pd.read_pickle('data/cna_all.pkl')
print(cna)
mrna = pd.read_pickle('data/mRNA_all.pkl')
cnacnt = cna.isna().mean(axis = 1)
cna = cna.loc[cnacnt<=0.2,:]
print(cna)
mrna = mrna.reindex(index = cna.index)
print(mrna)
cna.to_pickle('data/cna_all_na_filter.pkl')
mrna.to_pickle('data/mrna_all_na_filter.pkl')
#'''

'''
cna = pd.read_pickle('data/cna_all_na_filter.pkl')
mrna = pd.read_pickle('data/mrna_all_na_filter.pkl')
from scipy import stats
result = pd.DataFrame(np.nan, index=cna.index, columns = ['spearman_r','spearman_p','pearson_r','pearson_p'])
for i in cna.index:
    a = cna.loc[[i],:].T
    b = mrna.loc[[i],:].T
    tmp = pd.concat([a,b], axis = 1)
    tmp = tmp.dropna( how = 'any')
    if len(tmp) > 0:
        result.loc[i,'spearman_r'], result.loc[i,'spearman_p'] = stats.spearmanr( tmp.iloc[:,0], tmp.iloc[:,1])
        result.loc[i,'pearson_r'], result.loc[i,'pearson_p'] = stats.pearsonr( tmp.iloc[:,0], tmp.iloc[:,1])
print(result)
result.to_pickle('data/cor_cna_mrna.pkl')
#'''

'''
import pandas as pd
import numpy as np
from scipy import stats

# Load the data
cna = pd.read_pickle('data/cna_all_na_filter.pkl')
mrna = pd.read_pickle('data/mrna_all_na_filter.pkl')

# Initialize result DataFrame
result = pd.DataFrame(np.nan, index=cna.index, columns=['wilcoxon_stat', 'wilcoxon_p'])
result = pd.DataFrame(np.nan, index=cna.index, columns=['wilcoxon_stat', 'wilcoxon_p', 'wilcoxon_fdr'])

# Iterate through each gene (index = gene name)
for i in cna.index:
    a = cna.loc[[i], :].T  # CNA values for gene i
    b = mrna.loc[[i], :].T  # mRNA values for gene i

    tmp = pd.concat([a, b], axis=1)  # Combine into one DataFrame
    tmp.columns = ['cna', 'mrna']
    tmp = tmp.dropna(how='any')  # Drop rows with NA

    # Split into two groups: CNA == 2 vs CNA != 2
    group_2 = tmp[tmp['cna'] == 2]['mrna']
    group_other = tmp[tmp['cna'] != 2]['mrna']

    # Only run the test if both groups have enough samples
    if len(group_2) > 0 and len(group_other) > 0:
        stat, p = stats.ranksums(group_2, group_other)  # Wilcoxon rank-sum (Mann-Whitney U test)
        result.loc[i, 'wilcoxon_stat'] = stat
        result.loc[i, 'wilcoxon_p'] = p

# Save results
from statsmodels.stats.multitest import multipletests
valid_pvals = result['wilcoxon_p'].dropna()
rejected, pvals_corrected, _, _ = multipletests(valid_pvals, method='fdr_bh')
result.loc[valid_pvals.index, 'wilcoxon_fdr'] = pvals_corrected
print(result)
result.to_pickle('data/wilcoxon_cna2_vs_other_mrna.pkl')
#'''

'''cna = pd.read_pickle('data/cna_all_na_filter.pkl')
cna = cna.T
cnacnt = (cna>=2).mean()
cnacnt = pd.DataFrame(cnacnt.values, index = cnacnt.index, columns = ['amp_freq'])
cnacnt.loc[:,'del_freq'] = (cna<=-2).mean()
cnacnt.index.name = 'gene'
print(cnacnt)
cnacnt.to_csv('R_data/CNA_freq.csv')
#'''

#clinical data cleaning
'''study = pd.read_csv('data/study_breast.tsv',sep = '\t', index_col = 0)
study = study.loc[ study.loc[:,'select'] == 1,:]
print(study)
sample = []
for i in study.loc[:,'studyId']:
    print(i)
    clin = pd.read_csv('data/'+i+'_clinical_data.tsv', sep = '\t', index_col = 'Sample ID')
    print(clin)
    if i == 'brca_igr_2015':
        tmp = clin.loc[:,['Study ID','Patient ID']]
        tmp.loc[:,'ER'] = np.nan
        tmp.loc[ clin.loc[:,'IHC-HER2'].str.contains('HR+'),'ER'] = 'Positive'
        tmp.loc[ clin.loc[:,'IHC-HER2'].str.contains('HR-'),'ER'] = 'Negative'
        tmp.loc[:,'HER2'] = np.nan
        tmp.loc[ clin.loc[:,'IHC-HER2'].str.contains('HER+'),'HER2'] = 'Positive'
        tmp.loc[ clin.loc[:,'IHC-HER2'].str.contains('HER-'),'HER2'] = 'Negative'
        tmp.loc[:,'TMB'] = clin.loc[:,'TMB (nonsynonymous)']
        #tmp.loc[:,'Metastasis'] = True
        sample.append(tmp)
    if i == 'brca_cptac_2020':
        tmp = clin.loc[:,['Study ID','Patient ID','Age']]
        tmp.loc[:,'ER'] = clin.loc[:,'ER Updated Clinical Status']
        tmp.loc[:,'HER2'] = clin.loc[:,'ERBB2 Updated Clinical Status']
        tmp.loc[ clin.loc[:,'TNBC Updated Clinical Status'] == 'Positive', 'ER'] = 'Negative'
        tmp.loc[ clin.loc[:,'TNBC Updated Clinical Status'] == 'Positive', 'HER2'] = 'Negative'
        tmp.loc[:,'TMB'] = clin.loc[:,'TMB (nonsynonymous)']
        tmp.loc[:,'stage'] = clin.loc[:,'Tumor Stage']
        tmp.loc[tmp.loc[:,'stage'].str.contains('Stage III').fillna(False),'stage'] = 3
        tmp.loc[tmp.loc[:,'stage'].str.contains('Stage II').fillna(False),'stage'] = 2
        tmp.loc[tmp.loc[:,'stage'].str.contains('Stage I').fillna(False),'stage'] = 1
        tmp.loc[clin.loc[:,'Tumor Stage'].isna(),'stage'] = np.nan
        sample.append(tmp)
    if i == 'brca_mbcproject_2022':
        tmp = clin.loc[:,['Study ID','Patient ID']]
        age = clin.loc[:,'MedR Age at Diagnosis'].str.split('-',expand = True)
        age.iloc[age.iloc[:,0].str.contains('<').fillna(False),0] = '21'
        age = age.astype(float)
        tmp.loc[:,'Age'] = age.iloc[:,0]+2
        if False:
            check = clin.loc[['RP-1156_MBCProject_07t2uAFn_BLOOD_P_v2_Exome'],:]
            j = 0
            while j < len(check.columns):
                print(check.iloc[:,j:j+50].T)
                j = j + 50
        tmp.loc[:,'ER'] = clin.loc[:,'PATH Estrogen Receptor Status']
        need_sub = (tmp.loc[:,'ER'] != 'POSITIVE') & (tmp.loc[:,'ER'] != 'NEGATIVE')
        tmp.loc[ need_sub ,'ER'] = clin.loc[ need_sub, 'MedR Diagnostic ER Status']
        need_sub = (tmp.loc[:,'ER'] != 'POSITIVE') & (tmp.loc[:,'ER'] != 'NEGATIVE')
        tmp.loc[:,'ER'] = tmp.loc[:,'ER'].replace('POSITIVE','Positive')
        tmp.loc[:,'ER'] = tmp.loc[:,'ER'].replace('NEGATIVE','Negative')
        tmp.loc[need_sub,'ER'] = np.nan
        tmp.loc[:,'HER2'] = clin.loc[:,'PATH HER2 Status']
        tmp.loc[ tmp.loc[:,'HER2'] == 'EQUIVOCAL', 'HER2'] = 'NEGATIVE'
        need_sub = (tmp.loc[:,'HER2'] != 'POSITIVE') & (tmp.loc[:,'HER2'] != 'NEGATIVE')
        tmp.loc[ need_sub ,'HER2'] = clin.loc[ need_sub, 'MedR Diagnostic HER2 Status']
        tmp.loc[ tmp.loc[:,'HER2'] == 'EQUIVOCAL', 'HER2'] = 'NEGATIVE'
        need_sub = (tmp.loc[:,'HER2'] != 'POSITIVE') & (tmp.loc[:,'HER2'] != 'NEGATIVE')
        tmp.loc[:,'HER2'] = tmp.loc[:,'HER2'].replace('POSITIVE','Positive')
        tmp.loc[:,'HER2'] = tmp.loc[:,'HER2'].replace('NEGATIVE','Negative')
        tmp.loc[need_sub,'HER2'] = np.nan
        tmp.loc[:,'Grade'] = clin.loc[:,'PATH Sample Grade']
        tmp.loc[tmp.loc[:,'Grade'].str.contains('III').fillna(False),'Grade'] = 3
        tmp.loc[tmp.loc[:,'Grade'].str.contains('II').fillna(False),'Grade'] = 2
        tmp.loc[tmp.loc[:,'Grade'].str.contains('I').fillna(False),'Grade'] = 1
        need_sub = ~tmp.loc[:,'Grade'].isin([1,2,3])
        tmp.loc[ need_sub,'Grade'] = clin.loc[need_sub, 'MedR Diagnostic Grade'] 
        tmp.loc[tmp.loc[:,'Grade'].str.contains('III').fillna(False),'Grade'] = 3
        tmp.loc[tmp.loc[:,'Grade'].str.contains('II').fillna(False),'Grade'] = 2
        tmp.loc[tmp.loc[:,'Grade'].str.contains('I').fillna(False),'Grade'] = 1
        need_sub = ~tmp.loc[:,'Grade'].isin([1,2,3])
        tmp.loc[need_sub,'Grade'] = np.nan
        tmp.loc[:,'TMB'] = clin.loc[:,'TMB (nonsynonymous)']
        tmp.loc[:,'stage'] = clin.loc[:,'MedR Stage at Diagnosis']
        tmp.loc[tmp.loc[:,'stage'].str.contains('3'),'stage'] = 3
        tmp.loc[tmp.loc[:,'stage'].str.contains('2').fillna(False),'stage'] = 2
        tmp.loc[tmp.loc[:,'stage'].str.contains('1').fillna(False),'stage'] = 1
        tmp.loc[tmp.loc[:,'stage'].str.contains('4').fillna(False),'stage'] = 4
        tmp.loc[~tmp.loc[:,'stage'].isin([1,2,3,4]),'stage'] = np.nan
        tmp.loc[:,'Histology'] = clin.loc[:,'PATH Sample Histology']
        need_sub = tmp.loc[:,'Histology'] == 'N/A_BLOOD_SAMPLE'
        tmp.loc[need_sub,'Histology'] = clin.loc[need_sub,'MedR Diagnostic Histology']
        tmp.loc[ (tmp.loc[:,'Histology'] == 'ABSTRACTION_PENDING') | (tmp.loc[:,'Histology'] == 'NOT_FOUND_IN_RECORD'),'Histology'] = np.nan 
        print(tmp.loc[:,'Histology'].value_counts())
        sample.append(tmp)
    if i == 'brca_tcga':
        tmp = clin.loc[:,['Study ID','Patient ID']]
        tmp.loc[:,'Age'] = clin.loc[:,'Diagnosis Age']
        tmp.loc[:,'ER'] = clin.loc[:,'ER Status By IHC']
        tmp.loc[:,'HER2'] = clin.loc[:,'HER2 fish status']
        need_sub = (tmp.loc[:,'HER2'] != 'Positive') & (tmp.loc[:,'HER2'] != 'Negative')
        tmp.loc[ need_sub, 'HER2'] = clin.loc[ need_sub, 'HER2 ihc score']
        tmp.loc[ tmp.loc[:,'HER2'] == 3, 'HER2'] = 'Positive'
        tmp.loc[ (tmp.loc[:,'HER2'] == 2) | (tmp.loc[:,'HER2'] == 1), 'HER2'] = 'Negative'
        need_sub = (tmp.loc[:,'HER2'] != 'Positive') & (tmp.loc[:,'HER2'] != 'Negative')
        tmp.loc[ need_sub, 'HER2'] = clin.loc[ need_sub, 'IHC-HER2']
        tmp.loc[ tmp.loc[:,'HER2'] == 'EQUIVOCAL','HER2'] = 'Negative'
        print(tmp.loc[:,'HER2'].value_counts())
        tmp.loc[:,'stage'] = clin.loc[:,'Neoplasm Disease Stage American Joint Committee on Cancer Code']
        tmp.loc[tmp.loc[:,'stage'].str.contains('Stage IV').fillna(False),'stage'] = 4
        tmp.loc[tmp.loc[:,'stage'].str.contains('Stage III').fillna(False),'stage'] = 3
        tmp.loc[tmp.loc[:,'stage'].str.contains('Stage II').fillna(False),'stage'] = 2
        tmp.loc[tmp.loc[:,'stage'].str.contains('Stage I').fillna(False),'stage'] = 1
        tmp.loc[~tmp.loc[:,'stage'].isin([1,2,3,4]),'stage'] = np.nan
        tmp.loc[:,'TMB'] = clin.loc[:,'TMB (nonsynonymous)']
        tmp.loc[:,'Histology'] = clin.loc[:,'Neoplasm Histologic Type Name']
        tmp.loc[:,'Histology'] = tmp.loc[:,'Histology'].str.replace('Infiltrating Lobular Carcinoma','ILC')
        tmp.loc[:,'Histology'] = tmp.loc[:,'Histology'].str.replace('Infiltrating Ductal Carcinoma','IDC')
        print(tmp.loc[:,'Histology'].value_counts())
        sample.append(tmp)
    if i == 'brca_metabric':
        tmp = clin.loc[:,['Study ID','Patient ID']]
        tmp.loc[:,'Age'] = clin.loc[:,'Age at Diagnosis']
        tmp.loc[:,'ER'] = clin.loc[:,'ER Status']
        tmp.loc[:,'HER2'] = clin.loc[:,'HER2 Status']
        tmp.loc[:,'stage'] = clin.loc[:,'Tumor Stage']
        tmp.loc[~tmp.loc[:,'stage'].isin([1,2,3,4]),'stage'] = np.nan
        print(tmp.loc[:,'stage'].value_counts())
        tmp.loc[:,'TMB'] = clin.loc[:,'TMB (nonsynonymous)']
        tmp.loc[:,'Histology'] = clin.loc[:,'Oncotree Code']
        tmp.loc[tmp.loc[:,'Histology'] == 'BRCA','Histology'] = np.nan
        tmp.loc[:,'Grade'] = clin.loc[:,'Neoplasm Histologic Grade']
        sample.append(tmp)
#tcga = pd.read_csv('data/brca_tcga_pub_clinical_data.tsv', sep = '\t', index_col = 2)
sample = pd.concat(sample, axis = 0)
sample.loc[:,'molecular'] = np.nan
sample.loc[ sample.loc[:,'HER2'] == 'Positive', 'molecular'] = 'HER2+'
need_sub = sample.loc[:,'molecular'].isna()
sample.loc[ (sample.loc[:,'ER'] == 'Positive')&need_sub, 'molecular'] = 'ER+'
need_sub = sample.loc[:,'molecular'].isna()
sample.loc[ (sample.loc[:,'ER'] == 'Negative')&need_sub, 'molecular'] = 'TNBC'
need_sub = (sample.loc[:,'Histology'] != 'ILC')&(sample.loc[:,'Histology'] != 'IDC')&(~sample.loc[:,'Histology'].isna())
sample.loc[need_sub,'Histology'] = 'other'
print(sample.loc[:,'Histology'].value_counts())
print(sample.loc[:,'molecular'].value_counts())
print(sample)
sample.to_pickle('data/clinical_info_0623.pkl')
#'''

'''
sample = pd.read_pickle('data/clinical_info_0623.pkl')
print(sample)
sample.to_csv('data/clinical_info.txt',sep = '\t')
#'''

#subtype chi/fisher exact
'''from scipy import stats
cna = pd.read_pickle('data/cna_all_na_filter.pkl')
sample = pd.read_pickle('data/clinical_info_0623.pkl')
cna = cna.T
cna.index = cna.index.str.replace('.','-')
sample = sample.reindex(index = cna.index)
sample = sample.loc[ ~sample.loc[:,'molecular'].isna(),:]
cna = cna.reindex(index = sample.index)
stat = ['amp_freq','amp_p','del_freq','del_p']
#stat = ['amp_freq','amp_p']
col = []
for i in ['ER+','HER2+','TNBC']:
    for j in stat:
        col.append( i+'_'+j)
result = pd.DataFrame(np.nan, index = cna.columns, columns = col)
pos = []
pos_num = []
for i in ['ER+','HER2+','TNBC']:
    positive = sample.loc[:,'molecular'] == i
    pos.append(positive)
    pos_num.append(positive.sum())
print(pos_num)
molecular = ['ER+','HER2+','TNBC']
total_num = len(cna.index)
for i in cna.columns:
    amp = cna.loc[:,i] >= 2 
    de = cna.loc[:,i] <= -2
    amp_num = amp.sum()
    de_num = de.sum()
    for j in range(3):
        a = (amp&pos[j]).sum()
        result.loc[ i, molecular[j]+'_amp_freq'] = a/pos_num[j]
        b = pos_num[j]-a
        c = amp_num-a
        d = total_num-amp_num-pos_num[j]+a
        test_tab = [[a,b],[c,d]]
        if min(a,b,c,d) <= 5:
            res = stats.fisher_exact(test_tab)
            result.loc[i,molecular[j]+'_amp_p'] = res.pvalue
        else:
            res = stats.chi2_contingency(test_tab)
            result.loc[i,molecular[j]+'_amp_p'] = res.pvalue
        a = (de&pos[j]).sum()
        result.loc[ i, molecular[j]+'_del_freq'] = a/pos_num[j]
        b = pos_num[j]-a
        c = de_num-a
        d = total_num-de_num-pos_num[j]+a
        test_tab = [[a,b],[c,d]]
        if min(a,b,c,d) <= 5:
            res = stats.fisher_exact(test_tab)
            result.loc[i,molecular[j]+'_del_p'] = res.pvalue
        else:
            res = stats.chi2_contingency(test_tab)
            result.loc[i,molecular[j]+'_del_p'] = res.pvalue
print(result)
result.to_pickle('data/cna_subtype_freq_stat.pkl')
#'''
#correct for FDR
'''
from scipy import stats
cna = pd.read_pickle('data/cna_all_na_filter.pkl')
sample = pd.read_pickle('data/clinical_info_0623.pkl')
cna = cna.T
cna.index = cna.index.str.replace('.','-')
sample = sample.reindex(index = cna.index)
sample = sample.loc[ ~sample.loc[:,'molecular'].isna(),:]
cna = cna.reindex(index = sample.index)
result = pd.read_pickle('data/cna_subtype_freq_stat.pkl')
cna_anno = ['amp','del']
molecular =  ['ER+','HER2+','TNBC']
import statsmodels.stats.multitest as stat
for i in ['ER+','HER2+','TNBC']:
    for j in cna_anno:
        tmp = result.loc[:,i+'_'+j+'_p']
        result.loc[:,i+'_'+j+'_FDR'] = np.nan
        reject, correct, _ , _ = stat.multipletests(tmp, method = 'fdr_bh')
        result.loc[:,i+'_'+j+'_FDR'] = correct
print(result)
result.to_pickle('data/cna_subtype_freq_stat_FDR.pkl')
#'''

'''
from scipy import stats

cnaall = pd.read_csv('R_data/CNA_freq.csv', index_col = 0)
cnaall.columns = ['All_amp_freq','All_del_freq']
cna = pd.read_pickle('data/cna_all_na_filter.pkl')
sample = pd.read_pickle('data/clinical_info_0623.pkl')
cna = cna.T
cna.index = cna.index.str.replace('.','-')
sample = sample.reindex(index = cna.index)
sample = sample.loc[ ~sample.loc[:,'molecular'].isna(),:]
cna = cna.reindex(index = sample.index)
result = pd.read_pickle('data/cna_subtype_freq_stat.pkl')
cna_anno = ['amp','del']
molecular =  ['ER+','HER2+','TNBC']
result = pd.read_pickle('data/cna_subtype_freq_stat_FDR.pkl')
result.index.name = 'gene'
result = pd.concat([result, cnaall], axis = 1)
print(result)
result.to_csv('R_data/CNA_freq_subtype.csv')
for i in molecular:
    for j in cna_anno:
        tmp = result.loc[:,i+'_'+j+'_FDR']
        print(tmp.loc[tmp<0.0001])
        print(result.loc[ tmp < 0.0001, i+'_'+j+'_freq'])
        #tmp = result.loc[:,i+'_'+j+'_freq']
        #print(result.loc[ tmp > 0.2, i+'_'+j+'_FDR'])
#'''

#membrane protein database
'''
opm = pd.read_csv('data/proteins-2025-06-27.csv')
print(opm)
print(opm.loc[:,'membrane_id'].value_counts())
opm = opm.loc[ opm.loc[:,'species_name_cache'] == 'Homo sapiens',:]
opm = opm.loc[ opm.loc[:,'membrane_id'] == 4,:]
print(opm)
import requests

def pdb_to_uniprot(pdb_id):
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    r = requests.get(url)
    if r.status_code == 200:
        data = r.json()
        uniprot_ids = list(data[pdb_id.lower()]["UniProt"].keys())
        return uniprot_ids
    return []
opm.loc[:,'uniport'] = np.nan
# Example
uni = []
for i in opm.index:
    print(i)
    tmp = opm.loc[i,'pdbid'].replace("=",'').replace("\"",'')
    tmp = pdb_to_uniprot(tmp) 
    print(tmp)
    if len(tmp) == 0:
        continue
    uni.extend(tmp)
    if len(tmp) > 1:
        tmp = ' '.join(tmp)
    opm.loc[i,'uniport'] = tmp
print(opm)
opm.to_pickle('data/human_opm.pkl')
print(uni)
save_obj(uni, 'data/uniprot_list.pkl')
#'''

'''uni = load_obj('data/uniprot_list.pkl')
uni = pd.DataFrame(uni, columns = ['gene'])
uni.to_csv('data/uniprot_list.tsv', sep = '\t')
#'''

'''
freq = pd.read_csv('R_data/CNA_freq_subtype_chromosome.tsv', index_col = 1, sep = '\t')
print(freq)
tmp = freq.loc[ freq.loc[:,'chromosome_name'] == '11',:]
print(tmp.loc[:,'All_amp_freq'].nlargest(50))
tmp = freq.loc[ freq.loc[:,'chromosome_name'] == '8',:]
print(tmp.loc[:,'All_amp_freq'].nlargest(50))
tmp = freq.loc[ freq.loc[:,'chromosome_name'] == '1',:]
print(tmp.loc[:,'All_amp_freq'].nlargest(50))
#'''

'''cor = pd.read_pickle('data/cor_cna_mrna.pkl')
print(cor.loc['CD34',:])
#'''

'''
freq = pd.read_csv('R_data/CNA_freq_subtype_chromosome.tsv', index_col = 1, sep = '\t')
print(freq)
cor = pd.read_pickle('data/cor_cna_mrna.pkl')
cor = cor.reindex(index = freq.index)
freq = pd.concat([freq,cor], axis = 1)
print(freq.loc['PVRL2',:])
molecular = ['All','ER+','HER2+','TNBC']
#thres_freq = freq.loc['ERBB2','All_amp_freq']
thres_freq = 0.15
thres_cor = 0.2
freq.loc[:,'sel'] = 0
for i in molecular:
    if i == 'All':
        sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) & (freq.loc[:,'spearman_r'] > thres_cor),:]
        print(sel)
        freq.loc[sel.index,'sel'] = 1
        print(freq.loc[:,'sel'].value_counts())
        freq.loc[:,i] = 0
        freq.loc[sel.index,i] = 1
    else:
        sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) & (freq.loc[:,'spearman_r'] > thres_cor) & (freq.loc[:,i+'_amp_FDR'] < 0.0001)&(freq.loc[:,i+'_amp_freq']>freq.loc[:,'All_amp_freq']) ,:]
        print(sel)
        freq.loc[sel.index,'sel'] = 1
        freq.loc[:,i] = 0
        freq.loc[:,i] = 0
        freq.loc[sel.index,i] = 1
freq = freq.loc[ freq.loc[:,'sel'] == 1,:]
print(freq)
print(freq.loc[:,'chromosome_name'].value_counts())
print(freq.loc[:,'ER+'].value_counts())
print(freq.loc[:,'HER2+'].value_counts())
print(freq.loc[:,'TNBC'].value_counts())
#'''

#combining other omics annotation
#freq = pd.read_csv('R_data/CNA_freq_subtype_chromosome.csv', index_col = 0)
'''freq = pd.read_csv('R_data/CNA_freq_subtype_chromosome.tsv', index_col = 1, sep = '\t')
print(freq)
cor = pd.read_pickle('data/cor_cna_mrna.pkl')
cor = cor.reindex(index = freq.index)
freq = pd.concat([freq,cor], axis = 1)
membrane = pd.read_csv('data/subcell_location_Plasma.tsv', index_col = 0, sep = '\t')
wilcox = pd.read_pickle('data/wilcoxon_cna2_vs_other_mrna.pkl')
wilcox = wilcox.loc[~wilcox.index.duplicated(),:]
freq = freq.loc[~freq.index.duplicated(),:]
print(wilcox)
freq = pd.concat([freq, wilcox], axis = 1)
#membrane = pd.read_csv('data/uniprot_gprofiler.csv', index_col = 0, sep = ',')
if True:
    tmp_1 = pd.read_csv('data/PDB/rcsb_pdb_0001-2500.csv', sep = ',')
    tmp_2 = pd.read_csv('data/PDB/rcsb_pdb_2501-5000.csv', sep = ',')
    tmp_3 = pd.read_csv('data/PDB/rcsb_pdb_5001-5383.csv', sep = ',')
    mem = [tmp_1, tmp_2, tmp_3]
    gene = []
    for i in mem:
        tmp = i.loc[:,['Polymer EntityData']]
        tmp = tmp.dropna( how = 'all')
        for j in tmp.index:
            ad = tmp.loc[j,'Polymer EntityData'].split(', ')
            gene.extend(ad)
    #print(gene)
#print(membrane)
#freq = freq.loc[freq.index.isin(membrane.loc[:,'name']),:]
#freq = freq.loc[freq.index.isin(membrane.index),:]
freq.loc[:,'memebrane_protein'] = 'No'
freq.loc[:,'plasma_protein'] = 'No'
freq.loc[freq.index.isin(gene),'membrane_protein'] = 'Yes'
freq.loc[freq.index.isin(membrane.index),'plasma_protein'] = 'Yes'
print(freq)
if False:
    molecular = ['All','ER+','HER2+','TNBC']
    #thres_freq = freq.loc['ERBB2','All_amp_freq']
    thres_freq = 0.15
    thres_cor = 0.2
    freq.loc[:,'sel'] = 0
    for i in molecular:
        if i == 'All':
            sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) & (freq.loc[:,'spearman_r'] > thres_cor),:]
            print(sel)
            freq.loc[sel.index,'sel'] = 1
            freq.loc[:,i] = 0
            freq.loc[sel.index,i] = 1
        else:
            sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) & (freq.loc[:,'spearman_r'] > thres_cor) & (freq.loc[:,i+'_amp_FDR'] < 0.0001)&(freq.loc[:,i+'_amp_freq']>freq.loc[:,'All_amp_freq']) ,:]
            print(sel)
            freq.loc[sel.index,'sel'] = 1
            freq.loc[:,i] = 0
            freq.loc[:,i] = 0
            freq.loc[sel.index,i] = 1
print(freq)
if True:
    freq.loc[:,'CommonEss'] = 0
    control = pd.read_csv('../Prioritization/data/23Q4/common_ess_2021.csv', sep = ',',  index_col = 0)
    freq.loc[ freq.index.isin(control.index),'CommonEss'] = 1
    control = pd.read_csv('../Prioritization/data/23Q4/CRISPRInferredCommonEssentials.csv', sep = ' ',  index_col = 0)
    freq.loc[ freq.index.isin(control.index),'CommonEss'] = 1
    control = pd.read_csv('../Prioritization/data/23Q4/AchillesCommonEssentialControls.csv', sep = ' ', index_col = 0)
    freq.loc[ freq.index.isin(control.index),'CommonEss'] = 1
if False:
    sample = pd.read_pickle('../Prioritization/data/sample.pkl')
    dep = pd.read_pickle('../Prioritization/data/ceres_dep.pkl')
    sample = sample.reindex( index = dep.index)
    sample = sample.loc[ sample.loc[:,'OncotreeLineage'] == 'Breast',:]
    sample = sample.loc[:, ['StrippedCellLineName','LegacyMolecularSubtype','LegacySubSubtype']]
    sample.loc['ACH-001065','LegacySubSubtype'] = 'HER2_pos'
    sample.loc['ACH-001419','LegacySubSubtype'] = 'HER2_pos'
    sample.loc['ACH-002179','LegacySubSubtype'] = 'HER2_pos'
    sample.loc['ACH-002399','LegacySubSubtype'] = 'HER2_pos'
    sample.loc['ACH-001820','LegacySubSubtype'] = 'TNBC'
    sample = sample.loc[ ~sample.loc[:,'LegacySubSubtype'].isna(),:]
    sample.loc[ sample.loc[:,'LegacySubSubtype'].str.contains('HER2pos'),'LegacySubSubtype'] = 'HER2_pos'
    sample.loc[ sample.loc[:,'LegacySubSubtype'] == 'ERneg_HER2neg','LegacySubSubtype'] = 'TNBC'
    sample.loc[ sample.loc[:,'LegacySubSubtype'] == 'ERpos_HER2neg','LegacySubSubtype'] = 'ER_pos'
    eff = pd.read_pickle('../Prioritization/data/ceres_dep.pkl')
    eff = eff.loc[:,eff.isna().sum() < len(eff.index) * 0.2]
    eff = eff.reindex(index = sample.index, columns = freq.index)
    print(eff)
    er = eff.loc[ sample.loc[:,'LegacySubSubtype'] == 'ER_pos',:]
    her = eff.loc[ sample.loc[:,'LegacySubSubtype'] == 'HER2_pos',:]
    tnbc = eff.loc[ sample.loc[:,'LegacySubSubtype'] == 'TNBC',:]
    freq.loc[:,'ER+_dep'] = er.mean()
    freq.loc[:,'HER2+_dep'] = her.mean()
    freq.loc[:,'TNBC_dep'] = tnbc.mean()
    freq.loc[:,'All_dep'] = eff.mean()
freq = freq.drop(columns = ['gene'])
print(freq)
freq.to_pickle('data/summarization_large_table.pkl') 
#'''

#clinical data summary
'''sample = pd.read_pickle('data/clinical_info_0623.pkl')
cna = pd.read_pickle('data/cna_all_na_filter.pkl')
cna = cna.T
cna.index = cna.index.str.replace('.','-')
sample = sample.reindex( index = cna.index)
print(sample)
sample = sample.drop( columns = ['Patient ID'])
for i in sample.columns:
    if i != 'Age' and i != 'TMB':
        sample.loc[ sample.loc[:,i].isna(), i] = 'NA'
    if i == 'HER2' or i == 'ER':
        sample.loc[ ~sample.loc[:,'HER2'].isin(['Positive','Negative']), i] = 'NA'
for i in sample.columns:
    print(sample.loc[:,i].value_counts())
    print(sample.loc[:,i].value_counts(normalize = True))
    if i == 'Age':
        print(sample.loc[:,i].mean())
        print(sample.loc[:,i].std())
        print(sample.loc[:,i].median())
        print(sample.loc[:,i].quantile(0.25))
        print(sample.loc[:,i].quantile(0.75))
    if i == 'TMB':
        print( (sample.loc[:,i]>=10).sum() )
        print( (sample.loc[:,i]<10).sum() )
        print( (sample.loc[:,i].isna()).sum() )
        print( (sample.loc[:,i]>=10).sum()/len(sample.index) )
        print( (sample.loc[:,i]<10).sum()/len(sample.index) )
        print( (sample.loc[:,i].isna()).sum()/len(sample.index) )
#'''

#GTEx reading
#extract gene name and remove HGNC
'''
import pandas as pd
from collections import defaultdict
import numpy as np
import pickle

# Input and output paths
input_file = "R_data/GTEx_Analysis_v10_RSEMv1.3.3_transcripts_expected_count.txt"
output_pickle = 'gene_expected_counts.pkl'

# Step 1: Read header to get sample names
with open(input_file, 'r') as f:
    header = f.readline().rstrip('\n').split('\t')
    transcript_col, gene_col, sample_names = header[0], header[1], header[2:]

num_samples = len(sample_names)
gene_counts = defaultdict(lambda: np.zeros(num_samples, dtype=np.float32))

# Step 2: Stream the file line by line and aggregate counts by gene_id
with open(input_file, 'r') as f:
    next(f)
    for line_num, line in enumerate(f, start=2):
        fields = line.rstrip('\n').split('\t')
        gene_id = fields[1]
        try:
            counts = np.array(fields[2:], dtype=np.float32)
            gene_counts[gene_id] += counts
        except Exception as e:
            print(f"Skipping line {line_num} due to error: {e}")
            continue

# Step 3: Convert to DataFrame
gene_df = pd.DataFrame.from_dict(gene_counts, orient='index', columns=sample_names)

print(gene_df)
gene_df.to_pickle('data/gtex_mrna.pkl')
#'''

#resolve naming issue
'''
exp = pd.read_pickle('data/gtex_mrna.pkl')
sample = pd.read_csv('/Users/hao-kuen/yale/amp/R_data/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt', sep = '\t')
sample = sample.loc[ sample.loc[:,'SAMPID'].isin(exp.columns),:]
sample.index = sample.loc[:,'SAMPID']
sample = sample.reindex(index = exp.columns)
#print(sample.loc[:,['SMTS','SMTSD']])
print(sample.loc[:,'SMTS'].value_counts())
sample.index = sample.index.str.replace('-','_')
exp.columns = exp.columns.str.replace('-','_')
sample.loc[:,'SMTS'] = sample.loc[:,'SMTS'].str.replace(' ','_')
sample = sample.loc[:,['SMTS']]
sample.columns = ['tissue']
sample.loc[:,'tumor'] = 'normal'
import mygene

# Assume you already have your gene-level counts loaded into gene_df
# Index: ENSG IDs, Columns: sample names
# Example: gene_df = pd.read_pickle('gene_expected_counts.pkl')
gene_df = exp

mg = mygene.MyGeneInfo()

# Strip version numbers if present (e.g., ENSG00000141510.12 → ENSG00000141510)
gene_df.index = gene_df.index.str.split('.').str[0]

# Query mygene for mapping to HGNC symbol
ensembl_ids = gene_df.index.unique().tolist()
info = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human', as_dataframe=True)
print(info)
print(info['notfound'])
# Drop entries with no HGNC symbol
info = info[ (~info['symbol'].isna())]

# Build mapping ENSG → HGNC
ensg_to_symbol = info['symbol'].to_dict()

# Filter and remap the DataFrame
mapped_genes = gene_df.index.isin(ensg_to_symbol.keys())
gene_df = gene_df[mapped_genes]
gene_df.index = gene_df.index.map(ensg_to_symbol)

# Sum duplicate HGNC symbols
gene_df = gene_df.groupby(gene_df.index).sum()

# Now gene_df is HGNC_symbol x samples
# Save it
print(gene_df)
gene_df.to_pickle("data/gtx_hgnc.pkl")
sample.to_pickle('data/gtx_sample.pkl')
#'''

#tcga rsem
'''
exp = pd.read_csv('R_data/brca_tcga_mRNA_rsem.tsv', sep = '\t', index_col = 0)
exp = exp.drop(columns = ['group','group_name'])
exp.columns = exp.columns.str.replace('.','_')
print(exp)
exp.to_pickle('data/brca_tcga_mRNA_rsem.pkl')
#'''

'''
exp = pd.read_pickle('data/brca_tcga_mRNA_rsem.pkl')
print(exp)
gdc = pd.read_pickle("data/gdc_brca_tcga_rsem_hgnc.pkl")
print(gdc)
#'''

#gdc brca mrna
'''exp = pd.read_csv('R_data/gdc_tcga_brca_rsem.tsv', sep = '\t', index_col = 0)
exp.columns = exp.columns.str.replace('-','_')
import mygene

# Assume you already have your gene-level counts loaded into gene_df
# Index: ENSG IDs, Columns: sample names
# Example: gene_df = pd.read_pickle('gene_expected_counts.pkl')
gene_df = exp

mg = mygene.MyGeneInfo()

# Strip version numbers if present (e.g., ENSG00000141510.12 → ENSG00000141510)
gene_df.index = gene_df.index.str.split('.').str[0]

# Query mygene for mapping to HGNC symbol
ensembl_ids = gene_df.index.unique().tolist()
info = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human', as_dataframe=True)
print(info)
print(info['notfound'])
# Drop entries with no HGNC symbol
info = info[ (~info['symbol'].isna())]

# Build mapping ENSG → HGNC
ensg_to_symbol = info['symbol'].to_dict()

# Filter and remap the DataFrame
mapped_genes = gene_df.index.isin(ensg_to_symbol.keys())
gene_df = gene_df[mapped_genes]
gene_df.index = gene_df.index.map(ensg_to_symbol)

# Sum duplicate HGNC symbols
gene_df = gene_df.groupby(gene_df.index).sum()

# Now gene_df is HGNC_symbol x samples
# Save it
print(gene_df)
gene_df.to_pickle("data/gdc_brca_tcga_rsem_hgnc.pkl")
#'''

#combine gtx and tcga
'''exp = pd.read_csv('R_data/brca_tcga_mRNA_rsem.tsv', sep = '\t', index_col = 0)
exp = exp.drop(columns = ['group','group_name'])
exp.columns = exp.columns.str.replace('.','_')
print(exp)
sample = pd.read_pickle('data/gtx_sample.pkl')
gtx = pd.read_pickle('data/gtx_hgnc.pkl')
gtx = gtx.reindex(index = exp.index)
gtx = gtx.dropna( how ='all')
exp = exp.reindex(gtx.index)
exp = pd.read_pickle("data/gdc_brca_tcga_rsem_hgnc.pkl")
exp = exp.reindex(index = gtx.index)
exp = exp.dropna( how ='all')
gtx = gtx.reindex(exp.index)
tcga = pd.DataFrame(np.nan, index = exp.columns, columns = ['tissue','tumor'])
tcga.loc[:,'tissue'] = 'Breast_Cancer'
tcga.loc[:,'tumor'] = 'cancer'
sample = pd.concat([sample,tcga], axis = 0)
exp = pd.concat([gtx,exp], axis = 1)
exp = exp.astype('int')
print(exp)
print(sample)
exp.index.name = 'gene'
exp.to_csv('R_data/normal_tissue_gdc.tsv', sep = '\t')
sample.index.name = 'sample_id'
sample.to_csv('R_data/normal_tissue_info_gdc.tsv', sep = '\t')
#'''

'''exp = pd.read_csv('R_data/brca_tcga_mRNA_rsem.tsv', sep = '\t', index_col = 0)
exp = exp.drop(columns = ['group','group_name'])
exp.columns = exp.columns.str.replace('.','_')
print(exp)
sample = pd.read_pickle('data/gtx_sample.pkl')
gtx = pd.read_pickle('data/gtx_hgnc.pkl')
gtx = gtx.reindex(index = exp.index)
gtx = gtx.dropna( how ='all')
exp = exp.reindex(gtx.index)
tcga = pd.DataFrame(np.nan, index = exp.columns, columns = ['tissue','tumor'])
tcga.loc[:,'tissue'] = 'Breast_Cancer'
tcga.loc[:,'tumor'] = 'cancer'
sample = pd.concat([sample,tcga], axis = 0)
exp = pd.concat([gtx,exp], axis = 1)
exp = exp.astype('int')
print(exp)
print(sample)
exp.index.name = 'gene'
exp.to_csv('R_data/normal_tissue.tsv', sep = '\t')
sample.index.name = 'sample_id'
sample.to_csv('R_data/normal_tissue_info.tsv', sep = '\t')
#'''

#select normal tissue to compare (more vital organs)
'''
tissue_sel = ['Lung','Bladder','Blood_Vessel','Skin','Esophagus','Small_Intestine','Thyroid','Kidney','Muscle','Spleen','Liver','Stomach','Pituitary','Brain','Pancreas','Blood','Nerve','Prostate','Colon','Adrenal_gland','Heart']
t = pd.DataFrame(np.nan, index = sorted(tissue_sel), columns = ['tmp'])
t.to_csv('data/tissue_selected.csv')
sorted(tissue_sel)
#'''

#retreive differential expression results
'''
import os
import pandas as pd

sample = pd.read_csv('R_data/normal_tissue_info_gdc.tsv', sep = '\t', index_col = 0)
# === Parameters ===
deseq_dir = 'R_data/deseq'
threshold = 1.0  # Modify this threshold as needed

# === Load all DE files ===
combined_df_list = []
output = []
tissue_sel = ['Lung','Bladder','Blood_Vessel','Skin','Esophagus','Small_Intestine','Thyroid','Kidney','Muscle','Spleen','Liver','Stomach','Pituitary','Brain','Pancreas','Blood','Nerve','Prostate','Colon','Adrenal_gland','Heart']
for filename in os.listdir(deseq_dir):
    if filename.startswith("DE_") and filename.endswith("_vs_Breast_Cancer_gdc.csv"):
        tissue_name = filename.replace("DE_", "").replace("_vs_Breast_Cancer_gdc.csv", "")
        if tissue_name in tissue_sel:
            file_path = os.path.join(deseq_dir, filename)
            
            df = pd.read_csv(file_path, index_col = 0)
            df = df[['log2FoldChange', 'padj']].copy()
            df['gene'] = df.index  # assume gene names are row index
            df['log2FoldChange'] = -df['log2FoldChange']  # Flip the sign
            df['tissue'] = tissue_name
            
            combined_df_list.append(df)
            df_tmp = df.loc[:,['log2FoldChange','padj']]
            df_tmp.columns = [i+'_'+tissue_name for i in df_tmp.columns]
            output.append(df_tmp)

# Combine all data into one dataframe
combined_df = pd.concat(combined_df_list, ignore_index=True)
output = pd.concat(output, axis = 1)
print(output)

# === Pivot to per-gene summary ===
# Average log2FC per gene
avg_fc = combined_df.groupby('gene')['log2FoldChange'].mean().reset_index(name='avg_log2FoldChange')

# Count of tissues with significant FC above threshold
def count_sig_fc(df, threshold):
    return ((df['padj'] < 0.05) & (abs(df['log2FoldChange']) > threshold)).groupby(df['gene']).sum().rename('n_tissues_sig')

sig_counts = count_sig_fc(combined_df, threshold).reset_index()

# Merge summary
summary_df = avg_fc.merge(sig_counts, on='gene', how='left')
summary_df['n_tissues_sig'] = summary_df['n_tissues_sig'].fillna(0).astype(int)
summary_df.index = summary_df.loc[:,'gene']
summary_df = summary_df.drop(columns = ['gene'])
print(summary_df)
summary_df = pd.concat([summary_df, output], axis = 1)

# === Output ===
summary_df = summary_df.sort_values(by='avg_log2FoldChange', ascending=False)
print(summary_df)
summary_df.to_csv('data/DE_summary_by_gene_gdc.csv')

print("Summary saved to 'DE_summary_by_gene_gdc.csv'")
#'''

'''
import os
import pandas as pd

sample = pd.read_csv('R_data/normal_tissue_info_gdc.tsv', sep = '\t', index_col = 0)
# === Parameters ===
deseq_dir = 'R_data/deseq'
threshold = 1.0  # Modify this threshold as needed

# === Load all DE files ===
combined_df_list = []
output_f = []
output_p = []
tissue_sel = ['Lung','Bladder','Blood_Vessel','Skin','Esophagus','Small_Intestine','Thyroid','Kidney','Muscle','Spleen','Liver','Stomach','Pituitary','Brain','Pancreas','Blood','Nerve','Prostate','Colon','Adrenal_gland','Heart']
for filename in os.listdir(deseq_dir):
    if filename.startswith("DE_") and filename.endswith("_vs_Breast_Cancer_gdc.csv"):
        tissue_name = filename.replace("DE_", "").replace("_vs_Breast_Cancer_gdc.csv", "")
        if tissue_name in tissue_sel:
            file_path = os.path.join(deseq_dir, filename)
            
            df = pd.read_csv(file_path, index_col = 0)
            df = df[['log2FoldChange', 'padj']].copy()
            df['gene'] = df.index  # assume gene names are row index
            df['log2FoldChange'] = -df['log2FoldChange']  # Flip the sign
            df['tissue'] = tissue_name
            
            combined_df_list.append(df)
            df_tmp = df.loc[:,['log2FoldChange']]
            df_tmp.columns = [tissue_name]
            output_f.append(df_tmp)
            df_tmp = df.loc[:,['padj']]
            df_tmp.columns = [tissue_name]
            output_p.append(df_tmp)

# Combine all data into one dataframe
output_f = pd.concat(output_f, axis = 1)
output_p = pd.concat(output_p, axis = 1)
print(output_f)
print(output_p)
output_f.to_pickle('data/DE_result_gdc_f.pkl')
output_p.to_pickle('data/DE_result_gdc_p.pkl')
#'''

#resolving name issue for the previous summarization table
'''
import pandas as pd
from gprofiler import GProfiler
# Load your data
df = pd.read_pickle('data/summarization_large_table.pkl')
print(df)
# Initialize g:Profiler client
gp = GProfiler(return_dataframe=True)

# Step 1: Convert HGNC -> ENSG (and latest HGNC)
mapping = gp.convert(
    organism="hsapiens",
    query=list(df.index),
    target_namespace="HGNC"   # we want Ensembl Gene IDs
)
print(mapping)
save_obj(mapping, 'data/hgnc.pkl')
#'''

'''
engc = load_obj('data/engc.pkl')
print(engc)
tmp = engc.loc[engc.loc[:,'converted'] == 'None',:]
tmp.index = tmp.loc[:,'incoming']
tmp.to_csv('data/engc_none.txt', sep = '\t')
#'''

#merge the table with differential expression results
'''exp = pd.read_pickle('data/brca_tcga_mRNA_rsem.pkl')
exp = np.log2(exp+1)
print(exp)
exp.loc[:,'mean'] = exp.T.mean()
print(exp.loc[:,'mean'].quantile(0.5))
exp.loc[:,'mean_quantile'] = pd.qcut(exp.loc[:,'mean'], q = 4, labels = False) + 1
print(exp)
freq = pd.read_pickle('data/summarization_large_table.pkl')
freq = pd.concat([freq, exp.loc[:,['mean','mean_quantile']]], axis = 1)
print(freq)
#freq.to_pickle('data/summarization_large_table_quantile.pkl')
#'''


'''
hgnc = pd.read_csv('data/hgnc-symbol-check.csv', index_col = 0)
print(hgnc)
df = pd.read_pickle('data/summarization_large_table_quantile.pkl')
dup = hgnc.index.duplicated()
dup = hgnc.index[dup]
for i in hgnc.index:
    g = hgnc.loc[i,'Approved symbol']
    index = df.index.tolist()
    if (i in dup):
        if i in index:
            for gene in g:
                print(gene)
                if gene in index:
                    df.loc[gene,df.loc[ gene,:].isna()] = df.loc[i,df.loc[gene,:].isna()]
                else:
                    df.loc[gene,:] = df.loc[i,:]
            df = df.drop( index = [i])
    else:
        if g in index:
            df.loc[g,df.loc[ g,:].isna()] = df.loc[i,df.loc[g,:].isna()]
        else:
            df.loc[g,:] = df.loc[i,:]
        df = df.drop( index = [i])
print(df)
print(df.loc['NECTIN4',:])
df.to_pickle('data/summarization_large_table_hgnc_quantile.pkl')
#'''

#Utilizing Gemini to do literature review on membrane protein information
'''
import pandas as pd
import requests
import time
from openai import OpenAI  # Gemini is OpenAI-compatible

# 1. Setup API Configuration
# Get your key at: https://aistudio.google.com/app/apikey
GEMINI_API_KEY = "YOUR_GEMINI_API_KEY"
client = OpenAI(
    api_key="AIzaSyDrwIwApSnCfMgBF1DdcSSW4uXoVpKKOlo",
    base_url="https://generativelanguage.googleapis.com/v1beta/openai/"
)

def query_gemini_for_gene(gene_symbol):
    """Queries Gemini to perform a literature-based screening for a gene."""
    prompt = f"""
    Search the published scientific literature and public protein data bases and gene annotations 
    to determine if the gene '{gene_symbol}' is a transmembrane protein localized to the 
    cell membrane or not. Additionally, determine if it has the ability for internalization 
    upon ligand binding. Provide a Yes/No for both and a brief supporting citation if possible.
    """
    prompt = f"""
    Search the published scientific literature and public protein data bases and gene annotations 
    to determine if the gene '{gene_symbol}' is a membrane protein localized to the 
    cell membrane or not. Additionally, determine if it has the ability for internalization 
    upon ligand binding. Provide a Yes/No for both and a brief supporting citation if possible.
    """

    try:
        response = client.chat.completions.create(
            model="gemini-3-flash-preview", # Use flash for high-volume screening
            messages=[{"role": "user", "content": prompt}]
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Error: {str(e)}"

def get_uniprot_data(gene_symbol):
    """Retrieves validated subcellular location from UniProt API."""
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene_symbol} AND reviewed:yes&format=json"
    try:
        res = requests.get(url)
        data = res.json()
        if data['results']:
            # Extract Subcellular Location comments
            comments = data['results'][0].get('comments', [])
            locations = [c for c in comments if c.get('commentType') == 'SUBCELLULAR_LOCATION']
            return str(locations)
        return "Not found in reviewed UniProt"
    except:
        return "API Error"

# 2. Execution Loop
# Load your list of 228 genes
gene_list = ["ERBB2", "MAL2", "SUSD4"] # Replace with your full list
tmp = pd.read_csv('data/top_sel_gene.txt', sep = '\t', index_col = 0)
tmp = tmp.dropna( how = 'all')
print(tmp)

results = []

print(f"Starting review for {len(gene_list)} genes...")

for gene in tmp.index:
    print(f"Processing {gene}...")
    
    # Get biological evidence
    #uniprot_info = get_uniprot_data(gene)
    
    # Get LLM Literature Synthesis
    lit_review = query_gemini_for_gene(gene)
    
    results.append({
        "Gene": gene,
        #"UniProt_Evidence": uniprot_info,
        "LLM_Literature_Review": lit_review
    })
    
    # Rate limiting (adjust based on your API tier)
    time.sleep(1) 

# 3. Save Results
df = pd.DataFrame(results)
df.to_csv("data/Gene_Literature_Review_Results_membrane.txt",sep = '\t', index=False)
print("Review complete. Results saved to Gene_Literature_Review_Results.xlsx")
#'''

'''
df = pd.read_csv('Untitled.txt', sep = ',', index_col = 0)
print(df)
#'''

#retreive amplication targets
'''
import pandas as pd
import numpy as np

# Load the dataset
# Assuming the file is named 'Amp_1702.txt' and is tab-delimited
df = pd.read_csv('Amp_1707.txt', sep='\t', index_col = 0)

# 1. Select genes where 'sel' is 1
#df_filtered = df[df['sel'] == 1].copy()
df_filtered = df.copy()

# 2. Fill 'Finalist' NA values with 0
df_filtered['Finalist'] = df_filtered['Finalist'].fillna(0)

# 3. Create annotation column based on subtype prioritization
# Criteria: If only one subtype has exp_prioritize == 1, name it. 
# If more than one has 1, label as 'more'.
def annotate_subtype(row):
    subtypes = {
        'ER+': row['ER+_exp_prioritize'],
        'HER2+': row['HER2+_exp_prioritize'],
        'TNBC': row['TNBC_exp_prioritize']
    }
    
    # Identify which subtypes are active (value == 1)
    active_subtypes = [k for k, v in subtypes.items() if v == 1]
    
    if len(active_subtypes) == 1:
        return active_subtypes[0]
    elif len(active_subtypes) > 1:
        return 'shared'
    elif row['All_exp_prioritize'] == 1:
        return 'shared' # Or any default if no subtype is prioritized
    else:
        return 'None'

df_filtered['annotation'] = df_filtered.apply(annotate_subtype, axis=1)

# 4. Select columns needed for R
# Columns based on file structure: 'mean', 'avg_log2FoldChange', 'annotation', 'Finalist', 'gene_symbol'

df_filtered.to_csv('data/amp1707.txt', sep = '\t')

output_cols = ['mean', 'avg_log2FoldChange', 'annotation', 'Finalist']
df_r_input = df_filtered[output_cols]
print(df_r_input)

# 5. Output to a txt file for R
df_r_input.to_csv('R_data/gene_data_for_r_1707.txt', sep='\t')

print("Python processing complete. Output saved to 'gene_data_for_r.txt'.")
#'''

#GO plot, calculate subtype specific amplification
'''
from gprofiler import GProfiler
gp = GProfiler( return_dataframe = True)
def identify_recurrent_amplicons(filepath, subtype, frequency_threshold=0.15, gap_limit=1000000):
    """
    Identifies recurrent amplicons from GISTIC gene-level data.
    
    Parameters:
    - filepath: Path to 'all_thresholded.by_genes.txt'
    - frequency_threshold: Minimum % of samples required to consider an amplicon 'recurrent' (default 5%)
    - gap_limit: Max distance (bp) between genes to be considered the same amplicon.
    """
    
    # 1. Load Data
    # GISTIC files usually have 'Gene Symbol' and 'Locus ID' as the first two columns
    df = filepath.copy()
    
    # 2. Calculate Amplification Frequency
    # We look specifically for the value '2'
    
    # 3. Filter for Recurrent Genes
    df.loc[:,'gene_symbol'] = df.index
    recurrent_genes = df
    if subtype != 'All':
        if True:
            if subtype == 'ER+':
                recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'HER2+_freq_prioritize']==1)|(recurrent_genes.loc[:,'TNBC_freq_prioritize']==1))]
            if subtype == 'HER2+':
                recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'ER+_freq_prioritize']==1)|(recurrent_genes.loc[:,'TNBC_freq_prioritize']==1))]
            if subtype == 'TNBC':
                recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'HER2+_freq_prioritize']==1)|(recurrent_genes.loc[:,'ER+_freq_prioritize']==1))]
    recurrent_genes = recurrent_genes[recurrent_genes.loc[:,subtype+'_freq_prioritize'] == 1]
    subtype_col = f"{subtype}_amp_freq"    
    print(recurrent_genes)
    if recurrent_genes.empty:
        print("No genes met the frequency threshold.")
        return None

    # 4. Parse Genomic Coordinates
    # Sort by genomic position
    recurrent_genes = recurrent_genes.sort_values(['chromosome_name', 'start_position'])

    # 2. Calculate the gap (Current Start minus Previous End)
    # .shift() moves the 'End' column down by one row
    gap_distances = recurrent_genes['start_position'] - recurrent_genes['end_position'].shift()

    # 3. Define the break: 
    # New amplicon if chromosome changes OR if gap > 1MB
    recurrent_genes['Amplicon_ID'] = (
        (recurrent_genes['chromosome_name'] != recurrent_genes['chromosome_name'].shift()) | 
        (gap_distances > 1000000)
    ).cumsum()
    
    # 6. Identify Peak Genes (highest frequency in each cluster)
    results = []
    for amp_id, group in recurrent_genes.groupby('Amplicon_ID'):
        peak_gene = group.loc[group[subtype+'_amp_freq'].idxmax()]
        
        results.append({
            #'Amplicon_ID': amp_id,
            'Chrom': group['chromosome_name'].iloc[0],
            #'Start_Pos': group['start_position'].min(),
            #'End_Pos': group['end_position'].max(),
            'Peak_Gene': peak_gene['gene_symbol'],
            'Peak_Frequency': peak_gene[subtype_col],
            'Gene_Count': len(group),
            'All_Genes': ", ".join(group['gene_symbol'].tolist())
        })
        if False:
            print(", ".join(group['gene_symbol'].tolist()))
            result = gp.profile( organism = 'hsapiens', query = group['gene_symbol'].tolist(), sources = ['GO:MF','GO:CC','GO:BP'],  no_evidences = False, no_iea = True)
            leng = 0
            while len(result.index) > leng + 50:
                print(result.loc[ leng:leng+50,['name','p_value','intersections']])
                leng = leng + 50
            print(result.loc[ leng:leng+50,['name','p_value','intersections']])

    result = gp.profile( organism = 'hsapiens', query = recurrent_genes['gene_symbol'].tolist(), sources = ['GO:MF','GO:CC','GO:BP'],  no_evidences = False, no_iea = True)
    leng = 0
    while len(result.index) > leng + 50:
        print(result.loc[ leng:leng+50,['name','p_value','intersections']])
        leng = leng + 50
    print(result.loc[ leng:leng+50,['name','p_value','intersections']])
    df = result.copy()
    # keep top N terms
    TOP_N = 8   # reviewer sweet spot: 5–10
    #df = df.iloc[ [1,5,6,12,18,20],:]
    df = df.head(TOP_N)
    # add -log10(p)
    df['log10_p'] = -np.log10(df['p_value'])
    fig = plt.figure(figsize=(7, 5))
    import textwrap
    def wrap_go_terms(term, width=35):
        return "\n".join(textwrap.wrap(term, width=width))
    df['name'] = df['name'].apply(wrap_go_terms)
    sns.scatterplot(
        data=df,
        x='log10_p',
        y='name',
        size='intersection_size',
        #sizes=(50, 350),
        color='firebrick',
        edgecolor='black',
        alpha=0.85
    )
    plt.ylabel('')
    plt.legend(title = 'number of genes')
    plt.title('GO:BP Enrichment (g:Profiler)', fontsize=12)
    plt.grid(axis='x', linestyle='--', alpha=0.3)
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()
    return pd.DataFrame(results)

# Usage
freq = pd.read_pickle('data/summarization_large_table_hgnc_quantile.pkl')
freq.loc[:,'gene_symbol'] = freq.index
freq.loc[freq.loc[:,'membrane_protein'] != 'Yes', 'membrane_protein'] = ''
freq.loc[freq.loc[:,'plasma_protein'] != 'Yes', 'plasma_protein'] = ''
freq.loc[:,'membrane_protein_combined'] = ''
freq.loc[ (freq.loc[:,'membrane_protein'] == 'Yes'),'membrane_protein_combined'] = 'Yes'
freq.loc[ (freq.loc[:,'plasma_protein'] == 'Yes'),'membrane_protein_combined'] = 'Yes'
freq = freq.drop(columns = ['memebrane_protein'])
freq = freq.loc[:, ~freq.columns.str.contains('del_')]
freq = freq.loc[:, ~freq.columns.str.contains('amp_p')]
de = pd.read_csv('data/DE_summary_by_gene_gdc.csv', index_col = 0)
freq = pd.concat([freq,de], axis = 1)
freq.to_csv('data/sum_gdc.csv')
#print(freq.loc[:,'avg_log2FoldChange'].mean(), freq.loc[:,'avg_log2FoldChange'].std(), freq.loc[:,'avg_log2FoldChange'].sem())
freq = freq.loc[:,['gene_symbol','ER+_amp_freq', 'HER2+_amp_freq', 'TNBC_amp_freq', 'ER+_amp_FDR', 'HER2+_amp_FDR', 'TNBC_amp_FDR', 'All_amp_freq', 'chromosome_name', 'start_position', 'end_position', 'band', 'chr', 'spearman_r', 'spearman_p', 'wilcoxon_stat', 'wilcoxon_fdr', 'plasma_protein', 'membrane_protein', 'membrane_protein_combined', 'CommonEss', 'avg_log2FoldChange', 'n_tissues_sig','mean','mean_quantile']]
freq = freq.rename(columns = {'plasma_protein':'plasma_membrane_protein'})
if False:
    for i in ['ER+','HER2+','TNBC','All']:
        print(i)
        results = identify_recurrent_amplicons(freq,i)
        print(results)

import matplotlib.backends.backend_pdf
import seaborn as sns
pdf = matplotlib.backends.backend_pdf.PdfPages('figure/go_plot.pdf')
if True:
    molecular = ['All','ER+','HER2+','TNBC']
    thres_freq = 0.15
    thres_cor = 0.2
    thres_wilcox = 0.001
    freq.loc[:,'sel'] = 0
    for i in molecular:
        if i == 'All':
            print(i)
            sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq),:]
            freq.loc[:,i+'_freq_prioritize'] = 0
            freq.loc[sel.index,i+'_freq_prioritize'] = 1
            print(freq.loc[:,i+'_freq_prioritize'].sum())
            #sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) & (freq.loc[:,'spearman_r'] > thres_cor),:]
            sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) & (freq.loc[:,'wilcoxon_stat'] > 0) &  (freq.loc[:,'wilcoxon_fdr'] < thres_wilcox)&(freq.loc[:,'mean_quantile']>=3),:]
            freq.loc[sel.index,'sel'] = 1
            freq.loc[:,i+'_exp_prioritize'] = 0
            freq.loc[sel.index,i+'_exp_prioritize'] = 1
            print( (freq.loc[:,i+'_exp_prioritize'] == 1).sum())
        else:
            print(i)
            sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) & (freq.loc[:,i+'_amp_FDR'] < 0.0001)&(freq.loc[:,i+'_amp_freq']>freq.loc[:,'All_amp_freq']) ,:]
            freq.loc[:,i+'_freq_prioritize'] = 0
            freq.loc[sel.index,i+'_freq_prioritize'] = 1
            print(freq.loc[:,i+'_freq_prioritize'].sum())
            sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) &(freq.loc[:,'wilcoxon_stat'] > 0) &  (freq.loc[:,'wilcoxon_fdr'] < thres_wilcox) & (freq.loc[:,i+'_amp_FDR'] < 0.0001)&(freq.loc[:,i+'_amp_freq']>freq.loc[:,'All_amp_freq'])&(freq.loc[:,'mean_quantile']>=3),:]
            freq.loc[sel.index,'sel'] = 1
            freq.loc[:,i+'_exp_prioritize'] = 0
            freq.loc[sel.index,i+'_exp_prioritize'] = 1
            print( (freq.loc[:,i+'_exp_prioritize'] == 1).sum())
    print(freq.loc[(freq.loc[:,freq.columns.str.contains('freq_prioritize')].sum(axis = 1) >= 1),:])
    for i in molecular:
        print(i)
        results = identify_recurrent_amplicons(freq,i)
        print(results)
        if i != 'All':
            recurrent_genes = freq.copy()
            if i == 'ER+':
                recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'HER2+_exp_prioritize']==1)|(recurrent_genes.loc[:,'TNBC_exp_prioritize']==1))]
            if i == 'HER2+':
                recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'ER+_exp_prioritize']==1)|(recurrent_genes.loc[:,'TNBC_exp_prioritize']==1))]
            if i == 'TNBC':
                recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'HER2+_exp_prioritize']==1)|(recurrent_genes.loc[:,'ER+_exp_prioritize']==1))]
            recurrent_genes = recurrent_genes[recurrent_genes.loc[:,i+'_exp_prioritize'] == 1]
            print(recurrent_genes)

freq = freq.loc[(freq.loc[:,freq.columns.str.contains('freq_prioritize')].sum(axis = 1) >= 1),:]
ref = pd.read_csv('data/lit_rev_ref.txt', sep = '\t', index_col = 0)
print(ref)
ref.index = ref.loc[:,'gene']
freq = pd.concat([ref, freq], axis = 1)
print(freq)
pdf.close()
#freq.to_csv('data/new_final.txt', sep = '\t')
if False:
    cna_df = pd.read_pickle('data/cna_all_na_filter.pkl')
    mrna_df = pd.read_pickle('data/mrna_all_na_filter.pkl')
    print(cna_df)
    cna_df = cna_df.T
    print(mrna_df)
    mrna_df = mrna_df.T
    print(freq)
    gene_list = freq.index.tolist()
    gene_list.append('CD34')

    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from scipy.stats import spearmanr

    def plot_cna_vs_mrna_to_pdf(cna_df, mrna_df, gene_list, pdf_path='figure/cna_mrna_plots.pdf'):
        # Define CNA categories and palette with string keys
        cna_categories = ['-2', '-1', '0', '1', '2']
        cna_palette = {
            '-2': '#d73027',
            '-1': '#fc8d59',
             '0': '#ffffbf',
             '1': '#91bfdb',
             '2': '#4575b4'
        }

        # Ensure same patients
        common_patients = cna_df.index.intersection(mrna_df.index)
        cna_df = cna_df.loc[common_patients]
        mrna_df = mrna_df.loc[common_patients]

        with PdfPages(pdf_path) as pdf:
            for gene in gene_list:
                if gene not in cna_df.columns or gene not in mrna_df.columns:
                    print(f"Skipping {gene} (not found in both datasets)")
                    continue

                df = pd.DataFrame({
                    'CNA': cna_df[gene],
                    'mRNA': mrna_df[gene]
                }).dropna()

                # Convert CNA to string categorical to match palette keys
                df['CNA'] = df['CNA'].astype(int).astype(str)
                df['CNA'] = pd.Categorical(df['CNA'], categories=cna_categories, ordered=True)
                try:
                    rho, pval = spearmanr(df['CNA'], df['mRNA'])
                    stat_text = f"Spearman r = {rho:.2f}\np = {pval:.1e}"
                except Exception as e:
                    stat_text = "Spearman r = NA\np = NA"
                    print(f"⚠️  Spearman failed for {gene}: {e}")

                plt.figure(figsize=(6, 4))
                sns.boxplot(x='CNA', y='mRNA', data=df, order=cna_categories,
                            palette=cna_palette, showfliers=False)
                sns.stripplot(x='CNA', y='mRNA', data=df, order=cna_categories,
                              color='black', size=3, jitter=True, alpha=0.6)
                plt.text(x=-0.4, y=9.5, s=stat_text,
                         ha='left', va='top', fontsize=10,
                         bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.3'))

                plt.ylim(-4, 10)
                plt.title(f'{gene}: CNA vs mRNA expression')
                plt.xlabel('CNA value')
                plt.ylabel('mRNA expression (Z-score)')
                plt.tight_layout()

                pdf.savefig()
                plt.close()

        print(f"✅ Plots saved to: {pdf_path}")
    plot_cna_vs_mrna_to_pdf(cna_df, mrna_df, gene_list)
#'''


#add quantile selection
'''
from gprofiler import GProfiler
gp = GProfiler( return_dataframe = True)
def identify_recurrent_amplicons(filepath, subtype, frequency_threshold=0.15, gap_limit=1000000):
    """
    Identifies recurrent amplicons from GISTIC gene-level data.
    
    Parameters:
    - filepath: Path to 'all_thresholded.by_genes.txt'
    - frequency_threshold: Minimum % of samples required to consider an amplicon 'recurrent' (default 5%)
    - gap_limit: Max distance (bp) between genes to be considered the same amplicon.
    """
    
    # 1. Load Data
    # GISTIC files usually have 'Gene Symbol' and 'Locus ID' as the first two columns
    df = filepath.copy()
    
    # 2. Calculate Amplification Frequency
    # We look specifically for the value '2'
    
    # 3. Filter for Recurrent Genes
    df.loc[:,'gene_symbol'] = df.index
    recurrent_genes = df[df[subtype+'_amp_freq'] > frequency_threshold].copy()
    if subtype != 'All':
        if True:
            if subtype == 'ER+':
                recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'HER2+_amp_freq']>frequency_threshold)&(recurrent_genes.loc[:,'HER2+_amp_FDR']<0.0001))]
                recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'TNBC_amp_freq']>frequency_threshold)&(recurrent_genes.loc[:,'TNBC_amp_FDR']<0.0001))]
            if subtype == 'HER2+':
                recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'ER+_amp_freq']>frequency_threshold)&(recurrent_genes.loc[:,'ER+_amp_FDR']<0.0001))]
                recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'TNBC_amp_freq']>frequency_threshold)&(recurrent_genes.loc[:,'TNBC_amp_FDR']<0.0001))]
            if subtype == 'TNBC':
                recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'HER2+_amp_freq']>frequency_threshold)&(recurrent_genes.loc[:,'HER2+_amp_FDR']<0.0001))]
                recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'ER+_amp_freq']>frequency_threshold)&(recurrent_genes.loc[:,'ER+_amp_FDR']<0.0001))]
        recurrent_genes = recurrent_genes[recurrent_genes.loc[:,subtype+'_amp_FDR'] < 0.0001]
        recurrent_genes = recurrent_genes[recurrent_genes.loc[:,subtype+'_amp_freq'] > recurrent_genes.loc[:,'All_amp_freq']]
    if False:
        freq.loc[sel.index,'sel'] = 1
        freq.loc[:,i] = 0
        freq.loc[sel.index,i] = 1
    subtype_col = f"{subtype}_amp_freq"    
    print(len(recurrent_genes))
    if recurrent_genes.empty:
        print("No genes met the frequency threshold.")
        return None

    # 4. Parse Genomic Coordinates
    # Sort by genomic position
    recurrent_genes = recurrent_genes.sort_values(['chromosome_name', 'start_position'])

    # 2. Calculate the gap (Current Start minus Previous End)
    # .shift() moves the 'End' column down by one row
    gap_distances = recurrent_genes['start_position'] - recurrent_genes['end_position'].shift()

    # 3. Define the break: 
    # New amplicon if chromosome changes OR if gap > 1MB
    recurrent_genes['Amplicon_ID'] = (
        (recurrent_genes['chromosome_name'] != recurrent_genes['chromosome_name'].shift()) | 
        (gap_distances > 1000000)
    ).cumsum()
    
    # 6. Identify Peak Genes (highest frequency in each cluster)
    results = []
    for amp_id, group in recurrent_genes.groupby('Amplicon_ID'):
        peak_gene = group.loc[group[subtype+'_amp_freq'].idxmax()]
        
        results.append({
            #'Amplicon_ID': amp_id,
            'Chrom': group['chromosome_name'].iloc[0],
            #'Start_Pos': group['start_position'].min(),
            #'End_Pos': group['end_position'].max(),
            'Peak_Gene': peak_gene['gene_symbol'],
            'Peak_Frequency': peak_gene[subtype_col],
            'Gene_Count': len(group),
            'All_Genes': ", ".join(group['gene_symbol'].tolist())
        })
        if False:
            print(", ".join(group['gene_symbol'].tolist()))
            result = gp.profile( organism = 'hsapiens', query = group['gene_symbol'].tolist(), sources = ['GO:MF','GO:CC','GO:BP'],  no_evidences = False, no_iea = True)
            leng = 0
            while len(result.index) > leng + 50:
                print(result.loc[ leng:leng+50,['name','p_value','intersections']])
                leng = leng + 50
            print(result.loc[ leng:leng+50,['name','p_value','intersections']])

    result = gp.profile( organism = 'hsapiens', query = recurrent_genes['gene_symbol'].tolist(), sources = ['GO:MF','GO:CC','GO:BP'],  no_evidences = False, no_iea = True)
    leng = 0
    while len(result.index) > leng + 50:
        print(result.loc[ leng:leng+50,['name','p_value','intersections']])
        leng = leng + 50
    print(result.loc[ leng:leng+50,['name','p_value','intersections']])


    return pd.DataFrame(results)

# Usage
freq = pd.read_pickle('data/summarization_large_table_hgnc_quantile.pkl')
freq.loc[freq.loc[:,'membrane_protein'] != 'Yes', 'membrane_protein'] = ''
freq.loc[freq.loc[:,'plasma_protein'] != 'Yes', 'plasma_protein'] = ''
freq.loc[:,'membrane_protein_combined'] = ''
freq.loc[ (freq.loc[:,'membrane_protein'] == 'Yes'),'membrane_protein_combined'] = 'Yes'
freq.loc[ (freq.loc[:,'plasma_protein'] == 'Yes'),'membrane_protein_combined'] = 'Yes'
freq = freq.drop(columns = ['memebrane_protein'])
freq = freq.loc[:, ~freq.columns.str.contains('del_')]
freq = freq.loc[:, ~freq.columns.str.contains('amp_p')]
de = pd.read_csv('data/DE_summary_by_gene_gdc.csv', index_col = 0)
freq = pd.concat([freq,de], axis = 1)

freq = pd.read_pickle('data/summarization_large_table_hgnc_quantile.pkl')
freq.loc[freq.loc[:,'membrane_protein'] != 'Yes', 'membrane_protein'] = ''
freq.loc[freq.loc[:,'plasma_protein'] != 'Yes', 'plasma_protein'] = ''
freq.loc[:,'membrane_protein_combined'] = ''
#freq.loc[ (freq.loc[:,'plasma_protein'] == 'Yes')|(freq.loc[:,'membrane_protein'] == 'Yes'),'membrane_protein_combined'] = 'Yes'
freq.loc[ (freq.loc[:,'membrane_protein'] == 'Yes'),'membrane_protein_combined'] = 'Yes'
freq.loc[ (freq.loc[:,'plasma_protein'] == 'Yes'),'membrane_protein_combined'] = 'Yes'
freq = freq.drop(columns = ['memebrane_protein'])
freq = freq.loc[:, ~freq.columns.str.contains('del_')]
freq = freq.loc[:, ~freq.columns.str.contains('amp_p')]
de = pd.read_csv('data/DE_summary_by_gene_gdc.csv', index_col = 0)
freq = pd.concat([freq,de], axis = 1)
freq.to_csv('data/sum_gdc.csv')
print(freq.loc[:,'avg_log2FoldChange'].mean(), freq.loc[:,'avg_log2FoldChange'].std(), freq.loc[:,'avg_log2FoldChange'].sem())
freq = freq.loc[:,['ER+_amp_freq', 'HER2+_amp_freq', 'TNBC_amp_freq', 'ER+_amp_FDR', 'HER2+_amp_FDR', 'TNBC_amp_FDR', 'All_amp_freq', 'chromosome_name', 'band', 'chr', 'spearman_r', 'spearman_p', 'wilcoxon_stat', 'wilcoxon_fdr', 'plasma_protein', 'membrane_protein', 'membrane_protein_combined', 'CommonEss', 'avg_log2FoldChange', 'n_tissues_sig','mean','mean_quantile']]
freq = freq.loc[ freq.loc[:,'mean_quantile'] >= 3,:]
freq = freq.rename(columns = {'plasma_protein':'plasma_membrane_protein'})
print(freq.loc['ERBB2'])
print(freq.loc['NECTIN4'])

if False:
    for i in ['ER+','HER2+','TNBC','All']:
        print(i)
        results = identify_recurrent_amplicons(freq,i)
        print(results)
if True:
    molecular = ['All','ER+','HER2+','TNBC']
    #thres_freq = freq.loc['ERBB2','All_amp_freq']
    thres_freq = 0.15
    thres_cor = 0.2
    thres_wilcox = 0.001
    freq.loc[:,'sel'] = 0
    for i in molecular:
        if i == 'All':
            #sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) & (freq.loc[:,'spearman_r'] > thres_cor),:]
            sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) & (freq.loc[:,'wilcoxon_stat'] > 0) &  (freq.loc[:,'wilcoxon_fdr'] < thres_wilcox),:]
            freq.loc[sel.index,'sel'] = 1
            freq.loc[:,i] = 0
            freq.loc[sel.index,i] = 1
            print(freq.loc[:,'sel'].sum())
        else:
            sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) &(freq.loc[:,'wilcoxon_stat'] > 0) &  (freq.loc[:,'wilcoxon_fdr'] < thres_wilcox) & (freq.loc[:,i+'_amp_FDR'] < 0.0001)&(freq.loc[:,i+'_amp_freq']>freq.loc[:,'All_amp_freq']) ,:]
            freq.loc[sel.index,'sel'] = 1
            freq.loc[:,i] = 0
            freq.loc[sel.index,i] = 1
            print( (freq.loc[:,i] == 1).sum())
    for i in molecular:
        print(i)
        print(freq.loc[(freq.loc[:,freq.columns.isin(molecular)].sum(axis = 1) == 1)&(freq.loc[:,i] == 1),i].sum())
    print(freq.loc[:,'sel'].sum())
print( freq.loc[ (freq.loc[:,'ER+_amp_freq'] >= 0.15)&(freq.loc[:,'HER2+_amp_freq'] >= 0.15)&(freq.loc[:,'TNBC_amp_freq'] >= 0.15),:])
#print( freq.loc[ (freq.loc[:,'ER+'] == 1)&(freq.loc[:,'HER2+'] == 1)&(freq.loc[:,'TNBC'] == 1),:])
tmp = freq.loc[freq.loc[:,'sel'] == 1,(~freq.columns.str.contains('sel'))]
tmp.to_csv('data/sum_sel_gdc.csv')
if True:
    freq.loc[:,'diff_sel'] = 0
    #freq.loc[ (freq.loc[:,'avg_log2FoldChange'] >= 1.5) & (freq.loc[:,'n_tissues_sig'] >= 12), 'diff_sel'] = 1
    freq.loc[ (freq.loc[:,'avg_log2FoldChange'] >= 1.5), 'diff_sel'] = 1
    print(freq.loc[:,'diff_sel'].sum())
freq = freq.loc[ (freq.loc[:,'diff_sel'] == 1)&(freq.loc[:,'sel'] == 1),:]
#freq = freq.loc[ (freq.loc[:,'membrane_protein'] == 'Yes')|(freq.loc[:,'plasma_protein']=='Yes'),:]
freq = freq.loc[ (freq.loc[:,'membrane_protein_combined'] == 'Yes'),:]
#tmp = freq.loc[:,~freq.columns.str.contains('sel')]
tmp = freq
tmp.to_csv('data/sum_sel_sel_2_gdc.csv')
print(freq)
if False:
    cna_df = pd.read_pickle('data/cna_all_na_filter.pkl')
    mrna_df = pd.read_pickle('data/mrna_all_na_filter.pkl')
    print(cna_df)
    cna_df = cna_df.T
    print(mrna_df)
    mrna_df = mrna_df.T
    print(freq)
    gene_list = freq.index.tolist()
    gene_list.append('CD34')

    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from scipy.stats import spearmanr

    def plot_cna_vs_mrna_to_pdf(cna_df, mrna_df, gene_list, pdf_path='figure/cna_mrna_plots.pdf'):
        # Define CNA categories and palette with string keys
        cna_categories = ['-2', '-1', '0', '1', '2']
        cna_palette = {
            '-2': '#d73027',
            '-1': '#fc8d59',
             '0': '#ffffbf',
             '1': '#91bfdb',
             '2': '#4575b4'
        }

        # Ensure same patients
        common_patients = cna_df.index.intersection(mrna_df.index)
        cna_df = cna_df.loc[common_patients]
        mrna_df = mrna_df.loc[common_patients]

        with PdfPages(pdf_path) as pdf:
            for gene in gene_list:
                if gene not in cna_df.columns or gene not in mrna_df.columns:
                    print(f"Skipping {gene} (not found in both datasets)")
                    continue

                df = pd.DataFrame({
                    'CNA': cna_df[gene],
                    'mRNA': mrna_df[gene]
                }).dropna()

                # Convert CNA to string categorical to match palette keys
                df['CNA'] = df['CNA'].astype(int).astype(str)
                df['CNA'] = pd.Categorical(df['CNA'], categories=cna_categories, ordered=True)
                try:
                    rho, pval = spearmanr(df['CNA'], df['mRNA'])
                    stat_text = f"Spearman r = {rho:.2f}\np = {pval:.1e}"
                except Exception as e:
                    stat_text = "Spearman r = NA\np = NA"
                    print(f"⚠️  Spearman failed for {gene}: {e}")

                plt.figure(figsize=(6, 4))
                sns.boxplot(x='CNA', y='mRNA', data=df, order=cna_categories,
                            palette=cna_palette, showfliers=False)
                sns.stripplot(x='CNA', y='mRNA', data=df, order=cna_categories,
                              color='black', size=3, jitter=True, alpha=0.6)
                plt.text(x=-0.4, y=9.5, s=stat_text,
                         ha='left', va='top', fontsize=10,
                         bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.3'))

                plt.ylim(-4, 10)
                plt.title(f'{gene}: CNA vs mRNA expression')
                plt.xlabel('CNA value')
                plt.ylabel('mRNA expression (Z-score)')
                plt.tight_layout()

                pdf.savefig()
                plt.close()

        print(f"✅ Plots saved to: {pdf_path}")
    plot_cna_vs_mrna_to_pdf(cna_df, mrna_df, gene_list)
#'''

#chromosome investigation
'''
freq = pd.read_pickle('data/summarization_large_table_hgnc_quantile.pkl')
freq.loc[freq.loc[:,'membrane_protein'] != 'Yes', 'membrane_protein'] = ''
freq.loc[freq.loc[:,'plasma_protein'] != 'Yes', 'plasma_protein'] = ''
freq.loc[:,'membrane_protein_combined'] = ''
freq.loc[ (freq.loc[:,'membrane_protein'] == 'Yes'),'membrane_protein_combined'] = 'Yes'
freq.loc[ (freq.loc[:,'plasma_protein'] == 'Yes'),'membrane_protein_combined'] = 'Yes'
freq = freq.drop(columns = ['memebrane_protein'])
freq = freq.loc[:, ~freq.columns.str.contains('del_')]
freq = freq.loc[:, ~freq.columns.str.contains('amp_p')]
de = pd.read_csv('data/DE_summary_by_gene_gdc.csv', index_col = 0)
freq = pd.concat([freq,de], axis = 1)
freq = freq.loc[:,['ER+_amp_freq', 'HER2+_amp_freq', 'TNBC_amp_freq', 'ER+_amp_FDR', 'HER2+_amp_FDR', 'TNBC_amp_FDR', 'All_amp_freq', 'chromosome_name', 'band', 'chr', 'spearman_r', 'spearman_p', 'wilcoxon_stat', 'wilcoxon_fdr', 'plasma_protein', 'membrane_protein', 'membrane_protein_combined', 'CommonEss', 'avg_log2FoldChange', 'n_tissues_sig','mean','mean_quantile']]
#freq = freq.loc[ freq.loc[:,'mean_quantile'] >= 3,:]
freq = freq.rename(columns = {'plasma_protein':'plasma_membrane_protein'})
from gprofiler import GProfiler
gp = GProfiler( return_dataframe = True)

msk = pd.read_csv('data/cancerGeneList.tsv', sep = '\t')
print(msk.columns)
onco = msk.loc[msk.loc[:,'Is Oncogene'] == 'Yes',:]

if True:
    molecular = ['All','ER+','HER2+','TNBC']
    #thres_freq = freq.loc['ERBB2','All_amp_freq']
    thres_freq = 0.15
    thres_wilcox = 0.001
    freq.loc[:,'sel'] = 0
    for i in molecular:
        if i == 'All':
            #sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) & (freq.loc[:,'spearman_r'] > thres_cor),:]
            #sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) & (freq.loc[:,'wilcoxon_stat'] > 0) &  (freq.loc[:,'wilcoxon_fdr'] < thres_wilcox),:]
            sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq),:]
            freq.loc[sel.index,'sel'] = 1
            freq.loc[:,i] = 0
            freq.loc[sel.index,i] = 1
            print(freq.loc[:,'sel'].sum())
            for j in ['1','8']:
                print(j)
                ingene = freq.loc[freq.loc[:,'All'] == 1,:]
                ingene = ingene.loc[ingene.loc[:,'chromosome_name'] == j,:].index.tolist()
                print( onco.loc[ onco.loc[:,'Hugo Symbol'].isin(ingene),:])
                #result = gp.profile( organism = 'hsapiens', query = ingene, sources = ['GO:MF'], user_threshold = 0.0001, no_evidences = False, no_iea = True)
                result = gp.profile( organism = 'hsapiens', query = ingene, sources = ['GO:MF','GO:CC','GO:BP'],  no_evidences = False, no_iea = True)
                leng = 0
                while len(result.index) > leng + 50:
                    print(result.loc[ leng:leng+50,['name','p_value','intersections']])
                    leng = leng + 50
                print(result.loc[ leng:leng+50,['name','p_value','intersections']])
        else:
            #sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) &(freq.loc[:,'wilcoxon_stat'] > 0) &  (freq.loc[:,'wilcoxon_fdr'] < thres_wilcox) & (freq.loc[:,i+'_amp_FDR'] < 0.0001)&(freq.loc[:,i+'_amp_freq']>freq.loc[:,'All_amp_freq']) ,:]
            sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) & (freq.loc[:,i+'_amp_FDR'] < 0.0001)&(freq.loc[:,i+'_amp_freq']>freq.loc[:,'All_amp_freq']) ,:]
            freq.loc[sel.index,'sel'] = 1
            freq.loc[:,i] = 0
            freq.loc[sel.index,i] = 1
            print( (freq.loc[:,i] == 1).sum())
            print(i)
            for j in ['1','8','11','17']:
                print(j)
                ingene = freq.loc[freq.loc[:,i] == 1,:]
                ingene = ingene.loc[ingene.loc[:,'chromosome_name'] == j,:].index.tolist()
                if i == 'ER+' and j == '1':
                    eringene = set(ingene)
                print( onco.loc[ onco.loc[:,'Hugo Symbol'].isin(ingene),:])
                if False:
                    if len(ingene) == 0:
                        continue
                    #result = gp.profile( organism = 'hsapiens', query = ingene, sources = ['GO:MF'], user_threshold = 0.0001, no_evidences = False, no_iea = True)
                    result = gp.profile( organism = 'hsapiens', query = ingene, sources = ['GO:MF','GO:CC','GO:BP'],  no_evidences = False, no_iea = True)
                    leng = 0
                    while len(result.index) > leng + 50:
                        print(result.loc[ leng:leng+50,['name','p_value','intersections']])
                        leng = leng + 50
                    print(result.loc[ leng:leng+50,['name','p_value','intersections']])
    for i in molecular:
        print(i)
        print(freq.loc[(freq.loc[:,freq.columns.isin(molecular)].sum(axis = 1) == 1)&(freq.loc[:,i] == 1),i].sum())

    print(freq.loc[:,'sel'].sum())
if False:
    print( freq.loc[ (freq.loc[:,'ER+_amp_freq'] >= 0.15)&(freq.loc[:,'HER2+_amp_freq'] >= 0.15)&(freq.loc[:,'TNBC_amp_freq'] >= 0.15),:])
    #print( freq.loc[ (freq.loc[:,'ER+'] == 1)&(freq.loc[:,'HER2+'] == 1)&(freq.loc[:,'TNBC'] == 1),:])
    tmp = freq.loc[freq.loc[:,'sel'] == 1,(~freq.columns.str.contains('sel'))]
    tmp.to_csv('data/sum_sel_gdc.csv')
    if True:
        freq.loc[:,'diff_sel'] = 0
        #freq.loc[ (freq.loc[:,'avg_log2FoldChange'] >= 1.5) & (freq.loc[:,'n_tissues_sig'] >= 12), 'diff_sel'] = 1
        freq.loc[ (freq.loc[:,'avg_log2FoldChange'] >= 1.5), 'diff_sel'] = 1
        print(freq.loc[:,'diff_sel'].sum())
    freq = freq.loc[ (freq.loc[:,'diff_sel'] == 1)&(freq.loc[:,'sel'] == 1),:]
    #freq = freq.loc[ (freq.loc[:,'membrane_protein'] == 'Yes')|(freq.loc[:,'plasma_protein']=='Yes'),:]
    freq = freq.loc[ (freq.loc[:,'membrane_protein_combined'] == 'Yes'),:]
    #tmp = freq.loc[:,~freq.columns.str.contains('sel')]
    tmp = freq
    print(freq)
#'''



'''
freq = pd.read_csv('data/sum_sel_gdc.csv', index_col = 0)
ref = pd.read_csv('data/ref_update.txt', index_col = 0, sep = '\t')
print(freq)
ref = ref.loc[:,['location','classification','Reference','internalization','reference','ADC reported']]
print(ref)
freq = pd.concat([freq, ref], axis = 1)
print(freq)
freq.to_csv('data/sum_sel_gdc_ref.csv')
#'''

#plot generate
'''
freq = pd.read_csv('data/final_anno.csv', index_col = 0)
freq = pd.read_csv('Amp_1702.txt', sep='\t', index_col = 0)
print(freq)
freq.loc[:,'internalization'] = freq.loc[:,'Finalist']
freq = freq.loc[:,['internalization','mean','avg_log2FoldChange']]
freq = freq.loc[ freq.loc[:,'internalization'] == 1,:]
freq = freq.sort_values(by = 'avg_log2FoldChange', ascending = False)
output_f = pd.read_pickle('data/DE_result_gdc_f.pkl')
output_p = pd.read_pickle('data/DE_result_gdc_p.pkl')
output_f = output_f.reindex(freq.index)
output_p = output_p.reindex(freq.index)
print(output_f)
print(output_p)
output_f.to_csv('R_data/adc_target_log2.csv')
output_p.to_csv('R_data/adc_target_p.csv')
#'''

#fisher
'''
import pandas as pd
from gprofiler import GProfiler

gp = GProfiler(return_dataframe=True)

def identify_recurrent_amplicons(filepath, subtype, frequency_threshold=0.15, gap_limit=1000000):
    """
    Identifies recurrent amplicons from GISTIC gene-level data.
    """
    # 1. Load Data
    df = filepath.copy()

    # 2. Calculate Amplification Frequency
    df.loc[:, 'gene_symbol'] = df.index
    recurrent_genes = df[df[subtype+'_amp_freq'] > frequency_threshold].copy()

    if subtype != 'All':
        if subtype == 'ER+':
            recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'HER2+_amp_freq']>frequency_threshold)&(recurrent_genes.loc[:,'HER2+_amp_FDR']<0.0001))]
            recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'TNBC_amp_freq']>frequency_threshold)&(recurrent_genes.loc[:,'TNBC_amp_FDR']<0.0001))]
        if subtype == 'HER2+':
            recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'ER+_amp_freq']>frequency_threshold)&(recurrent_genes.loc[:,'ER+_amp_FDR']<0.0001))]
            recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'TNBC_amp_freq']>frequency_threshold)&(recurrent_genes.loc[:,'TNBC_amp_FDR']<0.0001))]
        if subtype == 'TNBC':
            recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'HER2+_amp_freq']>frequency_threshold)&(recurrent_genes.loc[:,'HER2+_amp_FDR']<0.0001))]
            recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'ER+_amp_freq']>frequency_threshold)&(recurrent_genes.loc[:,'ER+_amp_FDR']<0.0001))]

        recurrent_genes = recurrent_genes[recurrent_genes.loc[:,subtype+'_amp_FDR'] < 0.0001]
        recurrent_genes = recurrent_genes[recurrent_genes.loc[:,subtype+'_amp_freq'] > recurrent_genes.loc[:,'All_amp_freq']]

    subtype_col = f"{subtype}_amp_freq"
    print(f"Number of recurrent genes for {subtype}: {len(recurrent_genes)}")
    if recurrent_genes.empty:
        print("No genes met the frequency threshold.")
        return None

    # 4. Parse Genomic Coordinates
    recurrent_genes = recurrent_genes.sort_values(['chromosome_name', 'start_position'])

    # Calculate the gap (Current Start minus Previous End)
    gap_distances = recurrent_genes['start_position'] - recurrent_genes['end_position'].shift()

    # Define the break
    recurrent_genes['Amplicon_ID'] = (
        (recurrent_genes['chromosome_name'] != recurrent_genes['chromosome_name'].shift()) |
        (gap_distances > gap_limit)
    ).cumsum()

    # 6. Identify Peak Genes
    results = []
    for amp_id, group in recurrent_genes.groupby('Amplicon_ID'):
        peak_gene = group.loc[group[subtype+'_amp_freq'].idxmax()]

        results.append({
            'Chrom': group['chromosome_name'].iloc[0],
            'Peak_Gene': peak_gene['gene_symbol'],
            'Peak_Frequency': peak_gene[subtype_col],
            'Gene_Count': len(group),
            'All_Genes': ", ".join(group['gene_symbol'].tolist())
        })

    # Functional enrichment analysis via gprofiler
    result = gp.profile(organism='hsapiens', query=recurrent_genes['gene_symbol'].tolist(), sources=['GO:MF','GO:CC','GO:BP'], no_evidences=False, no_iea=True)
    leng = 0
    while len(result.index) > leng + 50:
        print(result.loc[leng:leng+50, ['name', 'p_value', 'intersections']])
        leng = leng + 50
    print(result.loc[leng:leng+50, ['name', 'p_value', 'intersections']])

    return pd.DataFrame(results)


# --- DATA PROCESSING MAIN PIPELINE ---

# 1. Load baseline genomic and clinical table
freq = pd.read_pickle('data/summarization_large_table_hgnc_quantile.pkl')
freq.loc[freq.loc[:,'membrane_protein'] != 'Yes', 'membrane_protein'] = ''
freq.loc[freq.loc[:,'plasma_protein'] != 'Yes', 'plasma_protein'] = ''
freq.loc[:,'membrane_protein_combined'] = ''
freq.loc[(freq.loc[:,'membrane_protein'] == 'Yes'), 'membrane_protein_combined'] = 'Yes'
freq.loc[(freq.loc[:,'plasma_protein'] == 'Yes'), 'membrane_protein_combined'] = 'Yes'

if 'memebrane_protein' in freq.columns:
    freq = freq.drop(columns=['memebrane_protein'])
freq = freq.loc[:, ~freq.columns.str.contains('del_')]
freq = freq.loc[:, ~freq.columns.str.contains('amp_p')]

# Load Expression Summary Data
de = pd.read_csv('data/DE_summary_by_gene_gdc.csv', index_col=0)
freq = pd.concat([freq, de], axis=1)

# 2. OVERWRITE WITH FISHER EXACT TEST FDR DATA
fisher_df = pd.read_pickle('data/cna_subtype_freq_stat_FDR_fisher.pkl')
fisher_fdr_cols = ['ER+_amp_FDR', 'HER2+_amp_FDR', 'TNBC_amp_FDR']

# Ensure old columns are dropped or overwritten precisely by Fisher's metrics
freq = freq.drop(columns=[col for col in fisher_fdr_cols if col in freq.columns])
freq = freq.join(fisher_df[fisher_fdr_cols], how='left')

# Save intermediate matrix with the fisher suffix
freq.to_csv('data/sum_gdc_fisher.csv')
print("Log2 Fold Change Mean, Std, SEM:", freq.loc[:,'avg_log2FoldChange'].mean(), freq.loc[:,'avg_log2FoldChange'].std(), freq.loc[:,'avg_log2FoldChange'].sem())

# Filter down to specific analysis subset (Including positional data needed for amplicon calculation)
keep_cols = [
    'ER+_amp_freq', 'HER2+_amp_freq', 'TNBC_amp_freq',
    'ER+_amp_FDR', 'HER2+_amp_FDR', 'TNBC_amp_FDR', 'All_amp_freq',
    'chromosome_name', 'start_position', 'end_position', 'band', 'chr',
    'spearman_r', 'spearman_p', 'wilcoxon_stat', 'wilcoxon_fdr',
    'plasma_protein', 'membrane_protein', 'membrane_protein_combined',
    'CommonEss', 'avg_log2FoldChange', 'n_tissues_sig', 'mean', 'mean_quantile'
]
freq = freq.loc[:, [c for c in keep_cols if c in freq.columns]]
freq = freq.loc[freq.loc[:,'mean_quantile'] >= 3, :]
freq = freq.rename(columns={'plasma_protein': 'plasma_membrane_protein'})

# Print check control targets
if 'ERBB2' in freq.index: print(freq.loc['ERBB2'])
if 'NECTIN4' in freq.index: print(freq.loc['NECTIN4'])

# Evaluate recurrent amplicons with updated Fisher stats
for i in ['ER+', 'HER2+', 'TNBC', 'All']:
    print(f"\n--- Processing Subtype: {i} ---")
    results = identify_recurrent_amplicons(freq, i)
    print(results)

# 3. SELECT CANDIDATES BASED ON FISHER SIGNIFICANCE
molecular = ['All', 'ER+', 'HER2+', 'TNBC']
thres_freq = 0.15
thres_cor = 0.2
thres_wilcox = 0.001
freq.loc[:, 'sel'] = 0

for i in molecular:
    if i == 'All':
        sel = freq.loc[(freq.loc[:, i+'_amp_freq'] > thres_freq) & (freq.loc[:, 'wilcoxon_stat'] > 0) & (freq.loc[:, 'wilcoxon_fdr'] < thres_wilcox), :]
        freq.loc[sel.index, 'sel'] = 1
        freq.loc[:, i] = 0
        freq.loc[sel.index, i] = 1
        print("All Subtypes Cumulative Selection Count:", freq.loc[:, 'sel'].sum())
    else:
        # Filtering with the fresh Fisher '_amp_FDR' columns
        sel = freq.loc[(freq.loc[:, i+'_amp_freq'] > thres_freq) & (freq.loc[:, 'wilcoxon_stat'] > 0) & (freq.loc[:, 'wilcoxon_fdr'] < thres_wilcox) & (freq.loc[:, i+'_amp_FDR'] < 0.0001) & (freq.loc[:, i+'_amp_freq'] > freq.loc[:, 'All_amp_freq']), :]
        freq.loc[sel.index, 'sel'] = 1
        freq.loc[:, i] = 0
        freq.loc[sel.index, i] = 1
        print(f"Subtype {i} Independent Count:", (freq.loc[:, i] == 1).sum())

for i in molecular:
    print(f"Exclusive to {i}:")
    print(freq.loc[(freq.loc[:, freq.columns.isin(molecular)].sum(axis=1) == 1) & (freq.loc[:, i] == 1), i].sum())
print("Total Unique Selection Count:", freq.loc[:, 'sel'].sum())

# Export selection hits to unique Fisher file
tmp = freq.loc[freq.loc[:, 'sel'] == 1, (~freq.columns.str.contains('sel'))]
tmp.to_csv('data/sum_sel_gdc_fisher.csv')

# 4. MEMBRANE-TARGET STRINGENT FILTRATION
freq.loc[:, 'diff_sel'] = 0
freq.loc[(freq.loc[:, 'avg_log2FoldChange'] >= 1), 'diff_sel'] = 1
print("Differentially Expressed selection total:", freq.loc[:, 'diff_sel'].sum())

freq = freq.loc[(freq.loc[:, 'diff_sel'] == 1) & (freq.loc[:, 'sel'] == 1), :]
freq = freq.loc[(freq.loc[:, 'membrane_protein_combined'] == 'Yes'), :]

# Export clean targets with the fisher suffix
freq.to_csv('data/sum_sel_sel_2_gdc_fisher.csv')
print("\nFinal Selected ADC Candidate Targets Utilizing Fisher Statistics:")
print(freq)
#'''

#fisher compare
'''
import pandas as pd

def compare_selections(chisq_path, fisher_path, stage_name):
    """Compares gene selections between Chi-square and Fisher pipelines."""
    try:
        df_chisq = pd.read_csv(chisq_path, index_col=0)
        df_fisher = pd.read_csv(fisher_path, index_col=0)

        genes_chisq = set(df_chisq.index)
        genes_fisher = set(df_fisher.index)

        shared = genes_chisq.intersection(genes_fisher)
        dropped_by_fisher = genes_chisq - genes_fisher
        added_by_fisher = genes_fisher - genes_chisq

        print(f"=== Comparison for {stage_name} ===")
        print(f"Chi-square unique genes: {len(genes_chisq)}")
        print(f"Fisher exact unique genes: {len(genes_fisher)}")
        print(f"Shared genes (both): {len(shared)}")
        print(f"Dropped by switching to Fisher (n={len(dropped_by_fisher)}): {sorted(list(dropped_by_fisher))}")
        print(f"Newly added by Fisher (n={len(added_by_fisher)}): {sorted(list(added_by_fisher))}\n")

        return {
            'shared': shared,
            'dropped': dropped_by_fisher,
            'added': added_by_fisher
        }
    except FileNotFoundError as e:
        print(f"Skipping {stage_name}: One of the files was not found. ({e.filename})")
        return None

# 1. Compare the broad FDR + Wilcoxon selections
fdr_comparison = compare_selections(
    chisq_path='data/sum_sel_gdc.csv',
    fisher_path='data/sum_sel_gdc_fisher.csv',
    stage_name='Initial FDR Selection (sum_sel_gdc)'
)

# 2. Compare the final filtered ADC candidates (Membrane bound + Expression filter)
final_comparison = compare_selections(
    chisq_path='data/sum_sel_sel_2_gdc.csv',
    fisher_path='data/sum_sel_sel_2_gdc_fisher.csv',
    stage_name='Final ADC Target Selection (sum_sel_sel_2_gdc)'
)
#'''

#fisher merge
'''
import pandas as pd
import numpy as np

# 1. Load the core dataframes
df = pd.read_pickle('data/summarization_large_table_quantile.pkl')
fisher_df = pd.read_pickle('data/cna_subtype_freq_stat_FDR_fisher.pkl')

# 2. Add Fisher FDR metrics into the base dataframe before symbol mapping
fisher_cols = ['ER+_amp_FDR', 'HER2+_amp_FDR', 'TNBC_amp_FDR']
# Drop the columns from df if they already exist to prevent duplicate suffix errors
df = df.drop(columns=[col for col in fisher_cols if col in df.columns], errors='ignore')
df = df.join(fisher_df[fisher_cols], how='left')

# 3. Load and clean the HGNC Symbol Mapping file
hgnc = pd.read_csv('data/hgnc-symbol-check.csv')
# Explicitly drop any rows that lack a valid raw input or approved target symbol
hgnc = hgnc.dropna(subset=[hgnc.columns[0], 'Approved symbol'])

# 4. Vectorized Mapping Engine
# Create a direct translation dictionary from the historical/input symbol to approved symbol
mapping_dict = dict(zip(hgnc.iloc[:, 0], hgnc['Approved symbol']))

# Reset index to treat the gene symbol as a regular column for manipulation
df = df.reset_index().rename(columns={'index': 'Input_Symbol'})

# Map to the new approved symbols. If a gene isn't in hgnc, it retains its original name.
df['Approved_Symbol'] = df['Input_Symbol'].map(mapping_dict).fillna(df['Input_Symbol'])

# 5. Resolve Overlaps / Fills Vectorially
# Group by the new Approved Symbol. If multiple old genes map to the same approved gene,
# 'first' keeps the earliest valid record, and `.bfill().ffill()` recovers missing values across rows.
df_clean = (
    df.groupby('Approved_Symbol')
    .apply(lambda g: g.bfill().ffill().iloc[0] if len(g) > 1 else g.iloc[0])
)

# Clean up structural tracking columns
if 'Input_Symbol' in df_clean.columns:
    df_clean = df_clean.drop(columns=['Input_Symbol'])

# 6. Verification & Export
print("--- Sample verification for NECTIN4 ---")
if 'NECTIN4' in df_clean.index:
    print(df_clean.loc['NECTIN4', :])
else:
    print("NECTIN4 not found in the processed table.")

# Save directly to the requested fisher target file
df_clean.to_pickle('data/summarization_large_table_hgnc_quantile_fisher.pkl')
print("\nSuccessfully updated symbols and saved to: data/summarization_large_table_hgnc_quantile_fisher.pkl")
#'''

#fisher code
'''
import pandas as pd
from gprofiler import GProfiler

gp = GProfiler(return_dataframe=True)

def identify_recurrent_amplicons(filepath, subtype, frequency_threshold=0.15, gap_limit=1000000):
    """
    Identifies recurrent amplicons from GISTIC gene-level data.
    """
    df = filepath.copy()
    df.loc[:, 'gene_symbol'] = df.index
    recurrent_genes = df[df[subtype+'_amp_freq'] > frequency_threshold].copy()

    if subtype != 'All':
        recurrent_genes = recurrent_genes[recurrent_genes.loc[:,subtype+'_amp_FDR'] < 0.0001]
        recurrent_genes = recurrent_genes[recurrent_genes.loc[:,subtype+'_amp_freq'] > recurrent_genes.loc[:,'All_amp_freq']]
        print(recurrent_genes)
        if subtype == 'ER+':
            recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'HER2+_amp_freq']>frequency_threshold)&(recurrent_genes.loc[:,'HER2+_amp_FDR']<0.0001))]
            recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'TNBC_amp_freq']>frequency_threshold)&(recurrent_genes.loc[:,'TNBC_amp_FDR']<0.0001))]
        if subtype == 'HER2+':
            recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'ER+_amp_freq']>frequency_threshold)&(recurrent_genes.loc[:,'ER+_amp_FDR']<0.0001))]
            recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'TNBC_amp_freq']>frequency_threshold)&(recurrent_genes.loc[:,'TNBC_amp_FDR']<0.0001))]
        if subtype == 'TNBC':
            recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'HER2+_amp_freq']>frequency_threshold)&(recurrent_genes.loc[:,'HER2+_amp_FDR']<0.0001))]
            recurrent_genes = recurrent_genes[~((recurrent_genes.loc[:,'ER+_amp_freq']>frequency_threshold)&(recurrent_genes.loc[:,'ER+_amp_FDR']<0.0001))]

        recurrent_genes = recurrent_genes[recurrent_genes.loc[:,subtype+'_amp_FDR'] < 0.0001]
        recurrent_genes = recurrent_genes[recurrent_genes.loc[:,subtype+'_amp_freq'] > recurrent_genes.loc[:,'All_amp_freq']]

    subtype_col = f"{subtype}_amp_freq"
    print(f"Number of recurrent genes for {subtype}: {len(recurrent_genes)}")
    if recurrent_genes.empty:
        print("No genes met the frequency threshold.")
        return None

    recurrent_genes = recurrent_genes.sort_values(['chromosome_name', 'start_position'])
    gap_distances = recurrent_genes['start_position'] - recurrent_genes['end_position'].shift()

    recurrent_genes['Amplicon_ID'] = (
        (recurrent_genes['chromosome_name'] != recurrent_genes['chromosome_name'].shift()) |
        (gap_distances > gap_limit)
    ).cumsum()

    results = []
    for amp_id, group in recurrent_genes.groupby('Amplicon_ID'):
        peak_gene = group.loc[group[subtype+'_amp_freq'].idxmax()]
        results.append({
            'Chrom': group['chromosome_name'].iloc[0],
            'Peak_Gene': peak_gene['gene_symbol'],
            'Peak_Frequency': peak_gene[subtype_col],
            'Gene_Count': len(group),
            'All_Genes': ", ".join(group['gene_symbol'].tolist())
        })

    result = gp.profile(organism='hsapiens', query=recurrent_genes['gene_symbol'].tolist(), sources=['GO:MF','GO:CC','GO:BP'], no_evidences=False, no_iea=True)
    leng = 0
    while len(result.index) > leng + 50:
        print(result.loc[leng:leng+50, ['name', 'p_value', 'intersections']])
        leng = leng + 50
    print(result.loc[leng:leng+50, ['name', 'p_value', 'intersections']])

    return pd.DataFrame(results)


# --- OPTIMIZED PIPELINE VIA PRE-MERGED FISHER DATAFRAME ---

# 1. Load the pre-merged table containing clean HGNC symbols and Fisher FDRs
freq = pd.read_pickle('data/summarization_large_table_hgnc_quantile_fisher.pkl')

# Normalize the protein localization annotations
freq.loc[freq['membrane_protein'] != 'Yes', 'membrane_protein'] = ''
freq.loc[freq['plasma_protein'] != 'Yes', 'plasma_protein'] = ''
freq['membrane_protein_combined'] = ''
freq.loc[(freq['membrane_protein'] == 'Yes'), 'membrane_protein_combined'] = 'Yes'
freq.loc[(freq['plasma_protein'] == 'Yes'), 'membrane_protein_combined'] = 'Yes'

if 'memebrane_protein' in freq.columns:
    freq = freq.drop(columns=['memebrane_protein'])
freq = freq.loc[:, ~freq.columns.str.contains('del_')]
freq = freq.loc[:, ~freq.columns.str.contains('amp_p')]

# 2. Add Differential Expression metrics
de = pd.read_csv('data/DE_summary_by_gene_gdc.csv', index_col=0)
# Use combine_first to preserve existing columns in freq and safely pull missing ones from de
freq = freq.combine_first(de)

# Save the updated complete table
freq.to_csv('data/sum_gdc_fisher.csv')
print("Log2 Fold Change Mean, Std, SEM:", freq['avg_log2FoldChange'].mean(), freq['avg_log2FoldChange'].std(), freq['avg_log2FoldChange'].sem())

# Filter down columns and rows to the analytical population
keep_cols = [
    'ER+_amp_freq', 'HER2+_amp_freq', 'TNBC_amp_freq',
    'ER+_amp_FDR', 'HER2+_amp_FDR', 'TNBC_amp_FDR', 'All_amp_freq',
    'chromosome_name', 'start_position', 'end_position', 'band', 'chr',
    'spearman_r', 'spearman_p', 'wilcoxon_stat', 'wilcoxon_fdr',
    'plasma_protein', 'membrane_protein', 'membrane_protein_combined',
    'CommonEss', 'avg_log2FoldChange', 'n_tissues_sig', 'mean', 'mean_quantile'
]
freq = freq.loc[:, [c for c in keep_cols if c in freq.columns]]
#freq = freq.loc[freq['mean_quantile'] >= 3, :]
freq = freq.rename(columns={'plasma_protein': 'plasma_membrane_protein'})

# Verification controls check
if 'ERBB2' in freq.index: print(freq.loc['ERBB2'])
if 'NECTIN4' in freq.index: print(freq.loc['NECTIN4'])

# 3. Evaluate Recurrent Amplicons
for i in ['ER+', 'HER2+', 'TNBC', 'All']:
    print(f"\n--- Processing Subtype: {i} ---")
    results = identify_recurrent_amplicons(freq, i)
    print(results)

# 4. Statistical Selection
molecular = ['All', 'ER+', 'HER2+', 'TNBC']
thres_freq, thres_cor, thres_wilcox = 0.15, 0.2, 0.001
freq['sel'] = 0

for i in molecular:
    if i == 'All':
        sel = freq.loc[(freq[i+'_amp_freq'] > thres_freq) & (freq['wilcoxon_stat'] > 0) & (freq['wilcoxon_fdr'] < thres_wilcox), :]
        freq.loc[sel.index, 'sel'] = 1
        freq.loc[:, i] = 0
        freq.loc[sel.index, i] = 1
        print("All Subtypes Cumulative Selection Count:", freq['sel'].sum())
    else:
        sel = freq.loc[(freq[i+'_amp_freq'] > thres_freq) & (freq['wilcoxon_stat'] > 0) & (freq['wilcoxon_fdr'] < thres_wilcox) & (freq[i+'_amp_FDR'] < 0.0001) & (freq[i+'_amp_freq'] > freq['All_amp_freq']), :]
        freq.loc[sel.index, 'sel'] = 1
        freq.loc[:, i] = 0
        freq.loc[sel.index, i] = 1
        print(f"Subtype {i} Independent Count:", (freq[i] == 1).sum())
# --- Replace your current exclusive loop block with this ---
for i in molecular:
    print(f"Exclusive to {i}:")
    # 1. Isolate just the molecular columns that actually exist in the dataframe
    existing_molecular_cols = [col for col in molecular if col in freq.columns]
    
    # 2. Count how many subtypes are active (equal to 1) for each gene row
    subtype_counts = (freq[existing_molecular_cols] == 1).sum(axis=1)
    
    # 3. Filter for rows where total count is exactly 1 AND the specific subtype is active
    is_exclusive = (subtype_counts == 1) & (freq[i] == 1)
    
    # 4. Sum up the occurrences safely
    print(is_exclusive.sum())
# Intermediate export for selected genes
tmp = freq.loc[freq['sel'] == 1, (~freq.columns.str.contains('sel'))]
tmp.to_csv('data/sum_sel_gdc_fisher.csv')

# 5. Membrane & Expression Filter (using log2FoldChange >= 1 threshold)
freq['diff_sel'] = 0
freq.loc[(freq['avg_log2FoldChange'] >= 1), 'diff_sel'] = 1
print("Differentially Expressed selection total:", freq['diff_sel'].sum())

freq = freq.loc[(freq['diff_sel'] == 1) & (freq['sel'] == 1), :]
freq = freq.loc[(freq['membrane_protein_combined'] == 'Yes'), :]

# Final target export
freq.to_csv('data/sum_sel_sel_2_gdc_fisher.csv')
print("\nFinal Selected ADC Candidate Targets Utilizing Fisher Statistics:")
print(freq)
#'''

#chromosome investigation
'''
freq = pd.read_pickle('data/summarization_large_table_hgnc_quantile_fisher.pkl')
freq.loc[freq.loc[:,'membrane_protein'] != 'Yes', 'membrane_protein'] = ''
freq.loc[freq.loc[:,'plasma_protein'] != 'Yes', 'plasma_protein'] = ''
freq.loc[:,'membrane_protein_combined'] = ''
freq.loc[ (freq.loc[:,'membrane_protein'] == 'Yes'),'membrane_protein_combined'] = 'Yes'
freq.loc[ (freq.loc[:,'plasma_protein'] == 'Yes'),'membrane_protein_combined'] = 'Yes'
freq = freq.drop(columns = ['memebrane_protein'])
freq = freq.loc[:, ~freq.columns.str.contains('del_')]
freq = freq.loc[:, ~freq.columns.str.contains('amp_p')]
de = pd.read_csv('data/DE_summary_by_gene_gdc.csv', index_col = 0)
freq = pd.concat([freq,de], axis = 1)
freq = freq.loc[:,['ER+_amp_freq', 'HER2+_amp_freq', 'TNBC_amp_freq', 'ER+_amp_FDR', 'HER2+_amp_FDR', 'TNBC_amp_FDR', 'All_amp_freq', 'chromosome_name', 'band', 'chr', 'spearman_r', 'spearman_p', 'wilcoxon_stat', 'wilcoxon_fdr', 'plasma_protein', 'membrane_protein', 'membrane_protein_combined', 'CommonEss', 'avg_log2FoldChange', 'n_tissues_sig','mean','mean_quantile']]
#freq = freq.loc[ freq.loc[:,'mean_quantile'] >= 3,:]
freq = freq.rename(columns = {'plasma_protein':'plasma_membrane_protein'})
from gprofiler import GProfiler
gp = GProfiler( return_dataframe = True)

msk = pd.read_csv('data/cancerGeneList.tsv', sep = '\t')
print(msk.columns)
onco = msk.loc[msk.loc[:,'Is Oncogene'] == 'Yes',:]

if True:
    molecular = ['All','ER+','HER2+','TNBC']
    #thres_freq = freq.loc['ERBB2','All_amp_freq']
    thres_freq = 0.15
    thres_wilcox = 0.001
    freq.loc[:,'sel'] = 0
    for i in molecular:
        if i == 'All':
            #sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) & (freq.loc[:,'spearman_r'] > thres_cor),:]
            #sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) & (freq.loc[:,'wilcoxon_stat'] > 0) &  (freq.loc[:,'wilcoxon_fdr'] < thres_wilcox),:]
            sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq),:]
            freq.loc[sel.index,'sel'] = 1
            freq.loc[:,i] = 0
            freq.loc[sel.index,i] = 1
            print(freq.loc[:,'sel'].sum())
            for j in ['1','8']:
                print(j)
                ingene = freq.loc[freq.loc[:,'All'] == 1,:]
                ingene = ingene.loc[ingene.loc[:,'chromosome_name'] == j,:].index.tolist()
                print( onco.loc[ onco.loc[:,'Hugo Symbol'].isin(ingene),:])
                #result = gp.profile( organism = 'hsapiens', query = ingene, sources = ['GO:MF'], user_threshold = 0.0001, no_evidences = False, no_iea = True)
                result = gp.profile( organism = 'hsapiens', query = ingene, sources = ['GO:MF','GO:CC','GO:BP'],  no_evidences = False, no_iea = True)
                leng = 0
                while len(result.index) > leng + 50:
                    print(result.loc[ leng:leng+50,['name','p_value','intersections']])
                    leng = leng + 50
                print(result.loc[ leng:leng+50,['name','p_value','intersections']])
        else:
            #sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) &(freq.loc[:,'wilcoxon_stat'] > 0) &  (freq.loc[:,'wilcoxon_fdr'] < thres_wilcox) & (freq.loc[:,i+'_amp_FDR'] < 0.0001)&(freq.loc[:,i+'_amp_freq']>freq.loc[:,'All_amp_freq']) ,:]
            sel = freq.loc[ (freq.loc[:,i+'_amp_freq'] > thres_freq) & (freq.loc[:,i+'_amp_FDR'] < 0.0001)&(freq.loc[:,i+'_amp_freq']>freq.loc[:,'All_amp_freq']) ,:]
            freq.loc[sel.index,'sel'] = 1
            freq.loc[:,i] = 0
            freq.loc[sel.index,i] = 1
            print( (freq.loc[:,i] == 1).sum())
            print(i)
            for j in ['1','8']:
                print(j)
                ingene = freq.loc[freq.loc[:,i] == 1,:]
                ingene = ingene.loc[ingene.loc[:,'chromosome_name'] == j,:].index.tolist()
                if i == 'ER+' and j == '1':
                    eringene = set(ingene)
                print( onco.loc[ onco.loc[:,'Hugo Symbol'].isin(ingene),:])
                if len(ingene) == 0:
                    continue
                #result = gp.profile( organism = 'hsapiens', query = ingene, sources = ['GO:MF'], user_threshold = 0.0001, no_evidences = False, no_iea = True)
                result = gp.profile( organism = 'hsapiens', query = ingene, sources = ['GO:MF','GO:CC','GO:BP'],  no_evidences = False, no_iea = True)
                leng = 0
                while len(result.index) > leng + 50:
                    print(result.loc[ leng:leng+50,['name','p_value','intersections']])
                    leng = leng + 50
                print(result.loc[ leng:leng+50,['name','p_value','intersections']])
            if i == 'TNBC':
                sel = freq.loc[ (freq.loc[:,i+'_amp_FDR'] < 0.0001)&(freq.loc[:,'chromosome_name'] == '1') ,:]
                print(sel)
                ingene = sel.index.tolist()
                if len(ingene) == 0:
                    continue
                #result = gp.profile( organism = 'hsapiens', query = ingene, sources = ['GO:MF'], user_threshold = 0.0001, no_evidences = False, no_iea = True)
                tningene = set(ingene)
                print('er',len(eringene))
                print('tnbc',len(tningene))
                print('and',len( eringene&tningene))
                result = gp.profile( organism = 'hsapiens', query = ingene, sources = ['GO:MF','GO:CC','GO:BP'],  no_evidences = False, no_iea = True)
                leng = 0
                while len(result.index) > leng + 50:
                    print(result.loc[ leng:leng+50,['name','p_value','intersections']])
                    leng = leng + 50
                print(result.loc[ leng:leng+50,['name','p_value','intersections']])
    for i in molecular:
        print(i)
        print(freq.loc[(freq.loc[:,freq.columns.isin(molecular)].sum(axis = 1) == 1)&(freq.loc[:,i] == 1),i].sum())

    print(freq.loc[:,'sel'].sum())
if False:
    print( freq.loc[ (freq.loc[:,'ER+_amp_freq'] >= 0.15)&(freq.loc[:,'HER2+_amp_freq'] >= 0.15)&(freq.loc[:,'TNBC_amp_freq'] >= 0.15),:])
    #print( freq.loc[ (freq.loc[:,'ER+'] == 1)&(freq.loc[:,'HER2+'] == 1)&(freq.loc[:,'TNBC'] == 1),:])
    tmp = freq.loc[freq.loc[:,'sel'] == 1,(~freq.columns.str.contains('sel'))]
    tmp.to_csv('data/sum_sel_gdc.csv')
    if True:
        freq.loc[:,'diff_sel'] = 0
        #freq.loc[ (freq.loc[:,'avg_log2FoldChange'] >= 1.5) & (freq.loc[:,'n_tissues_sig'] >= 12), 'diff_sel'] = 1
        freq.loc[ (freq.loc[:,'avg_log2FoldChange'] >= 1), 'diff_sel'] = 1
        print(freq.loc[:,'diff_sel'].sum())
    freq = freq.loc[ (freq.loc[:,'diff_sel'] == 1)&(freq.loc[:,'sel'] == 1),:]
    #freq = freq.loc[ (freq.loc[:,'membrane_protein'] == 'Yes')|(freq.loc[:,'plasma_protein']=='Yes'),:]
    freq = freq.loc[ (freq.loc[:,'membrane_protein_combined'] == 'Yes'),:]
    #tmp = freq.loc[:,~freq.columns.str.contains('sel')]
    tmp = freq
    print(freq)
#'''
'''
import pandas as pd

# 1. Load the two files
df_1702 = pd.read_csv('Amp_1702.txt', sep='\t', index_col=0)
df_fisher = pd.read_csv('data/new_final_fisher.txt', sep='\t', index_col=0)

# 2. Identify the new genes
genes_1702 = set(df_1702.index)
genes_fisher = set(df_fisher.index)
new_genes = genes_fisher - genes_1702

print("=== New Genes Identified ===")
print(f"Found {len(new_genes)} new genes: {sorted(list(new_genes))}\n")

# 3. Target columns to explicitly overwrite from the fisher file
overwrite_cols = [
    'sel', 'All_freq_prioritize', 'All_exp_prioritize',
    'ER+_freq_prioritize', 'ER+_exp_prioritize',
    'HER2+_freq_prioritize', 'HER2+_exp_prioritize',
    'TNBC_freq_prioritize', 'TNBC_exp_prioritize',
    'ER+_amp_FDR', 'HER2+_amp_FDR', 'TNBC_amp_FDR'
]
overwrite_cols = [col for col in overwrite_cols if col in df_fisher.columns]

# 4. Clear old column versions from the baseline file to prevent conflicts
df_1702_clean = df_1702.drop(columns=[col for col in overwrite_cols if col in df_1702.columns], errors='ignore')

# 5. Isolate the updated metrics matrix
df_fisher_metrics = df_fisher[overwrite_cols]

# 6. Merge baseline with the new metrics
# This creates rows for the new genes, leaving non-metric columns as NaN initially
df_1707 = df_1702_clean.join(df_fisher_metrics, how='outer')

# 7. FILL IN THE NEW GENES' VALUES FROM THE FISHER FILE
# combine_first looks at any remaining NaN values (like chromosome, position, etc. for the new genes)
# and patches them using the complete information available inside df_fisher
df_1707 = df_1707.combine_first(df_fisher)

# Reorder columns to match original Amp_1702 structure as closely as possible
final_cols = [col for col in df_1702.columns if col in df_1707.columns] + [col for col in overwrite_cols if col not in df_1702.columns]
df_1707 = df_1707.reindex(columns=final_cols)

# 8. Export to the updated tracking sheet
df_1707.to_csv('Amp_1707.txt', sep='\t')
print(f"Successfully generated 'Amp_1707.txt'. Total shape: {df_1707.shape}")
#'''

'''
import pandas as pd

# 1. Load the files
df_1702 = pd.read_csv('Amp_1702.txt', sep='\t', index_col=0)
df_fisher = pd.read_csv('data/new_final_fisher.txt', sep='\t', index_col=0)

# 2. Identify which genes are being removed from the historical list
# These are genes that were in the old file but are missing from the new Fisher file
genes_1702 = set(df_1702.index)
genes_fisher = set(df_fisher.index)
removed_genes = genes_1702 - genes_fisher

print("=== REMOVED GENES REPORT ===")
print(f"Number of historical genes removed (not present in Fisher final file): {len(removed_genes)}")
print("\n--- List of Removed Genes ---")
print(sorted(list(removed_genes)))

# 3. Clean and prepare columns to be overwritten from the fisher data
overwrite_cols = [
    'sel', 'All_freq_prioritize', 'All_exp_prioritize',
    'ER+_freq_prioritize', 'ER+_exp_prioritize',
    'HER2+_freq_prioritize', 'HER2+_exp_prioritize',
    'TNBC_freq_prioritize', 'TNBC_exp_prioritize',
    'ER+_amp_FDR', 'HER2+_amp_FDR', 'TNBC_amp_FDR'
]
overwrite_cols = [col for col in overwrite_cols if col in df_fisher.columns]

# Drop the columns from the baseline dataframe to prevent overlap artifacts
df_1702_clean = df_1702.drop(columns=[col for col in overwrite_cols if col in df_1702.columns], errors='ignore')

# 4. Filter the baseline dataframe to ONLY keep genes that are in the new Fisher file
# (This acts as a 'right join' constraint)
df_1702_filtered = df_1702_clean.loc[df_1702_clean.index.isin(genes_fisher)]

# 5. Merge the datasets safely
df_1707 = df_1702_filtered.join(df_fisher[overwrite_cols], how='left')

# 6. Backfill any missing annotation details for newly introduced genes using combine_first
df_1707 = df_1707.combine_first(df_fisher)

# 7. Ensure column alignment matches the historical formatting layout
final_cols = [col for col in df_1702.columns if col in df_1707.columns] + [col for col in overwrite_cols if col not in df_1702.columns]
df_1707 = df_1707.reindex(columns=final_cols)

# 8. Save clean tracking table
df_1707.to_csv('Amp_1707.txt', sep='\t')
print(f"\nSuccessfully generated clean 'Amp_1707.txt'. Total shape: {df_1707.shape}")
print(f"The number of genes matches your final fisher matrix exactly: {len(df_1707)}")
#'''

'''
import pandas as pd

# 1. Load the two files
df_1702 = pd.read_csv('Amp_1702_final.txt', sep='\t', index_col=0)
df_1707 = pd.read_csv('data/amp1707.txt', sep='\t', index_col=0)

fdr_cols = ['ER+_amp_FDR', 'HER2+_amp_FDR', 'TNBC_amp_FDR']
annotation_col = 'subtype_amplification_with_increased_expression'

# 2. Track down the gene anomalies
genes_1702 = set(df_1702.index)
genes_1707 = set(df_1707.index)

removed_genes = genes_1702 - genes_1707
added_genes = genes_1707 - genes_1702

print("=== Gene Population Discrepancy ===")
print(f"Gene being removed from 1702 (absent in 1707): {sorted(list(removed_genes))}")
print(f"Additional genes being added and backfilled from 1707: {sorted(list(added_genes))}\n")

# 3. Step 1: Remove the missing gene from the 1702 baseline dataframe
df_1702_filtered = df_1702.drop(index=list(removed_genes), errors='ignore')

# 4. Map the annotation header from the new run file
if 'annotation' in df_1707.columns:
    df_1707_mapped = df_1707.rename(columns={'annotation': annotation_col})
else:
    df_1707_mapped = df_1707.copy()

# 5. Drop target overwrite columns to prevent duplication suffix artifacts (_x, _y)
overwrite_cols = fdr_cols + [annotation_col]
overwrite_cols = [col for col in overwrite_cols if col in df_1707_mapped.columns]
df_1702_clean = df_1702_filtered.drop(columns=[col for col in overwrite_cols if col in df_1702_filtered.columns], errors='ignore')

# 6. Step 2: Use an outer join so additional genes from 1707 are added as new rows
df_merged = df_1702_clean.join(df_1707_mapped[overwrite_cols], how='outer')

# 7. Step 3: Backfill annotations/coordinates for the additional genes using combine_first
# For the new genes, any NaN column from 1702 will be populated with data directly from amp1707
df_final = df_merged.combine_first(df_1707_mapped)

# 8. Clean up column alignment to preserve original 1702 layout
final_cols = [col for col in df_1702.columns if col in df_final.columns] + [col for col in overwrite_cols if col not in df_1702.columns]
df_final = df_final.reindex(columns=final_cols)

# 9. Save final integrated baseline tracker
df_final.to_csv('Amp_1707_merged_final.txt', sep='\t')

print("=== Complete Validation ===")
print(f"Original Amp_1702_final Genes: {len(df_1702)}")
print(f"New data/amp1707 Genes:        {len(df_1707)}")
print(f"Final Cleaned & Merged Genes:  {len(df_final)}")
print("Output matrix successfully saved to: 'Amp_1707_merged.txt'")
#'''

'''
import pandas as pd

# 1. Load your newly generated file
# (Adjust file path/name if you saved it as 'Amp_1707_merged.txt')
df = pd.read_csv('Amp_1707.txt', sep='\t', index_col=0)

# Target prioritization columns
prioritize_cols = [
    'ER+_freq_prioritize',
    'HER2+_freq_prioritize',
    'TNBC_freq_prioritize',
    'All_freq_prioritize'
]

print("=== Expression Prioritization Counts (== 1) ===")
for col in prioritize_cols:
    if col in df.columns:
        # pd.to_numeric cleanly strips any hidden text prefixes like "'"
        # and converts the data into math-ready integers/floats
        numeric_series = pd.to_numeric(df[col], errors='coerce').fillna(0)

        # Count how many genes equal 1
        active_count = (numeric_series == 1).sum()
        print(f"{col}: {active_count} genes")
    else:
        print(f"{col}: Column missing from this file metadata layout.")
#'''

import pandas as pd

# 1. Load the source file containing prioritization columns
df_source = pd.read_csv("Amp_1707.txt", sep="\t", index_col=0)

# Map prioritization columns to clean subtype labels
prioritize_map = {
    "ER+_freq_prioritize": "ER+",
    "HER2+_freq_prioritize": "HER2+",
    "TNBC_freq_prioritize": "TNBC",
    "All_freq_prioritize": "All",
}

# Clean numeric values
for col in prioritize_map:
    if col in df_source.columns:
        df_source[col] = (
            pd.to_numeric(df_source[col], errors="coerce").fillna(0).astype(int)
        )


# 2. Build the freq_prioritization annotation per gene
# Exclude "All" if any specific subtype is active
def get_prioritized_subtypes(row):
    active = [
        label
        for col, label in prioritize_map.items()
        if col in row and row[col] == 1
    ]

    # If any specific subtype is present, drop "All"
    specific_subtypes = [s for s in active if s != "All"]
    if specific_subtypes:
        active = specific_subtypes

    return ";".join(active) if active else "None"


df_source["freq_prioritization"] = df_source.apply(
    get_prioritized_subtypes, axis=1
)

# 3. Print breakdown per chromosome for each subtype
print("=== Prioritized Gene Counts per Chromosome ===")
if "chromosome_name" in df_source.columns:
    for col, label in prioritize_map.items():
        if col in df_source.columns:
            subset = df_source[df_source[col] == 1]
            chrom_counts = (
                subset["chromosome_name"]
                .astype(str)
                .value_counts()
                .sort_index()
            )
            print(f"\n--- Subtype: {label} (Total: {len(subset)}) ---")
            print(chrom_counts.to_string())
else:
    print("Warning: 'chromosome_name' column not found in Amp_1707.txt.")

# 4. Load target file and merge the annotation column
df_target = pd.read_csv("Amp_1707_merged_final.txt", sep="\t", index_col=0)

# Map freq_prioritization using index/gene identifiers
df_target["freq_prioritization"] = df_target.index.map(
    df_source["freq_prioritization"]
).fillna("None")

# 5. Save output
output_file = "Amp_1707_20260831.txt"
df_target.to_csv(output_file, sep="\t")
print(f"\nSuccessfully saved annotated dataset to {output_file}")
#'''
