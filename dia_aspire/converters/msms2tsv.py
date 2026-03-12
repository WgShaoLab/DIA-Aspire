import pandas as pd
import numpy as np
import re

def removeline(txt):
    return re.sub('_','',txt)

def submod(text):
    t1 = re.sub('\(ox\)', '(UniMod:35)', text)
    t2 = re.sub('\(ac\)', '(UniMod:1)', t1)
    t3 = re.sub('\(ph\)', '(Unimod:21)',t2)
    t4 = re.sub('\(de\)', '(Unimod:7)',t3)
    return t4 

def convertdata1(dax):
    
    damass = dax['Masses'].str.split(';',expand=True)
    col = [] 
    for i in range(damass.shape[1]):
        coli = 'm' + str(i)
        col.append(coli)
    damass.columns = col
    damass['ions'] = dax['ions']
    damass
    dame = []
    for i in range(len(col)):
        col_s = ['ions']
        col_s.append(col[i])
    #     print(col_s)
        das = damass[col_s]
        das.columns = ['ions','FragmentMz']
        dame.append(das)
    damer1 = pd.concat(dame)
    damer1 = damer1[~damer1['FragmentMz'].isnull()]
    damer1.index = np.arange(damer1.shape[0])
    return damer1   

def convertdata2(dax):
    
    damass = dax['Intensities'].str.split(';',expand=True)
    col = [] 
    for i in range(damass.shape[1]):
        coli = 'm' + str(i)
        col.append(coli)
    damass.columns = col
    damass['ions'] = dax['ions']
    damass
    dame = []
    for i in range(len(col)):
        col_s = ['ions']
        col_s.append(col[i])
    #     print(col_s)
        das = damass[col_s]
        das.columns = ['ions','Intensity']
        dame.append(das)
    damer1 = pd.concat(dame)
    damer1 = damer1[~damer1['Intensity'].isnull()]
    damer1.index = np.arange(damer1.shape[0])
    return damer1  

def convertdata3(dax):
    
    damass = dax['Matches'].str.split(';',expand=True)
    col = [] 
    for i in range(damass.shape[1]):
        coli = 'm' + str(i)
        col.append(coli)
    damass.columns = col
    damass['ions'] = dax['ions']
    damass
    dame = []
    for i in range(len(col)):
        col_s = ['ions']
        col_s.append(col[i])
    #     print(col_s)
        das = damass[col_s]
        das.columns = ['ions','FragmentType']
        dame.append(das)
    damer1 = pd.concat(dame)
    damer1 = damer1[~damer1['FragmentType'].isnull()]
    damer1['FragmentTypes'] = damer1['FragmentType'].astype(str)
    damer1.index = np.arange(damer1.shape[0])
    return damer1  

def extract_charge(text):
    pattern = r'\((.*?)\)'  # 匹配()之间的内容
    result = re.findall(pattern, text)
    if result:
        results = int(result[0][0])
    else:
        results = 1
    return results

def remove_charge(text):
    return re.sub('\(.*?\)', '', text)

def get_final(da,n):
    dax = pd.concat([group.head(n) for _, group in da.groupby(da['ions'])])
    return dax

def convert_msms2tsv(inp):
    num1 = 12
    
    da = pd.read_csv(inp,sep='\t')
    da = da[da['Reverse']!='+']
    col = ['Sequence','Modified sequence','Protein Names','Charge','m/z','Retention time','Score',
           'Intensities','Matches','Masses']
    da1 = da[col]

    da2 = da1.sort_values(by=['Modified sequence','Charge','Score'],ascending=[False,True,False])

    da3 = da2.drop_duplicates(subset=['Modified sequence','Charge'],keep='first')

    da4 = da3.copy()
    da4['Modified sequence'] = da4['Modified sequence'].apply(submod)
    da4['modifiedpep'] = da4['Modified sequence'].apply(removeline)
    da4['ions'] = da4['modifiedpep'] + da4['Charge'].astype(str)

    da4_tmp1 = da4[da4['modifiedpep'].str.contains('UniMod')]
    da4_tmp2 = da4[da4['modifiedpep'].str.contains('\(')==False]
    da4 = pd.concat([da4_tmp2,da4_tmp1])

    damerge1 = convertdata1(da4)

    damerge2 = convertdata2(da4)

    damerge3 = convertdata3(da4)

    damerge4 = damerge3.copy()
    damerge4['FragmentCharge'] = damerge4['FragmentType'].apply(extract_charge)
    damerge4['Fragmentx'] = damerge4['FragmentType'].apply(remove_charge)

    byions = damerge4['Fragmentx'].str.split('-',expand=True)
    byions.columns = ['ion2','loss']
    byions['FragmentNumber'] = byions['ion2'].str[1:]

    damerge5 = damerge4.copy()
    damerge5['FragmentNumber'] = pd.to_numeric(byions5['FragmentNumber'], errors='coerce')
    damerge5 = damerge5.drop(['Fragmentx'],axis=1)

    damerge6 = damerge5.copy()
    damerge6['ProductMz'] = damerge1['FragmentMz']
    damerge6['LibraryIntensity'] = damerge2['Intensity']

    da4.columns
    da5 = da4[['Sequence', 'Protein Names', 'Charge', 'm/z',
           'Retention time','modifiedpep', 'ions']]
    da5 = da5.copy()
    da5.columns = ['StrippedPeptide','Protein_name','PrecursorCharge','PrecursorMz','iRT','ModifiedPeptide','ions']
    da5['uniprot_id'] = da5['Protein_name']
    da5['Tr_recalibrated'] = da5['iRT']
    da5['decoy'] = 'FALSE'
    da5['shared'] = 'TRUE'
    da5.loc[da5['Protein_name'].str.contains(';')==True,'shared'] = 'FALSE'

    dataf = pd.merge(da5,damerge6,on='ions',how='outer')

    dataf1 = dataf.sort_values(by=['ModifiedPeptide','PrecursorCharge','LibraryIntensity'],ascending=[True,True,False])
    dataf1 = dataf1[dataf1['FragmentNumber']!=1]

    datafinal = get_final(dataf1,num1)
    col = ['PrecursorMz','ProductMz','LibraryIntensity','iRT','Protein_name','ModifiedPeptide',\
               'StrippedPeptide','FragmentType','FragmentNumber','PrecursorCharge','FragmentCharge',\
               'uniprot_id','Tr_recalibrated','shared', 'decoy']
    datafinal = datafinal[col]
    return datafinal
