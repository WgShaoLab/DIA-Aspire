import csv
import os,sys
import re
import numpy as np
import pandas as pd

#print("please input 'argv1: inputname of sptxt','argv2: number of top fragments' ")

def remove_parent1(text):
    return re.sub('\[.*?\]', '', text)

def submod(text):
    t1 = re.sub('M\[147\]', 'M(UniMod:35)', text)
    t2 = re.sub('S\[167\]', 'S(UniMod:21)', t1)
    t3 = re.sub('T\[181\]', 'T(UniMod:21)', t2)
    t4 = re.sub('Y\[243\]', 'Y(UniMod:21)', t3)
    t5 = re.sub('N\[115\]', 'N(UniMod:7)', t4)
    t6 = re.sub('Q\[129\]', 'Q(UniMod:7)', t5)
    t7 = re.sub('S\[129\]', 'Q(UniMod:1)', t6)
    t8 = re.sub('n\[43\]', '(UniMod:1)', t7)
    return t8 

def extract(text):
    pattern = r'\((.*?)\)'  # 匹配()之间的内容
    result = re.findall(pattern, text)
    return result

def remove_parent(text):
    return re.sub('\(.*?\)', '', text)

def getdata(da,flg):

    da3 = da.copy()
    if flg==1:
        da3 = da3[da3['FragmentType'].str.contains('[\?]')==False]
    elif flg==0:
        da3 = da3[da3['FragmentType'].str.contains('[\?]')==True]
    else:
        da3 = da.copy()

    da3 = da3[da3['ModifiedPeptide'].str.contains('\[')==False]
    da3['IonMobility'] = 0
    da3['FragmentMZ'] = da3['FragmentMZ'].astype(float)
    da3['RelativeIntensity'] = da3['RelativeIntensity'].astype(float)
    da3['PrecursorMZ'] = da3['PrecursorMZ'].astype(float)
    da3['PrecursorCharge'] = da3['PrecursorCharge'].astype(int)
    da3['FragmentCharge'] = da3['FragmentCharge'].astype(int)
    da3 = da3.rename(columns={'rt':'iRT'})

    da3['ions'] = da3['ModifiedPeptide'] + da3['PrecursorCharge'].astype(str)
    da4 = da3.sort_values(by=['ModifiedPeptide','PrecursorCharge','RelativeIntensity'],ascending=[True,True,False])
    return da4

def get_final(da,n):
    dax = pd.concat([group.head(n) for _, group in da.groupby(da['ions'])])
    
    col2 = ['PrecursorMZ','FragmentMZ','RelativeIntensity','iRT','Protein_name','ModifiedPeptide',\
           'StrippedPeptide','FragmentType','FragmentNumber','PrecursorCharge','FragmentCharge',\
           'uniprot_id','Tr_recalibrated','shared', 'decoy']
    dax2 = dax.loc[:,col2]
                   
    dax2.columns = ['PrecursorMz','ProductMz','LibraryIntensity','iRT','Protein_name','ModifiedPeptide',\
           'StrippedPeptide','FragmentType','FragmentNumber','PrecursorCharge','FragmentCharge',\
           'uniprot_id','Tr_recalibrated','shared', 'decoy']
    return dax2

def get_final2(da):
    dax = da.copy()
    col2 = ['PrecursorMZ','FragmentMZ','RelativeIntensity','iRT','Protein_name','ModifiedPeptide',\
           'StrippedPeptide','FragmentType','FragmentNumber','PrecursorCharge','FragmentCharge',\
           'uniprot_id','Tr_recalibrated','shared', 'decoy']
    dax2 = dax.loc[:,col2]
                   
    dax2.columns = ['PrecursorMz','ProductMz','LibraryIntensity','iRT','Protein_name','ModifiedPeptide',\
           'StrippedPeptide','FragmentType','FragmentNumber','PrecursorCharge','FragmentCharge',\
           'uniprot_id','Tr_recalibrated','shared', 'decoy']
    return dax2


#inp = sys.argv[1] 
def convert_sptxt2tsv(inp):
    num1 = 12

    f = open(inp,'r')
    spts = f.readlines()
    f.close()

    ionsname=[]
    preMZ = []
    bys = []
    retime = []
    npeaks = []
    prot = []

    for spt in spts:
        if re.match('^Name',spt):
            ion = re.sub('Name: ','',spt)
            ionsname.append(ion[0:-1])
        elif re.match('^PrecursorMZ',spt):
            pmz = re.sub('PrecursorMZ: ','',spt)
            preMZ.append(pmz[0:-1])
        elif re.match('^Comment',spt):
            regx = '(?<=RetentionTime=).[0-9|,|.]*'
            a = re.findall(regx,spt)
            retime.append(a)
            
            match = re.search(r'Protein=(.*?)\s', spt)
            prot.append(match.group(1))
            
        elif re.match('^NumPeaks',spt):
            npk = re.sub('NumPeaks: ','',spt)
            npeaks.append(npk[0:-1])
        elif re.match('^\d',spt):
            bys.append(spt[0:-1])

    peaks = pd.DataFrame(npeaks)

    RT = pd.DataFrame(retime)

    df0 = pd.DataFrame({'peptide':ionsname,'PrecursorMZ':preMZ,'Protein_name':prot})

    byx = pd.DataFrame(bys)
    byions = byx.iloc[:,0].str.split('\t',expand=True)
    byions1 = byions.iloc[:,[0,1,2]]
    byions1.columns = ['FragmentMZ','RelativeIntensity','Fragment']
    byions2 = byions1['Fragment'].str.split(',',expand=True)

    RT.columns = ['rt']
    RT1 = RT['rt'].str.split(',',expand=True)

    for col in RT1.columns:
        RT1[col] = RT1[col].astype(float)  

    df = df0.copy()
    df['iRT'] = RT1.median(axis=1)

    nump = peaks.iloc[:,0].tolist()
    df_r = df.reindex(df.index.repeat(nump)).reset_index(drop=True)

    da0 = pd.concat([byions1,df_r],axis=1)

    da = da0.copy()
    da['Fragment'] = byions2.iloc[:,0]
    da['Fragment'] = da['Fragment'].apply(remove_parent1)
    da['modpep'] = da['peptide'].apply(submod) #更改mod格式
    da = da[da['Fragment'].str.contains('[i|p|\+]')==False]

    byions3 = da['Fragment'].str.split('/',expand=True)
    byions3.columns = ['ion','error']
    byions3['FragmentType'] = byions3['ion'].str[0]

    byions4 = byions3['ion'].str.split('^',expand=True)
    byions4.columns = ['ion1','charge']
    byions4['FragmentCharge'] = byions4['charge']
    byions4['FragmentType'] = byions4['ion1']
    byions4['FragmentCharge'].fillna(1,inplace=True)

    byions5 = byions4['ion1'].str.split('-',expand=True)
    byions5.columns = ['ion2','loss']
    byions5['FragmentNumber'] = byions5['ion2'].str[1:]
    byions5.loc[byions5['ion2'].str.contains('[I|m|\?]')==True,'FragmentNumber'] = np.nan

    da1 = da.copy()
    da1['ModifiedPeptide'] = da1['modpep'].str[0:-2]

    da1['PrecursorCharge'] = da1['peptide'].str[-1]
    #da1['LabeledPeptide'] = da1['ModifiedPeptide'] 
    da1['StrippedPeptide'] = da1['peptide'].str[0:-2].apply(remove_parent1)

    da1['FragmentNumber'] = byions5['FragmentNumber']
    da1['FragmentType'] =  byions4['FragmentType']
    da1['FragmentCharge'] = byions4['FragmentCharge']

    da2 = da1.drop('Fragment',axis=1)
    da2 = da2.drop('peptide',axis=1)
    da2 = da2.drop('modpep',axis=1)

    da2x = da2.copy()
    da2x = da2x[da2x['Protein_name'].str.contains('^\d\/DECOY')==False]
    da2x = da2x[da2x['Protein_name'].str.contains('^\d\/rev')==False]
    da2x['shared'] = 'TRUE'
    da2x['decoy']='FALSE'
    da2x.loc[da2x['Protein_name'].str.contains('^1/')==True,'shared'] = 'FALSE'
    da2x['uniprot_id'] = da2x['Protein_name']
    da2x['Tr_recalibrated'] = da2x['iRT']

    da5 = getdata(da2x,1)

    daout5x = da5[da5['FragmentType'].str.contains('[y|b|\-|a|m]')]
    daout5 = get_final(daout5x,num1)

    # outname5 = 'top12_bynam_library.tsv'
    return daout5
