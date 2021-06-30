#!/usr/bin/env python
# coding: utf-8

# In[1309]:


#ASIST program for phenotype based on Antibiotics profile
# create a profile based on selected antibiotics only
# rakesh4osdd@gmail.com, 14-June-2021


# In[1]:


import pandas as pd
import sys
import os
from collections import Counter


# In[176]:


input_file=sys.argv[1]
output_file=sys.argv[2]
#input_file='test-data/asist_input.csv'
#output_file='test-data/asist_output.csv'


# In[177]:


# strain_profile to phenotype condition
def s_phen(sus,res,intm,na,pb_sus):
    if (sus>0 and res==0 and na>=0):
        #print('Possible Susceptible')
        phen='Possible Susceptible'
    elif (sus>=0 and 3<=res<7 and na>=0 and pb_sus==0):
        #print('Possible MDR')
        phen='Possible MDR'
    elif (sus>=0 and 7<=res<9 and na>=0 and pb_sus==0):
        #print('Possible XDR')
        phen='Possible XDR'
    #special cases
    elif (sus>=1 and res>0 and na>=0 and pb_sus==1):
        #print('Possible XDR')
        phen='Possible XDR'
    #special cases
    elif (sus>0 and res==9 and na>=0):
        #print('Possible XDR')
        phen='Possible XDR'
    elif (sus==0 and res==9 and na>=0):
        #print('Possible TDR')
        phen='Possible TDR'
    else:
        #print('Strain could not be classified')
        phen='Strain could not be classified ('+ str(intm)+' | ' + str(na) +')'
    return(phen)

#print(s_phen(1,9,0,0))


# In[178]:


# define Antibiotic groups as per antibiotic of CLSI breakpoints MIC
#Aminoglycoside
cat1=['Amikacin','Tobramycin','Gentamycin','Netilmicin']
#Beta-lactams- Carbapenems
cat2=['Imipenem','Meropenam','Doripenem']
#Fluoroquinolone
cat3=['Ciprofloxacin','Levofloxacin']
#Beta-lactam inhibitor
cat4=['Piperacillin/tazobactam','Ticarcillin/clavulanicacid']
#Cephalosporin
cat5=['Cefotaxime','Ceftriaxone','Ceftazidime','Cefepime']
#Sulfonamides
cat6=['Trimethoprim/sulfamethoxazole']
#Penicillins/beta-lactamase
cat7=['Ampicillin/sulbactam']
#Polymyxins
cat8=['Colistin','Polymyxinb']
#Tetracycline
cat9=['Tetracycline','Doxicycline','Minocycline']

def s_profiler(pd_series):
    #print(type(pd_series),'\n', pd_series)
    #create a dictionary of dataframe series
    cats={'s1':cat1,'s2':cat2,'s3':cat3,'s4':cat4,'s5':cat5,'s6':cat6,'s7':cat7,'s8':cat8,'s9':cat9}
    # find the antibiotics name in input series
    for cat in cats:
        #print(cats[cat])
        cats[cat]=pd_series.filter(cats[cat])
        #print(cats[cat])
    #define res,sus,intm,na,pb_sus
    res=0
    sus=0
    intm=0
    na=0
    pb_sus=0
    # special case of 'Polymyxin b' for its value
    if 'Polymyxinb' in pd_series:
        ctp=cats['s8']['Polymyxinb'].strip().lower()
        if ctp == 'susceptible':
            pb_sus=1
        #print((ctp,p_sus))
    # check all categories
    for cat in cats:
        #ctp=cats['s8'].iloc[i:i+1].stack().value_counts().to_dict()
        #print(ctp)
        # Pandas series
        ct=cats[cat].value_counts().to_dict()
        #print(ct)
        # remove whitespace and convert to lowercase words
        ct =  {k.strip().lower(): v for k, v in ct.items()}
        #print(ct)
        k=Counter(ct)
        #j=Counter(ct)+Counter(j)
        #print(j)
        # category wise marking
        if k['resistant']>=1:
            res=res+1
        if k['susceptible']>=1:
            sus=sus+1
        if k['intermediate']>=1:
            intm=intm+1
        if k['na']>=1:
            na=na+1
    #print(sus,res,intm,na,pb_sus)
    #print(s_phen(sus,res,intm,na,pb_sus))
    return(s_phen(sus,res,intm,na,pb_sus))


# In[179]:


#input_file='input2.csv_table.csv'
#output_file=input_file+'_output.txt'
strain_profile=pd.read_csv(input_file, sep=',',na_filter=False,skipinitialspace = True)


# In[180]:


old_strain_name=strain_profile.columns[0]
new_strain_name=old_strain_name.capitalize().strip().replace(' ', '')


# In[181]:


# make header capitalization, remove leading,lagging, and multiple whitespace for comparision
strain_profile.columns=strain_profile.columns.str.capitalize().str.strip().str.replace('\s+', '', regex=True)
#print(strain_profile.columns)
#strain_profile.head()
#strain_profile.columns


# In[182]:


# add new column in dataframe on second position
strain_profile.insert(1, 'Strain phenotype','')
#strain_profile.head()


# In[183]:


strain_profile['Strain phenotype'] = strain_profile.apply(lambda x: (s_profiler(x)), axis=1)


# In[184]:


#strain_profile.head()


# In[185]:


#rename headers for old name
strain_profile=strain_profile.rename(columns = {new_strain_name:old_strain_name, 'Ticarcillin/clavulanicacid':'Ticarcillin/ clavulanic acid','Piperacillin/tazobactam':'Piperacillin/ tazobactam','Trimethoprim/sulfamethoxazole': 'Trimethoprim/ sulfamethoxazole','Ampicillin/sulbactam':'Ampicillin/ sulbactam', 'Polymyxinb': 'Polymyxin B'} )


# In[186]:


#strain_profile.columns


# In[187]:


#strain_profile


# In[188]:


strain_profile.to_csv(output_file,na_rep='NA',index=False)


# In[189]:


# Open a file with access mode 'a'
with open(output_file, "a") as file_object:
    # Append 'hello' at the end of file
    file_object.write("Note: \n1. 'MDR': Multidrug-resistant, 'XDR': Extensively drug-resistant, 'TDR':totally drug resistant, NA': Data Not Available.\n2. 'Strain could not be classified' numbers follow the format as ('Number of antibiotics categories count as Intermediate' | 'Number of antibiotics categories count as NA')")


# In[ ]:




