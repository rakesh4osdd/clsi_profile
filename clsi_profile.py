#!/usr/bin/env python
# coding: utf-8

# In[206]:


# ASIST module2 | map AST result to the CLSI breakporints with combination antibiotics
# By rakesh4osdd@gmail.com, 06-Jun-2021
import pandas as pd
import re
import sys


# In[207]:


#print(pd.__version__, re.__version__)


# In[208]:


# compare two MIC value strings
def check_mic(mic1,mic2,mic_type):
    #print(mic1,mic2,mic_type)
    try:
        if '/' in mic1:
            m1a = mic1.split('/')[0]
            m1b = mic1.split('/')[1]
            if float(m1a)==0 or float(m1b)==0:
                strain_type='Strain could not be classified'
                return(strain_type)          
        elif '/' in mic2:
            m1a = mic1
            if float(m1a)==0:
                strain_type='Strain could not be classified'
                return(strain_type)            
            m1b = '1'
        elif float(mic1)==0:
            strain_type='Strain could not be classified'
            return(strain_type)
        else:
            m1a = mic1
            
        if '-' in mic2:
            m2a = mic2.split('-')[0]
            m2b = mic2.split('-')[1]           
         
    except ValueError:
        strain_type='Strain could not be classified' 
        return(strain_type)
    try:
        if '-' in mic2 and mic_type == 'i':   # for intermediate only
            if '/' in mic2:
                m2a = mic2.split('-')[0].split('/')[0]
                m2b = mic2.split('-')[0].split('/')[1]
                m2aa = mic2.split('-')[1].split('/')[0]
                m2bb = mic2.split('-')[1].split('/')[1]
                if (float(m2aa)>=float(m1a)>=float(m2a) and float(m2bb)>=float(m1b)>=float(m2b)):
                    #print('intermediate')
                    m_type='Intermediate'
                else:
                    #print('not define')
                    m_type='Strain could not be classified'
            else:
                m2a = mic2.split('-')[0]
                m2b = mic2.split('-')[1] 
                if (float(m2b)>=float(m1a)>=float(m2a)):
                    #print('intermediate')
                    m_type='Intermediate'
                else:
                    #print('not define')
                    m_type='Strain could not be classified'                
            #print (m1a,m1b,m2a,m2b,m2aa,m2bb)
        elif '/' in mic2:
            m2a = mic2.split('/')[0]
            m2b = mic2.split('/')[1]
            #print(m1a,m1b,m2a,m2b,mic_type)
            if (mic_type=='s' and (float(m1a)<=float(m2a) and float(m1b)<=float(m2b))):
                m_type='Susceptible'
            elif (mic_type=='r' and (float(m1a)>=float(m2a) and float(m1b)>=float(m2b))):
                m_type='Resistant'
            elif (mic_type=='i' and (float(m1a)==float(m2a) and float(m1b)==float(m2b))):
                m_type='Intermediate'
            else:
                m_type='Strain could not be classified'
        elif '-' in mic2:
                m_type='Strain could not be classified'
        else:
            m2a=mic2
            if (mic_type=='s' and (float(m1a)<=float(m2a))):
                m_type='Susceptible'
            elif (mic_type=='r' and (float(m1a)>=float(m2a))):
                m_type='Resistant'
            elif (mic_type=='i' and (float(m1a)==float(m2a))):
                m_type='Intermediate'
            else:
                m_type='Strain could not be classified-1'        
    except IndexError:
        strain_type='Strain could not be classified-2' 
        return(strain_type)
    
    return(m_type)

#check_mic('65','32-64','i')


# In[209]:


# compare MIC value in pandas list
def sus_res_int(mic):
    #print(mic)
    o_mic = mic[0].replace(' ', '')
    s_mic = mic[1].replace(' ', '')
    r_mic = mic[2].replace(' ', '')
    i_mic = mic[3].replace(' ', '')
    try:
        if check_mic(o_mic,s_mic,'s')=='Susceptible':
            strain_type='Susceptible'
        elif check_mic(o_mic,r_mic,'r')=='Resistant':
            strain_type='Resistant'
        elif check_mic(o_mic,i_mic,'i')=='Intermediate':
            strain_type='Intermediate'                    
        else:
            strain_type='Strain could not be classified'
    except ValueError:
        strain_type='Strain could not be classified'            
    return(strain_type)

#mic=['128','16/4','128/4','32/4-64/4']
#sus_res_int(mic)


# In[210]:


# for input argument
input_user = sys.argv[1]
input_clsi = sys.argv[2]
output_table = sys.argv[3]


# In[211]:

"""
input_user='input.csv'
input_clsi='clsi.csv'
output_profile=input_user+'_profile.csv'
output_table=input_user+'_table.csv'
"""

# In[212]:


# read user AST data with selected 3 columns
strain_mic=pd.read_csv(input_user, sep=',', usecols =['Strain name', 'Antibiotics', 'MIC'],na_filter=False)


# In[213]:


clsi_bp=pd.read_csv(input_clsi,sep=',')


# In[214]:


#clsi_bp
#strain_mic


# In[215]:


# convert MIC to numbers sMIC, rMIC
clsi_bp['s_mic'] =clsi_bp[['Susceptible']].applymap(lambda x: (re.sub(r'[^0-9.\/-]', '', x)))
clsi_bp['r_mic'] =clsi_bp[['Resistant']].applymap(lambda x: (re.sub(r'[^0-9.\/-]', '', x)))
clsi_bp['i_mic'] = clsi_bp[['Intermediate']].applymap(lambda x: (re.sub(r'[^0-9.\/-]', '', x)))


# In[216]:


#clsi_bp['i_mic'] = clsi_bp[['Intermediate']].applymap(lambda x: (re.sub(r'[^0-9.\/-]', '', x)))


# In[217]:


# Read only numbers in MIC values
#try:
strain_mic['o_mic']=strain_mic[['MIC']].applymap(lambda x: (re.sub(r'[^0-9.\/]','', x)))
#except TypeError:
#    print('Waring: Error in MIC value')


# In[218]:


#strain_mic


# In[219]:


# capitalize each Antibiotic Name for comparision with removing whitespace
strain_mic['Antibiotics']=strain_mic['Antibiotics'].str.capitalize().str.replace(" ","")
clsi_bp['Antibiotics']=clsi_bp['Antibiotics'].str.capitalize().str.replace(" ","")


# In[220]:


#compare CLSI Antibiotics only
#result=pd.merge(strain_mic, clsi_bp, on='Antibiotics',how='inner',  indicator=True)[['Strain name','Antibiotics', 'MIC', 'o_mic', 's_mic', 'r_mic','_merge']]
try:
    result=pd.merge(strain_mic, clsi_bp, on='Antibiotics',how='inner')[['Strain name','Antibiotics', 'MIC', 'o_mic', 's_mic', 'r_mic','i_mic']]
except KeyError:
    print('Waring: Error in input Values')


# In[221]:


#compare MIC values and assign Susceptible and Resistant to Strain
#try:
result[['CLSI_profile']] = result[['o_mic','s_mic','r_mic','i_mic']].apply(sus_res_int,axis = 1)
#except ValueError:
#    print('Waring: Error in input MIC value')


# In[222]:


#result


# In[223]:


#result[['Strain name', 'Antibiotics', 'MIC','s_mic','r_mic','CLSI_profile']].to_csv(output_profile,sep=',', index=False, encoding='utf-8-sig')


# In[224]:


#create a pivot table for ASIST
table=result[['Strain name', 'Antibiotics','CLSI_profile']].drop_duplicates()
result_table=pd.pivot_table(table, values ='CLSI_profile', index =['Strain name'],columns =['Antibiotics'], aggfunc = lambda x: ' '.join(x))


# In[225]:


#result_table


# In[226]:


#result_table.to_csv(output_table,na_rep='NA')


# In[227]:


# reorder the Antibiotics for ASIST
clsi_ab=['Amikacin','Tobramycin','Gentamycin','Imipenem','Meropenem','Doripenem','Ciprofloxacin','Levofloxacin',
         'Piperacillin/tazobactam','Ticarcillin/clavulanicacid','Cefotaxime','Ceftriaxone','Ceftazidime','Cefepime',
         'Trimethoprim/sulfamethoxazole','Ampicillin/sulbactam','Colistin','Polymyxinb','Tetracycline','Doxicycline ',
         'Minocycline']
result_selected=result_table.filter(clsi_ab)


# In[228]:


#print(result_selected.shape, result_table.shape)


# In[229]:


result_selected.insert(0,'Resistance_phenotype','')


# In[230]:


#rename headers
result_selected=result_selected.rename(columns = {'Ticarcillin/clavulanicacid':'Ticarcillin/clavulanic acid','Piperacillin/tazobactam':'Piperacillin/ tazobactam','Trimethoprim/sulfamethoxazole': 'Trimethoprim/ sulfamethoxazole','Ampicillin/sulbactam':'Ampicillin/ sulbactam', 'Polymyxinb': 'Polymyxin B'} )


# In[231]:


#result_selected


# In[232]:


result_selected.to_csv(output_table,na_rep='NA')

