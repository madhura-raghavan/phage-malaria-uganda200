"""MADHURA RAGHAVAN 
   DERISI LAB, UCSF 
   2020 DEC"""  


### This is the first part of the pipeline described in Fig6a 
### Pipeline to identify the list of significantly enriched motifs in the seroreactive peptides over background 


# Instructions 

## Prior to running this pipeline, one has to run the SliMFinder pipeline to identify the list of motifs shared by at least two peptides in a set of seroreactive peptides. In this case, this was run with the set of 9927 seroreactive peptides identified in this study  

## python SLiMSuite-master/tools/slimfinder.py seqin=HITS_Round2_Pfonly_250kfil_techclean_3zscorefil_5patientsfil_allhits_sequences_9227.fa teiresias=T efilter=F blastf=F masking=F ftmask=F imask=F compmask=F metmask=F slimlen=7 absmin=2 absminamb=2  slimchance=F maxwild=1 maxseq=10000 walltime=2400 minocc=0.0002 ambocc=0.0002 wildvar=False equiv=conservative_subs_equivalency.txt resdir=. 

## conservative_subs_equivalency.txt has the following [AG,TS,DE,NQ,RHK,LVI,YFW]


## A motif file is output from SlimFinder. Each line in the SlimFinder output file contains the motif, number of occurences, number of peptides contaning the motif, list of positions where it occurs. This is the input 'motif file' in the pipline below. 

## The set of seroreactive peptides input in the SliMsearch pipeline is the 'peptide file' (HITS_Round2_Pfonly_250kfil_techclean_3zscorefil_5patientsfil_allhits_sequences_9227.fa)

## wildcard yes/no tells the script whethere you want to look for motifs with or without wildcards 

## minmotiflength=<6/7/8/9> tells what the length of motif you are looking for 

## identity=<identity5aa/allconssub> tells whether you are looking for motifs with at least 5 identical residue positions or conservative substitutions allowed at all positions. If identitical residues are required, at least 5 have to be identical and rest can be wildcard or sub depending on motif length ## 1 wildcard allowed for a 6/7-mer and upto 2 wildcards for 8/9-mer. If conservative subs are allowed at all positions - wildcard wise 1 wildcard allowed for 6 and 7aa, upto 2 wildcards for 8 and 9aa motifs 

## command to run for this study - python Fig6a_Pipeline_Part1_enriched_motifs_over_background.py <motif file> HITS_Round2_Pfonly_250kfil_techclean_3zscorefil_5patientsfil_allhits_sequences_9227.fa wildcard=<> minmotiflength=<> identity=<>


import re
import sys 
from collections import defaultdict 
import os 

def help():
    print ("Usage is \n")
    print ("python Fig6a_Pipeline_Part1_enriched_motifs_over_background.py <motif file> <peptide file> wildcard=<yes/no> minmotiflength=<6/7/8/9> identity=<identity5aa/allconssub>\n Without spaces between argument and parameters \n")
 

### import input file (slim finder hits in any sample) and the corresponding peptide file\\

try:
    motif_file = sys.argv[1]
    peptide_file = sys.argv[2]
    wildcard = sys.argv[3]
    minmotiflen = str(sys.argv[4].split("=")[1])
    identity = sys.argv[5].split("=")[1]
    
except:
    help() 


folder = os.mkdir(motif_file.split("/")[-1].split(".out")[0] + "_minmotiflen_" + str(minmotiflen)+"_"+identity+'_'+'wildcard'+wildcard.split("=")[1] + '/')

    
    
filename = motif_file.split("/")[-1].split(".out")[0] + "_minmotiflen_" + str(minmotiflen)+"_"+identity+'_'+'wildcard'+wildcard.split("=")[1] +'/Motif'

### get the character set to look for in the motif to filter them 
if wildcard == 'wildcard=yes':
    chars = set('^${')
    filename = filename+'_withWILDCARD' 
    

elif wildcard == 'wildcard=no':
    chars = set('^${.')
    filename = filename+'_NOWILDCARD' 
    
else:
    raise ValueError('Use wildcard=yes or wildcard=no without spaces for the wildcard argument')
  
### creating a logfile 

log = open(filename+'.log','w')

#### 1 - Reading through the motif file to gather peptide clusters 

dict_clusters = defaultdict(list)


fh = open(motif_file,'r')

for i in range(0,5):
    fh.readline()  ### Read through first five lines because that is the header seq 


for line in fh:
    
    item = line.strip("\n").split("\t")
    motif = item[2].split(" ")[0]

    pos = item[2].split(" ")[1:]
    cluster = set(pos[0::2])
    dict_clusters[motif] = list(cluster)
    
fh.close()

log.write("The total number of motifs is " + str(len(dict_clusters.keys())) + '\n')

### 2 - collect motifs which satisfy the length and identity criteria 
### Also expand the motifs into the patterns they recognize in the dataset so you can be more specific about residues accepted in the wildcard position and get rid of redundant motifs reepresenting the same pattern 


### collect the sequence content of the list of peptides analyzed (hit peptides or random peptides) 

f = open(peptide_file,'r')
total_seq = ''
for line in f:
    total_seq = total_seq + next(f)
f.close()

        
### collect motifs which have atleast 5 identical items in the motif


dict_minlength_cluster = {}

i = 0

for motif in dict_clusters:
    
    if i% 100 == 0:
        print (i) 
    i = i+1
    
    
    if any((c in chars) for c in motif): ### no flexible wildcards or position fixed motifs and if condition chosen, no other wildcards too 
        continue
    elif (('EEEEE' not in motif) & ('NNNNN' not in motif) & ('KKKKK' not in motif)):    ### consider only motifs that do not have homopolymeric stretches in the motif. These stretches were common in P. falciparum motifs 
        
#calculating the minumum length of the motif possible
        motif_2 = re.sub(r'\[.[^]]*?\]',"X",motif) ### assign the flexible positions (conservative subs) as a single letter X 

#        min_length_alphabets = sum(c.isalpha() for c in motif_2) ### get the number of non-wildcard positions in the motif 
        min_length_alphabets = len(motif_2)
        wildcard = motif.count('.') ## count number of wildcards in the motif 
        subs = motif_2.count('X')
   
        if identity == 'identity5aa': ### if identitical residues are required, at least 5 have to be identical and rest can be wildcard or sub depending on motif length ## 1 wildcard allowed for a 6/7-mer and upto wildcards for 8/9-mer
 
            if (minmotiflen == str(6)):

                if (((min_length_alphabets == 6) & (wildcard == 0) & (subs <= 1) ) | 
                ((min_length_alphabets == 6) & (wildcard <= 1) & (subs == 0))):   

                    pattern = sorted(set(re.findall(motif,total_seq)))

                    try:
                        motif_expanded = '|'.join(pattern) 
                    except:
                        motif_expanded = motif
                    dict_minlength_cluster[motif_expanded] = dict_clusters[motif]
                    log.write( motif_expanded + '\n')
                    
            elif (minmotiflen == str(7)):

                if (((min_length_alphabets == 7) & (wildcard <= 1) & (subs <= 1) )):   

                    pattern = sorted(set(re.findall(motif,total_seq)))

                    try:
                        motif_expanded = '|'.join(pattern) 
                    except:
                        motif_expanded = motif
                    dict_minlength_cluster[motif_expanded] = dict_clusters[motif]
                    log.write( motif_expanded + '\n')

            elif (minmotiflen == str(8)):

                if (

                ((min_length_alphabets == 8) & (wildcard <= 2) & (subs <= 1))|
                ((min_length_alphabets == 8) & (wildcard <= 1) & (subs <= 2))
               ):   

                    pattern = sorted(set(re.findall(motif,total_seq)))

                    try:
                        motif_expanded = '|'.join(pattern) 
                    except:
                        motif_expanded = motif
                    dict_minlength_cluster[motif_expanded] = dict_clusters[motif]
                    log.write( motif_expanded + '\n')

            elif (minmotiflen == str(9)):

                if (

                ((min_length_alphabets == 9) & (wildcard <= 2) & (subs <= 2))
               ):   

                    pattern = sorted(set(re.findall(motif,total_seq)))

                    try:
                        motif_expanded = '|'.join(pattern) 
                    except:
                        motif_expanded = motif
                    dict_minlength_cluster[motif_expanded] = dict_clusters[motif]
                    log.write( motif_expanded + '\n')
            
        elif identity == 'allconssub': ### conservative subs are allowed at all positions - wildcard wise 1 wildcard allowed for 6 and 7aa, upto 2 wildcards for 8 and 9aa motifs
 
            if (minmotiflen == str(6)):

                if (((min_length_alphabets == 6) & (wildcard <= 1))):   

                    pattern = sorted(set(re.findall(motif,total_seq)))

                    try:
                        motif_expanded = '|'.join(pattern) 
                    except:
                        motif_expanded = motif
                    dict_minlength_cluster[motif_expanded] = dict_clusters[motif]
                    log.write( motif_expanded + '\n')
                    
            elif (minmotiflen == str(7)):

                if (((min_length_alphabets == 7) & (wildcard <= 1) )):   

                    pattern = sorted(set(re.findall(motif,total_seq)))

                    try:
                        motif_expanded = '|'.join(pattern) 
                    except:
                        motif_expanded = motif
                    dict_minlength_cluster[motif_expanded] = dict_clusters[motif]
                    log.write( motif_expanded + '\n')

            elif (minmotiflen == str(8)):

                if (

                ((min_length_alphabets == 8) & (wildcard <= 2))):   

                    pattern = sorted(set(re.findall(motif,total_seq)))

                    try:
                        motif_expanded = '|'.join(pattern) 
                    except:
                        motif_expanded = motif
                    dict_minlength_cluster[motif_expanded] = dict_clusters[motif]
                    log.write( motif_expanded + '\n')

            elif (minmotiflen == str(9)):

                if (

                ((min_length_alphabets == 9) & (wildcard <= 2)
               )):   

                    pattern = sorted(set(re.findall(motif,total_seq)))

                    try:
                        motif_expanded = '|'.join(pattern) 
                    except:
                        motif_expanded = motif
                    dict_minlength_cluster[motif_expanded] = dict_clusters[motif]
                    log.write( motif_expanded + '\n')
            


log.write("The total number of motifs of min length" + str(minmotiflen) + " is " + str(len(dict_minlength_cluster.keys())) + '\n') 
### make the reverse dictionary for the motifs meeting mimimum length which has the motifs present as values for each peptide as key

#write the dict_min_length_clusterr into a file 

fh899 = open(filename+'_minlength_motif_clusters.txt','w')

fh899.write('motif pattern' + '\t' + 'hit_peptides_with_motif' + '\n')

for motif in dict_minlength_cluster.keys():
    fh899.write(motif + '\t' + str(dict_minlength_cluster[motif])+ '\n')

fh899.close() 

dict_peptide_clusters = defaultdict(dict) 

for motif,cluster in dict_minlength_cluster.items():
    for peptide in cluster:
        dict_peptide_clusters[peptide][motif]=dict_minlength_cluster[motif]
        

#### 3 - calculate the background frequency of the identified motifs 

###  Background p-value calculation 

## Randomly sample 9927 peptides from the background. Bootstrap this 1000 times. 
## Each sampling, loop through the set of shared motifs and collect the number of occurences of the motif in the randomly sampled set using grep 
## Given this background distribution, calculate a p-value for the number of times the motif is shared in the real dataset 
### calculate the background freq of the motif and p-value for the freq observed in the hit peptides 
import pandas as pd 
import random 

      
### Get the falciparome peptidee fasta file (If starting from the Whole_library_peptide_to_gene_mapping.csv file, change the code below to get peptide and sequence from that file into a df) 

fh2 = open('FinalPlasmodiumFalciparumPhagePeptides.v2.fa','r')

falc_seq = {}
for line in fh2: 
    name = line.strip("\n")[1:] ## don't include the '>' symbol at the beginning
    name = name.replace(",", "_") ## some peptidees have a , in the file, but _ in the dataset
    sequence = next(fh2)
    falc_seq[name] = sequence

fh2.close()

df_peptide_fasta = pd.DataFrame.from_dict(falc_seq, orient='index')
df_peptide_fasta.columns = ['sequence']


df_peptide_fasta_PfOnly = df_peptide_fasta.drop(list(df_peptide_fasta.filter(regex = r'(virus|toxin|sapiens)',axis=0).index))

df_peptide_fasta_PfOnly.reset_index(inplace=True)



### Collect all the pattern matches for each motif in the whole background 
dict_bckd = {}

df_bckd_freq = pd.DataFrame()
i = 0


for motif in dict_minlength_cluster.keys():
    print (i)
    i = i +1
    match = list(df_peptide_fasta_PfOnly[df_peptide_fasta_PfOnly['sequence'].str.contains(motif)].index)
    dict_bckd[motif] = match 


df_dict_bckd = pd.DataFrame.from_dict(dict_bckd, orient='index')

df_dict_bckd.to_csv(filename + 'df_dict_bckd_temp.csv')


for i in range(1,1001):
    print ('co' + str(i) )
    rand = random.sample(range(1, 230556), 9927)
    
    count = df_dict_bckd.isin(rand).sum(axis=1)
    df_bckd_freq[i]=count 


  ### Calculate p-value

### Mean random bckd freq 
    
mean_bckd = df_bckd_freq.mean(axis=1)
## The mean and variance of the sampels are more or less similar. So, Poisson distribution will be used to calculate p-value

### Use Poisson distribution to calculate p-value of occurence of motif in the dataset 
from scipy.stats import poisson
dict_largest_cluster_motif_freq_hits = {}
for motif in dict_minlength_cluster.keys():
    dict_largest_cluster_motif_freq_hits[motif] = len(dict_minlength_cluster[motif])
 
for motif in list(df_bckd_freq.index):
    df_bckd_freq.loc[motif, 'p-value'] = poisson.pmf( dict_largest_cluster_motif_freq_hits[motif], mu = mean_bckd.loc[motif])

df_bckd_freq.to_csv(filename + 'df_bckd_freq_temp.csv')

### Multiple hypothesis correction with an FDR of 0.1% (0.001) 

import statsmodels.stats.multitest as sm
pvalue_adj = sm.multipletests(df_bckd_freq['p-value'], alpha=0.001, method='fdr_bh', is_sorted=False)
df_bckd_freq['pvalue_adj'] = pvalue_adj[1]
df_bckd_freq['significant'] = pvalue_adj[0]


df_bckd_freq.to_csv(filename + 'df_bckd_freq_temp.csv')

#### Get the significant motifs after MH correction and also get motifs where the number of occurences in HIT set is higher than mean (since significance could be called for those that fall on the left tail of the background distribution too 


significant = set(df_bckd_freq[(df_bckd_freq['significant']== True)].index)
significant_largest = set([x for x in significant if dict_largest_cluster_motif_freq_hits[x] > df_bckd_freq.loc[x,df_bckd_freq.columns[0:1000]].mean()])

log.write("The total number of MH corrected significant min length motifs is " + str(len(significant_largest)) + '\n') 

 ## Save the significant motifs into a file along with what hiit peptides have them as well as the background set

fh30 = open(filename+'_MHcorrected_significant_motifs_FDR_E_3.txt','w')
fh30.write('motif' + '\t' + 'hit_peptides_with_motif' + '\t' + 'allPfOnlypeptides_with_motif' + '\t' + 'pvalue' + '\t' + 'pvalue_adj'+ '\n')
for motif in significant_largest:
    fh30.write(motif + '\t' + str(dict_minlength_cluster[motif]) + '\t' + str(list(df_dict_bckd.loc[motif]))   +'\t' + str(df_bckd_freq.loc[motif,'p-value'])  + '\t' + str(df_bckd_freq.loc[motif,'pvalue_adj']) + '\n')
fh30.close()


### 4 - THIS IS AN OPTIONAL STEP. This didn't change the result for this study. 
## collect the most shared significantly eenriched motifs in each peptide (>= mean number of shared peptides for the motifs in each peptide)


def Average(lst): 
    return sum(lst) / len(lst) 

dict_largest_cluster_motif = defaultdict(list)

for i in range(0,9927):  ### Going through peptide numbers and clusters each one belongs to 
    
    d = dict_peptide_clusters[str(i)]
    print (i)
    #### collect all the motifs with the largest clusters for each peptide 
             
    try:
        sharing = set([len(v) for k,v in d.items() if k in significant_largest]) ## the number of sharing peptides for the different motifs that are significantly enriched over background 
        mean = Average(sharing)
        terms = [key for key in d if ((len(d[key]) >= mean) & (key in significant_largest))] ### collect all motifs with >= mean number of sharing
        
        dict_largest_cluster_motif[i] = terms
  
    except:  ### if the peptide doesn't have any minlength motifs identified, this process is not done  
        continue

set_largest_cluster_motifs = set()

for key in dict_largest_cluster_motif:
    for motif in dict_largest_cluster_motif[key]:
        set_largest_cluster_motifs.add(motif)

log.write("The total number of largest shared min length motifs is " + str(len(set_largest_cluster_motifs)) + '\n') 


 ## Save the significantly largest motifs into a file along with what hiit peptides have them as well as the background set

fh30 = open(filename+'_MHcorrected_significant_motifs_largest50percentshared_FDR_E_3.txt','w')
fh30.write('motif' + '\t' + 'hit_peptides_with_motif' + '\t' + 'allPfOnlypeptides_with_motif' + '\t' + 'pvalue' + '\t' + 'pvalue_adj'+ '\n')
for motif in set_largest_cluster_motifs:
    fh30.write(motif + '\t' + str(dict_minlength_cluster[motif]) + '\t' + str(list(df_dict_bckd.loc[motif]))   +'\t' + str(df_bckd_freq.loc[motif,'p-value'])  + '\t' + str(df_bckd_freq.loc[motif,'pvalue_adj']) + '\n')
fh30.close()


