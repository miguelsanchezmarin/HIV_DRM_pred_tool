import pandas as pd
import re

def HIVDB_singlemut_annot(mut_list, dataset: str):
    '''Retrieves the HIVDB annotation for each mutation in a list of mutations.
    INPUT:
    mut_list: list of mutations i.e [M41L,L90M].
    dataset: drug class to be used. 'INI', 'NNRTI', 'NRTI', 'PI'.

    OUTPUT:
    mut_annot_df: dataframe with the mutations and their HIVDB annotations.
    '''
    #we read the comments files for the dataset
    comments_file = pd.read_csv(f'HIVDB_rules/{dataset}_comments_Stanford_HIVDB', sep=',')
    comments_file.columns = comments_file.columns.str.replace('\n','')

    mut_list = [mut[1:] for mut in mut_list] #we remove the first letter of the mutations
    mut_annot = []
    for mut in mut_list:
        position = re.sub(r'\D+', '', mut)
        mut_aa = re.sub(r'\d+', '', mut)
        comment_in = False
        for i, comm in enumerate(comments_file['Condition'].values):
            if re.match(r'^\d+[A-Z]+$', comm) and position==re.sub(r'\D+', '', comm): #we check if the mutation is in the comment
                comm_ambig = re.sub(r'\d+', '', comm)
                if mut_aa in comm_ambig:
                    mut_annot.append([mut, 
                                     comments_file['Comment/Mutation Type'].values[i], 
                                     comments_file['Comment'].values[i]])
                    comment_in = True
                    break
        
        if not comment_in:
            mut_annot.append([mut, "Unknown", "Unknown"])
                       
    mut_annot_df = pd.DataFrame(mut_annot, columns=['Mutation', 'Annotation', 'Comment'])
    return mut_annot_df

def HIVDB_table(mut_tsv, lower_cutoff: float = 0.015):
    '''Generates a table with the observed mutations and their HIVDB classifications for all drugs and drug class.
    INPUT:
    mut_tsv: tsv file output from annotate_vcf.py with the observed mutations. Columns are Position, Ref, Mut, Freq and Prot. 
    lower_cutoff: minimum frequency to be considered as a mutation.

    OUTPUT:
    mutations_df: dataframe with the mutations for all drugs and drug classes and their HIVDB annotations.
    '''
    mutations_df = pd.read_csv(mut_tsv, sep='\t')#we read the tsv file
    mutations_filter_freq = mutations_df[mutations_df['Freq'] > lower_cutoff]#we apply the frequency cutoff
    mutations_filter_freq = mutations_filter_freq[mutations_filter_freq['Ref'] != mutations_filter_freq['Mut']]##We filter out the rows where Ref == Mut
    PR_mut_freq = extract_prot_mut_freq(mutations_filter_freq, 'PR')
    RT_mut_freq = extract_prot_mut_freq(mutations_filter_freq, 'RT')
    IN_mut_freq = extract_prot_mut_freq(mutations_filter_freq, 'IN')
    INI_annot_df = HIVDB_singlemut_annot(list(IN_mut_freq.keys()), 'INI')
    NNRTI_annot_df = HIVDB_singlemut_annot(list(RT_mut_freq.keys()), 'NNRTI')
    NRTI_annot_df = HIVDB_singlemut_annot(list(RT_mut_freq.keys()), 'NRTI')
    PI_annot_df = HIVDB_singlemut_annot(list(PR_mut_freq.keys()), 'PI')
    INI_annot_df['Freq'] = IN_mut_freq.values()
    INI_annot_df['Prot'] = 'IN'
    NNRTI_annot_df['Freq'] = RT_mut_freq.values()
    NNRTI_annot_df['Prot'] = 'RT'
    NRTI_annot_df['Freq'] = RT_mut_freq.values()
    NRTI_annot_df['Prot'] = 'RT'
    PI_annot_df['Freq'] = PR_mut_freq.values()
    PI_annot_df['Prot'] = 'PR'
    
    merged_NNRTI_NRTI = pd.merge(NNRTI_annot_df, NRTI_annot_df, on=['Mutation', 'Freq', 'Prot'], how='outer', suffixes=('_NNRTI', '_NRTI'))##We unify NNRTI and NRTI dataframes
    merged_NNRTI_NRTI['Annotation'] = merged_NNRTI_NRTI.apply(lambda row: row['Annotation_NNRTI'] if row['Annotation_NNRTI'] != 'Unknown' else row['Annotation_NRTI'], axis=1) ##The NRTI and NNRTI commented positions do not overlap
    merged_NNRTI_NRTI['Comment'] = merged_NNRTI_NRTI.apply(lambda row: row['Comment_NNRTI'] if row['Comment_NNRTI'] != 'Unknown' else row['Comment_NRTI'], axis=1)
    merged_NNRTI_NRTI = merged_NNRTI_NRTI[['Mutation', 'Freq', 'Prot', 'Annotation', 'Comment']]
    merged_NNRTI_NRTI["Pos"] = merged_NNRTI_NRTI['Mutation'].apply(lambda x: int(re.sub(r'\D+', '', x)))
    merged_NNRTI_NRTI = merged_NNRTI_NRTI.sort_values(by=['Pos'], ascending = True).drop(columns=['Pos'], axis=1)


    mutations_df = pd.concat([INI_annot_df, merged_NNRTI_NRTI, PI_annot_df], axis=0, ignore_index=True)
    mutations_df = mutations_df[['Mutation', 'Freq', 'Prot', 'Annotation', 'Comment']]

    return mutations_df

def extract_prot_mut_freq(mutations_df, prot = ''):
    '''Outputs a list of dictionaries, one for each alternative list of mutations.
       Each dictionary presents the mutations for a given protein and their observed frequency, coming from a mutation_freq.tsv file.
    INPUT:
    mutations_df: tsv file output from annotate_vcf.py with the observed mutations. Columns are Position, Ref, Mut, Freq and Prot.
    prot: protein to filter the mutations. PR, RT or IN. If not one of those, it will return all mutations.
    
    OUTPUT:
    prot_mut_freq: dictionary with the mutations for the given proteinas keys and the frequencies as values.
    '''
    if prot in ['PR', 'RT', 'IN']:
        mutations_df = mutations_df[mutations_df['Prot'] == prot]
    
    prot_mut_freq = {} #we get a dictionary with the mutations as keys and the frequencies as values
    for i, row in mutations_df.iterrows():
        
        full_mut = mutations_df.loc[i, 'Ref'] + mutations_df.loc[i, 'Position'].astype(str) + mutations_df.loc[i, 'Mut']

        if full_mut not in prot_mut_freq:
            prot_mut_freq[full_mut] = row['Freq']
        else:
            prot_mut_freq[full_mut] += row['Freq']
    
    return prot_mut_freq
