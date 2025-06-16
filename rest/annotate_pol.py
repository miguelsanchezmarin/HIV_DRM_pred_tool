import pandas as pd
from Bio.Seq import Seq

###Functions required to handle the subtype information for annotation
def get_subtype_from_fasta(fasta_file):
    """
    Get the subtype from a FASTA file. Just gets the subtype from a reference FASTA where it is after >.
    
    Args:
        fasta_file (str): Path to the FASTA file.
        
    Returns:
        str: Subtype extracted from the FASTA header.
    """
    with open(fasta_file, 'r') as f:
        header = f.readline().strip()
        # Assuming the subtype is the first part of the header after '>'
        subtype = header.split('>')[1].split()[0] if '>' in header else header.split()[0]
    
    subtype = subtype.replace('CONSENSUS', '').strip().upper()

    if subtype.find('01_AE'):
        subtype = '01_AE'
    elif subtype.find('02_AG'):
        subtype = '02_AG'
    elif subtype.find('03_AB'):
        subtype = '03_AB'
    elif subtype.find('04_CPX'):
        subtype = '04_CPX'
    elif subtype.find('06_CPX'):
        subtype = '06_CPX'
    elif subtype.find('08_BC'):
        subtype = '08_BC'
    elif subtype.find('09_CPX'):
        subtype = '09_CPX'
    elif subtype.find('10_CD'):
        subtype = '10_CD'
    elif subtype.find('11_CPX'):
        subtype = '11_CPX'
    elif subtype.find('12_BF'):
        subtype = '12_BF'
    elif subtype.find('14_BG'):
        subtype = '14_BG'
    elif subtype.find('A1'):
        subtype = 'A1'
    elif subtype.find('A2'):
        subtype = 'A2'
    elif subtype.find('F1'):
        subtype = 'F1'
    elif subtype.find('F2'):
        subtype = 'F2'
    elif subtype.find('G'):
        subtype = 'G'
    elif subtype.find('H'):
        subtype = 'H'
    elif subtype.find('D'):
        subtype = 'D'
    elif subtype.find('B'):
        subtype = 'B'
    elif subtype.find('C'):
        subtype = 'C'
    
    return subtype

def load_vcf_dict(vcf_file):
    """
    Load a VCF file into a dictionary.
    
    INPUT:
        vcf_file (str): Path to the VCF file.
        
    OUTPUT:
        dict: Dictionary with positions as keys and variant information as values.
    """
    vcf_dict = {}
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            pos = int(parts[1])  # Position in the VCF
            ref = parts[3]       # Reference allele
            alt = parts[4]       # Alternative allele
            info = parts[7]     # INFO field
            af = float(info.split('AF=')[1].split(';')[0]) if 'AF=' in info else None # Allele frequency
            vcf_dict[pos] = {'ref': ref, 'alt': alt, 'freq': af}
    return vcf_dict

def load_coverage_file(coverage_file):
    """
    Load a coverage file into a dictionnary.
    
    INPUT:
        coverage_file (str): Path to the coverage file.
        
    OUTPUT:
        coverage_dict: dictionary with positions as keys and coverage information as values.
    """
    cov_dict = {}
    with open(coverage_file, 'r') as f:
        for i, line in enumerate(f):  
            if i == 0: #we skip the first line as it is a header
                continue
            parts = line.strip().split('\t')
            pos = int(parts[1])
            coverage = float(parts[2])  # Assuming coverage is in the third column
            cov_dict[pos] = coverage
    return cov_dict

def get_ref_codons(subtype):
    """
    Get the reference codons for a given subtype from the reference FASTAs.
    
    INPUT:
        subtype (str): The subtype to get codons for.
        
    OUTPUT:
        list: a list with a reference codon for each position in the subtype.
    """
    fasta_ref = f"ref/pol_refs/pol_ref_CONSENSUS_{subtype}.fasta"
    with open(fasta_ref, 'r') as f:
        lines = f.readlines()
        ref_codons = []
        for line in lines:
            if line.startswith('>'):
                continue
            seq = line.strip()
            # Convert the sequence to codons
            for i in range(0, len(seq), 3):
                codon = seq[i:i+3]
                if len(codon) == 3:
                    ref_codons.append(codon)
    return ref_codons

def annotate_vcf(vcf_dict, ref_codons, subtype):
    """
    Annotates from the nn vcf the aminoacid changes and their associated allele frequency.
    INPUT:
        vcf_dict (dict): Dictionary with positions as keys and variant information as values.
        ref_codons (list): List of reference codons for the subtype.
        subtype (str): Subtype of the virus.
    OUTPUT:
        annotated_vcf_dict: Dictionary with positions as keys and ref aa, alt aa and allele frequency as values."""

    #first we get the start and end aa positions for the subtype for each protein
    indexes = pd.read_csv("ref/pol_refs/pol_indexes.tsv", sep='\t')
    subtype_indexes = indexes[indexes['Subtype'] == "CONSENSUS_" + subtype]
    PR_start, PR_end = subtype_indexes['PR_start'].values[0], subtype_indexes['PR_end'].values[0]
    RT_start, RT_end = subtype_indexes['RT_start'].values[0], subtype_indexes['RT_end'].values[0]
    IN_start, IN_end = subtype_indexes['IN_start'].values[0], subtype_indexes['IN_end'].values[0]
    
    annotated_vcf_dict = {}
    for pos in vcf_dict.keys():
        codon_index = (pos-1) // 3
        pos_in_codon = (pos-1) % 3

        if codon_index < len(ref_codons):
            ref_codon = ref_codons[codon_index]
            ref_aa = Seq(ref_codon).translate()
            alt_codon = ref_codon[:pos_in_codon] + vcf_dict[pos]['alt'] + ref_codon[pos_in_codon+1:]
            alt_aa = Seq(alt_codon).translate() if len(alt_codon) == 3 else None

            #we annotate the PR, RT and IN proteins
            if codon_index >= PR_start and codon_index <= PR_end:
                protein = 'PR'
                prot_pos = codon_index + 1 - PR_start
            elif codon_index >= RT_start and codon_index <= RT_end:
                protein = 'RT'
                prot_pos = codon_index + 1 - RT_start
            elif codon_index >= IN_start and codon_index <= IN_end:
                protein = 'IN'
                prot_pos = codon_index + 1 - IN_start
            
            print(f"Position: {prot_pos}, Ref Codon: {ref_codon}, Alt Codon: {alt_codon}, Ref AA: {ref_aa}, Alt AA: {alt_aa}, Frequency: {vcf_dict[pos]['freq']}, Protein: {protein}")

            annotated_vcf_dict[prot_pos] = {
                'ref_aa': str(ref_aa),
                'alt_aa': str(alt_aa) if alt_aa else None,
                'freq': vcf_dict[pos]['freq'],
                'protein': protein
            }

    return annotated_vcf_dict

def annotate_coverage(coverage_dict, subtype):
    """
    Annotates coverage information for a given subtype.
    
    INPUT:
        coverage_dict (str): dictionary with positions as keys and coverage information as values.
        subtype (str): Subtype of the virus.
        
    OUTPUT:
        coverage_df: DataFrame with annotated coverage information.
    """
    #first we get the start and end aa positions for the subtype for each protein
    indexes = pd.read_csv("ref/pol_refs/pol_indexes.tsv", sep='\t')
    subtype_indexes = indexes[indexes['Subtype'] == "CONSENSUS_" + subtype]
    PR_start, PR_end = subtype_indexes['PR_start'].values[0], subtype_indexes['PR_end'].values[0]
    RT_start, RT_end = subtype_indexes['RT_start'].values[0], subtype_indexes['RT_end'].values[0]
    IN_start, IN_end = subtype_indexes['IN_start'].values[0], subtype_indexes['IN_end'].values[0] 
    
    annotated_coverage = []
    # Annotate coverage positions
    covered_positions = sorted(coverage_dict.keys()) 
    first_position = covered_positions[0]
    if first_position % 3 == 0: #then the first covered nucleotide is the third of a codon
        coverage_codon_positions = [list(covered_positions[1:][i:i + 3]) for i in range(0, len(covered_positions[1:]), 3)]  # Group positions into codons
        coverage_codon_positions.insert(0, [first_position]) #add the first remaining codon
    elif first_position % 3 == 1: #then the first covered nucleotide is the first of a codon
        coverage_codon_positions = [list(covered_positions[i:i + 3]) for i in range(0, len(covered_positions), 3)]
    elif first_position % 3 == 2: #then the first covered nucleotide is the second of a codon
        coverage_codon_positions = [list(covered_positions[2:][i:i + 3]) for i in range(0, len(covered_positions[2:]), 3)] 
        coverage_codon_positions.insert(0, [first_position, covered_positions[1]]) #add the first two remaining codons
    
        
    for codon in coverage_codon_positions:
        if len(codon) == 3:
            if coverage_dict[codon[0]]==0 or coverage_dict[codon[1]]==0 or coverage_dict[codon[2]]==0: ##if not all positions are covered we considered the codon uncovered.
                pos_coverage = 0
            else:
                pos_coverage = (coverage_dict[codon[0]] + coverage_dict[codon[1]] + coverage_dict[codon[2]]) / 3 ##if all positions are covered we take the mean
            
            codon_index = (codon[0]-1) // 3
            if codon_index >= PR_start and codon_index <= PR_end:
                protein = 'PR'
                prot_pos = codon_index + 1 - PR_start
            elif codon_index >= RT_start and codon_index <= RT_end:
                protein = 'RT'
                prot_pos = codon_index + 1 - RT_start
            elif codon_index >= IN_start and codon_index <= IN_end:
                protein = 'IN'
                prot_pos = codon_index + 1 - IN_start
            
            annotated_coverage.append([subtype, prot_pos, pos_coverage, protein])
    coverage_df = pd.DataFrame(annotated_coverage, columns=['Subtype', 'CodonPosition', 'Coverage', 'Protein'])
            
    return coverage_df


def vcf_to_df(vcf_dict):
    """
    Converts annotated VCF dictionary to a DataFrame, adding PR, RT and IN protein names.
    
    INPUT:
        vcf_dict (dict): Dictionary with positions as keys and variant information as values.
    OUTPUT:
        annotated_vcf_df: DataFrame with annotated VCF information."""
    
    annotated_vcf_df = pd.DataFrame.from_dict(vcf_dict, orient='index')
    annotated_vcf_df.reset_index(inplace=True)
    annotated_vcf_df.rename(columns={'index': 'Position', 'ref_aa': 'Ref', 'alt_aa': 'Mut', 'freq': 'Freq', 'protein':'Prot'}, inplace=True)
    
    return annotated_vcf_df


def main(snvs_file, coverage_file, fasta_ref, snvs_output, coverage_output):
    """
    Main function to extract subtype from a FASTA file and print it.
    Args:
        snvs_file (str): Path to the SNVs file.
        coverage_file (str): Path to the coverage file.
        fasta_ref (str): Path to the reference FASTA file used for variant calling.
    """

    subtype = get_subtype_from_fasta(fasta_ref)
    print(f"Extracted subtype: {subtype.replace('CONSENSUS_', '')}")

    snvs_dic = load_vcf_dict(snvs_file)

    coverage_dic = load_coverage_file(coverage_file)

    ref_codons = get_ref_codons(subtype)

    annotated_vcf = annotate_vcf(snvs_dic, ref_codons, subtype)
    annotated_vcf_df = vcf_to_df(annotated_vcf)

    annotated_coverage = annotate_coverage(coverage_dic, subtype)

    annotated_vcf_df.to_csv(snvs_output, sep='\t', index=False)
    annotated_coverage.to_csv(coverage_output, sep='\t', index=False)
    

    

if __name__ == "__main__":
    main(
         snakemake.input.fname_vcf,
         snakemake.input.fname_cov,
         snakemake.input.ref,
         snakemake.output.fname_vcf,
         snakemake.output.fname_cov
         )
