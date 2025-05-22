#!/usr/bin/env python3
"""
CORRECTED:
Basically copied this file:
https://github.com/rpetit3/vcf-annotator/blob/master/vcf-annotator.py
We copied it from:
https://github.com/LaraFuhrmann/Scan-for-mutations-of-interest-NGS-samples/blob/main/workflow/scripts/annotate_vcf.py 

Annotate mutation calls using https://github.com/rpetit3/vcf-annotator.
Position of mutation calls are with respect to the modified parental_stock_consensus
and therefore must be mapped back to EB reference space.

"""
from fuc import pyvcf
import subprocess
import collections
import pandas as pd
import os

# ---------------
# from here it is copied

from Bio import SeqIO
from Bio.Seq import Seq
import vcf


class Annotator(object):
    """Annotate a given VCF file and coverage.tsv file according to the reference GenBank.""" #modified by MSM

    def __init__(self, gb_file=False, vcf_file=False, coverage_file=False):
        """Initialize variables."""
        self.__annotated_features = ["CDS", "tRNA", "rRNA", "ncRNA",
                                     "misc_feature"]
        self.__gb = GenBank(gb_file)
        self.__vcf = VCFTools(vcf_file)
        self.__coverage_file = coverage_file
        self.add_annotation_info()

    def add_annotation_info(self):
        """Add custom VCF info fields."""
        self.__vcf.add_information_fields([
            ['RefCodon', None, 'String', 'Reference codon'],
            ['AltCodon', None, 'String', 'Alternate codon'],
            ['RefAminoAcid', None, 'String', 'Reference amino acid'],
            ['AltAminoAcid', None, 'String', 'Alternate amino acid'],
            ['CodonPosition', '1', 'Integer', 'Codon position in the gene'],
            ['SNPCodonPosition', '1', 'Integer', 'SNP position in the codon'],
            ['AminoAcidChange', None, 'String', 'Amino acid change'],
            ['IsSynonymous', '1', 'Integer',
             '0:nonsynonymous, 1:synonymous, 9:N/A or Unknown'],
            ['IsTransition', '1', 'Integer',
             '0:transversion, 1:transition, 9:N/A or Unknown'],
            ['IsGenic', '1', 'Integer', '0:intergenic, 1:genic'],
            ['IsPseudo', '1', 'Integer', '0:not pseudo, 1:pseudo gene'],
            ['LocusTag', None, 'String', 'Locus tag associated with gene'],
            ['Gene', None, 'String', 'Name of gene'],
            ['Note', None, 'String', 'Note associated with gene'],
            ['Inference', None, 'String', 'Inference of feature.'],
            ['Product', None, 'String', 'Description of gene'],
            ['ProteinID', None, 'String', 'Protein ID of gene'],
            ['Comments', None, 'String', 'Example: Negative strand: T->C'],
            ['VariantType', None, 'String', 'Indel, SNP, Ambiguous_SNP'],
            ['FeatureType', None, 'String', 'The feature type of variant.'],
        ])

    def annotate_vcf_records(self):
        """Annotate each record in the VCF acording to the input GenBank."""
        # mutation_freq_list = []
        for record in self.__vcf.records:
            self.__gb.accession = record.CHROM
            self.__gb.version = record.CHROM
            self.__gb.index = record.POS

            # Set defaults
            record.INFO['RefCodon'] = '.'
            record.INFO['AltCodon'] = '.'
            record.INFO['RefAminoAcid'] = '.'
            record.INFO['AltAminoAcid'] = '.'
            record.INFO['CodonPosition'] = '.'
            record.INFO['SNPCodonPosition'] = '.'
            record.INFO['AminoAcidChange'] = '.'
            record.INFO['IsSynonymous'] = 9
            record.INFO['IsTransition'] = 9
            record.INFO['Comments'] = '.'
            record.INFO['IsGenic'] = '0'
            record.INFO['IsPseudo'] = '0'
            record.INFO['LocusTag'] = '.'
            record.INFO['Gene'] = '.'
            record.INFO['Note'] = '.'
            record.INFO['Inference'] = '.'
            record.INFO['Product'] = '.'
            record.INFO['ProteinID'] = '.'
            record.INFO['FeatureType'] = 'inter_genic'

            # Get annotation info
            if self.__gb.feature_exists:
                record.INFO['FeatureType'] = self.__gb.feature.type
                if self.__gb.feature.type in self.__annotated_features:
                    feature = self.__gb.feature
                    if feature.type == "CDS":
                        record.INFO['IsGenic'] = '1'

                    qualifiers = {
                        'Note': 'note', 'LocusTag': 'locus_tag',
                        'Gene': 'gene', 'Product': 'product',
                        'ProteinID': 'protein_id',
                        'Inference': 'inference'
                    }

                    if feature.type == "tRNA":
                        qualifiers['Note'] = 'anticodon'
                    for k, v in qualifiers.items():
                        if v in feature.qualifiers:
                            # Spell out semi-colons, commas and spaces
                            record.INFO[k] = feature.qualifiers[v][0].replace(
                                ';', '[semi-colon]'
                            ).replace(
                                ',', '[comma]'
                            ).replace(
                                ' ', '[space]'
                            )
                            if v == 'anticodon':
                                record.INFO[k] = 'anticodon{0}'.format(
                                    record.INFO[k]
                                )

                    if 'pseudo' in feature.qualifiers:
                        record.INFO['IsPseudo'] = '1'

            # Determine variant type
            if record.is_indel:
                if record.is_deletion:
                    record.INFO['VariantType'] = 'Deletion'
                else:
                    record.INFO['VariantType'] = 'Insertion'
            elif ("-" in record.ALT) or("-" in record.REF) :
                    record.INFO['VariantType'] = 'Indel'
            else:
                if len(record.ALT) > 1:
                    record.ALT = self.__gb.determine_iupac_base(record.ALT)
                    record.INFO['VariantType'] = 'Ambiguous_SNP'
                else:
                    if record.is_transition:
                        record.INFO['IsTransition'] = 1
                    else:
                        record.INFO['IsTransition'] = 0
                    record.INFO['VariantType'] = 'SNP'

                if int(record.INFO['IsGenic']):
                    alt_base = str(record.ALT[0])

                    # Determine codon information
                    codon = self.__gb.codon_by_position(record.POS)  ###Modified by MSM to fit the Genbank .gbff file
                    # print(record.POS, codon)
                    record.INFO['RefCodon'] = ''.join(list(codon[0]))
                    record.INFO['SNPCodonPosition'] = codon[1]
                    record.INFO['CodonPosition'] = codon[2]
                    # # MODIFICATION  ###SILENCED BY MSM
                    # # Instead of base from ref sequence use base from vcf file
                    # ref_base = str(record.REF[0])
                    # record.INFO['RefCodon'] = list(record.INFO['RefCodon'])
                    # record.INFO['RefCodon'][
                    #     record.INFO['SNPCodonPosition']
                    # ] = ref_base
                    # record.INFO['RefCodon'] = ''.join(record.INFO['RefCodon'])
                    # # END MODIFICATION

                    # Adjust for ambiguous base and negative strand.
                    if feature.strand == -1:
                        alt_base = str(
                            Seq(alt_base).complement()
                        )

                        record.INFO['Comments'] = 'Negative:{0}->{1}'.format(
                            Seq(record.REF).complement(),
                            alt_base
                        )

                    # Determine alternates
                    record.INFO['AltCodon'] = list(record.INFO['RefCodon'])
                    record.INFO['AltCodon'][
                        record.INFO['SNPCodonPosition']
                    ] = alt_base
                    record.INFO['AltCodon'] = ''.join(record.INFO['AltCodon'])
                    record.INFO['RefAminoAcid'] = Seq(
                        record.INFO['RefCodon']
                    ).translate()
                    record.INFO['AltAminoAcid'] = Seq(
                        record.INFO['AltCodon']
                    ).translate()
                    record.INFO['AminoAcidChange'] = '{0}{1}{2}'.format(
                        str(record.INFO['RefAminoAcid']),
                        record.INFO['CodonPosition'],
                        str(record.INFO['AltAminoAcid'])
                    )

                    if record.INFO['VariantType'] != 'Ambiguous_SNP':
                        ref = str(record.INFO['RefAminoAcid'])
                        alt = str(record.INFO['AltAminoAcid'])
                        if ref == alt:
                            record.INFO['IsSynonymous'] = 1
                        else:
                            record.INFO['IsSynonymous'] = 0
    
    def annotate_coverage(self):
        """Annotate the coverage file."""
        if self.__coverage_file:
            self.__coverage_file = pd.read_csv(
                self.__coverage_file, sep='\t'
            )
            self.__coverage_file.columns = ["Chromosome", "Position", "Coverage"]
            self.__coverage_file['CodonPosition'] = '.'
            # print(self.__coverage_file)
            for i in range(self.__coverage_file.shape[0]):
                self._accession = self.__coverage_file.iloc[i]['Chromosome']
                self._index = self.__coverage_file.iloc[i]['Position']
                # print(self.__coverage_file.iloc[i]['Position'], i)
                # print(self.__gb.codon_by_position(self.__coverage_file.iloc[i]['Position'])[2])
                ## we set CodonPosition to the value of the position in the coverage file
                self.__coverage_file.at[i, 'CodonPosition'] = self.__gb.codon_by_position(
                    self.__coverage_file.iloc[i]['Position'])[2]
                # print("Pos:", self.__coverage_file['Position'][i], self.__coverage_file['CodonPosition'][i])
                if i == 4308:
                    break



    def write_vcf(self, output='/dev/stdout'):
        """Write the VCF to the specified output."""
        self.__vcf.write_vcf(output)
    
    def write_tsv(self, output='/dev/stdout'): ###ADDED BY MSM TO GET THE TSV FILE WE WANT
        """Write the VCF mutations to a given output file in .tsv format."""
        self.__vcf.write_tsv(output)

    def write_coverage(self, output='/dev/stdout'):
        """Write the coverage file to the specified output by proteins."""
   
        chr = self.__coverage_file['Chromosome'][0]
        self.__coverage_file = self.__coverage_file.drop(['Chromosome', 'Position'], axis=1)
        self.__coverage_file = self.__coverage_file.groupby('CodonPosition', as_index = False).agg(lambda x: 0 if 0 in x.values else x.mean())
        self.__coverage_file['Subtype'] = chr
        self.__coverage_file = self.__coverage_file[['Subtype', 'CodonPosition', 'Coverage']]
        self.__coverage_file = self.__coverage_file[self.__coverage_file['CodonPosition'] != '.']
        self.__coverage_file['CodonPosition'] = self.__coverage_file['CodonPosition'].astype(int)
        self.__coverage_file["Protein"] = "Unknown"
        for n in range(self.__coverage_file.shape[0]):
            if self.__coverage_file['CodonPosition'][n] > 488 and self.__coverage_file['CodonPosition'][n] < 588: ##Protease aminoacids
                self.__coverage_file.at[n, 'Protein'] = "PR"
                self.__coverage_file.at[n, 'CodonPosition'] = self.__coverage_file['CodonPosition'][n] - 488
            elif self.__coverage_file['CodonPosition'][n] > 587 and self.__coverage_file['CodonPosition'][n] < 1028: ##RT aminoacids
                self.__coverage_file.at[n, 'Protein'] = "RT"
                self.__coverage_file.at[n, 'CodonPosition'] = self.__coverage_file['CodonPosition'][n] - 587
            elif self.__coverage_file['CodonPosition'][n] > 1147 and self.__coverage_file['CodonPosition'][n] < 1436: ##IN aminoacids
                self.__coverage_file.at[n, 'Protein'] = "IN"
                self.__coverage_file.at[n, 'CodonPosition'] = self.__coverage_file['CodonPosition'][n] - 1147

        self.__coverage_file.to_csv(output, sep='\t', index=False)


class GenBank(object):
    """A class for parsing GenBank files."""

    def __init__(self, gb=False):
        """Inititalize variables."""
        self.genbank_file = gb
        self.records = {}
        self.record_index = {}
        self.record_ids = {}
        self.__gb = None
        self._index = None
        self._accession = None
        self.__position_index = None
        self.feature = None
        self.features = ["CDS", "rRNA", "tRNA", "ncRNA", "repeat_region",
                         "misc_feature"]
        self.gene_codons = {}
        self.parse_genbank()

    @property
    def accession(self):
        """Accession for records."""
        return self._index

    @accession.setter
    def accession(self, value):

        if value not in self.records:
            # print("self.record_ids[value]", self.record_ids)
            value = self.record_ids[value]

        self._accession = value
        self.__gb = self.records[value]
        self.__position_index = self.record_index[value]

    @property
    def index(self):
        """Postion index for features."""
        return self._index

    @index.setter
    def index(self, value):
        try:
            self._index = self.__position_index[value - 1]
        except:
            self._index = None
        self.__set_feature()

    def parse_genbank(self):
        with open(self.genbank_file, 'r') as gb_fh:
            for record in SeqIO.parse(gb_fh, 'genbank'):
                self.records[record.name] = record
                self.gene_codons[record.name] = {}
                self.record_index[record.name] = [None] * len(record.seq)
                self.record_ids[record.id] = record.name
                for i in range(len(record.features)):
                    if record.features[i].type in self.features:
                        start = int(record.features[i].location.start)
                        end = int(record.features[i].location.end)
                        self.record_index[record.name][start:end] = [i] * (end - start)

    def __set_feature(self):
        if self._index is None:
            self.feature_exists = False
            self.feature = None
        else:
            self.feature_exists = True
            self.feature = self.records[self._accession].features[self._index]

    def codon_by_position(self, pos):
        """Retrieve the codon given a postion of a CDS feature."""
        # if self._index not in self.gene_codons[self._accession]: ###COMMENTED OUT BY MSM  
        #     self.split_into_codons()
        self.split_into_codons()
        gene_position = self.position_in_gene(pos) 
        codon_position = gene_position // 3

        # translation = []
        # for cod in self.gene_codons[self._accession][self._index]:
        #     if len(cod) == 3:
        #         translation.append(Seq(cod).translate())
        #     else:
        #         print("Codon length is not 3: ", cod)
        # print(translation)

        # print("Length:", len(self.gene_codons[self._accession][self._index]))

        return [self.gene_codons[self._accession][self._index][codon_position],  
                gene_position % 3,
                codon_position + 1]

    def split_into_codons(self):
        """Split the complete CDS feature in to a list of codons."""
        start = self.feature.location.start
        end = self.feature.location.end
        seq = ''.join(list(self.__gb.seq[start:end]))
        ###MODIFIED BY MSM TO FIX DIFFERENCE BETWEEN ONE REFERENCE AND THE OTHER
        ##we add in this seq an 'a' on the position 1302
        seq = seq[:1302] + 'A' + seq[1302:]
        ###END MODIFICATION MSM

        if self.feature.strand == -1:
            seq = Seq(seq).reverse_complement()

        self.gene_codons[self._accession][self._index] = [
            seq[i:i + 3] for i in range(0, len(seq), 3)
        ]

    def position_in_gene(self, pos):
        """Return a codon postion in a gene."""
        if self.feature.strand == 1:
            # return pos - self.feature.location.start - 1 ###MODIFIED BY MSM AS THIS IS ALREADY ACCOUNTED FOR IN THE REFERENCE
            return pos -2
        else:
            return self.feature.location.end - pos

    def base_by_pos(self, pos):
        """Print the base by position."""
        print(self.__gb.seq[pos - 1])

    def determine_iupac_base(self, bases):
        """
        Determine the IUPAC symbol for a list of nucleotides.
        Source: https://en.wikipedia.org/wiki/Nucleic_acid_notation
        List elements are in this order: [A,C,G,T]
        """
        if len(bases) > 1:
            iupac_notation = {
                'W': [True, False, False, True],
                'S': [False, True, True, False],
                'M': [True, True, False, False],
                'K': [False, False, True, True],
                'R': [True, False, True, False],
                'Y': [False, True, False, True],
                'B': [False, True, True, True],
                'D': [True, False, True, True],
                'H': [True, True, False, True],
                'V': [True, True, True, False],
                'N': [False, False, False, False]
            }

            base_condition = [base in bases for base in ['A', 'C', 'G', 'T']]
            for symbol, iupac_condition in iupac_notation.items():
                if iupac_condition == base_condition:
                    return symbol

    def is_transition(self, ref_base, alt_base):
        """
        Identify SNP as being a transition or not.
        1: Transition, 0:Transversion
        """
        substitution = ref_base + alt_base
        transition = ['AG', 'GA', 'CT', 'TC']

        if substitution in transition:
            return 1

        return 0


class VCFTools(object):
    """A class for parsing VCF formatted files."""

    def __init__(self, vcf_file):
        """Initialize variables."""
        self.reader = vcf.Reader(open(vcf_file, 'r'))
        self.records = [record for record in self.reader]

    def add_information_fields(self, info_list):
        """Add a given list of information fields to the VCF."""
        for i in info_list:
            id, num, type, desc = i
            self.__add_information_field(id, num, type, desc)

    def __add_information_field(self, id, num, type, desc):
        """Add a given information field to the VCF."""
        _Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc'])
        if id:
            self.reader.infos[id] = _Info(id, num, type, desc)

    def write_vcf(self, output='/dev/stdout'):
        """Write the VCF to a given output file."""
        vcf_writer = vcf.Writer(open(output, 'w'), self.reader)
        for record in self.records:
            vcf_writer.write_record(record)
    
    def write_tsv(self, output='/dev/stdout'): ###ADDED BY MSM TO GET THE TSV FILE WE WANT
        """Write the VCF mutations to a given output file in .tsv format."""
        df_list = []
        for record in self.records:
            if record.INFO['CodonPosition'] > 488 and record.INFO['CodonPosition'] < 588: ##Protease aminoacids
                df_list.append([record.INFO['CodonPosition']-488, record.INFO['RefAminoAcid'], record.INFO['AltAminoAcid'], record.INFO['AF'], "PR"])
            elif record.INFO['CodonPosition'] > 587 and record.INFO['CodonPosition'] < 1028: ##RT aminoacids ###827 is the last aminoacid from our ref sequence, but on the .gb the final aa is 1027.
                df_list.append([record.INFO['CodonPosition']-587, record.INFO['RefAminoAcid'], record.INFO['AltAminoAcid'], record.INFO['AF'], "RT"])
            elif record.INFO['CodonPosition'] > 1147 and record.INFO['CodonPosition'] < 1436: ##IN aminoacids
                df_list.append([record.INFO['CodonPosition']-1147, record.INFO['RefAminoAcid'], record.INFO['AltAminoAcid'], record.INFO['AF'], "IN"])

        df = pd.DataFrame(df_list, columns=['Position', 'Ref', 'Mut', 'Freq', 'Prot'])
        df.to_csv(output, index=False, sep='\t')

def update_vcf_chrom(in_vcf, out_vcf, chrom_name):
    vf = pyvcf.VcfFrame.from_file(in_vcf)
    vf.df['CHROM'] = chrom_name
    #vf.df = vf.df[vf.df['ALT']!='-']
    vf.to_file(out_vcf)

def update_coverage_chrom(in_cov, out_cov, chrom_name):
    cov = pd.read_csv(in_cov, sep="\t")
    #we remove the ref column and add a new one in first position
    # cov = cov.drop(columns=['ref'])
    cov['ref'] = chrom_name

    # cov_columns = cov.columns.tolist()
    # cov = cov[[cov_columns[-1]] + cov_columns[:-1]]
    # print(cov)
    cov.to_csv(out_cov, sep="\t", index=False)

def main(fname_snv_in, fname_cov_in, fname_genbank_file, fname_snv_out, fname_cov_out):

    # chrom_name = 'NC_045512.2'
    chrom_name = 'NC_001802.1' ##CHANGED MSM

    # sample = str(fname_snv_out).split("/variants")[0].split("/")[-4] SILENCED MSM
    fname_snv_temp = str(fname_snv_in).split('.vcf')[0]+'.temp.vcf'
    fname_cov_temp = str(fname_cov_in).split('.tsv')[0]+'.temp.tsv'

    # update_vcf_chrom(fname_snv_in, fname_snv_temp, chrom_name)
    update_vcf_chrom(fname_snv_in, fname_snv_temp, chrom_name)
    update_coverage_chrom(fname_cov_in, fname_cov_temp, chrom_name)

    # for vcf_file, coverage_file, out_vcf_file, out_coverage_file in zip(fname_snv_in, fname_cov_in, fname_snv_out, fname_cov_out):
    #     # annotator = Annotator(gb_file=fname_genbank_file, vcf_file=fname_snv_in + '/' + sample_id + '/mix_12_variants_chromchange.vcf', coverage_file=fname_snv_in + '/' + sample_id + '/coverage.tsv')
    #     annotator = Annotator(gb_file=fname_genbank_file, vcf_file=vcf_file, coverage_file=coverage_file)
    #     annotator.annotate_vcf_records()
    #     annotator.annotate_coverage()
    #     # annotator.write_vcf(fname_snv_out)
    #     if not os.path.exists(fname_snv_out):
    #         os.makedirs(fname_snv_out)
    #     # os.makedirs(fname_snv_out + '/' + sample_id, exist_ok=True)
    #     # annotator.write_tsv(fname_snv_out + '/' + sample_id + '/mut_freq.tsv')
    #     # annotator.write_coverage(fname_snv_out + '/' + sample_id + '/coverage_annotated.tsv')
    #     annotator.write_tsv(out_vcf_file)
    #     annotator.write_coverage(out_coverage_file)
    annotator = Annotator(gb_file=fname_genbank_file, vcf_file=fname_snv_temp, coverage_file=fname_cov_temp)
    annotator.annotate_vcf_records()
    annotator.annotate_coverage()
    annotator.write_tsv(fname_snv_out)
    annotator.write_coverage(fname_cov_out)

    os.remove(fname_snv_temp), os.remove(fname_cov_temp) #we remove the temp files created



if __name__ == "__main__":
    main(
        snakemake.input.fname_vcf,
        snakemake.input.fname_cov,
        snakemake.input.ref,
        # snakemake.input.fname_output
        # snakemake.input.fname_cov,
        snakemake.output.fname_vcf,
        snakemake.output.fname_cov
        # "./example_files/mix_12_variants_chromchange.vcf",
        # "./example_files/CAP257/week_54/alignments/coverage.tsv",
        # "./ref/reference_gagpol_only.gb",
        # "./example_files/mutation_freq.tsv",
        # "./example_files/coverage_annotated.tsv"
    )