# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 13:27:36 2019

@author: Yuxi
"""

import pybiomart as bm
import pandas as pd
import numpy as np
import requests, sys
import csv

class Codebook(object):
    def __init__(self, gene_list_file, codebook_name='merfish', version='1.0',
                 readout_file=r".\AllReadouts.fasta", readout_offset=0,
                 bulk_seq_file=None, bulk_seq_cutoff=500,
                 code_file=None, is_shuffle_code=True, verbose=True):
        
        # Read gene list
        self.gene_list = pd.read_csv(gene_list_file, header=0, sep='\t', encoding='utf-8')

        self.codebook_name = codebook_name
        self.version = version
        
        # Readout file to determine bit names
        self.n_bit = 0
        self.bit_names = []

        with open(readout_file, 'r') as infile:
            for line in infile:
                if line.strip().startswith('>'):
                    self.bit_names.append(line.strip()[1:])
        self.bit_names = self.bit_names[readout_offset:]
        
        # Abundance information to assign smELT
        self.bulk_seq_file = bulk_seq_file
        self.bulk_seq_cutoff = bulk_seq_cutoff
            
        # Hamming code file
        if code_file == None:
            self.code_file = \
                r'.\16bit.mhd4.txt'
        else:
            self.code_file = code_file
        self.is_shuffle_code = is_shuffle_code
            
        self.verbose = verbose
    
    def generate(self, out_file=None):
        """
        Call everything and write to out_file.
        """
        self._get_ensembl_id()
        self._get_attributes()
        if not (self.bulk_seq_file == None):
            self._get_abundance()
        else:
            self.gene_list['abundance'] = 1 # Arbitrary default

        self._assign_smELT()
        codebook_df = self._assign_barcode()
        
        self._write_codebook_csv(codebook_df, out_file)

    def _get_ensembl_id(self):
        """
        Convert gene names to ensembl gene id using the REST API.
        TODO: Check gene lost.
        """
        # Create a new column
        self.gene_list['ensembl_gene_id'] = ""
        
        # Server config
        server = "http://rest.ensembl.org"

        # Loop over each gene
        for gene_ind, gene_symbol in enumerate(self.gene_list['mgi_symbol'].tolist()):
            
            ensembl_id = None
            
            # Request by gene symbol
            ext = "/xrefs/symbol/mus_musculus/" + gene_symbol.upper() \
                    + "?object_type=gene"
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
            
            if not r.ok:
                r.raise_for_status()
                sys.exit()
            
            if self.verbose:
                print("\n===Ensembl ID query results for %s===\n" % gene_symbol)

            # Check gene name backwards to prevent hitting other genes with 
            # same acronyms
            for result in r.json():
                # Request by gene stable id
                ext = "/xrefs/id/" + result['id'] + "?"
 
                r_id = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
                 
                if not r_id.ok:
                  r_id.raise_for_status()
                  sys.exit()
                
                # Grab MGI symbol
                for seg in r_id.json():
                    if seg['dbname'] == 'MGI':
                        this_gene_symbol = seg['display_id']
                        break
                    
                if self.verbose:
                    print(result)
                    print(this_gene_symbol)
                    
                # If the ID corresponds to the correct MGI symbol
                if gene_symbol == this_gene_symbol:
                    ensembl_id = result['id']
                    break
                
            # If failed to find matching ID
            if ensembl_id == None:
                print('\n%s\nFailed to find ID for %s, check MGI symbol spelling!\n'
                      % (30*'!', gene_symbol))
            else:    
                # Write this id to the table
                self.gene_list.at[gene_ind, 'ensembl_gene_id'] = ensembl_id
    
    def _get_attributes(self, attributes=None,
                       dataset_name='mmusculus_gene_ensembl'):
        """
        Get gene attributes and find principal transcripts.
        Called after ensembl_gene_id query. Dependent on pybiomart package.
        """
        # Set the dataset. Default to mouse genes.
        dataset = bm.Dataset(name=dataset_name,
                     host='http://www.ensembl.org')
        
        # Set the attributes and filters for query.
        # Some temporary hard-coding here.
        if attributes == None:
            attributes = ['mgi_symbol', 'ensembl_gene_id', 'ensembl_gene_id_version', \
                          'ensembl_transcript_id', 'ensembl_transcript_id_version', \
                          'transcript_appris', 'transcript_length', \
                          'gene_biotype', 'transcript_count']
        filters = {'link_ensembl_gene_id': self.gene_list['ensembl_gene_id'].tolist()}
        
        # Retrieve information
        query_result = dataset.query(attributes=attributes, filters = filters)
        
        #####################################################
        ##### Find the transcript to use for each gene. #####
        #####################################################
        
        # Create a new column for the chosen transcript
        self.gene_list['ensembl_transcript_id_version'] = ''
        
        # For lncRNA, choose the longest transcript.
        lnc_qr = query_result[query_result['Gene type'] == 'lncRNA']
        lnc_ind = lnc_qr.groupby(['Gene stable ID']) \
                ['Transcript length (including UTRs and CDS)'].idxmax()
        lnc_qr = query_result.loc[lnc_ind]
        for ind,row in lnc_qr.iterrows():
            self.gene_list.at[self.gene_list['ensembl_gene_id'] == row['Gene stable ID'],
                              'ensembl_transcript_id_version'] = row.loc['Transcript stable ID version']
            if self.verbose:
                print('\n%s' % self.gene_list[self.gene_list['ensembl_gene_id'] == row['Gene stable ID']])        
        
        # For protein coding genes, select the one with smallest APPRIS annotation
        prot_qr = query_result[query_result['Gene type'] == 'protein_coding']
        prot_qr['appris_rank'] = prot_qr['APPRIS annotation'].apply(self._appris2rank)
        # lowest APPRIS rank
        prot_ind = prot_qr.groupby(['Gene stable ID'])['appris_rank'].idxmin()
        prot_qr = query_result.loc[prot_ind]
        # longest
        prot_ind = prot_qr.groupby(['Gene stable ID'])['Transcript length (including UTRs and CDS)'].idxmax()
        prot_qr = query_result.loc[prot_ind]
        # Write the selected transcript
        for ind,row in prot_qr.iterrows():
            self.gene_list.at[self.gene_list['ensembl_gene_id'] == row['Gene stable ID'],
                              'ensembl_transcript_id_version'] = row.loc['Transcript stable ID version']
            self.gene_list.at[self.gene_list['ensembl_gene_id'] == row['Gene stable ID'],
                              'appris_rank'] = row.loc['APPRIS annotation']
            self.gene_list.at[self.gene_list['ensembl_gene_id'] == row['Gene stable ID'],
                              'transcript_length'] = row.loc['Transcript length (including UTRs and CDS)']

    def _get_abundance(self):
        """
        Temporarily takes in a specified format. I/O hardcoded.
        """
        
        # Warning: specific format
        abundance = pd.read_csv(self.bulk_seq_file, header=4, sep='\t', encoding='utf-8')
        abundance = abundance[['Gene ID', 'Gene Name', 'postnatal day 63, forebrain']]
        abundance.columns = ['Gene ID', 'Gene Name', 'tpm']
        
        # Map to Ensembl gene ID, add to the attribute
        self.gene_list['abundance'] = 0
        for ind,row in self.gene_list.iterrows():
            self.gene_list.at[ind, 'abundance'] = abundance[abundance['Gene ID'] == row['ensembl_gene_id']]['tpm'].values[0]
        
    def _assign_smELT(self):
        # Create a new column by copying
        is_smELT = self.gene_list['is_force_smELT']
        
        # Assign genes that are expressed above a given value to smELT
        is_smELT[self.gene_list['abundance'] >= self.bulk_seq_cutoff] = 1
        self.gene_list['is_smELT'] = is_smELT
        print(self.gene_list[['mgi_symbol', 'ensembl_gene_id', 'abundance', 'is_smELT']])
        
    def _assign_barcode(self):
        """
        Assign Hamming code or sequential code. Add Blank codes. Sort genes.
        Returns:
            out_list - a dataframe with gene symbols, transcript_id and barcode
        """
        #####################################################
        ############# Prepare Hamming codebook ##############
        #####################################################

        hamming_code = pd.read_csv(self.code_file, sep='\t',
                                    header=None, index_col=False)
        # Shuffle
        if self.is_shuffle_code:
            np.random.seed(2)
            hamming_code = hamming_code.sample(frac=1).reset_index(drop=True)
        
        # Convert to a list
        hamming_code = hamming_code.values.squeeze().tolist()
        self.n_bit = len(hamming_code[0].replace(' ',''))
        
        #####################################################
        ##### Generate sorted gene list for the codebook ####
        #####################################################
        
        # Barcoded genes
        # Sort gene names alphabetically to dataframe out_list
        barcoded_list = self.gene_list[self.gene_list['is_smELT'] == 0].sort_values('mgi_symbol') \
                        [['mgi_symbol','ensembl_transcript_id_version']]
        # Check gene number
        n_barcoded_gene = len(barcoded_list)
        print('\nNumber of barcoded genes: %s.\n' % n_barcoded_gene)
        if n_barcoded_gene > 140:
            print('Too many barcoded genes!')
        # Blank barcode names
        blank_name = [('Blank-%d' % (i+1)) for i in range(len(hamming_code) - n_barcoded_gene)]
        blank_id = ['' for i in range(len(hamming_code) - n_barcoded_gene)]
        
        # smELT
        smelt_list = self.gene_list[self.gene_list['is_smELT'] == 1].sort_values('mgi_symbol') \
                     [['mgi_symbol','ensembl_transcript_id_version']]
        n_smelt_gene = len(smelt_list)
        
        # Concatenate
        list_data = {'mgi_symbol':
                barcoded_list['mgi_symbol'].tolist() + blank_name + smelt_list['mgi_symbol'].tolist(),
                'ensembl_transcript_id_version':
                barcoded_list['ensembl_transcript_id_version'].tolist() + blank_id + smelt_list['ensembl_transcript_id_version'].tolist()}
        out_list = pd.DataFrame(data = list_data)
        
        #####################################################
        ################# Assign all barcodes ###############
        #####################################################
        # Assembly barcode - a list
        smelt_code = np.eye(n_smelt_gene,dtype=int).astype(str).tolist()
        smelt_code = [''.join(row) for row in smelt_code]
        barcode = [(row.replace(' ', '') + n_smelt_gene * '0') for row in hamming_code] \
                  + [(self.n_bit * '0' + row) for row in smelt_code]
                  
        # Truncate bit_names to the proper length
        self.bit_names = self.bit_names[:len(barcode[0])]
        
        # Add barcode to out_list
        out_list['barcode'] = barcode
        
        # Rename columns to match the final format
        out_list.columns = ['name', 'id', 'barcode']

        return out_list
        
    def _write_codebook_csv(self, codebook_df, out_file):
        # Default out_file name
        if out_file == None:
            out_file = r'.\codebook_Musmusculus_' + self.codebook_name + '_v' \
                       + self.version + '.csv'
        
        # Write
        with open(out_file, 'w', newline='') as cb_file:
            cb_writer = csv.writer(cb_file, delimiter=',')
            cb_writer.writerow(['version', self.version])
            cb_writer.writerow(['codebook_name', self.codebook_name])
            cb_writer.writerow(['bit_names'] + self.bit_names)
            cb_writer.writerow(['name', 'id', 'barcode'])
            for _,row in codebook_df.iterrows():
                cb_writer.writerow(row.tolist())
                
        if self.verbose:
            print('\nFinished writing codebook to %s.\n' % out_file)
        
    def _appris2rank(self, anno):
        """
        Helper function that converts APPRIS annotation to a number for comparison.
        The smaller, the better.
        e.g. 'principal1' -> 1, 'alternative2' -> 11
        Returns:
            r - int.
        """
        if not isinstance(anno, str):
            r = int(99)
        elif 'principal' in anno:
            r = int(anno[-1])
        elif 'alternative' in anno:
            r = int(anno[-1]) + 9
        else:
            r = int(99)
        return r

class Codebook2hot(Codebook):
    def _assign_barcode(self):
        """
        Assign 2-hot barcodes neglecting is_smELT. Sort genes.
        Returns:
            out_list - a dataframe with gene symbols, transcript_id and barcode
        """
        
        #####################################################
        ##### Generate sorted gene list for the codebook ####
        #####################################################
        
        # Barcoded genes
        # Sort gene names alphabetically to dataframe out_list
        out_list = self.gene_list.sort_values('mgi_symbol') \
                        [['mgi_symbol','ensembl_transcript_id_version']]
        # Check gene number
        n_barcoded_gene = len(out_list)
        print('\nNumber of barcoded genes: %s.\n' % n_barcoded_gene)
        if n_barcoded_gene > 100:
            print('Too many barcoded genes!')
                
        # Generate 2-hot code as np.array
        barcode = np.zeros((n_barcoded_gene, 2*n_barcoded_gene), dtype=int)
        for i in range(n_barcoded_gene):
            barcode[i, 2*i:2*i+2] = 1
        # Assembly barcode - a list
        barcode = barcode.astype(str).tolist()
        barcode = [''.join(row) for row in barcode]

        #####################################################
        ################# Assign all barcodes ###############
        #####################################################
                  
        # Truncate bit_names to the proper length
        self.bit_names = self.bit_names[:len(barcode[0])]
        
        # Add barcode to out_list
        out_list['barcode'] = barcode
        
        # Rename columns to match the final format
        out_list.columns = ['name', 'id', 'barcode']

        return out_list


if __name__ == "__main__":
    gene_list_file = r".\lib01_merfish.txt"
    bulk_seq_file = r".\E-MTAB-6798-query-results.tpms.tsv"
    codebook_merfish = Codebook(gene_list_file, "lib01_merfish", bulk_seq_file=bulk_seq_file, verbose=False)
    codebook_merfish.generate()

    gene_list_file = r".\lib01_2hot.txt"
    codebook_2hot = Codebook2hot(gene_list_file, "lib01_2hot", readout_offset=25, verbose=False)
    codebook_2hot.generate()