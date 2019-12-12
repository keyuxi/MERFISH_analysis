# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 13:27:36 2019

@author: Yuxi
"""

import pybiomart as bm
import pandas as pd
import numpy as np
import requests, sys
import json
import tqdm

class Codebook(object):
    def __init__(self, gene_list_file, codebook_name='',
                 readout_file=None, readout_offset=0,
                 bulk_seq_file=None, bulk_seq_unit='fkpm',
                 code_file=None, shuffle_code=True, verbose=True):
        
        ### Read gene list ###
        self.gene_list = pd.read_csv(gene_list_file, header=0, sep='\t', encoding='utf-8')

        if readout_file == None:
            self.readout_file = r"C:\Users\Yuxi\workspace\genomeData\CommonOligos"
        self.verbose = verbose
    
    def generate(self, out_file=None):
        pass

    def get_ensembl_id(self):
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
            ext = "/xrefs/symbol/mus_musculus/" + gene_symbol.upper() \
                    + "?object_type=gene"
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
            
            if not r.ok:
                r.raise_for_status()
                sys.exit()
               
            # A dictionary of the result
            if self.verbose:
                print("\n===Ensembl ID query results for %s===\n" % gene_symbol)
                print(r.json()[0])
            ensembl_id = r.json()[0]['id']
            
            # Write to the table
            self.gene_list.at[gene_ind, 'ensembl_gene_id'] = ensembl_id
            if self.verbose:
                print(self.gene_list.iloc[gene_ind])

    
    def get_attributes(self, attributes=None,
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
        attr = ['mgi_symbol', 'ensembl_gene_id', 'ensembl_gene_id_version', \
                'ensembl_transcript_id', 'ensembl_transcript_id_version', \
                'transcript_appris', 'transcript_length', \
                'gene_biotype', 'transcript_count']
        filters = {'link_ensembl_gene_id': self.gene_list['ensembl_gene_id'].tolist()}
        
        # Retrieve information
        query_result = dataset.query(attributes=attr, filters = filters)
        
        # Filter out protein coding entries without an Appris annotation
        query_result[]
        
        return query_result


    
    def assign_smELT(self):
        # create a new column by copying
        self.gene_list['is_smELT'] = self.gene_list['is_force_smELT']
    
    def assign_barcode(self):
        pass
    

if __name__ == "__main__":
    #%%
    gene_list_file = r"C:\Users\Yuxi\workspace\probedesign\gene_list\FISH_markers.txt"
    # readout_file = r"C:\Users\Yuxi\workspace\genomeData\CommonOligos\ReadoutSeqs.fasta"
    out_file = r"C:\Users\Yuxi\workspace\probedesign\gene_list\codebook_BLA.csv"
    codebook = Codebook(gene_list_file, "BLA")
    codebook.get_ensembl_id()
    query_result = codebook.get_attributes()
    print(query_result)
    #codebook.generate(out_file)