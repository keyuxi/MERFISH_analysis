# Codebook construction

## Dependencies

See `requirements.txt`.

## Example and key attributes

Example of using the `Codebook` class:

```python
gene_list_file = r".\gene_list_example.tsv"
bulk_seq_file = r".\E-MTAB-6798-query-results.tpms.tsv"
codebook_merfish = Codebook(gene_list_file, "lib01_merfish",bulk_seq_file=bulk_seq_file, verbose=False)

codebook_merfish.generate()
```

- `gene_list_file` is a tsv file with columns:
  - `mgi_symbol`: MGI symbols of the gene names. Will be checked by the script. If not the standard MGI symbol, manually confirm and change to the desired gene to avoid confusion.
  - `is_forced_smELT`: 0/1 logical, specify which genes are singled out for sequential barcodes.
- `bulk_seq_file` contains tpm abundance information for the genes. The script is temporarily hard-coded to read expression level at the adult forebrain. Genes with abundance >= `bulk_seq_cutoff` will be s

## Template sequence files

### AllReadouts.fasta

Barcode sequences whose reverse compliment are to be concatenated to the primary probes. Same as the part on the adapters.

Derived from `AllAdapters.fasta`. Readout 053, 054 are replaced with 188, 189 in the original file. barcode188, barcode189 are manually deleted from the file.

Generated with the following Matlab snippet:

```matlab
adt = fastaread('AllAdapters_adpt2.fasta');
for i =1:length(adt)
    strparts = strsplit(adt(i).Header, '_');
    adt(i).Header = ['readout' strparts{1}(end-2:end)];
    adt(i).Sequence = adt(i).Sequence(1:20);
end
fastawrite('AllReadouts.fasta', adt);
```

