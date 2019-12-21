## Template sequence files

### AllBarcodes.fasta

Barcode sequences whose reverse compliment are to be concatenated to the primary probes. Same as the part on the adapters.

Derived from `AllAdapters.fasta`. Readout 053, 054 are replaced with 188, 189 in the original file. barcode188, barcode189 are manually deleted from the file.

Generated with the following Matlab snippet:

```matlab
adt = fastaread('AllAdapters_adpt2.fasta');
for i =1:length(adt)
    strparts = strsplit(adt(i).Header, '_');
    adt(i).Header = ['barcode' strparts{1}(end-2:end)];
    adt(i).Sequence = adt(i).Sequence(1:20);
end
fastawrite('AllBarcodes.fasta', adt);
```

