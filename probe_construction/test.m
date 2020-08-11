clear;clc;
regionProps = zeros(6,8);
regionProps(1,:) = [1,5,27,35,57,78,82,104];
regionProps(2,:) = 20;
padLength = 0;
selectedRegionData = TRDesigner.TileSplitRegions(regionProps, padLength);
expected = [5,27,82,104];
fprintf('Expected v.s. returned:\n')
disp(expected)
disp(selectedRegionData(1,:))

%%
bridgeFile = 'G:\My Drive\S_lab\reactions\literature_supp\2020_split_FISH\splitFISH_bridge.xlsx';
bridge = readtable(bridgeFile);

fastawrite('G:\My Drive\S_lab\reactions\literature_supp\2020_split_FISH\AllBridges.fasta', ...
    bridge.Bit, bridge.BridgeSequence);

%%
bcseq = [];
func = @(seq) seqrcomplement(seq(isstrprop(seq, 'lower')));
for i = 1:26
    bcseq{i} = func(bridge.BridgeSequence{i});
end

fastawrite('G:\My Drive\S_lab\reactions\literature_supp\2020_split_FISH\All6ntBarcodes.fasta', ...
    bridge.Bit, bcseq);

%%
codebook = LoadCodebook('C:\Users\Yuxi\workspace\transcriptomic_analysis\codebook\out\codebook_mm_mer.csv')
