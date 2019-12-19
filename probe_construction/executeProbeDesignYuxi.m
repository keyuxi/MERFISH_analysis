%% Run from scratch
clear;clc;
codebookPath = 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\codebook_Musmusculus_BLA_v2.0.csv';
pd = probeDesign('BLA_lib_v2', 'mouse', codebookPath);

set(pd, 'MERFISHAnalysisPath', 'C:\Users\Yuxi\workspace\genomeData');
%set(pd, 'basePath', 'C:\Users\Yuxi\workspace\probedesign\MERFISH_out_');
set(pd, 'basePath', 'C:\Users\Yuxi\workspace\genomeData');
set(pd, 'fpkmPath', 'Mus_musculus_proxy.fpkm_tracking');
set(pd, 'readoutPath', 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\ReadoutSeqs.fasta');

set(pd, 'regionGC', [0.43, 0.63], 'regionTm', [66,76], 'isoSpecificity', [0, 1], 'specificity', [0.75, 1]);
set(pd, 'FPKMabundanceThreshold', 0, 'numProbesPerGene', 92);
set(pd, 'probeSpacing', -20, 'doubleHeadedsmELT', true);
pd.buildLibrary();

%% Run from existing transcriptome objects
clear;clc;
addpath(genpath('C:\Users\Yuxi\workspace\MERFISH_analysis'));
addpath(genpath('C:\Users\Yuxi\workspace\genomeData'));
codebookPath = 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\codebook_Musmusculus_BLA_v2.0.csv';
pd = probeDesign('BLA_lib_v2', 'mouse', codebookPath);

set(pd, 'MERFISHAnalysisPath', 'C:\Users\Yuxi\workspace\genomeData');
set(pd, 'basePath', 'C:\Users\Yuxi\workspace\genomeData');
set(pd, 'fpkmPath', 'Mus_musculus_proxy.fpkm_tracking');
set(pd, 'rRNAtRNAPath', 'C:\Users\Yuxi\workspace\genomeDatahypothalamusLibrary\rRNAtRNA.fa');
set(pd, 'transcriptomePath', 'C:\Users\Yuxi\workspace\genomeDatahypothalamusLibrary\transcriptomeObj');
set(pd, 'readoutPath', 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\ReadoutSeqs.fasta');

set(pd, 'regionGC', [0.43, 0.63], 'regionTm', [66,76], 'isoSpecificity', [0, 1], 'specificity', [0.75, 1]);
set(pd, 'FPKMabundanceThreshold', 0, 'numProbesPerGene', 92);
set(pd, 'probeSpacing', -20, 'doubleHeadedsmELT', true);
pd.buildLibrary();