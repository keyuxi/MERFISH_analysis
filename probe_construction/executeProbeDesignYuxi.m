%% lib01_merfish: Run from scratch
clear;clc;
codebookPath = 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\codebook_Musmusculus_lib01_merfish_v1.0.csv';
pd = probeDesign('lib01_merfish', 'mouse', codebookPath);

set(pd, 'MERFISHAnalysisPath', 'C:\Users\Yuxi\workspace\genomeData');
set(pd, 'basePath', 'C:\Users\Yuxi\workspace\genomeData');
set(pd, 'fpkmPath', 'Mus_musculus_proxy.fpkm_tracking');
set(pd, 'readoutPath', 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\AllReadouts.fasta');

set(pd, 'regionGC', [0.43, 0.63], 'regionTm', [66,76], 'isoSpecificity', [0, 1], 'specificity', [0.75, 1]);
set(pd, 'FPKMabundanceThreshold', 0, 'numProbesPerGene', 92);
set(pd, 'probeSpacing', -20, 'tripleHeadedsmELT', true);

set(pd, 'isPredesignedPrimer', true, 'primerID', [1 1]);

pd.buildLibrary();

%% lib01_merfish: Run from existing transcriptome objects
clear;clc;
addpath(genpath('C:\Users\Yuxi\workspace\MERFISH_analysis'));
addpath(genpath('C:\Users\Yuxi\workspace\genomeData'));
codebookPath = 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\codebook_Musmusculus_lib01_merfish_v1.0.csv';
pd = probeDesign('lib01_merfish', 'mouse', codebookPath);

set(pd, 'MERFISHAnalysisPath', 'C:\Users\Yuxi\workspace\genomeData');
set(pd, 'basePath', 'C:\Users\Yuxi\workspace\genomeData');
set(pd, 'fpkmPath', 'Mus_musculus_proxy.fpkm_tracking');
set(pd, 'rRNAtRNAPath', 'C:\Users\Yuxi\workspace\genomeDatahypothalamusLibrary\rRNAtRNA.fa');
set(pd, 'transcriptomePath', 'C:\Users\Yuxi\workspace\genomeDatahypothalamusLibrary\transcriptomeObj');
set(pd, 'readoutPath', 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\AllReadouts.fasta');

set(pd, 'regionGC', [0.43, 0.63], 'regionTm', [66,76], 'isoSpecificity', [0, 1], 'specificity', [0.75, 1]);
set(pd, 'FPKMabundanceThreshold', 0, 'numProbesPerGene', 92);
set(pd, 'probeSpacing', -20, 'tripleHeadedsmELT', true);

set(pd, 'isPredesignedPrimer', true, 'primerID', [1 1]);

pd.buildLibrary();

%% lib01_2hot: Run from scratch
clear;clc;
addpath(genpath('C:\Users\Yuxi\workspace\MERFISH_analysis'));
addpath(genpath('C:\Users\Yuxi\workspace\genomeData'));

codebookPath = 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\codebook_Musmusculus_lib01_2hot_v1.0.csv';
pd = probeDesign('lib01_2hot', 'mouse', codebookPath);

set(pd, 'MERFISHAnalysisPath', 'C:\Users\Yuxi\workspace\genomeData');
set(pd, 'basePath', 'C:\Users\Yuxi\workspace\genomeData');
set(pd, 'fpkmPath', 'Mus_musculus_proxy.fpkm_tracking');
set(pd, 'readoutPath', 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\AllReadouts.fasta');

set(pd, 'regionGC', [0.43, 0.63], 'regionTm', [66,76], 'isoSpecificity', [0, 1], 'specificity', [0.75, 1]);
set(pd, 'FPKMabundanceThreshold', 0, 'numProbesPerGene', 48);
set(pd, 'probeSpacing', -20, 'tripleHeadedsmELT', true);

set(pd, 'isPredesignedPrimer', true, 'primerID', [2 2]);

pd.buildLibrary();