%% lib01_merfish
clear;clc;
lib_name = 'lib01_merfish';
codebookPath = 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\codebook_Musmusculus_lib01_merfish_v1.0.csv';
pd = probeDesign('lib01_merfish', 'mouse', codebookPath);

set(pd, 'MERFISHAnalysisPath', 'C:\Users\Yuxi\workspace\genomeData');
set(pd, 'basePath', 'C:\Users\Yuxi\workspace\genomeData');
set(pd, 'fpkmPath', 'Mus_musculus_proxy.fpkm_tracking');
set(pd, 'readoutPath', ['C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\used_readouts_' lib_name '.fasta']);

set(pd, 'regionGC', [0.43, 0.63], 'regionTm', [66,76], 'isoSpecificity', [0, 1], 'specificity', [0.75, 1]);
set(pd, 'FPKMabundanceThreshold', 0, 'numProbesPerGene', 92);
set(pd, 'probeSpacing', -20, 'tripleHeadedsmELT', true);

set(pd, 'isPredesignedPrimer', true, 'primerID', [1 1]);

pd.buildLibrary();

%% lib01_2hot
clear;clc;
lib_name = 'lib01_2hot';
addpath(genpath('C:\Users\Yuxi\workspace\MERFISH_analysis'));
addpath(genpath('C:\Users\Yuxi\workspace\genomeData'));

codebookPath = 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\codebook_Musmusculus_lib01_2hot_v1.0.csv';
pd = probeDesign('lib01_2hot', 'mouse', codebookPath);

set(pd, 'MERFISHAnalysisPath', 'C:\Users\Yuxi\workspace\genomeData');
set(pd, 'basePath', 'C:\Users\Yuxi\workspace\genomeData');
set(pd, 'fpkmPath', 'Mus_musculus_proxy.fpkm_tracking');
set(pd, 'readoutPath', ['C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\used_readouts_' lib_name '.fasta']);

set(pd, 'regionGC', [0.43, 0.63], 'regionTm', [66,76], 'isoSpecificity', [0, 1], 'specificity', [0.75, 1]);
set(pd, 'FPKMabundanceThreshold', 0, 'numProbesPerGene', 48);
set(pd, 'probeSpacing', -20, 'tripleHeadedsmELT', true);

set(pd, 'isPredesignedPrimer', true, 'primerID', [2 2]);

pd.buildLibrary();

%% Fill into order template
fn1 = 'C:\Users\Yuxi\workspace\genomeDatalib01_merfish\lib01_merfish_oligos.fasta';
fn2 = 'C:\Users\Yuxi\workspace\genomeDatalib01_2hot\lib01_2hot_oligos.fasta';
template = 'C:\Users\Yuxi\workspace\MERFISH_analysis\OutputForLIMS\Genscript_Schnitzer_lab_v0.xlsx';
out = 'C:\Users\Yuxi\workspace\MERFISH_analysis\OutputForLIMS\Genscript_Schnitzer_lab_v1.xlsx';

oligos = [fastaread(fn1); fastaread(fn2)];
oligos = struct2cell(oligos);
oligos = oligos(2,:)';
oligos = cellfun(@(x) x(~isspace(x)), oligos, 'UniformOutput', false);

status = copyfile(template, out);
xlswrite(out, oligos, 'Oligo pool sequence form', 'A3');
xlswrite(out, length(oligos), 'Quatition request form', 'I14');