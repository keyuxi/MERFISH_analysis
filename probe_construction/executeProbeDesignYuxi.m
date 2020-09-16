addpath(genpath('C:\Users\Yuxi\workspace\MERFISH_analysis'))

% clear;clc;
% designCodebook('mer', 'conventional');
% clear;clc;
% designCodebook('trihead', 'conventional');
% clear;clc;
% designCodebook('opioid', 'conventional');
% clear;clc;
% designCodebook('codeword', 'conventional');
clear;clc;
designCodebook('split', 'split');



function designCodebook(cbname, splitType)
    codebookDir = 'C:\Users\Yuxi\workspace\transcriptomic_analysis\codebook\out\';

    codebookPath = fullfile(codebookDir, ['codebook_mm_lib01_' cbname '.csv']);
    pd = probeDesign(['lib01_' cbname], 'mouse', codebookPath);

    set(pd, 'MERFISHAnalysisPath', 'C:\Users\Yuxi\workspace\probeOutput\lib01_');
    set(pd, 'basePath', 'C:\Users\Yuxi\workspace\genomeData');
    set(pd, 'fpkmPath', 'Mus_musculus_proxy.fpkm_tracking');
    set(pd, 'readoutPath', [codebookDir '\used_readouts_' cbname '.fasta']);

    set(pd, 'FPKMabundanceThreshold', 0);
    
    if strcmp(splitType, 'conventional')
        set(pd, 'probeSpacing', -20, 'tripleHeadedsmELT', true);
        set(pd, 'forbiddenSeqs', {'AAAA','TTTT','CCCC','GGGG'});
        set(pd, 'numProbesPerGene', 192);
        set(pd, 'regionGC', [0.43, 0.63], 'regionTm', [66,76]);
        set(pd, 'isoSpecificity', [0, 1], 'specificity', [0.75, 1]);

    elseif strcmp(splitType, 'split')
        set(pd, 'probeSpacing', 0, 'splitType', 'split');
        set(pd, 'regionLength', 25);
        set(pd, 'numProbesPerGene', 144);
        set(pd, 'regionGC', [0.2, 0.8], 'regionTm', [56, 72]);
        set(pd, 'isoSpecificity', [0, 1], 'specificity', [0.2, 1]);
    end

    set(pd, 'isPredesignedPrimer', true, 'primerID', [1 1]);

    pd.buildLibrary();

end