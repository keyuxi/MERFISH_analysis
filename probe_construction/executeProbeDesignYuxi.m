clear;clc;
addpath(genpath('C:\Users\Yuxi\workspace\MERFISH_analysis'))

%% lib01_merfish

% designCodebook('mer', 'conventional');
% designCodebook('trihead', 'conventional');
%designCodebook('opioid', 'conventional');
%clear;clc;
designCodebook('split', 'split');

%{
%% Fill into order template
n_oligos = zeros(7,1);
oligo_file_nm = {'merfish', 'merfish', 'EI', 'valence', 'seq', 'seq', 'seq'};
gene_ranges = {...
    {'Adarb2', 'Upk1b'},...
    {'Camk2a', 'Vip'},...
    {'Slc17a7', 'Gad1'},...
    {'Rspo2', 'Ppp1r1b'},...
    {'Cnr1', 'Stmn1'},...
    {'Fos', 'Bdnf'},...
    {'Oprk1', 'Oprm1'}
    };
template = 'C:\Users\Yuxi\workspace\MERFISH_analysis\OutputForLIMS\Genscript_Schnitzer_lab_v0.xlsx';

fn = 'C:\Users\Yuxi\workspace\genomeDatalib01_seq_0306\lib01_seq_0306_oligos.fasta';
out = 'C:\Users\Yuxi\workspace\MERFISH_analysis\OutputForLIMS\Genscript_Schnitzer_lab_lib01_test_v4.xlsx';

% loop over subpools
for i = 1:7
    fprintf('\n\n\n======subpool%02d======\n', i);
    fn = sprintf('C:\\Users\\Yuxi\\workspace\\genomeDatalib01_%s_0306\\lib01_%s_0306_oligos.fasta',...
        oligo_file_nm{i}, oligo_file_nm{i});
    out = sprintf('C:\\Users\\Yuxi\\workspace\\MERFISH_analysis\\OutputForLIMS\\lib01_v4\\Genscript_Schnitzer_lab_lib01_subpool%02d_v4.xlsx', i);
    disp(fn);
    disp(out)
    n_oligos(i) = genescript_form_new_primer(fn, out, template, gene_ranges{i}, [i i]);
end

%QC
fn = sprintf('C:\\Users\\Yuxi\\workspace\\genomeDatalib01_%s_0306\\lib01_%s_0306_oligos.fasta',...
    'EI', 'EI');
merfish_oligos = fastaread(fn);
out = sprintf('C:\\Users\\Yuxi\\workspace\\MERFISH_analysis\\OutputForLIMS\\lib01_v4\\Genscript_Schnitzer_lab_lib01_subpool%02d_v4.xlsx', 8);

fwdPrimerPath = 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\FwdPrimers.fasta';
revPrimerPath = 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\RevPrimers.fasta';
fwd = fastaread(fwdPrimerPath);
rev = fastaread(revPrimerPath);

qc = struct2cell(merfish_oligos([1, 48, 96]));
qc = qc(2,:)';
qc_parts = cellfun(@(x) strsplit(x), qc, 'UniformOutput', false);
for i = 1:3
    qc_parts{i} = [fwd(i+7).Sequence, qc_parts{i}(2:end-1), seqrcomplement(rev(i+7).Sequence)];
end
qc = cellfun(@(x) strjoin(x), qc_parts, 'UniformOutput', false);
qc = reshape(repmat(qc',3,1),9,1);
qc = cellfun(@(x) x(~isspace(x)), qc, 'UniformOutput', false);
status = copyfile(template, out);
xlswrite(out, qc, 'Oligo pool sequence form', 'A3');
xlswrite(out, length(qc), 'Quatition request form', 'I14');

% all primers
primer_out = 'C:\Users\Yuxi\workspace\MERFISH_analysis\OutputForLIMS\lib01_v4\lib01_primers.fasta';
revPrimerT7Path = 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\RevPrimersT7.fasta';
revT7 = fastaread(revPrimerT7Path);
used_primers = reshape([fwd(1:10), revT7(1:10)]', 20, 1);
fastawrite(primer_out, used_primers);
%}

function n_oligos = genescript_form_new_primer(fn, out, template, gene_range, primer_id)
% input: 
% generated oligo fasta file, generated filled-out order forms, 
% order form template from genescript,
% gene_range: {first gene name, last gene name}. 
% primer_id: e.g. [1 1], fwd and rev

% read primer seq -> two fasta structs
fwdPrimerPath = 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\FwdPrimers.fasta';
revPrimerPath = 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\RevPrimers.fasta';
fwd = fastaread(fwdPrimerPath);
rev = fastaread(revPrimerPath);
fwd = fwd(primer_id(1));
rev = rev(primer_id(2));

oligos = fastaread(fn);
oligos = struct2cell(oligos);

% find the gene range
oligo_nm = oligos(1,:)';
oligo_ind_range = zeros(2,1);
% find the start
for i = 1:length(oligo_nm)
   is_found = strfind(oligo_nm{i}, gene_range{1});
   if ~isempty(is_found)
      oligo_ind_range(1) = i;
      break;
   end
end
% find the end
flag = 0;
for i = 1:length(oligo_nm)
   is_found = strfind(oligo_nm{i}, gene_range{2});
   if ~isempty(is_found) && flag==0
       flag = 1;
   elseif isempty(is_found) && flag==1    
      oligo_ind_range(2) = i-1;
      break;
   end
end
if flag==1 && oligo_ind_range(2)==0
   oligo_ind_range(2) = length(oligo_nm); 
end


disp(oligo_ind_range);
disp(oligo_nm{oligo_ind_range(1)});
disp(oligo_nm{oligo_ind_range(2)});

% process the sequence
oligos = oligos(2,:)';
oligos = oligos(oligo_ind_range(1):oligo_ind_range(2));
oligo_parts = cellfun(@(x) strsplit(x), oligos, 'UniformOutput', false);
oligo_parts = cellfun(@(x) [fwd.Sequence, x(2:end-1), seqrcomplement(rev.Sequence)], oligo_parts, 'UniformOutput', false);
oligos = cellfun(@(x) strjoin(x), oligo_parts, 'UniformOutput', false);
oligos = cellfun(@(x) x(~isspace(x)), oligos, 'UniformOutput', false);

% write `oligos` to xls form
status = copyfile(template, out);
xlswrite(out, oligos, 'Oligo pool sequence form', 'A3');
xlswrite(out, length(oligos), 'Quatition request form', 'I14');
fprintf('#oligos: %d\n', length(oligos));
n_oligos = length(oligos);
end

function genescript_form(fn, out, template)
% generated oligo fasta file, generated filled-out order forms, order form
% template from genescript
oligos = fastaread(fn);
oligos = struct2cell(oligos);
oligos = oligos(2,:)';
%oligo_parts = cellfun(@(x) strsplit(x), oligos, 'UniformOutput', false);
%oligo_parts = cellfun(@(x) x(2:end-1), oligo_parts, 'UniformOutput', false);
%oligos = cellfun(@(x) strjoin(x), oligo_parts, 'UniformOutput', false);
oligos = cellfun(@(x) x(~isspace(x)), oligos, 'UniformOutput', false);
n_oligos = length(oligos);

status = copyfile(template, out);
xlswrite(out, oligos, 'Oligo pool sequence form', 'A3');
xlswrite(out, length(oligos), 'Quatition request form', 'I14');
fprintf('#oligos: %d\n', n_oligos);

end

function genescript_form_primerless(fn, out, template)
% generated oligo fasta file, generated filled-out order forms, order form
% template from genescript
oligos = fastaread(fn);
oligos = struct2cell(oligos);
oligos = oligos(2,:)';
oligo_parts = cellfun(@(x) strsplit(x), oligos, 'UniformOutput', false);
oligo_parts = cellfun(@(x) x(2:end-1), oligo_parts, 'UniformOutput', false);
oligos = cellfun(@(x) strjoin(x), oligo_parts, 'UniformOutput', false);
oligos = cellfun(@(x) x(~isspace(x)), oligos, 'UniformOutput', false);

status = copyfile(template, out);
xlswrite(out, oligos, 'Oligo pool sequence form', 'A3');
xlswrite(out, length(oligos), 'Quatition request form', 'I14');
fprintf('#oligos: %d\n', length(oligos));

end

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
        set(pd, 'numProbesPerGene', 96);
        set(pd, 'regionGC', [0.43, 0.63], 'regionTm', [66,76]);
        set(pd, 'isoSpecificity', [0, 1], 'specificity', [0.75, 1]);

    elseif strcmp(splitType, 'split')
        set(pd, 'probeSpacing', 0, 'splitType', 'split');
        set(pd, 'regionLength', 25);
        set(pd, 'numProbesPerGene', 144);
        set(pd, 'regionGC', [0.2, 0.8], 'regionTm', [66,76]);
        set(pd, 'isoSpecificity', [0, 1], 'specificity', [0.2, 1]);
    end

    %set(pd, 'isPredesignedPrimer', true, 'primerID', [1 1]);

    pd.buildLibrary();

end