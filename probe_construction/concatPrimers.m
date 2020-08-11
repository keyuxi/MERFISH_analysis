%% Fill into order template
genescript_template = 'C:\Users\Yuxi\workspace\MERFISH_analysis\OutputForLIMS\Genscript_Schnitzer_lab_v0.xlsx';
chipconfig_fn = 'C:\Users\Yuxi\workspace\transcriptomic_analysis\codebook\out\chip_sumary_lib01.xlsx';
% designed oligos, input to this script
fasta_base_dir = 'C:\Users\Yuxi\workspace\genomeData';
genescript_out = 'C:\Users\Yuxi\workspace\MERFISH_analysis\OutputForLIMS\Genscript_Schnitzer_lab_lib01_v5.xlsx';

% read chip config
chipconfig = readtable(chipconfig_fn);
n_subpool = sum(~isnan(chipconfig.fwdPrimer));
n_oligos = zeros(n_subpool, 1);

% loop over subpools
for i = 1:n_subpool
    fprintf('\n\n\n======subpool %02d======\n', i);
    codebook = chipconfig.codebook{i};
    subpool_nm = chipconfig.label{i};
    fn = sprintf('C:\\Users\\Yuxi\\workspace\\genomeDatalib01_%s\\lib01_%s_oligos.fasta',...
        codebook, codebook);
    out = sprintf('G:\\Shared drives\\RNA_imaging\\probes\\lib01_v5\\genescript_forms\\Genscript_Schnitzer_lab_lib01_subpool%02d_%s_v5.xlsx', i, subpool_nm);
    fasta_out = sprintf('G:\\Shared drives\\RNA_imaging\\probes\\lib01_v5\\lib01_oligo_fasta\\lib01_subpool%02d_%s_v5.fasta', i, subpool_nm);
    disp(fn);
    disp(out)
    primer_inds = [chipconfig.fwdPrimer(i), chipconfig.revPrimer(i)];
    n_oligos(i) = genescript_form_new_primer(fn, out, fasta_out, genescript_template, subpool_nm, [i i]);
end

% all primers
primer_out = 'G:\Shared drives\RNA_imaging\probes\lib01_v5\genescript_forms\lib01_primers.fasta';
fwdPrimerPath = 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\FwdPrimers.fasta';
fwd = fastaread(fwdPrimerPath);
revPrimerT7Path = 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\RevPrimersT7.fasta';
revT7 = fastaread(revPrimerT7Path);
used_primers = reshape([fwd(1:n_subpool), revT7(1:n_subpool)]', 2*n_subpool, 1);
if isfile(primer_out)
   delete(primer_out); 
end
fastawrite(primer_out, used_primers);

% update and save chipconfig
chipconfig.numActualOligo(1:n_subpool) = n_oligos;
chipconfig.numActualOligo(end) = sum(n_oligos);
writetable(chipconfig, chipconfig_fn, 'sheet', 'chipconfig');
chipconfig_copy = 'G:\Shared drives\RNA_imaging\probes\lib01_v5\genescript_forms\chipconfig.xlsx';
copyfile(chipconfig_fn, chipconfig_copy);
writetable(chipconfig, chipconfig_copy, 'sheet', 'chipconfig');

function n_oligos = genescript_form_new_primer(fn, out, fasta_out, template, subpool_nm, primer_id)
    % input: 
    % generated oligo fasta file, generated filled-out order forms, 
    % order form template from genescript,
    % subpool_nm: str. 
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

    % filter oligos that belong to the subpool
    headers = oligos(1,:)';
    header_parts = cellfun(@(x) strsplit(x), headers, 'UniformOutput', false);
    % `header_subpool` - cell array of subpool names
    header_subpool = cellfun(@(x) strsplit(x{1}, '__'), header_parts, 'UniformOutput', false);
    header_subpool = cellfun(@(x) x{2}, header_subpool, 'UniformOutput', false);
    % inpool_inds - array, whether an oligo is in the pool
    inpool_inds = strcmp(header_subpool, subpool_nm);
    % add primer info to headers
    header_final = cellfun(@(x) strjoin([x(1) fwd.Header x(3:end-1) rev.Header]), header_parts, 'UniformOutput', false);
    header_final = header_final(inpool_inds);

    % process the sequence
    oligos = oligos(2,:)';
    % filter those that belong
    oligos = oligos(inpool_inds);
    oligo_parts = cellfun(@(x) strsplit(x), oligos, 'UniformOutput', false);
    % concat primers
    oligo_parts = cellfun(@(x) [fwd.Sequence, x(2:end-1), seqrcomplement(rev.Sequence)], oligo_parts, 'UniformOutput', false);
    oligos = cellfun(@(x) strjoin(x), oligo_parts, 'UniformOutput', false);
    % squeeze whitespace
    oligos_compact = cellfun(@(x) x(~isspace(x)), oligos, 'UniformOutput', false);

    % write `oligos` to xls form
    status = copyfile(template, out);
    xlswrite(out, oligos_compact, 'Oligo pool sequence form', 'A3');
    xlswrite(out, length(oligos_compact), 'Quatition request form', 'I14');

    %write `header_final` and `oligos` to fasta
    fprintf('Writing to fasta %s\n', fasta_out);
    oligo_fasta = cell2struct([header_final'; oligos'], {'Header','Sequence'});
    if isfile(fasta_out)
       delete(fasta_out); 
    end
    fastawrite(fasta_out, oligo_fasta);

    fprintf('#oligos: %d\n', length(oligos_compact));
    n_oligos = length(oligos_compact);

end
