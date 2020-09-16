%% Fill into order template
genescript_template = 'C:\Users\Yuxi\workspace\MERFISH_analysis\OutputForLIMS\Genscript_Schnitzer_lab_v0.xlsx';
IDTtemplate = 'C:\Users\Yuxi\workspace\MERFISH_analysis\OutputForLIMS\IDT_plate_template.xlsx';
chipconfig_fn = 'C:\Users\Yuxi\workspace\transcriptomic_analysis\codebook\out\chip_summary_lib01.xlsx';
% designed oligos, input to this script
fasta_base_dir = 'C:\Users\Yuxi\workspace\genomeData';
genescript_out = 'C:\Users\Yuxi\workspace\MERFISH_analysis\OutputForLIMS\Genscript_Schnitzer_lab_lib01_v5.xlsx';

% read chip config
chipconfig = readtable(chipconfig_fn);
n_subpool = sum(~isnan(chipconfig.fwdPrimer));
n_oligos = zeros(n_subpool, 1);

% loop over subpools
all_oligos = cell(0,1);
for i = 1:n_subpool
    fprintf('\n\n\n======subpool %02d======\n', i);
    codebook = chipconfig.codebook{i};
    subpool_nm = chipconfig.label{i};
    fn = sprintf('C:\\Users\\Yuxi\\workspace\\genomeDatalib01_%s\\lib01_%s_oligos.fasta',...
        codebook, codebook);
    out = sprintf('G:\\Shared drives\\RNA_imaging\\probes\\lib01_v5\\genescript_forms\\Genscript_Schnitzer_lab_lib01_subpool%02d_v5.xlsx', i);
    fasta_out = sprintf('G:\\Shared drives\\RNA_imaging\\probes\\lib01_v5\\lib01_oligo_fasta\\lib01_subpool%02d_%s_v5.fasta', i, subpool_nm);
    disp(fn);
    disp(out)
    primer_inds = [chipconfig.fwdPrimer(i), chipconfig.revPrimer(i)];
    [n_oligos(i), oligos_compact] = genescript_form_new_primer(fn, out, fasta_out, genescript_template, subpool_nm, [i i]);
    all_oligos = [all_oligos; oligos_compact];
end

% write to merged genescript form
outall = sprintf('G:\\Shared drives\\RNA_imaging\\probes\\lib01_v5\\genescript_forms\\Genscript_Schnitzer_lab_lib01_v5.xlsx');
status = copyfile(genescript_template, outall);
xlswrite(outall, all_oligos, 'Oligo pool sequence form', 'A3');
xlswrite(outall, length(all_oligos), 'Quatition request form', 'I14');
fprintf('Finished writing to final Genescript form:\n%s', outall);

%% all primers
primer_out = 'G:\Shared drives\RNA_imaging\probes\lib01_v5\primers\lib01PCRprimers.fasta';
primer_out_IDT = 'G:\Shared drives\RNA_imaging\probes\lib01_v5\primers\lib01PCRprimers.xlsx';
fwdPrimerPath = 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\FwdPrimers.fasta';
fwd = fastaread(fwdPrimerPath);
revPrimerT7Path = 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\RevPrimersT7.fasta';
revT7 = fastaread(revPrimerT7Path);
used_primers = reshape([fwd(1:n_subpool), revT7(1:n_subpool)]', 2*n_subpool, 1);
% fasta
if isfile(primer_out)
   delete(primer_out); 
end
fastawrite(primer_out, used_primers);
%xlsx IDT
used_primers_xlsx = struct2cell(reshape([fwd(1:n_subpool), revT7(1:n_subpool)], 2*n_subpool, 1));
status = copyfile(IDTtemplate, primer_out_IDT);
xlswrite(primer_out_IDT, used_primers_xlsx', 1, 'B2');

%% digestion primers
split_inds = [15,16,17];
cut_primer_out = 'G:\Shared drives\RNA_imaging\probes\lib01_v5\primers\lib01CutPrimers.fasta';

rev = fastaread('C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\RevPrimers.fasta');

cutPrimers = [];
for ind = split_inds
   fwdCutSeq = ['NN' ...
       'GGTACC' ... %  KpnI
       seqrcomplement(fwd(ind).Sequence(end-13:end))];
   fwdCutHeader = ['Cut' fwd(ind).Header];
   fwdCut = struct('Header', fwdCutHeader, 'Sequence', fwdCutSeq);
   
   revCutSeq = [rev(ind).Sequence(end-13:end) ...
       'GAATTC' ... % EcoRI
       'NN'];
   revCutHeader = ['Cut' rev(ind).Header];
   revCut = struct('Header', revCutHeader, 'Sequence', revCutSeq);
   
   cutPrimers = [cutPrimers; fwdCut; revCut];
end

fastawrite(cut_primer_out, cutPrimers);

%% Bridge sequences
bridge_out = 'G:\Shared drives\RNA_imaging\probes\lib01_v5\readouts\lib01_split_cortex_bridges.fasta';

[~,bridge_template,~] = xlsread('G:\Shared drives\RNA_imaging\probes\lib01_v5\readouts\splitPaperBridgeSeq.xlsx');
bridge_template = bridge_template(3:end,:);

adt = {};
adt{1} = 'TGGGACGGTTCCAATCGGATC';
adt{2} = 'ACCTCCGTTAGACCCGTCAG';
adt{3} = 'CTCACCTGCACCTCCAACCG';
adt{4} = 'CCACCCATTCTGGGAGTACG';

bridge_inds = [17,18,19];
bridge_fasta = [];
for i = 1:3
    adt_ind = i;
    bridge_parts = strsplit(bridge_template{bridge_inds(i), 2});
    bridge_core = bridge_parts{2};
    bridge_seq = [seqrcomplement(adt{adt_ind}) ' '...
        bridge_core ' ' ...
        seqrcomplement(adt{adt_ind})];
    if adt_ind == 1
        bridge_seq = bridge_seq(1:end-1);% to fit into 60 nt
    end
    bridge_header = [bridge_template{bridge_inds(i), 1} '_adpt' num2str(adt_ind)];
    bridge_fasta = [bridge_fasta; struct('Header', bridge_header, 'Sequence', bridge_seq)];
end

fastawrite(bridge_out, bridge_fasta);

%% update and save chipconfig
chipconfig.numActualOligo(1:n_subpool) = n_oligos;
chipconfig.numActualOligo(end) = sum(n_oligos);
writetable(chipconfig, chipconfig_fn, 'sheet', 'chipconfig');
chipconfig_copy = 'G:\Shared drives\RNA_imaging\probes\lib01_v5\genescript_forms\chipconfig.xlsx';
copyfile(chipconfig_fn, chipconfig_copy);
writetable(chipconfig, chipconfig_copy, 'sheet', 'chipconfig');

fprintf('\n\nCompleted!\n');

function [n_oligos, oligos_compact] = genescript_form_new_primer(fn, out, fasta_out, template, subpool_nm, primer_id)
    % input: 
    %   fn - source, generated oligo fasta file
    %   out - generated filled-out order forms
    %   fasta_out - generated final oligos in fasta format
    %   template - order form template from genescript,
    %   subpool_nm: str. 
    %   primer_id: e.g. [1 1], fwd and rev
    % returns:
    %   n_oligos - int
    %   oligos_compact - cell array, oligo sequences

    % read primer seq -> two fasta structs
    fwdPrimerPath = 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\FwdPrimers.fasta';
    revPrimerPath = 'C:\Users\Yuxi\workspace\MERFISH_analysis\codebook_construction\RevPrimers.fasta';
    fwd = fastaread(fwdPrimerPath);
    rev = fastaread(revPrimerPath);
    fwd = fwd(primer_id(1));
    rev = rev(primer_id(2));

    % readout path
    readoutPath = 'C:\Users\Yuxi\workspace\transcriptomic_analysis\codebook\ref\AllReadouts.fasta';
    readouts = fastaread(readoutPath);
    extra_readouts = readouts(65:67);
    extra_genes = {'Cux2', 'Grm2', 'Foxp2'};
    
    % read designed oligos
    oligos = fastaread(fn);
    oligos = struct2cell(oligos);
    
    % preprocess header information
    headers = oligos(1,:)';
    header_parts = cellfun(@(x) strsplit(x), headers, 'UniformOutput', false);
    % `header_subpool` - cell array of subpool names
    header_subpool = cellfun(@(x) strsplit(x{1}, '__'), header_parts, 'UniformOutput', false);
    header_subpool = cellfun(@(x) x{2}, header_subpool, 'UniformOutput', false);
    
    % filter oligos that belong to the subpool by inpool_inds & add primer info to headers
    if endsWith(subpool_nm, '_25mer')
        % handle 25mer control subpool which uses the same TR as another subpool
        ref_subpool_nm = strsplit(subpool_nm, '_25mer');
        ref_subpool_nm = ['split_' ref_subpool_nm{1}];
        inpool_inds = strcmp(header_subpool, ref_subpool_nm);
        % add primer info to headers
        header_final = cellfun(@(x) [['lib01_split_' subpool_nm] fwd.Header x(3:end-1) rev.Header], header_parts, 'UniformOutput', false);
    else
        % inpool_inds - array, whether an oligo is in the pool
        inpool_inds = strcmp(header_subpool, subpool_nm);
        % add primer info to headers
        header_final = cellfun(@(x) strjoin([x(1) fwd.Header x(3:end-1) rev.Header]), header_parts, 'UniformOutput', false);
    end
    % filter headers
    header_final = header_final(inpool_inds); 

    % process the sequence
    oligos = oligos(2,:)';
    % filter those that belong
    oligos = oligos(inpool_inds);
    oligo_parts = cellfun(@(x) strsplit(x), oligos, 'UniformOutput', false);

    if endsWith(subpool_nm, '_25mer')
        % deal with readout swap. Hardcoded.
        gene_fasta = cell(0,2);
        for idg = 1:3
            gene = extra_genes{idg};
            gene_mask = cellfun(@(x) contains(strjoin(x), gene), header_final);
            gene_headers = header_final(gene_mask);
            gene_oligos = oligo_parts(gene_mask);
            % function, arg: cell array - header parts, returns: int - index of TRegion wiithin an oligo
            getTR = @(parts) find(cellfun(@(part) contains(part, gene), parts), 1); 
            tr_inds = cellfun(getTR, gene_headers);
            rdout = extra_readouts(idg);
            rdout.Sequence = seqrcomplement(rdout.Sequence);
            % loop over each oligo
            for ido = 1:length(gene_headers)
                rndint = randi([0 1]);
                if rndint == 1
                    gene_headers{ido} = strjoin([gene_headers{ido}(1:2) rdout.Header rdout.Header gene_headers{ido}(tr_inds(ido)) rdout.Header gene_headers{ido}(end)]);
                    gene_oligos{ido} = strjoin([fwd.Sequence rdout.Sequence rdout.Sequence 'A' gene_oligos{ido}(tr_inds(ido)) 'A' rdout.Sequence rev.Sequence]);
                elseif rndint == 0
                    gene_headers{ido} = strjoin([gene_headers{ido}(1:2) rdout.Header gene_headers{ido}(tr_inds(ido)) rdout.Header rdout.Header gene_headers{ido}(end)]);
                    gene_oligos{ido} = strjoin([fwd.Sequence rdout.Sequence 'A' gene_oligos{ido}(tr_inds(ido)) 'A' rdout.Sequence rdout.Sequence rev.Sequence]);
                end
            end % end loop over each oligo
            gene_fasta = [gene_fasta; [gene_headers, gene_oligos]];
        end % end loop over each gene
        oligos = gene_fasta(:,2);
        oligo_fasta = cell2struct(gene_fasta', {'Header','Sequence'});
    else
        % for those with original readouts, just concat primers
        oligo_parts = cellfun(@(x) [fwd.Sequence, x(2:end-1), seqrcomplement(rev.Sequence)], oligo_parts, 'UniformOutput', false);
        oligos = cellfun(@(x) strjoin(x), oligo_parts, 'UniformOutput', false);
        oligo_fasta = cell2struct([header_final'; oligos'], {'Header','Sequence'});
    end

    % squeeze whitespace
    oligos_compact = cellfun(@(x) x(~isspace(x)), oligos, 'UniformOutput', false);

    % write `oligos` to Genescript xls form
%     status = copyfile(template, out);
%     xlswrite(out, oligos_compact, 'Oligo pool sequence form', 'A3');
%     xlswrite(out, length(oligos_compact), 'Quatition request form', 'I14');

    %write `header_final` and `oligos` to fasta
    fprintf('Writing to fasta %s\n', fasta_out);
    
    if isfile(fasta_out)
       delete(fasta_out); 
    end
    fastawrite(fasta_out, oligo_fasta);

    fprintf('#oligos: %d\n', length(oligos_compact));
    n_oligos = length(oligos_compact);

end
