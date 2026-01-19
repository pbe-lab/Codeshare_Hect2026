 %% setup
addpath(genpath('C:\Users\jasmi\OneDrive - University of Pittsburgh\gsp'))

gcs=loadSettings();
[path, seeg, nStim, vl_setup, gsp_setup, classLabels, classId] = gsp_setup_short({'vl','gsp'});
dataDir = 'C:\Users\jasmi\OneDrive - University of Pittsburgh\gsp\save\process_sEEG\';

% P23N010 no 30k audio
% P22N002 no consent
% P24N004 passive and partial
% P24N004 two sessions, only trials 1-288, montage change
% P24N016 no DAQ (presby)
% P25N001 no DAQ (presby)
nats_pts = listPtsExpts('naturalSounds');

which_data = 'erp';

vl_pts = listPtsExpts('voiceLocalizer');
vl_pts = vl_pts(~contains(vl_pts,'P22N002')&~contains(vl_pts,'P24N004')&~contains(vl_pts,'P25N013'));
opts = struct('mode',which_data, ...
    'fsNr',1000,'fsDs',1000, ...
    'avgRepeats',false, ...
    'hp_cutoff',0.2,'hp_order',2, ...
    'erpFilterPhase','zero', ...
    'erpGDComp',false,'refMode','CAR', ...
    'baselineNorm',false,...
    'saveDirRoot','','dataDirRoot','C:\Users\jasmi\OneDrive - University of Pittsburgh\gsp\data\');
for iPt = 1:numel(vl_pts)
    if ~isfile(fullfile(dataDir, which_data, 'voiceLocalizer', sprintf('%s.mat',vl_pts{iPt})))
        process_sEEG(vl_pts{iPt},'voiceLocalizer', opts)
    end
end

gsp_pts = listPtsExpts('gsp');
opts = struct('mode',which_data, ...
    'fsNr',1000,'fsDs',1000, ...
    'avgRepeats',false, ...
    'hp_cutoff',0.2,'hp_order',2, ...
    'erpFilterPhase','zero', ...
    'erpGDComp',false,'refMode','CAR', ...
    'baselineNorm',false,...
    'saveDirRoot','','dataDirRoot','C:\Users\jasmi\OneDrive - University of Pittsburgh\gsp\data\');
for iPt = 1:numel(gsp_pts)
    if ~isfile(fullfile(dataDir, which_data, 'gsp', sprintf('%s.mat',gsp_pts{iPt})))
        process_sEEG(gsp_pts{iPt},'gsp', opts)
    end
end
gsp_pts = gsp_pts(~contains(gsp_pts,'P24N011'));
all_tasks = ismember(gsp_pts, nats_pts) & ismember(gsp_pts, vl_pts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load ratings
gorilla = importGorilla;
x = [400, 250, 144];  %gorilla data organized as tsp, gsp, vl **
idxOff = cumsum(x);
idxOn = [1, idxOff(1:2)];
idx_rate(1,:) = [idxOn(1,2)+1 idxOff(1,2)]; %gsp
idx_rate(2,:) = [idxOn(1,3)+1 idxOff(1,3)]; %vl

ratings_gsp_all = gorilla.response(idx_rate(1,1):idx_rate(1,2),:);
% exclude listeners with sd = 0
% use log10(voicelikeness)
drop_gorillas = std(ratings_gsp_all,'omitnan') <= 0.2;
sum(drop_gorillas)
ratings_gsp_all = ratings_gsp_all(:, ~drop_gorillas);
ratings_vl_all = gorilla.response(idx_rate(2,1):idx_rate(2,2),:);

ratings_gsp = mean(ratings_gsp_all, 2,'omitnan');
ratings_vl = mean(gorilla.response(idx_rate(2,1) : idx_rate(2,2),:), 2, 'omitnan');
[~, sort_idx_vl_grp] = sort(ratings_vl,'ascend');
[~, sort_idx_gsp_grp] = sort(ratings_gsp,'ascend');

ratings_gsp_sort = ratings_gsp(sort_idx_gsp_grp);
ratings_vl_sort = ratings_vl(sort_idx_vl_grp);

clr_map_roi = brighten(cbrewer2('Set3', 9),-0.5);
clr_map_roi(2,:) = clr_map_roi(end,:);
clr_map_task =  [0.00,0.45,0.74;0.85,0.33,0.10];


% % define quartiles for gsps
cutoffs = [0 0.25 0.5 0.75 1];
idx_gsp_grp = nan(size(ratings_gsp_sort));
idx_gsp_grp(ratings_gsp_sort < quantile(ratings_gsp_sort, cutoffs(2))) = 1; % least voice-like
idx_gsp_grp(ratings_gsp_sort >= quantile(ratings_gsp_sort, cutoffs(2))...
    & ratings_gsp_sort < quantile(ratings_gsp_sort, cutoffs(3))) = 2;
idx_gsp_grp(ratings_gsp_sort >= quantile(ratings_gsp_sort, cutoffs(3)) ...
    & ratings_gsp_sort < quantile(ratings_gsp_sort, cutoffs(4))) = 3;
idx_gsp_grp(ratings_gsp_sort >= quantile(ratings_gsp_sort, cutoffs(4))) = 4; % most voice-like
plot_voice_gsp = idx_gsp_grp;

% % define quartiles for vl
idx_vl_grp = nan(size(ratings_vl_sort));
idx_vl_grp(ratings_vl_sort < quantile(ratings_vl_sort, cutoffs(2))) = 1; % least voice-like
idx_vl_grp(ratings_vl_sort >= quantile(ratings_vl_sort, cutoffs(2)) ...
    & ratings_vl_sort < quantile(ratings_vl_sort, cutoffs(3))) = 2;
idx_vl_grp(ratings_vl_sort >= quantile(ratings_vl_sort, cutoffs(3)) ...
    & ratings_vl_sort < quantile(ratings_vl_sort, cutoffs(4))) = 3;
idx_vl_grp(ratings_vl_sort >= quantile(ratings_vl_sort, cutoffs(4))) = 4; % most voice-like
plot_voice_vl = idx_vl_grp;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get yamnet
rmpath('C:\Users\jasmi\OneDrive - University of Pittsburgh\gsp\code\external\Gaussian_Sound_Code_for_Distribution')
%conflicting hann function

yamout_gsp = runYamClass({'gsp'});
yamout_voice_gsp = yamout_gsp.topCat{1}(1,:)';
yamout_voice_gsp_sort = zscore(yamout_voice_gsp(sort_idx_gsp_grp));
yamout_animal_gsp = yamout_gsp.topCat{1}(2,:)';
yamout_animal_gsp_sort = zscore(yamout_animal_gsp(sort_idx_gsp_grp));

yamout_vl = runYamClass({'vl'});
yamout_voice_vl = yamout_vl.topCat{1}(1,:)';
yamout_voice_vl_sort = zscore(yamout_voice_vl(sort_idx_vl_grp));

yamout_animal_vl = yamout_vl.topCat{1}(2,:)';
yamout_animal_vl_sort = zscore(yamout_animal_vl(sort_idx_vl_grp));

[max_yam_out_vl, max_yam_id_vl] = max(yamout_vl.topCat{1}, [], 1);
max_yam_out_vl_sort = max_yam_out_vl(sort_idx_vl_grp)';
max_yam_id_vl_sort = max_yam_id_vl(sort_idx_vl_grp)';
yamout_vl_id_sort = yamout_vl.topCat{1}(:, sort_idx_vl_grp)';
[max_yam_out_gsp, max_yam_id_gsp] = max(yamout_gsp.topCat{1}, [], 1);
max_yam_out_gsp_sort = max_yam_out_gsp(sort_idx_gsp_grp)';
max_yam_id_gsp_sort = max_yam_id_gsp(sort_idx_gsp_grp)';
yamout_gsp_id_sort = yamout_gsp.topCat{1}(:, sort_idx_gsp_grp)';

cutoffs = [0 0.25 0.5 0.75 1];
% % define quartiles for gsps
idx_gsp_grp_yam = nan(size(yamout_voice_gsp));
idx_gsp_grp_yam(yamout_voice_gsp < quantile(yamout_voice_gsp, cutoffs(2))) = 1; % least voice-like
idx_gsp_grp_yam(yamout_voice_gsp >= quantile(yamout_voice_gsp, cutoffs(2))...
    & yamout_voice_gsp < quantile(yamout_voice_gsp, cutoffs(3))) = 2;
idx_gsp_grp_yam(yamout_voice_gsp >= quantile(yamout_voice_gsp, cutoffs(3)) ...
    & yamout_voice_gsp < quantile(yamout_voice_gsp, cutoffs(4))) = 3;
idx_gsp_grp_yam(yamout_voice_gsp >= quantile(yamout_voice_gsp, cutoffs(4))) = 4; % most voice-like
plot_voice_gsp_yam = idx_gsp_grp_yam(sort_idx_gsp_grp);

% % define quartiles for vl
idx_vl_grp_yam = nan(size(yamout_voice_vl));
idx_vl_grp_yam(yamout_voice_vl < quantile(yamout_voice_vl, cutoffs(2))) = 1; % least voice-like
idx_vl_grp_yam(yamout_voice_vl >= quantile(yamout_voice_vl, cutoffs(2)) ...
    & yamout_voice_vl < quantile(yamout_voice_vl, cutoffs(3))) = 2;
idx_vl_grp_yam(yamout_voice_vl >= quantile(yamout_voice_vl, cutoffs(3)) ...
    & yamout_voice_vl < quantile(yamout_voice_vl, cutoffs(4))) = 3;
idx_vl_grp_yam(yamout_voice_vl >= quantile(yamout_voice_vl, cutoffs(4))) = 4; % most voice-like
plot_voice_vl_yam = idx_vl_grp_yam(sort_idx_vl_grp);

% make concatenated stuff
cat_rate = cat(1, ratings_vl_sort, ratings_gsp_sort);
cat_yam = cat(1, yamout_voice_vl_sort, yamout_voice_gsp_sort);
cat_clr = cat(1, plot_voice_vl, plot_voice_gsp);
cat_clr_yam = cat(1, plot_voice_vl_yam, plot_voice_gsp_yam);
is_nat = logical(zeros(numel(cat_clr), 1));
is_nat(1:size(plot_voice_vl, 1), 1) = true;
cat_yamout_id = cat(1, yamout_vl_id_sort, yamout_gsp_id_sort);
bin_cat = ones(size(cat_clr, 1), 1);
bin_cat(cat_clr == 3 | cat_clr ==4) = 2;
bin_cat_yam = ones(size(cat_clr_yam, 1), 1);
bin_cat_yam(cat_clr_yam == 3 | cat_clr_yam ==4) = 2;

% make color maps
clr_map_enc = brighten(cmocean('balance',numel(unique(plot_voice_gsp)) + 4), 0.2);
clr_map_enc = clr_map_enc(2:end - 1, :);
clr_map_nr = brighten(cbrewer2('Set3', numel(unique(plot_voice_gsp)) + 2), -0.6);
clr_map = cmocean('rain', numel(unique(plot_voice_gsp)) + 2);
clr_map = clr_map(2:end, :);
clr_map_roi = cbrewer2('Set1', 8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------------get stimulus acoustic features-------------------------
if isfile(fullfile(gcs.fOneDrive,'../gsp/stimuli/stimAcFts.mat','vl','gsp'))
    load(fullfile(gcs.fOneDrive,'../gsp/stimuli/stimAcFts.mat'),'vl','gsp')
else
    [vl, gsp] = getStimAcFts();
    save(fullfile(gcs.fOneDrive,'../gsp/stimuli/stimAcFts.mat'),'vl','gsp')
end
% sort on listener ratings
gsp.gen_params = gsp.gen_params (sort_idx_gsp_grp, :);
gsp.tcoeff = gsp.tcoeff(sort_idx_gsp_grp, :);
gsp.fcoeff = gsp.fcoeff(sort_idx_gsp_grp, :);
gsp.acfts = gsp.acfts(:, :, sort_idx_gsp_grp);
gsp.t = 0:10:1500-50;
vl.tcoeff = vl.tcoeff(sort_idx_vl_grp, :);
vl.fcoeff = vl.fcoeff(sort_idx_vl_grp, :);
vl.acfts = vl.acfts(:, :, sort_idx_vl_grp);
vl.t = 0:10:550-50;
kp_stim_vl = bin_cat(is_nat) > 0;
kp_stim_gsp = bin_cat(~is_nat) > 0;

corr_acfts_vl = nan(numel(vl.t), 88);
corr_acfts_vl_pval = nan(numel(vl.t), 88);
corr_acfts_gsp = nan(numel(gsp.t), 88);
corr_acfts_gsp_pval = nan(numel(gsp.t), 88);

for iFt = 1:88
    [corr_acfts_vl(:, iFt), corr_acfts_vl_pval(:, iFt)] = corr(permute(vl.acfts(:, iFt, kp_stim_vl), [3,1,2]), ...
        ratings_vl_sort(kp_stim_vl), 'rows','pairwise','type','spearman');

    [corr_acfts_gsp(:,iFt), corr_acfts_gsp_pval(:,iFt)] = corr(permute(gsp.acfts(:,iFt,kp_stim_gsp), [3,1,2]), ...
        ratings_gsp_sort(kp_stim_gsp), 'rows','pairwise','type','spearman');
end

% at least 8 windows of significance, less than 10% missing time points
kp_fts = find(sum(corr_acfts_vl_pval < 0.01, 1)' >= 10 &...
    (sum(isnan(corr_acfts_vl_pval),1)/size(corr_acfts_vl_pval,1) < 0.05)' &...
    ~contains(table2cell(gsp.acft_names(:,2)),'voiced','IgnoreCase',true));

acfts_gsp = gsp.acfts;
acfts_vl = vl.acfts;
% 
% makeBar(cat(1, median(corr_acfts_vl(:, kp_fts),1,'omitnan'), ...
%     median(corr_acfts_gsp(:, kp_fts),1,'omitnan'))', ...
%     table2cell(gsp.acft_names(kp_fts, 2)), [],...
%     {'natural','synthetic'}, 'rho', 'Ratings (M)')
% ylim([-0.5 0.5]); set(gca, 'fontsize',12)
% grid on
% 
% steps = linspace(min(ratings_vl_sort),max(ratings_vl_sort),300);
% for i=1:numel(kp_fts)
%     figure('color','w','Position',[100+(240*i),590,240,180]); hold on
%     for iQ = 1:numel(steps)-1
%         scatter(mean(ratings_vl_sort(ratings_vl_sort>=steps(iQ) & ratings_vl_sort<=steps(iQ+1))), ...
%             mean(mean(acfts_vl(1:26, kp_fts(i), ratings_vl_sort>=steps(iQ) & ratings_vl_sort<=steps(iQ+1)),3,'omitnan'),1,'omitnan'),...
%             'filled','SizeData', 20,'markerfacealpha',1,'markerfacecolor',clr_map_task(1,:))
%         scatter(mean(ratings_gsp_sort(ratings_gsp_sort>=steps(iQ) & ratings_gsp_sort<=steps(iQ+1))), ...
%             mean(mean(acfts_gsp(1:26, kp_fts(i), ratings_gsp_sort>=steps(iQ) & ratings_gsp_sort<=steps(iQ+1)),3,'omitnan'),1,'omitnan'),...
%             'filled','sizedata', 20,'markerfacealpha',1,'markerfacecolor',clr_map_task(2,:))
%     end
%     set(gca, 'fontsize',12);xlabel('Rating');xticks([-1 1.8]);xlim([-1,1.8]);xticklabels({'NV','V'});yticks([])
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------curate data & save----------------------------
% which data to use?
which_data = 'erp';
which_filter = 'none'; % none or gsp
which_tasks = 'vl'; % all vl nats or gsp
nBs = 1000;

if ~isfile(fullfile(gcs.fOneDrive, '..\gsp\save\', sprintf('dataIn_%s_%s_%s_process_sEEG.mat', ...
        which_data, which_filter, which_tasks)))
    if strcmp(which_filter, 'none')
        if strcmp(which_tasks,'gsp')
            pts_rsa = gsp_pts;
        elseif strcmp(which_tasks, 'vl')
            pts_rsa = vl_pts;
        elseif strcmp(which_tasks, 'nats')
            pts_rsa = nats_pts;
        end
    else
        if strcmp(which_task,'all')
            pts_rsa = gsp_pts(all_tasks);
        elseif strcmp(which_task, 'vl')
            pts_rsa = gsp_pts(ismember(gsp_pts, vl_pts));
        end
    end

    pts_rsa = pts_rsa(~contains(pts_rsa,'P21N007') & ~contains(pts_rsa,'P23N010') & ~contains(pts_rsa,'P24N016'));
    clear hga_matrix  nr_matrix pt_matrix chan_roi_matrix
    hga_matrix = [];
    nr_matrix = [];
    chan_roi_matrix = {};
    pt_matrix = {};

    for iRsa = 1:length(pts_rsa)
        pt = pts_rsa{iRsa};
        fprintf('%s\n',pt)
        if strcmp(which_filter, 'gsp') || strcmp(which_tasks, 'gsp')
            ld_gsp = load(fullfile(path.save,'process_sEEG',which_data,'gsp',...
                sprintf('%s.mat',pt)));
            % drop 1-back repeats
            is_repeat = ld_gsp.info.isRepeat;
            kpt = (ld_gsp.t >= -500 & ld_gsp.t <= 2000);
            if isfield(ld_gsp,'erp')
                [~, mn, sd] = zscore(ld_gsp.erp(ld_gsp.t < -100, :, ~is_repeat),0,1);
                gsp_nr = (ld_gsp.erp(kpt, :, ~is_repeat) - mn)./sd;
                gsp_nr_noRpts = nan(size(gsp_nr,1),size(gsp_nr,2), 144);
            else
                [~, mn, sd] = zscore(ld_gsp.hga(ld_gsp.t < -100, :, ~is_repeat),0,1);
                gsp_nr = (ld_gsp.hga(kpt, :, ~is_repeat) - mn)./sd;
                gsp_nr_noRpts = nan(size(gsp_nr,1),size(gsp_nr,2), 144);
            end 
            stim_ids = unique(ld_gsp.info.trialStimIDs(~is_repeat));
            for i = 1:numel(stim_ids)
                iStim = stim_ids(i);
                gsp_nr_noRpts(:,:,iStim) = mean(gsp_nr(:,:,ld_gsp.info.trialStimIDs(~is_repeat)==iStim),3,'omitnan');
            end
        elseif  strcmp(which_tasks,'vl') || strcmp(which_tasks,'all')
            ld_vl = load(fullfile(path.save,'process_sEEG',which_data,'voiceLocalizer',sprintf('%s.mat',pt)));
            
            if ~strcmpi(pt, 'P22N006')
                is_repeat = ld_vl.info.isRepeat;
            end
            kpt = (ld_vl.t >= -500 & ld_vl.t <= 1000);
            if isfield(ld_vl,'erp')
                [~, mn, sd] = zscore(ld_vl.erp(ld_vl.t < -100, :, ~is_repeat),0,1);
                vl_nr = (ld_vl.erp(kpt, :, ~is_repeat) - mn)./sd;
                vl_nr_noRpts = nan(size(vl_nr,1),size(vl_nr,2), 144);
            else
                [~, mn, sd] = zscore(ld_vl.hga(ld_vl.t < -100, :, ~is_repeat),0,1);
                vl_nr = (ld_vl.hga(kpt, :, ~is_repeat) - mn)./sd;
                vl_nr_noRpts = nan(size(vl_nr,1),size(vl_nr,2), 144);
            end
            stim_ids = unique(ld_vl.info.trialStimIDs(~is_repeat));
            for i = 1:numel(stim_ids)
                iStim = stim_ids(i);
                vl_nr_noRpts(:,:, iStim) = mean(vl_nr(:,:,ld_vl.info.trialStimIDs(~is_repeat)==iStim),3,'omitnan');
            end
        end

        if strcmp(which_tasks,'all')
            % reconcile remaining differences in montages - use gsp montage to find
            % mutually inclusive chans in the vl and nats
            ld_nats = load(fullfile(path.save,'process_sEEG',which_data,'naturalSounds',sprintf('%s.mat',pt)));
            %             kp_stim = (natS.cat==2 | natS.cat==1); % exclude singing
            is_repeat = ld_nats.info.isRepeat;
            kpt = (ld_nats.t >= -500 & ld_nats.t <= 2000);
            [~,mn,sd] = zscore(ld_nats.erp(ld_nats.t<0, :, ~is_repeat),0,1);
            nats_erp = (ld_nats.erp(kpt, :, ~is_repeat) - mn)./sd;
            if ~all(ismember(ld_gsp.info.chanNames, ld_vl.info.chanNames)) ...
                    || ~all(ismember(ld_gsp.info.chanNames, ld_nats.info.chanNames)) ...
                    || size(gsp_nr,2) ~= size(vl_nr,2) || size(gsp_nr,2) ~= size(nat_erp,2)
                kp_gsp_chans = ismember(ld_gsp.info.chanNames,ld_vl.info.chanNames) ...
                    & ismember(ld_gsp.info.chanNames,ld_nats.info.chanNames);
                gsp_nr = gsp_nr(:,kp_gsp_chans,:);
                kp_chans = ismember(ld_nats.info.chanNames,ld_gsp.info.chanNames(kp_gsp_chans));
                nat_erp = nat_erp(:, kp_chans, :);
                kp_chans = ismember(ld_vl.info.chanNames,ld_gsp.info.chanNames(kp_gsp_chans));
                vl_nr = vl_nr(:, kp_chans, :);
            end
        end
        if strcmp(which_filter,'gsp') && strcmp(which_tasks,'all') && ...
                (~all(ismember(ld_gsp.info.chanNames, ...
                ld_vl.info.chanNames)) || size(gsp_nr,2) ~= size(vl_nr,2))
            kp_gsp_chans = ismember(ld_gsp.info.chanNames,ld_vl.info.chanNames);
            gsp_nr = gsp_nr(:, kp_gsp_chans, :);
            cat_syn_nr = gsp_nr_noRpts(:, kp_gsp_chans, sort_idx_gsp_grp);
            kp_vl_chans = ismember(ld_vl.info.chanNames,ld_gsp.info.chanNames(kp_gsp_chans));
            vl_nr = vl_nr(:, kp_vl_chans, :);
            cat_vl_nr = vl_nr_noRpts(:, kp_vl_chans, sort_idx_vl_grp);
            chans = ld_vl.info.chanNames(kp_vl_chans);
            pt_rep = repmat({pt}, size(ld_vl.info.chanNames(kp_vl_chans), 1), 1);
        elseif strcmp(which_tasks,'gsp')
            cat_syn_nr = gsp_nr_noRpts(:, :, sort_idx_gsp_grp);
            chans =  ld_gsp.info.chanNames;
            pt_rep = repmat({pt}, size(ld_gsp.info.chanNames,1),1);
            t = ld_gsp.t(kpt);
        else
            cat_vl_nr = vl_nr_noRpts(:, :, sort_idx_vl_grp);
            chans =  ld_vl.info.chanNames;
            pt_rep = repmat({pt}, size(ld_vl.info.chanNames,1),1);
            t = ld_vl.t(kpt);
        end

        if ~strcmp(which_filter,'none')
            [~,pval,ci,stats] = ttest(squeeze(mean(cat(3,vl_nr(t < -100,:,:),...
                gsp_nr(t < -100,:,:)),1))', squeeze(mean(cat(3, vl_nr(t > 0 & t < 550,:,:),...
                gsp_nr(t > 0 & t < 550,:,:)),1))','tail','both');
        elseif strcmp(which_tasks,'gsp')
            [~,pval,ci,stats] = ttest(squeeze(mean(gsp_nr(t < -100, :,:),1))',...
                squeeze(mean(gsp_nr(t > 0 & t < 550,:,:),1))','tail','both');
        else
            [~,pval,ci,stats] = ttest(squeeze(mean(vl_nr(t < -100, :,:), 1))',...
                squeeze(mean(vl_nr(t > 0 & t < 550,:,:), 1))','tail','both');
        end
        fdr = mafdr(pval,'BHFDR',true);

        % downselect to channels with auditory response p<0.01
        kp_chans = (fdr < 0.05)';
        if any(kp_chans)
            chan_labels = getChanLabels(pt, chans(kp_chans), 'useCvs');
            if strcmp(which_filter,'gsp')
                nr_matrix = cat(2, nr_matrix, cat(3, cat_vl_nr(:, kp_chans,:), cat_syn_nr(:, kp_chans, :)));
            elseif strcmp(which_tasks,'gsp')
                nr_matrix = cat(2, nr_matrix, cat_syn_nr(:, kp_chans,:));
            else
                nr_matrix = cat(2, nr_matrix, cat_vl_nr(:, kp_chans,:));
            end
            pt_matrix = cat(1, pt_matrix, pt_rep{kp_chans});
            chan_roi_matrix = cat(1, chan_roi_matrix, cat(2,{chan_labels.chans}',{chan_labels.HCPex}'));
        end
    end
    if strcmp(which_tasks,'gsp')
        kpt = (ld_gsp.t >= -500 & ld_gsp.t <= 2000);
        t = ld_gsp.t(kpt);
    else
        kpt = (ld_vl.t >= -500 & ld_vl.t <= 1000);
        t = ld_vl.t(kpt);
    end
    data_matrix = nr_matrix; clear nr_matrix

    nan_chans = any(any(isnan(data_matrix), 3), 1);
    data_matrix = data_matrix(:, ~nan_chans, :);
    pt_matrix = pt_matrix(~nan_chans);

    chan_roi_matrix = chan_roi_matrix(~nan_chans, :);

    tds = t(1):1000/200:t(end);
    if strcmp(which_data,'hga') 
        data_matrix_ds = data_matrix;
        tds = t;
    else
        data_matrix_ds = nan(size(tds, 2), size(data_matrix, 2), size(data_matrix, 3));
        for iChan = 1:size(data_matrix, 2)
            data_matrix_ds(:, iChan, :) = resample(data_matrix(:, iChan,:), 200, 1000);
        end
    end
    save(fullfile(gcs.fOneDrive, '..\gsp\save\', sprintf('dataIn_%s_%s_%s_process_sEEG.mat',...
        which_data, which_filter, which_tasks)), ...
        'data_matrix', 'pt_matrix', 'chan_roi_matrix', 'pts_rsa', ...
        't', 'tds', 'data_matrix_ds', '-v7.3')
else
    clear nr_matrix hga_matrix data_matrix data_matrix_ds cat_vl_nr cat_syn_nr
    clear all_stim_gsp all_stim_vl gsp_nr noRpts vl_nr chan_labels tds t chan_roi_matrix

    % ----------------------------load neural data-----------------------------
    ld_vl = load(fullfile(gcs.fOneDrive, '..\gsp\save\',sprintf('dataIn_%s_none_vl_process_sEEG.mat', which_data)));
    ld_vl.pt_chan = (string(ld_vl.pt_matrix) + ...
        repmat('_', size(ld_vl.data_matrix_ds, 2), 1) + string({ld_vl.chan_roi_matrix{:, 1}}'));
    [voice_rois_vl, is_L_vl, ~] = defineRois(table(cat(1, ld_vl.chan_roi_matrix(:, 1)),...
        cat(1, ld_vl.chan_roi_matrix(:, 2)), 'VariableNames', {'chans', 'HCPex'}));

    ld_gsp = load(fullfile(gcs.fOneDrive, '..\gsp\save\',sprintf('dataIn_%s_none_gsp_process_sEEG.mat', which_data)));
    ld_gsp.pt_chan = (string(ld_gsp.pt_matrix) + ...
        repmat('_', size(ld_gsp.data_matrix_ds, 2), 1) + string({ld_gsp.chan_roi_matrix{:, 1}}'));
    [voice_rois_gsp, is_L_gsp, ~] = defineRois(table(cat(1, ld_gsp.chan_roi_matrix(:, 1)),...
        cat(1, ld_gsp.chan_roi_matrix(:, 2)), 'VariableNames', {'chans', 'HCPex'}));
end
plot_auditory_roi = {'A1', 'Belt', 'Parabelt', 'STGS','vlPFC','OFC'}; %,'vSTS',...
%     'dlPFC', 'vlPrCG','dlPrCG','FOP'};
plotChansSurf('P25N002',{'Up_02'},'rois',{'lateralorbitofrontal L+lateralorbitofrontal R+medialorbitofrontal L+medialorbitofrontal R',...
    'parsopercularis L+parsopercularis R+parsorbitalis L+parsorbitalis R+parstriangularis L+parstriangularis R',...
    'superiortemporal L+superiortemporal R+bankssts L+bankssts R','transversetemporal L+transversetemporal R'},'roiColor',cbrewer2('Set1',8))
% 'L_medial_belt+R_medial_belt+L_lateral_belt+R_lateral_belt','L_parabelt+R_parabelt'

voice_rois_vl = replace(replace(voice_rois_vl,'dSTS','STG'),'STG','STGS');
voice_rois_gsp = replace(replace(voice_rois_gsp,'dSTS','STG'),'STG','STGS');
voice_rois_vl = replace(voice_rois_vl,'PBelt','Parabelt');
voice_rois_gsp = replace(voice_rois_gsp,'PBelt','Parabelt');

% if you want to down select only to channels with both gsp and vl
kp_gsp_chans = contains(ld_gsp.pt_chan, ld_vl.pt_chan);
kp_vl_chans = contains(ld_vl.pt_chan, ld_gsp.pt_chan);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------voice localizer------------------------------
kpt_acfts = vl.t >= 0 & vl.t <= 500;
t = vl.t(vl.t >= 0 & vl.t <= 500);
stepSize = 0.01;
frameSize = 200;
selT   = find(ld_vl.tds == t(1)) : stepSize*frameSize : find(ld_vl.tds == t(end));
lm_rois = voice_rois_vl(contains(voice_rois_vl, plot_auditory_roi));
nChans  = numel(lm_rois);
ftNames = replace(table2cell(gsp.acft_names(kp_fts,2)),'_',' ');
roi_map = brighten(cbrewer2('Set3',12),-0.4);
which_stim = cat_clr(is_nat) > 0;
nPerm  = 1000;        % number of permutations
usePar = true;        % set false if you don't want parfor
rng(1);               % reproducible

% -------------fit LME on acfts for all channels in an ROI-----------------
nTrials=sum(which_stim);
nTwin = numel(selT);
% X = (acfts_vl(kpt_acfts,:,which_stim) - mean(acfts_vl(kpt_acfts,:,which_stim),...
%     1, 'omitnan')) ./ std(acfts_vl(kpt_acfts,:,which_stim), 1, 'omitnan');
X = permute(acfts_vl(kpt_acfts,:,which_stim), [2,1,3]);
X = squeeze(permute(reshape(X(kp_fts,:,:), numel(kp_fts), nTrials*nTwin),[1,3,2]));
X(isnan(X)) = 0;
[coeff,score] = pca(X);

lm_data = permute(ld_vl.data_matrix_ds(selT,...
    contains(voice_rois_vl, plot_auditory_roi), which_stim), [2,1,3]);
clear beta tstat lo hi pvalue R2
ld_vl.roi_lme_R2 = nan(numel(plot_auditory_roi), 1);
ld_vl.roi_lme_R2_p = nan(numel(plot_auditory_roi), 1);
ld_vl.roi_lme_R2_null = nan(numel(plot_auditory_roi), 1);
ld_vl.roi_lme_beta   =  nan(3,numel(plot_auditory_roi));
ld_vl.roi_lme_tstat  = nan(3,numel(plot_auditory_roi));
ld_vl.roi_lme_lo     = nan(3,numel(plot_auditory_roi));
ld_vl.roi_lme_hi     = nan(3,numel(plot_auditory_roi));
ld_vl.roi_lme_pvalue = nan(3,numel(plot_auditory_roi));
ld_vl.fdr= nan(3,numel(plot_auditory_roi));
for iRoi = 1:numel(plot_auditory_roi)
    iChan = contains(lm_rois, plot_auditory_roi{iRoi});
    fprintf('roi %i/%i\n',iRoi,numel(plot_auditory_roi))
    pt_chan = reshape(repmat(ld_vl.pt_chan(iChan), nTwin,nTrials),[],1);
    y_full = permute(reshape(lm_data(iChan,:,:), sum(iChan), nTrials*nTwin),[3 1 2]);
    fts = reshape(repmat(coeff(:,1:3),1,sum(iChan),1), 3,[])';
    y1 = reshape(y_full, [],1);
    tbl = table(categorical(pt_chan),...
        y1,  fts(:,1),fts(:,2),fts(:,3),...
        'variablenames',{'pt_chan','erp','ft1','ft2','ft3'});
    out = fitlme(tbl,'erp ~ ft1 + ft2 + ft3 + (1|pt_chan)')
    ld_vl.roi_lme_beta(:,iRoi)   = out.Coefficients.Estimate(2:end);
    ld_vl.roi_lme_tstat(:,iRoi)  = out.Coefficients.tStat(2:end);
    ld_vl.roi_lme_lo(:,iRoi)     = out.Coefficients.Lower(2:end);
    ld_vl.roi_lme_hi(:,iRoi)     = out.Coefficients.Upper(2:end);
    ld_vl.roi_lme_pvalue(:,iRoi) = out.Coefficients.pValue(2:end);
    ld_vl.roi_lme_R2(iRoi)       = out.Rsquared.Ordinary;

    coefNames = out.CoefficientNames;                    % cellstr
    isFeat    = startsWith(coefNames,'ft');              % logical index
    H = eye(sum(isFeat));                                % one row per feature
    Hfull = zeros(sum(isFeat), numel(coefNames));
    Hfull(:, isFeat) = H;
    [p_joint, F_joint] = coefTest(out, Hfull);           % joint F-test
    ld_vl.roi_lme_F_joint(iRoi) = F_joint;
    ld_vl.roi_lme_p_joint(iRoi) = p_joint;
    ld_vl.fdr(:,iRoi) = mafdr(ld_vl.roi_lme_pvalue(:,iRoi),'BHFDR', true);
end

models = plot_auditory_roi;
xlabels = {'pc1','pc2','pc3'};
beta = ld_vl.roi_lme_beta;
lo = ld_vl.roi_lme_lo;
hi = ld_vl.roi_lme_hi;
pvals = ld_vl.fdr;
ylab = 'Coefficients';
y0 = 'erp natural';
makeLMEBar(models, xlabels, beta, lo, hi, pvals, ylab, y0)

% ----------------------------gsp-------------------------------------
kpt_acfts = gsp.t>=0 & gsp.t<=500;
t = gsp.t(gsp.t>=0 & gsp.t<=500);
stepSize = 0.01;
frameSize = 200;
selT   = find(ld_gsp.tds == t(1)): stepSize*frameSize : find(ld_gsp.tds == t(end));
lm_rois = voice_rois_gsp(contains(voice_rois_gsp, plot_auditory_roi));
nChans  = numel(lm_rois);
ftNames = replace(table2cell(gsp.acft_names(kp_fts,2)),'_',' ');
roi_map = brighten(cbrewer2('Set3',12),-0.4);
which_stim = cat_clr(~is_nat) > 0;

% ------fit LME on acfts for all channels in an ROI-----
nTrials = sum(which_stim);
nTwin  = numel(selT);
% X = (acfts_gsp(kpt_acfts,:,which_stim) - mean(acfts_gsp(kpt_acfts,:,which_stim),...
%     1, 'omitnan')) ./ std(acfts_gsp(kpt_acfts,:,which_stim), 1, 'omitnan');
% X = permute(X, [2,1,3]);
% X = squeeze(permute(reshape(X(kp_fts,:,:), numel(kp_fts), nTrials*nTwin),[1,3,2]));
X = permute(acfts_gsp(kpt_acfts,:,which_stim), [2,1,3]);
X = squeeze(permute(reshape(X(kp_fts,:,:), numel(kp_fts), nTrials*nTwin),[1,3,2]));
X(isnan(X)) = 0;
nTwin=numel(selT);
X(isnan(X)) = 0;
[coeff,score] = pca(X);

lm_data = permute(ld_gsp.data_matrix_ds(selT,...
    contains(voice_rois_gsp, plot_auditory_roi), which_stim), [2,1,3]);
clear beta tstat lo hi pvalue R2
ld_gsp.roi_lme_R2 = nan(numel(plot_auditory_roi), nTwin);
ld_gsp.roi_lme_beta   =  nan(3,numel(plot_auditory_roi));
ld_gsp.roi_lme_tstat  = nan(3,numel(plot_auditory_roi));
ld_gsp.roi_lme_lo     = nan(3,numel(plot_auditory_roi));
ld_gsp.roi_lme_hi     = nan(3,numel(plot_auditory_roi));
ld_gsp.roi_lme_pvalue = nan(3,numel(plot_auditory_roi));
ld_gsp.fdr = nan(3,numel(plot_auditory_roi));

for iRoi = 1:numel(plot_auditory_roi)
    iChan = contains(lm_rois, plot_auditory_roi{iRoi});
    fprintf('roi %i/%i\n',iRoi,numel(plot_auditory_roi))
    pt_chan = reshape(repmat(ld_gsp.pt_chan(iChan), nTwin,nTrials),[],1);
    y_full = permute(reshape(lm_data(iChan,:,:), sum(iChan), nTrials*nTwin),[3 1 2]);
    fts = reshape(repmat(coeff(:,1:3),1,sum(iChan),1), 3,[])';
    y1 = reshape(y_full, [],1);
    tbl = table(categorical(pt_chan),...
        y1,  fts(:,1),fts(:,2),fts(:,3),...
        'variablenames',{'pt_chan','erp','ft1','ft2','ft3'});
    out = fitlme(tbl,'erp ~ ft1 + ft2 + ft3 + (1|pt_chan)')
    ld_gsp.roi_lme_beta(:,iRoi)   = out.Coefficients.Estimate(2:end);
    ld_gsp.roi_lme_tstat(:,iRoi)   = out.Coefficients.tStat(2:end);
    ld_gsp.roi_lme_lo(:,iRoi)     = out.Coefficients.Lower(2:end);
    ld_gsp.roi_lme_hi(:,iRoi)     = out.Coefficients.Upper(2:end);
    ld_gsp.roi_lme_pvalue(:,iRoi) = out.Coefficients.pValue(2:end);
    ld_gsp.roi_lme_R2(iRoi) = out.Rsquared.Ordinary;
    ld_gsp.fdr(:,iRoi) = mafdr(ld_gsp.roi_lme_pvalue(:,iRoi),'BHFDR', true);

    coefNames = out.CoefficientNames;                    % cellstr
    isFeat    = startsWith(coefNames,'ft');              % logical index
    H = eye(sum(isFeat));                                % one row per feature
    Hfull = zeros(sum(isFeat), numel(coefNames));
    Hfull(:, isFeat) = H;
    [p_joint, F_joint] = coefTest(out, Hfull);           % joint F-test
    ld_gsp.roi_lme_F_joint(iRoi) = F_joint;
    ld_gsp.roi_lme_p_joint(iRoi) = p_joint;

end
models = plot_auditory_roi;
xlabels = {'pc1','pc2','pc3'};
beta = ld_gsp.roi_lme_beta;
lo = ld_gsp.roi_lme_lo;
hi = ld_gsp.roi_lme_hi;
pvals = ld_gsp.fdr;
ylab = 'Coefficients';
y0 = 'erp synthetic';
makeLMEBar(models, xlabels, beta, lo, hi, pvals, ylab, y0)

%--------------------------corr clustering------------------------------
clear corr_rate* corr_yam* clust_pcorr* clust_corr*
nBs = 1000;
which_tasks = 'gsp';
if isfile(fullfile(gcs.fOneDrive, '..\gsp\save\', ...
        sprintf('all_out_webb_%s_%s_%s_process_sEEG_nbs_%i_new.mat',...
        which_data, which_filter, which_tasks, nBs)))
    load(fullfile(gcs.fOneDrive, '..\gsp\save\', ...
        sprintf('all_out_webb_%s_%s_%s_process_sEEG_nbs_%i_new.mat',...
        which_data, which_filter, which_tasks, nBs)))
else
    if strcmp(which_tasks, 'gsp') && ~exist('corr_rate_syn','var')
        for iChan = 1:size(ld_gsp.data_matrix_ds,2)
            [corr_rate_syn(iChan,:), ~] = corr(permute(ld_gsp.data_matrix_ds(:, iChan, :),...
                [3 1 2]), cat_rate(~is_nat));
        end
        [clust_corr_rate_syn, sig_mask_corr_rate_syn] = permClustStat(ld_gsp.data_matrix_ds,...
            ld_gsp.tds, [], false, 'corr', cat_rate(~is_nat), [], 1000, 0.01);
        sig_mask_corr_rate_syn(sig_mask_corr_rate_syn==0) = nan;
        sig_mask_corr_rate_syn(:, all(sig_mask_corr_rate_syn==1,1)) = nan;
        save(fullfile(gcs.fOneDrive, '..\gsp\save\',  sprintf('all_out_webb_%s_%s_%s_process_sEEG_nbs_%i_new.mat',...
            which_data, which_filter, which_tasks, nBs)),...
            'corr_rate_syn','clust_corr_rate_syn', 'sig_mask_corr_rate_syn');
        %            'corr_yam_syn', 'clust_corr_yam_syn','sig_mask_corr_yam_syn');
    end
end

which_tasks = 'vl';
if isfile(fullfile(gcs.fOneDrive, '..\gsp\save\', sprintf('all_out_webb_%s_%s_%s_process_sEEG_nbs_%i_new.mat',...
        which_data, which_filter, which_tasks, nBs)))
    load(fullfile(gcs.fOneDrive, '..\gsp\save\', sprintf('all_out_webb_%s_%s_%s_process_sEEG_nbs_%i_new.mat',...
        which_data, which_filter, which_tasks, nBs)))
else
    if strcmp(which_tasks, 'vl') && ~exist('corr_rate_nat','var')
        for iChan = 1:size(ld_vl.data_matrix_ds,2)
            [corr_rate_nat(iChan,:), ~] = corr(permute(ld_vl.data_matrix_ds(:,iChan, :),...
                [3 1 2]), cat_rate(is_nat));
        end
        [clust_corr_rate_nat, sig_mask_corr_rate_nat] = permClustStat(ld_vl.data_matrix_ds,...
            ld_vl.tds, [], false, 'corr', cat_rate(is_nat), [], 1000, 0.01);
        sig_mask_corr_rate_nat(sig_mask_corr_rate_nat==0) = nan;
        sig_mask_corr_rate_nat(:, all(sig_mask_corr_rate_nat==1,1)) = nan;
        save(fullfile(gcs.fOneDrive, '..\gsp\save\',  sprintf('all_out_webb_%s_%s_%s_process_sEEG_nbs_%i_new.mat',...
            which_data, which_filter, which_tasks, nBs)),...
            'corr_rate_nat','clust_corr_rate_nat', 'sig_mask_corr_rate_nat');
        %             'corr_yam_nat', 'clust_corr_yam_nat','sig_mask_corr_yam_nat');
    end
end

% -----------------get trial wise onsets and peak characteristics---------
which_tasks = 'vl';
if ~isfield(ld_vl,'peaks') || ~isfield(ld_vl,'onsets')
    opts.num_bins = 12;
    opts.smooth_ms = 40;
    opts.snr_on_z = true;
    opts.winKp = [0 500];
    opts.sustain_ms = 50;
    opts.baseline = 'none';
    opts.p_thresh = 0.05;
    onsets = get_trial_onsets(ld_vl.data_matrix, ld_vl.t, opts);
    peaks = get_erp_peaks(ld_vl.data_matrix, ld_vl.t, opts);
    peaks.peak_amp(std(peaks.peak_amp,[],2)<0.5,:)= nan;
    peaks.peak_lat_ms(std(peaks.peak_amp,[],2)<0.5,:)= nan;
    save(fullfile(gcs.fOneDrive, '..\gsp\save\',  sprintf('all_out_webb_%s_%s_%s_process_sEEG_nbs_%i_new.mat',...
        which_data, which_filter, which_tasks, nBs)),'onsets','peaks','-append')
else
    load(fullfile(gcs.fOneDrive, '..\gsp\save\',  sprintf('all_out_webb_%s_%s_%s_process_sEEG_nbs_%i_new.mat',...
        which_data, which_filter, which_tasks, nBs)),'peaks','onsets')
end
    ld_vl.peaks = peaks;
    ld_vl.onsets = onsets;

which_tasks = 'gsp';
if ~isfield(ld_gsp,'peaks') || ~isfield(ld_gsp,'onsets')
    opts.num_bins = 12;
    opts.smooth_ms = 40;
    opts.snr_on_z = true;
    opts.winKp = [0 500];
    opts.sustain_ms = 50;
    opts.baseline = 'none';
    opts.p_thresh = 0.05;
    onsets = get_trial_onsets(ld_gsp.data_matrix, ld_gsp.t, opts);
    peaks = get_erp_peaks(ld_gsp.data_matrix, ld_gsp.t, opts);
    peaks.peak_amp(std(peaks.peak_amp, [], 2) < 0.5,:)= nan;
    peaks.peak_lat_ms(std(peaks.peak_amp, [], 2) < 0.5,:)= nan;
    save(fullfile(gcs.fOneDrive, '..\gsp\save\',  sprintf('all_out_webb_%s_%s_%s_process_sEEG_nbs_%i_new.mat',...
        which_data, which_filter, which_tasks, nBs)),'onsets','peaks','-append')
else
    load(fullfile(gcs.fOneDrive, '..\gsp\save\',  sprintf('all_out_webb_%s_%s_%s_process_sEEG_nbs_%i_new.mat',...
        which_data, which_filter, which_tasks, nBs)),'peaks','onsets')
end
    ld_gsp.peaks = peaks;
    ld_gsp.onsets = onsets;
% ---------------------------- make tables-------------------------------
pts = unique(ld_vl.pt_matrix);
demo = getDemographics(pts);
clear age_rep
for iPt = 1:numel(pts)
    age_rep{iPt} = repmat(demo.AgeAtSEEG_years_(strcmp(demo.sEEGID, pts{iPt})),...
        sum(contains(ld_vl.pt_matrix,pts{iPt})),1);
end
age_rep = cat(1,age_rep{:});
pts_rep_trials = repmat(categorical(ld_vl.pt_matrix), size(cat_rate(is_nat), 1), 1);
pt_chan_rep_trials = repmat(categorical(string(ld_vl.pt_matrix) + ...
    repmat('_', size(ld_vl.data_matrix_ds,2), 1) + string(cat(1,ld_vl.chan_roi_matrix(:,1)))), size(cat_rate(is_nat), 1), 1);
age_rep_trials = repmat(age_rep, size(cat_rate(is_nat),1), 1);
chan_rep_trials = repmat(cat(1,ld_vl.chan_roi_matrix(:,1)),size(cat_rate(is_nat), 1),1);
voice_rois_trials = repmat(replace(replace(voice_rois_vl','L_',''),'R_',''), size(cat_rate(is_nat), 1), 1);
is_L_trials = repmat(logical(contains(voice_rois_vl','L_')), size(cat_rate(is_nat), 1), 1);
data_trials_mean = reshape(mean(ld_vl.data_matrix_ds(ld_vl.tds>=100 & ld_vl.tds<=400, :, :)), [], 1);
data_trials_mean_early = reshape(mean(ld_vl.data_matrix_ds(ld_vl.tds>=0 & ld_vl.tds<=200, :, :)), [], 1);
data_trials_mean_late = reshape(mean(ld_vl.data_matrix_ds(ld_vl.tds>=200 & ld_vl.tds<=400, :, :)), [], 1);
cat_trial_onsets = reshape(ld_vl.onsets.detected_onsets,[],1);
cat_peak_lat_ms = reshape(ld_vl.peaks.peak_lat_ms,[],1);
cat_peak_amp = reshape(ld_vl.peaks.peak_amp, [], 1);
cat_rate_trials = reshape(repmat(rescale(cat_rate(is_nat),0,1), 1, size(ld_vl.data_matrix_ds, 2))', [], 1);
cat_clr_trials = reshape(repmat(cat_clr(is_nat), 1, size(ld_vl.data_matrix_ds, 2))', [], 1);
cat_yam_trials = reshape(repmat(cat_yam(is_nat), 1, size(ld_vl.data_matrix_ds, 2))', [], 1);
cat_clr_yam_trials = reshape(repmat(cat_clr_yam(is_nat), 1, size(ld_vl.data_matrix_ds, 2))', [], 1);
cat_tcoeff_trials = reshape(repmat(double(vl.tcoeff(:,1)), 1, size(ld_vl.data_matrix_ds, 2))', [], 1);
cat_fcoeff_trials = reshape(repmat(double(vl.fcoeff(:,1)), 1, size(ld_vl.data_matrix_ds, 2))', [], 1);
cat_task = categorical(reshape(repmat(is_nat(is_nat), 1, size(ld_vl.data_matrix_ds, 2))', [], 1));
cat_kp_vl_chan_trials = logical(repmat(kp_vl_chans, size(cat_rate(is_nat), 1), 1));
tbl_vl = table(pts_rep_trials, pt_chan_rep_trials, age_rep_trials, voice_rois_trials, ...
    chan_rep_trials, is_L_trials, ...
    data_trials_mean_early, data_trials_mean_late, data_trials_mean, ...
    cat_rate_trials, ...
    cat_clr_trials, cat_tcoeff_trials, cat_fcoeff_trials, cat_trial_onsets, cat_peak_lat_ms,...
    cat_peak_amp, double(cat_yam_trials), cat_clr_yam_trials, cat_task,cat_kp_vl_chan_trials,...
    'VariableNames',{'pt','pt_chan','age','roi','chan','isL','erp_early', 'erp_late', ...
    'erp', 'ratings',...
    'rate_quartile', 'tcoeff', 'fcoeff','onset','peak_latency','peak_amp',...
    'yamout', 'yamout_quartile', 'task','kp_chan'});
clear chan_rep_trials is_L_trials data_trials_mean_early data_trials_mean_late data_trials_mean...
    cat_rate_trials cat_clr_trials cat_tcoeff_trials cat_fcoeff_trials cat_trial_onsets cat_peak_lat_ms...
    cat_peak_amp cat_yam_trials cat_clr_yam_trials cat_task cat_kp_vl_chan_trials

% tbl_vl.isL = categorical(tbl_vl.isL, [false true], {'Right', 'Left'});
% tbl_vl.task = categorical(tbl_vl.task, [false true], {'Syn', 'Nat'});
kp_rois_vl = contains(tbl_vl.roi, plot_auditory_roi);

pts = unique(ld_gsp.pt_matrix);
demo = getDemographics(pts);
clear age_rep
for iPt = 1:numel(pts)
    age_rep{iPt} = repmat(demo.AgeAtSEEG_years_(strcmp(demo.sEEGID, pts{iPt})),...
        sum(contains(ld_gsp.pt_matrix,pts{iPt})),1);
end
age_rep = cat(1,age_rep{:});
pts_rep_trials = repmat(categorical(ld_gsp.pt_matrix), size(cat_rate(~is_nat), 1), 1);
pt_chan_rep_trials = repmat(categorical(string(ld_gsp.pt_matrix) + ...
    repmat('_', size(ld_gsp.data_matrix_ds,2), 1) + string(cat(1,ld_gsp.chan_roi_matrix(:,1)))), ...
    size(cat_rate(~is_nat), 1), 1);
age_rep_trials = repmat(age_rep, size(cat_rate(~is_nat),1), 1);
chan_rep_trials = repmat(cat(1,ld_gsp.chan_roi_matrix(:,1)),size(cat_rate(~is_nat), 1),1);
voice_rois_trials = repmat(replace(replace(voice_rois_gsp','L_',''),'R_',''), size(cat_rate(~is_nat), 1), 1);
is_L_trials = repmat(logical(contains(voice_rois_gsp','L_')), size(cat_rate(~is_nat), 1), 1);
data_trials_mean = reshape(mean(ld_gsp.data_matrix_ds(ld_gsp.tds>=100 & ld_gsp.tds<=400, :, :)), [], 1);
data_trials_mean_early = reshape(mean(ld_gsp.data_matrix_ds(ld_gsp.tds>=0 & ld_gsp.tds<=200, :, :)), [], 1);
data_trials_mean_late = reshape(mean(ld_gsp.data_matrix_ds(ld_gsp.tds>=200 & ld_gsp.tds<=400, :, :)), [], 1);
cat_trial_onsets = reshape(ld_gsp.onsets.detected_onsets,[],1);
cat_peak_lat_ms = reshape(ld_gsp.peaks.peak_lat_ms,[],1);
cat_peak_amp = reshape(ld_gsp.peaks.peak_amp, [], 1);
cat_rate_trials = reshape(repmat(rescale(cat_rate(~is_nat),0,1), 1, size(ld_gsp.data_matrix_ds, 2))', [], 1);
cat_clr_trials = reshape(repmat(cat_clr(~is_nat), 1, size(ld_gsp.data_matrix_ds, 2))', [], 1);
cat_yam_trials = reshape(repmat(cat_yam(~is_nat), 1, size(ld_gsp.data_matrix_ds, 2))', [], 1);
cat_clr_yam_trials = reshape(repmat(cat_clr_yam(~is_nat), 1, size(ld_gsp.data_matrix_ds, 2))', [], 1);
cat_tcoeff_trials = reshape(repmat(double(gsp.tcoeff(:,1)), 1, size(ld_gsp.data_matrix_ds, 2))', [], 1);
cat_fcoeff_trials = reshape(repmat(double(gsp.fcoeff(:,1)), 1, size(ld_gsp.data_matrix_ds, 2))', [], 1);
cat_task = categorical(reshape(repmat(is_nat(~is_nat), 1, size(ld_gsp.data_matrix_ds, 2))', [], 1));
cat_kp_gsp_chan_trials = logical(repmat(kp_gsp_chans, size(cat_rate(~is_nat), 1), 1));
tbl_gsp = table(pts_rep_trials, pt_chan_rep_trials, age_rep_trials, voice_rois_trials, ...
    chan_rep_trials, is_L_trials, ...
    data_trials_mean_early, data_trials_mean_late, data_trials_mean, cat_rate_trials, ...
    cat_clr_trials, cat_tcoeff_trials, cat_fcoeff_trials, cat_trial_onsets, cat_peak_lat_ms,...
    cat_peak_amp, double(cat_yam_trials), cat_clr_yam_trials, cat_task, cat_kp_gsp_chan_trials,...
    'VariableNames',{'pt','pt_chan','age','roi','chan','isL','erp_early', 'erp_late', ...
    'erp', 'ratings',...
    'rate_quartile', 'tcoeff', 'fcoeff','onset','peak_latency','peak_amp',...
    'yamout', 'yamout_quartile', 'task', 'kp_chan'});
% tbl_vl.isL = categorical(tbl_vl.isL, [false true], {'Right', 'Left'});
% tbl_vl.task = categorical(tbl_vl.task, [false true], {'Syn', 'Nat'});
kp_rois_gsp = contains(tbl_gsp.roi, plot_auditory_roi);
clear chan_rep_trials is_L_trials data_trials_mean_early data_trials_mean_late data_trials_mean...
    cat_rate_trials cat_clr_trials cat_tcoeff_trials cat_fcoeff_trials cat_trial_onsets cat_peak_lat_ms...
    cat_peak_amp cat_yam_trials cat_clr_yam_trials cat_task cat_kp_vl_chan_trials


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ----------------------- fit lmes and plot coefficients------------------
eqs = {'erp  ~  ratings + (1|pt_chan) + (ratings|pt_chan)'};%,...
%     'peak_latency  ~  ratings*isL + task + (1|task) + (1|pt) + (ratings|pt_chan)',...
%     'peak_amp  ~  ratings*isL + task + (1|task) + (1|pt) + (ratings|pt_chan)',...
%     'onset  ~  ratings*isL + task + (1|task) + (1|pt) + (ratings|pt_chan)'};

models = {{'ratings'}};%,{'isL*ratings'},{'task'}};...
%      {'isL'},{'ratings'},{'isL*ratings'},{'task'};...
%       {'isL'},{'ratings'},{'isL*ratings'},{'task'};...
%        {'isL'},{'ratings'},{'isL*ratings'},{'task'}};
% models = models(:,1:3);
ylab = 'Coefficients';
shift = [-0.1250, 0.1250];
y0 = {'erp'};%,'peak_latency','peak_amp','onset'};

big_tbl = false;

for iLme =1:numel(eqs)
    clear beta lo hi pvalue yerr
    for iTask = 1:2
        if big_tbl
            tbl  = cat(1, tbl_vl, tbl_gsp);
        elseif ~big_tbl && iTask == 1
            tbl = tbl_vl;
            %             tbl = tbl_vl(tbl_vl.kp_chan, :);
        else
            %             tbl = tbl_gsp(tbl_gsp.kp_chan, :);
            tbl = tbl_gsp;
        end
        for iRoi = 1:numel(plot_auditory_roi)
            lme = fitlme(tbl(contains(tbl.roi, plot_auditory_roi{iRoi}),:),...
                eqs{iLme});
            display(lme)
            beta(iRoi,:,iTask) = lme.Coefficients.Estimate(2:end);
            lo(iRoi,:,iTask)   = lme.Coefficients.Lower(2:end);
            hi(iRoi,:,iTask)   = lme.Coefficients.Upper(2:end);
            pvalue(iRoi,:,iTask)   = mafdr(lme.Coefficients.pValue(2:end),'BHFDR',true);
        end
    end
    xlabels = plot_auditory_roi;
    x = 1:numel(xlabels);
    for imdl = 1:numel(models(iLme,:))
        figure('Color','w','Position',[332,689,292,200]);
        if big_tbl
            b = bar(x, squeeze(beta(:,imdl,iTask)), ...
                'EdgeColor','none','barwidth',0.3);
        else
            b = bar(x, squeeze(beta(:,imdl,:)),'grouped',...
                'EdgeColor','none');
        end
        ttl = sprintf('%s on %s', string(models{iLme,imdl}), y0{iLme});
        hold on;

        for iTask = 1:2
            yerr = [beta(:,imdl,iTask) - lo(:,imdl,iTask), hi(:,imdl,iTask)- beta(:,imdl,iTask)];
            if big_tbl
                errorbar(x , squeeze(beta(:,imdl,iTask)), yerr(:,1), yerr(:,2), 'k', ...
                    'LineStyle','none', 'CapSize',4, 'LineWidth',1);
            else
                errorbar(x + shift(iTask) , squeeze(beta(:,imdl,iTask)), yerr(:,1), yerr(:,2), 'k', ...
                    'LineStyle','none', 'CapSize',4, 'LineWidth',1);
            end
            for i = 1:numel(xlabels)
                s = p2stars(pvalue(i,imdl,iTask));
                if ~isempty(s)
                    y = beta(i,imdl,iTask);
                    if y >= 0
                        ystar = y + yerr(i,1);
                        va = 'bottom';
                    else
                        ystar = y - yerr(i,2);
                        va = 'top';
                    end
                    if big_tbl
                        text(i, ystar, s, 'HorizontalAlignment','center', ...
                            'VerticalAlignment',va, 'FontSize',10, 'FontWeight','bold');
                    else
                        text(i + shift(iTask), ystar, s, 'HorizontalAlignment','center', ...
                            'VerticalAlignment',va, 'FontSize',10, 'FontWeight','bold');
                    end
                end

            end
        end
        title(ttl);
        ylim([min(beta(:,imdl,:) - yerr(:,1,:),[],'all') - max(yerr,[],'all')*1.5,...
            max(beta(:,imdl,:) + yerr(:,2,:),[],'all') + max(yerr,[],'all')*1.5])
        xlim([0.5 numel(xlabels)+0.5]);
        set(gca,'XTick',x,'XTickLabel',erase(xlabels,'_'),...
            'FontSize',12,'LineWidth',1);
        xlabel('ROI');
        ylabel(ylab);
        box on;
    end
end

clear ylab shift y0 eqs models ylab big_tbl iLme iRoi iChan beta lo hi pvalue xlabels x
% ------------------ make time varying stats plots------------------------
clr_map_roi = brighten(cbrewer2('Set3', 9),-0.5);
clr_map_roi(2,:) = clr_map_roi(end,:);
clr_map_task =  [0.00,0.45,0.74;0.85,0.33,0.10];
% % plot sig ROI channels for cluster over time
% sum_sig_mask_voice = nan(numel(plot_auditory_roi),numel(ld_vl.tds(ld_vl.tds>=-70 & ld_vl.tds<=650)));
% sum_sig_mask_onset = nan(numel(plot_auditory_roi),numel(ld_vl.tds(ld_vl.tds>=-70 & ld_vl.tds<=650)));
% for iRoi = 1:numel(plot_auditory_roi)
%     voice_rois = erase(erase(voice_rois_vl,'L_'),'R_');
%     iChan = contains(voice_rois, plot_auditory_roi{iRoi});
%     if any(any(sig_mask_corr_rate_nat(ld_vl.tds>=-70  & ld_vl.tds<=650, iChan)))
%         figure('color','w','Position',[100+230*iRoi,600,250,230]); hold on
%         sum_sig_mask_voice(iRoi,:) = smoothdata(sum(sig_mask_corr_rate_nat(ld_vl.tds>=-70  & ld_vl.tds<=650, iChan), 2,'omitnan')/...
%             sum(any(sig_mask_corr_rate_nat(ld_vl.tds>=-70  & ld_vl.tds<=650,iChan),1)),1,'gaussian',20)'*100;
%         plot(ld_vl.tds(ld_vl.tds>=-70  & ld_vl.tds<=650)+10, sum_sig_mask_voice(iRoi,:), 'color', clr_map_task(1,:),'linewidth',1.5)
%         voice_rois = erase(erase(voice_rois_gsp,'L_'),'R_');
%         iChan = contains(voice_rois, plot_auditory_roi{iRoi});
%         sum_sig_mask = smoothdata(sum(sig_mask_corr_rate_syn(ld_gsp.tds>=-70  & ld_gsp.tds<=650, iChan), 2,'omitnan')/...
%             sum(any(sig_mask_corr_rate_syn(ld_gsp.tds>=-70  & ld_gsp.tds<=650,iChan),1)),1,'gaussian',20)'*100;
%         plot(ld_gsp.tds(ld_gsp.tds>=-70  & ld_gsp.tds<=650), sum_sig_mask, 'color',  clr_map_task(2,:),'linewidth',1.5)
%         axis tight; xlim([-70 650]);
%         legend({'nat','syn'},'location','best')
%         title(sprintf('%s\nNat n= %i, Syn n= %i', plot_auditory_roi{iRoi},...
%             sum(contains(voice_rois_vl,  plot_auditory_roi{iRoi})), sum(contains(voice_rois_gsp,  plot_auditory_roi{iRoi}))))
%         ylabel('percent of ROI'); xlabel('time (ms)')
%         set(gca,'fontsize',12,'fontname','arial')
%         if iRoi <7
%             ylim([0 60])
%         else
%             ylim([0 20])
%         end
%     end
% end

% % plot sig ROI channels for cluster over time
% for iRoi = 4%:numel(plot_auditory_roi)
%     iChan = contains(voice_rois_vl, plot_auditory_roi{iRoi})';
%     plotChan = sum(sig_mask_corr_rate_nat(ld_vl.tds>=0  & ld_vl.tds<=650,:),1,'omitnan')'>8;
%     iChan = iChan & plotChan;
%     if any(any(sig_mask_corr_rate_nat(ld_vl.tds>=-70  & ld_vl.tds<=650, iChan)))
%         figure('color','w','Position',[100+250*iRoi,600,250,230]); hold on
%         %         sum_sig_mask_voice = smoothdata(sum(sig_mask_corr_rate_nat(ld_vl.tds>=-70  & ld_vl.tds<=650, :), 1,'omitnan')/...
%         %             sum(any(sig_mask_corr_rate_nat(ld_vl.tds>=-70  & ld_vl.tds<=650,:),2)),1,'gaussian',20)'*100;
%         y= smoothdata(median(corr_rate_nat(iChan, ld_vl.tds>=-70  & ld_vl.tds<=650),1,'omitnan'),1,'gaussian',3);
%         plot(ld_vl.tds(ld_vl.tds>=-70  & ld_vl.tds<=650), y,'color', clr_map_task(1,:),'linewidth', 3)
%
%         face_alpha= rescale(sum(sig_mask_corr_rate_nat(ld_vl.tds>=0  & ld_vl.tds<=650, iChan),2,'omitnan')/...
%             mean(sum(sig_mask_corr_rate_nat(ld_vl.tds>=0  & ld_vl.tds<=650, iChan), 1,'omitnan'),'all'),0,0.5);
%         plott = ld_vl.tds(ld_vl.tds>=0 & ld_vl.tds<=650);
%         markers = any(sig_mask_corr_rate_nat(ld_vl.tds>=0  & ld_vl.tds<=650,iChan),2)*-0.3;
%         markers(markers==0)=nan;
%         for iT = 1:numel(face_alpha)
%             scatter(plott(iT), markers(iT),30,'filled','MarkerFaceAlpha', face_alpha(iT),'markerfacecolor',clr_map_task(1,:))
%         end
%         vl_chans = sum(iChan);
%
%
%         iChan = contains(voice_rois_gsp, plot_auditory_roi{iRoi})';
%         plotChan = sum(sig_mask_corr_rate_syn(ld_gsp.tds>=0  & ld_gsp.tds<=650,:),1,'omitnan')'>8;
%         iChan = iChan & plotChan;
%         %         sum_sig_mask_voice = smoothdata(sum(sig_mask_corr_rate_syn(ld_gsp.tds>=-70  & ld_gsp.tds<=650, :), 1,'omitnan')/...
%         %             sum(any(sig_mask_corr_rate_syn(ld_gsp.tds>=-70  & ld_gsp.tds<=650,:),2)),1,'gaussian',20)'*100;
%
%         y=smoothdata(median(corr_rate_syn(iChan, ld_gsp.tds>=-70  & ld_gsp.tds<=650),1,'omitnan'),1,'gaussian',25);
%         plot(ld_gsp.tds(ld_gsp.tds>=-70  & ld_gsp.tds<=650),y,'color', clr_map_task(end,:),'linewidth', 3)
%         face_alpha= rescale(sum(sig_mask_corr_rate_syn(ld_gsp.tds>=0  & ld_gsp.tds<=650, iChan),2,'omitnan')/...
%             mean(sum(sig_mask_corr_rate_syn(ld_gsp.tds>=0  & ld_gsp.tds<=650, iChan), 1,'omitnan'),'all'),0,0.5);
%         plott = ld_gsp.tds(ld_gsp.tds>=0 & ld_gsp.tds<=650);
%         markers = any(sig_mask_corr_rate_syn(ld_gsp.tds>=0  & ld_gsp.tds<=650, iChan),2)*-0.32;
%         markers(markers==0)=nan;
%         for iT = 1:numel(face_alpha)
%             scatter(plott(iT), markers(iT),30,'filled','MarkerFaceAlpha', face_alpha(iT),'markerfacecolor',clr_map_task(end,:))
%         end
%
%         axis tight; xlim([-70 650]);
% %         legend({'nat','syn'},'location','best')
%         title(sprintf('%s\nNat n= %i, Syn n= %i', plot_auditory_roi{iRoi},...
%             vl_chans, sum(iChan)))
%         ylabel('rho'); xlabel('time (ms)')
%         set(gca,'fontsize',12,'fontname','arial')
%         ylim([-0.35 0.2])
%     end
% end

figure('color','w','Position',[100, 590, 240, 180]); hold on

for iRoi = 4%:numel(plot_auditory_roi)
    iChan = contains(voice_rois_vl, plot_auditory_roi{iRoi})';
    plotChan = sum(sig_mask_corr_rate_nat(ld_vl.tds >= 200 & ld_vl.tds <= 550,:),1,'omitnan')' > 8;
    iChan = find(iChan & plotChan);
    steps = linspace(min(ratings_vl_sort), max(ratings_vl_sort), 40);
    for i=1:numel(iChan)
        for iQ = 1:numel(steps)-1
            scatter(mean(ratings_vl_sort(ratings_vl_sort >= steps(iQ) & ratings_vl_sort <= steps(iQ+1))), ...
                squeeze(mean(ld_vl.data_matrix_ds(ld_vl.tds==300, iChan(i), (ratings_vl_sort >= steps(iQ) & ratings_vl_sort <= steps(iQ + 1))),3))',...
                'filled','SizeData', 20, 'markerfacealpha', 0.5, 'markerfacecolor', clr_map_task(1,:))
        end
    end

    steps = linspace(min(ratings_gsp_sort), max(ratings_gsp_sort), 20);
    iChan = contains(voice_rois_gsp, plot_auditory_roi{iRoi})';
    plotChan = sum(sig_mask_corr_rate_syn(ld_gsp.tds >= 200 & ld_gsp.tds <= 550,:),1,'omitnan')' > 8;
    iChan = find(iChan & plotChan);
    for i=1:numel(iChan)
        for iQ = 1:numel(steps)-1
            scatter(mean(ratings_gsp_sort(ratings_gsp_sort >= steps(iQ) & ratings_gsp_sort <= steps(iQ + 1))), ...
                squeeze(mean(ld_gsp.data_matrix_ds(ld_gsp.tds==200, iChan(i), (ratings_gsp_sort >= steps(iQ) & ratings_gsp_sort <= steps(iQ + 1))),3))',...
                'filled','sizedata', 20, 'markerfacealpha', 0.5, 'markerfacecolor', clr_map_task(2, :))
        end
        set(gca, 'fontsize',12); xlabel('Rating'); xticks([-1 1.8]); xlim([-1,1.8]); xticklabels({'NV','V'}); yticks([])
    end
end

%% rsa stuff
cat_clr_map = [cbrewer2('set1', 4); brighten(cbrewer2('set1', 4), 0.6)];
nPerm = 5000;
t0 = 0;
t1 = 500;
stepSize = 10;
frameSize = 200;
fs = 200;
rate_cutoff = 0;
nr_samps_per_frame = 200/ (1000/200);
% -----openSmile opts-----
% frameMode = fixed
% frameSize = 0.200
% frameStep = 0.01
% frameCenterSpecial = left
% noPostEOIprocessing = 0
% ----------------------------natural sounds ------------------------------
    metric = 'correlation';
if ~exist('rsa_ratings_vl','var')
    kpt_acfts = vl.t>=  t0 & vl.t<= t1;
    t = vl.t(kpt_acfts);
    kpT   = find(ld_vl.tds == t0):find(ld_vl.tds == t1);
    selT   = 1: stepSize*frameSize/1000 : find(ld_vl.tds == t(end));
    which_stim = cat_clr(is_nat) > rate_cutoff;
    D_rating_vl = squareform(pdist(cat_rate(is_nat & cat_clr > rate_cutoff),'euclidean'));
   
    corr_coeffs = cat(2,vl.tcoeff(which_stim,1),vl.fcoeff(which_stim,1));
    corr_coeffs(abs(zscore(corr_coeffs(:,1),[],'omitnan'))>6,1) = nan;
    corr_coeffs(abs(zscore(corr_coeffs(:,2),[],'omitnan'))>6,2) = nan;
    D_acoustics_vl = squareform(pdist(cat(2,zscore(corr_coeffs,[],'omitnan'),...
        permute(zscore(mean(vl.acfts(kpt_acfts, kp_fts,...
        which_stim),'omitnan'),[],3),[3,2,1])),'correlation')); % use correlation! zscore across stimuli after averaging across time.
%     D_acoustics_vl = squareform(pdist(permute(mean(zscore(vl.acfts(kpt_acfts, kp_fts,...
%         which_stim),[],3),'omitnan'),[3,2,1]),'correlation'));
    rsa_ac_vl = zscore(vl.acfts(kpt_acfts, kp_fts,which_stim),[],3); %zscore across stimuli
    rsa_data_vl = permute(ld_vl.data_matrix_ds(:, :, which_stim), [3 2 1]);
    for iRoi = 1:numel(plot_auditory_roi)
        iChan = contains(voice_rois_vl, plot_auditory_roi{iRoi});
        D_neural_vl(:,:,iRoi) = squareform(pdist(reshape(rsa_data_vl(:, iChan, kpT),sum(which_stim),[]),metric));
        [rho_vl(iRoi), p_vl(iRoi)] = mantel_spearman(D_neural_vl(:,:,iRoi), D_rating_vl, nPerm, 1);
        for iT = 1:numel(t)
            nr_win_idx = find(ld_vl.tds==t(iT)) : (find(ld_vl.tds==t(iT)) + nr_samps_per_frame/4);
            D_neural_iT= squareform(pdist(reshape(rsa_data_vl(:, iChan, nr_win_idx), sum(which_stim), []), metric));
            [rho_iT_vL(iRoi,iT), p_iT_vL(iRoi,iT)] = mantel_spearman(D_neural_iT, D_rating_vl, nPerm, 1);
        end
        [rho_ac_vl(iRoi), p_ac_vl(iRoi)] = mantel_spearman(D_neural_vl(:,:,iRoi), D_acoustics_vl, nPerm, 1);
        for iT = 1:numel(t)
            nr_win_idx = find(ld_vl.tds==t(iT)) : (find(ld_vl.tds==t(iT)) + nr_samps_per_frame/4);
            D_acoustics_vl_iT = squareform(pdist(squeeze(rsa_ac_vl(iT,:,:))','correlation'));
            D_neural_iT= squareform(pdist(reshape(rsa_data_vl(:, iChan, nr_win_idx), sum(which_stim), []), metric));
            [rho_ac_iT_vL(iRoi,iT), p_ac_iT_vL(iRoi,iT)] = mantel_spearman(D_neural_iT, D_acoustics_vl_iT, nPerm, 1);
        end
        [rho_pcorr_vl(iRoi), p_pcorr_vl(iRoi)] = mantel_spearman(D_neural_vl(:,:,iRoi), ...
            D_rating_vl, nPerm, 1,'pcorr', D_acoustics_vl);
        for iT = 1:numel(t)
            nr_win_idx = find(ld_vl.tds==t(iT)) : (find(ld_vl.tds==t(iT)) + nr_samps_per_frame/4);
            D_acoustics_vl_iT = squareform(pdist(squeeze(rsa_ac_vl(iT,:,:))','correlation'));
            D_neural_iT= squareform(pdist(reshape(rsa_data_vl(:, iChan, nr_win_idx), sum(which_stim), []), metric));
            [rho_pcorr_iT_vL(iRoi,iT), p_pcorr_iT_vL(iRoi,iT)] = mantel_spearman(D_neural_iT, ...
                D_rating_vl, nPerm, 1, 'pcorr', D_acoustics_vl_iT);
        end
    end
    rsa_ratings_vl = struct('rho',rho_vl,'p',p_vl, 'p_iT',p_iT_vL,'rho_iT',rho_iT_vL,...
        'rho_ac',rho_ac_vl,'p_ac',p_ac_vl, 'p_ac_iT',p_ac_iT_vL,'rho_ac_iT',rho_ac_iT_vL,...
        'rho_pcorr',rho_pcorr_vl,'p_pcorr',p_pcorr_vl, 'p_pcorr_iT',p_pcorr_iT_vL,...
        'rho_pcorr_iT',rho_pcorr_iT_vL,'nPerm',nPerm);
    save(fullfile(gcs.fOneDrive, '..\gsp\save\',  sprintf('all_out_webb_%s_%s_%s_process_sEEG_nbs_%i_new.mat',...
        which_data, which_filter, 'vl', nBs)),'rsa_ratings_vl','-append');
else
    load(fullfile(gcs.fOneDrive, '..\gsp\save\',  sprintf('all_out_webb_%s_%s_%s_process_sEEG_nbs_%i_new.mat',...
        which_data, which_filter, 'vl', nBs)), 'rsa_ratings_vl')
end
  
kpT   = find(ld_vl.tds == t0):find(ld_vl.tds == t1);
which_stim = cat_clr(is_nat) > rate_cutoff;
rsa_data_vl = permute(ld_vl.data_matrix_ds(kpT, :, which_stim), [3 2 1]);
rsa_vl_U = nan(4,4,numel(plot_auditory_roi));
rsa_vl_pval = nan(4,4,numel(plot_auditory_roi));
figure('color','w','position',[100 100 900 200]);
for iRoi = 1:numel(plot_auditory_roi)
    subplot(1,numel(plot_auditory_roi),iRoi); hold on
    iChan = contains(voice_rois_vl, plot_auditory_roi{iRoi});
    D_neural_plot = squareform(pdist(reshape(zscore(rsa_data_vl(:,iChan,:),[],3),...
        sum(which_stim), []),metric));
    for iQ = 1:4
        imask = cat_clr(is_nat) == iQ;
        within = reshape(triu(D_neural_plot(imask, imask)), [] ,1);
        for jQ = 1:4
            jmask = cat_clr(is_nat) == jQ;
            between = reshape(D_neural_plot(imask, jmask), [], 1);
                [rsa_vl_pval(iQ,jQ,iRoi), ~, stats] = ranksum(within, between);
                n1 = numel(within);
                n2 = numel(between);
                U = stats.ranksum - (n1*(n1+1))/2;  % convert ranksum to Mann–Whitney U
                rsa_vl_U(iQ, jQ,iRoi) = (2 * U) / (n1 * n2) - 1;
        end
    end
    imagesc(rsa_vl_U(:,:,iRoi));title(plot_auditory_roi{iRoi})
    colormap(flipud(cbrewer2('Reds')));colorbar()
    axis square xy tight; set(gca,'fontsize',12)
end

% figure('color','w'); boxplot(rho_ac_vl, erase(erase(voice_rois_kp,'R_'),'L_'),...
%     'boxstyle','filled','grouporder',plot_auditory_roi,'colors','k','symbol','.k')
% title('natural - acoustics'); ylim([-0.15 0.25]); set(gca,'fontsize',12,'fontname','arial')
% figure('color','w'); boxplot(rho_vl, erase(erase(voice_rois_kp,'R_'),'L_'),...
%     'boxstyle','filled','grouporder',plot_auditory_roi,'colors','k','symbol','.k')
% title('natural - ratings'); ylim([-0.15 0.25]); set(gca,'fontsize',12,'fontname','arial')
% figure('color','w'); boxplot(rho_pcorr_vl, erase(erase(voice_rois_kp,'R_'),'L_'),...
%     'boxstyle','filled','grouporder',plot_auditory_roi,'colors','k','symbol','.k')
% title('natural - ratings(acoustics)'); set(gca,'fontsize',12,'fontname','arial')
%
%
% tbl = table(ld_vl.pt_matrix, ld_vl.chan_roi_matrix(:, 1),...
%     'VariableNames',{'pt','chans'});
% plotChansSurf(tbl.pt(idx_voice_chans), tbl.chans(idx_voice_chans), 'color', ...
%     rho_vl,'map','amp','clim',[0 0.11],'size', p_vl < 0.05, 'szmap',[0.5 1.5])
% title('natural - ratings');  set(gca,'fontsize',12,'fontname','arial')
% plotChansSurf(tbl.pt(idx_voice_chans), tbl.chans(idx_voice_chans), 'color', ...
%     rho_ac_vl,'map','amp','clim',[0 0.11],'size', p_ac_vl < 0.05, 'szmap',[0.5 1.5])
% title('natural - acoustics'); set(gca,'fontsize',12,'fontname','arial')
% plotChansSurf(tbl.pt(idx_voice_chans), tbl.chans(idx_voice_chans), 'color', ...
%     rho_pcorr_vl,'map','amp','size',p_pcorr_vl < 0.05, 'szmap',[0.5 1.3])
% title('natural - ratings(acoustics)'); set(gca,'fontsize',12,'fontname','arial')

figure('color','w'); hold on
for iRoi = 1:numel(plot_auditory_roi)
    plot(t,rsa_ratings_vl.rho_iT(iRoi,:),'linewidth',3)
end
legend(plot_auditory_roi)
set(gca,'fontsize',12)
ylabel('Spearman rho');xlabel('time (ms)')
title('natural - ratings')

figure('color','w'); hold on
for iRoi = 1:numel(plot_auditory_roi)
    plot(t,rsa_ratings_vl.rho_pcorr_iT(iRoi,:),'linewidth',3)
end
legend(plot_auditory_roi)
set(gca,'fontsize',12)
ylabel('Spearman rho');xlabel('time (ms)')
title('natural - ratings(acoustics)')

figure('color','w'); hold on
for iRoi = 1:numel(plot_auditory_roi)
    plot(t,rsa_ratings_vl.rho_ac_iT(iRoi ,:),'linewidth',3)
end
legend(plot_auditory_roi)
set(gca,'fontsize',12)
ylabel('Spearman rho');xlabel('time (ms)')
title('natural - acoustics')


% --------------------------synthetic sounds ------------------------------
if ~exist('rsa_ratings_gsp','var')
    kpt_acfts = gsp.t>=  t0 & gsp.t<= t1;
    t = gsp.t(kpt_acfts);
    kpT   = find(ld_gsp.tds == t0):find(ld_gsp.tds == t1);
    selT   = 1: stepSize*frameSize/1000 : find(ld_gsp.tds == t(end));
    which_stim = cat_clr(~is_nat) > rate_cutoff;

    D_rating_gsp = squareform(pdist(cat_rate(~is_nat & cat_clr > rate_cutoff), 'euclidean'));

    corr_coeffs = cat(2,gsp.tcoeff(which_stim,1),gsp.fcoeff(which_stim,1));
    corr_coeffs(abs(zscore(corr_coeffs(:,1),[],'omitnan'))>6,1) = nan;
    corr_coeffs(abs(zscore(corr_coeffs(:,2),[],'omitnan'))>6,2) = nan;
    D_acoustics_gsp = squareform(pdist(cat(2,zscore(corr_coeffs,[],'omitnan'),...
        permute(zscore(mean(gsp.acfts(kpt_acfts, kp_fts,...
        which_stim),'omitnan'),[],3),[3,2,1])),'correlation')); % use correlation! zscore across stimuli after averaging across time.

    %     D_acoustics_gsp = squareform(pdist(reshape(permute(zscore(gsp.acfts(kpt_acfts, kp_fts,which_stim),...
    %         [],'omitnan'),[3,2,1]),sum(which_stim),[]),'correlation'));
    rsa_ac_gsp = zscore(gsp.acfts(kpt_acfts, kp_fts,which_stim),[],3); %zscore across stimuli
    rsa_data_gsp = permute(ld_gsp.data_matrix_ds(:, : , which_stim), [3 2 1]);
    for iRoi = 1:numel(plot_auditory_roi)
        iChan = contains(voice_rois_gsp, plot_auditory_roi{iRoi});
        D_neural_gsp(:,:,iRoi) = squareform(pdist(reshape(rsa_data_gsp(:, iChan, kpT),sum(which_stim),[]),metric));
        [rho_gsp(iRoi), p_gsp(iRoi)] = mantel_spearman(D_neural_gsp(:,:,iRoi), D_rating_gsp, nPerm, 1);
        for iT = 1:numel(t)
            nr_win_idx = find(ld_gsp.tds==t(iT)) : (find(ld_gsp.tds==t(iT)) + nr_samps_per_frame/4);
            D_neural_iT= squareform(pdist(reshape(rsa_data_gsp(:, iChan, nr_win_idx), sum(which_stim), []), metric));
            [rho_iT_gsp(iRoi,iT), p_iT_gsp(iRoi,iT)] = mantel_spearman(D_neural_iT, D_rating_gsp, nPerm, 1);
        end
        [rho_ac_gsp(iRoi), p_ac_gsp(iRoi)] = mantel_spearman(D_neural_gsp(:,:,iRoi), D_acoustics_gsp, nPerm, 1);
        for iT = 1:numel(t)
            nr_win_idx = find(ld_gsp.tds==t(iT)) : (find(ld_gsp.tds==t(iT)) + nr_samps_per_frame/4);
            D_neural_iT= squareform(pdist(reshape(rsa_data_gsp(:, iChan, nr_win_idx), sum(which_stim), []), metric));
            D_acoustics_gsp_iT = squareform(pdist(squeeze(rsa_ac_gsp(iT,:,:))','correlation'));
            [rho_ac_iT_gsp(iRoi,iT), p_ac_iT_gsp(iRoi,iT)] = mantel_spearman(D_neural_iT, D_acoustics_gsp_iT, nPerm, 1);
        end
        [rho_pcorr_gsp(iRoi), p_pcorr_gsp(iRoi)] = mantel_spearman(D_neural_gsp(:,:,iRoi), ...
            D_rating_gsp, nPerm, 1,'pcorr', D_acoustics_gsp);
        for iT = 1:numel(t)
            nr_win_idx = find(ld_gsp.tds==t(iT)) : (find(ld_gsp.tds==t(iT)) + nr_samps_per_frame/4);
            D_acoustics_gsp_iT = squareform(pdist(squeeze(rsa_ac_gsp(iT,:,:))','correlation'));
            D_neural_iT= squareform(pdist(reshape(rsa_data_gsp(:, iChan, nr_win_idx), sum(which_stim), []), metric));
            [rho_pcorr_iT_gsp(iRoi,iT), p_pcorr_iT_gsp(iRoi,iT)] = mantel_spearman(D_neural_iT, ...
                D_rating_gsp, nPerm, 1, 'pcorr', D_acoustics_gsp_iT);
        end
    end
    rsa_ratings_gsp = struct('rho',rho_gsp,'p',p_gsp, 'p_iT',p_iT_gsp,'rho_iT',rho_iT_gsp,...
        'rho_ac',rho_ac_gsp,'p_ac',p_ac_gsp, 'p_ac_iT',p_ac_iT_gsp,'rho_ac_iT',rho_ac_iT_gsp,...
        'rho_pcorr',rho_pcorr_gsp,'p_pcorr',p_pcorr_gsp, 'p_pcorr_iT',p_pcorr_iT_gsp,'rho_pcorr_iT',...
        rho_pcorr_iT_gsp,'nPerm',nPerm);
        save(fullfile(gcs.fOneDrive, '..\gsp\save\',  sprintf('all_out_webb_%s_%s_%s_process_sEEG_nbs_%i_new.mat',...
            which_data, which_filter, 'gsp', nBs)),'rsa_ratings_gsp','-append');
else
    load(fullfile(gcs.fOneDrive, '..\gsp\save\',  sprintf('all_out_webb_%s_%s_%s_process_sEEG_nbs_%i_new.mat',...
        which_data, which_filter, 'gsp', nBs)), 'rsa_ratings_gsp')
end
kpT   = find(ld_gsp.tds == t0):find(ld_gsp.tds == t1);
which_stim = cat_clr(~is_nat) > rate_cutoff;
rsa_gsp_U = nan(4,4,numel(plot_auditory_roi));
rsa_gsp_pval = nan(4,4,numel(plot_auditory_roi));
rsa_data_gsp = permute(ld_gsp.data_matrix_ds(kpT, :, which_stim), [3 2 1]);
figure('color','w','position',[100 100 900 200]);
for iRoi = 1:numel(plot_auditory_roi)
    subplot(1,numel(plot_auditory_roi),iRoi); hold on
    iChan = contains(voice_rois_gsp, plot_auditory_roi{iRoi});
    D_neural_plot = squareform(pdist(reshape(rsa_data_gsp(:,iChan,:),...
        sum(which_stim), []),metric));
    for iQ = 1:4
        imask = cat_clr(~is_nat) == iQ;
        within = reshape(triu(D_neural_plot(imask, imask)), [] ,1);
        for jQ = 1:4
            jmask = cat_clr(~is_nat) == jQ;
            between = reshape(D_neural_plot(imask, jmask), [], 1);
            [rsa_gsp_pval(iQ,jQ,iRoi), ~, stats] = ranksum(within, between);
            n1 = numel(within);
            n2 = numel(between);
            U = stats.ranksum - (n1*(n1+1))/2;  % convert ranksum to Mann–Whitney U
            rsa_gsp_U(iQ, jQ,iRoi) = (2 * U) / (n1 * n2) - 1;
        end
    end
    imagesc(rsa_gsp_U(:,:,iRoi));title(plot_auditory_roi{iRoi})
    colormap(flipud(cbrewer2('Reds')));colorbar()
    axis square xy tight; set(gca,'fontsize',12)
end

figure('color','w'); hold on
for iRoi = 1:numel(plot_auditory_roi)
    plot(t,rsa_ratings_gsp.rho_iT(iRoi,:),'linewidth',3)
end
legend(plot_auditory_roi)
set(gca,'fontsize',12)
ylabel('Spearman rho');xlabel('time (ms)')
title('Synthetic - ratings')

figure('color','w'); hold on
for iRoi = 1:numel(plot_auditory_roi)
    plot(t,rsa_ratings_gsp.rho_pcorr_iT(iRoi,:),'linewidth',3)
end
legend(plot_auditory_roi)
set(gca,'fontsize',12)
ylabel('Spearman rho');xlabel('time (ms)')
title('Synthetic - ratings(acoustics)')

figure('color','w'); hold on
for iRoi = 1:numel(plot_auditory_roi)
    plot(t,rsa_ratings_gsp.rho_ac_iT(iRoi ,:),'linewidth',3)
end
legend(plot_auditory_roi)
set(gca,'fontsize',12)
ylabel('Spearman rho');xlabel('time (ms)')
title('Synthetic - acoustics')

%% across both tasks

rsa_data_vl_plot = rsa_data_vl(idx_voice_chans,:,:);
rsa_vl_U = nan(8,8,numel(plot_auditory_roi));
rsa_vl_tpval = nan(8,8,numel(plot_auditory_roi));
figure('color','w','position',[100 100 900 200]);
for iRoi = 1:numel(plot_auditory_roi)
    subplot(1,numel(plot_auditory_roi),iRoi); hold on
    iChan = contains(find(contains(voice_rois_vl,plot_auditory_roi)), plot_auditory_roi{iRoi});
    D_neural_plot_vl = squareform(1 - abs(pdist(reshape(permute(rsa_data_vl_plot(voice_rois_vl(contains(voice_rois_vl,plot_auditory_roi)),:,:),...
    [2,3,1]), sum(cat_clr(is_nat) > rate_cutoff), []),'correlation') + 1)); 

    iChan = contains(find(contains(voice_rois_gsp, plot_auditory_roi)), plot_auditory_roi{iRoi});
    D_neural_plot_vl = squareform(1 - abs(pdist(reshape(permute(rsa_data_gsp_plot(voice_rois_gsp(contains(voice_rois_gsp,plot_auditory_roi)),:,:),...
    [2,3,1]), sum(cat_clr(~is_nat) > rate_cutoff), []),'correlation') + 1)); 
    for iQ = 1:4
        imask_vl = triu(true(size(D_neural_plot)), 1) & (cat_clr(is_nat) == iQ);
        for jQ = iQ:4
            jmask = triu(true(size(D_neural_plot)), 1) & (cat_clr(is_nat) == jQ);
            [~, rsa_vl_tpval(iQ, jQ, iRoi), ~, stat] = ttest2(reshape(D_neural_plot(imask_vl), [], 1),...
                reshape(D_neural_plot(jmask), [], 1),'tail','both');
            rsa_vl_U(iQ,jQ,iRoi) = stat.tstat;
        end
    end
    imagesc(rsa_vl_U(:,:,iRoi));title(plot_auditory_roi{iRoi})
    clim([max(abs(rsa_vl_U(:,:,iRoi)),[],'all')*-1 max(abs(rsa_vl_U(:,:,iRoi)),[],'all')]); 
    colormap(flipud(cbrewer2('RdBu')));colorbar()
    axis square xy tight; set(gca,'fontsize',12)
end
%% ------------------------- same channels only ------------------------
kp_gsp_chans = contains(ld_gsp.pt_chan, ld_vl.pt_chan)'; kp_vl_chans = contains(ld_vl.pt_chan, ld_gsp.pt_chan)';
rsa_data_cat = cat(2, rsa_data_vl(kp_vl_chans,:,:), rsa_data_gsp(kp_gsp_chans,:,:));
idx_voice_chans = find(contains(voice_rois_vl(kp_vl_chans),plot_auditory_roi));
rho_iT = nan(numel(idx_voice_chans), numel(selT));
p_iT = nan(numel(idx_voice_chans), numel(selT));
for iChan = 1:numel(idx_voice_chans)
    fprintf('%g/%g\n',iChan,numel(idx_voice_chans))
    D_neural_gsp= squareform(pdist(squeeze(rsa_data_cat(idx_voice_chans(iChan), :, :)),'euclidean'));
    [rho(iChan), rsa_vl_pval(iChan)] = mantel_spearman(D_neural_gsp, D_rating_vl, nPerm, 1);
    % Significant positive rho ⇒ larger |Δrating| ↔ larger neural distance
    for iT = 1:numel(selT)
        if iT + 9 < numel(selT)
            D_neural_iT= squareform(pdist(squeeze(rsa_data_cat(idx_voice_chans(iChan), :, iT:iT+9)),'euclidean'));
            [rho_iT(iChan,iT), p_iT(iChan,iT)] = mantel_spearman(D_neural_iT, D_rating_vl, nPerm, 1);
        end
    end
end

rsa_cat_tstat = nan(8,8,numel(plot_auditory_roi));
rsa_cat_tpval = nan(8,8,numel(plot_auditory_roi));
figure('color','w','position',[100 100 900 200]);
for iRoi = 1:numel(plot_auditory_roi)
    subplot(1,numel(plot_auditory_roi),iRoi); hold on
    iChan = contains(voice_rois_gsp, plot_auditory_roi(iRoi));
    D_neural_plot = reshape(permute(zscore(rsa_data_gsp(iChan,:,:),0,3),...
    [2,3,1]), sum(which_stim),[]); 
    for iQ = 1:4
        mask = triu(true(size(D_neural_plot)),1) & (cat_clr(~is_nat)==iQ);
        for jQ = iQ:4
            jmask = triu(true(size(D_neural_plot)),1) & (cat_clr(~is_nat)==jQ);
            [~,rsa_gsp_tpval(iQ,jQ,iRoi),~,stat] = ttest2(reshape(D_neural_plot(mask),[],1),...
                reshape(D_neural_plot(jmask),[],1),'tail','both');
            rsa_cat_tstat(iQ,jQ,iRoi) = stat.tstat;
        end
    end
    imagesc(rsa_gsp_tstat(:,:,iRoi));title(plot_auditory_roi{iRoi})
    clim([min(rsa_cat_tstat(:,:,iRoi),[],'all'), min(rsa_cat_tstat(:,:,iRoi),[],'all')*-1]); 
    colormap(flipud(cbrewer2('RdBu')))
    axis square xy tight;
end

%% load stimuli and make cochleograms
% load(fullfile(path.save,'P19N001_test.mat'),'processed_data');
% tHga = processed_data.t; clear processed_data
% tSample = 0:100:2000;
% gcs = loadSettings();
% [wav,fs] = arrayfun(@(x) audioread(fullfile(x.folder,x.name)),dir(fullfile(gcs.fTasks,'naturalSounds/stimuli_48000Hz/*.wav')),'UniformOutput',false);
% wav_names = dir(fullfile(gcs.fTasks,'naturalSounds/stimuli_48000Hz/*.wav'));
% nStim = length(wav);
% paras = [10 8 .1 0];
% global COCHBA;load('aud24')
% for ii = 1:nStim
%     thisWavDs = resample(wav{ii},16e3,fs{ii});
%     coch = log10(wav2aud(thisWavDs,paras) + 1e-3);
%     cochs(:,:,ii) = resample(coch,32,128,'Dimension',2);
% end
% clear coch wav thisWavDs
% cochs = zscore(cochs,0,3);
% z_cochs_vec = reshape(cochs,[],size(cochs,3)); clear cochs


function [rho,p,rho_perm] = mantel_spearman(D1,D2,nPerm,rngSeed,varargin)
if nargin<3||isempty(nPerm), nPerm=5000; end
if nargin>3&&~isempty(rngSeed), rng(rngSeed); end

mask = triu(true(size(D1)),1);
v1 = D1(mask); v2 = D2(mask);
if any(strcmp(varargin,'pcorr'))
    D3 = varargin{find(strcmp(varargin,'pcorr'))+1};
    v3 = D3(mask);
    rho = partialcorr(v1, v2, v3,'type','Spearman','rows','complete');
    n = size(D1,1); idx = 1:n; rho_perm = zeros(nPerm,1);
    for i=1:nPerm
        perm = idx(randperm(n));
        vp = D1(perm,perm); vp = vp(mask);
        rho_perm(i) = corr(vp, v2, 'type','Spearman','rows','complete');
    end
    p = (1 + sum(rho_perm >= rho)) / (nPerm + 1); % one-sided: expect positive

else
    rho = corr(v1, v2, 'type','Spearman','rows','complete');
    n = size(D1,1); idx = 1:n; rho_perm = zeros(nPerm,1);
    for i=1:nPerm
        perm = idx(randperm(n));
        vp = D1(perm,perm); vp = vp(mask);
        rho_perm(i) = corr(vp, v2, 'type','Spearman','rows','complete');
    end
    p = (1 + sum(rho_perm >= rho)) / (nPerm + 1); % one-sided: expect positive
end
end


function s = p2stars(p)
if p < 0.001
    s = '***';
elseif p < 0.01
    s = '**';
elseif p < 0.05
    s = '*';
else
    s = '';
end
end

function makeLMEBar(models, xlabels, beta, lo, hi, pvals, ylab, y0)
for imdl = 1:numel(models)
    figure('Color','w','Position',[290*imdl,200,272,254]);
    x = 1:numel(xlabels);
    % Bar
    b = bar(x, beta(:,imdl),'FaceColor',[0.30 0.55 0.95], ...
        'EdgeColor','none');
    ttl = sprintf('%s on %s', models{imdl}, y0);

    hold on;
    % Error bars from CI
    yerr = [beta(:,imdl)-lo(:,imdl), hi(:,imdl)-beta(:,imdl)];
    errorbar(x, beta(:,imdl), yerr(:,1), yerr(:,2), 'k', ...
        'LineStyle','none', 'CapSize',6, 'LineWidth',1);
    % Zero line
    yline(0,'k-','LineWidth',1);
    xlim([0.5 numel(xlabels)+0.5]);
    set(gca,'XTick',x,'XTickLabel',erase(xlabels,'_'),...
        'FontSize',8,'LineWidth',1);
    xlabel('ROI');
    ylabel(ylab);
    title(ttl);
    box on;
    ylim([min(beta(:,imdl)-yerr(:,1))-max(yerr,[],'all'),...
        max(beta(:,imdl)+yerr(:,2))+max(yerr,[],'all')])

    % Sig stars
    for i = 1:numel(xlabels)
        s = p2stars(pvals(i, imdl));
        if ~isempty(s)
            y = beta(i,imdl);
            % lift star a bit above the error bar/extreme of bar
            lift = max(yerr(i,2), yerr(i,1));
            if y >= 0
                ystar = y + lift;
                va = 'bottom';
            else
                ystar = y - lift;
                va = 'top';
            end
            text(i, ystar, s, 'HorizontalAlignment','center', ...
                'VerticalAlignment',va, 'FontSize',14, 'FontWeight','bold');
        end
    end
    hold off;


    % Save
    %     exportgraphics(gcf, [savebase '.png'], 'Resolution', 300);
    %     exportgraphics(gcf, [savebase '.svg']);
end

end


function makeBar(X, xlabels, pvals, bar_group, ylab, ttl)
figure('color','w','position', [680,627,1300,351]);
x = 1:numel(xlabels);
if numel(bar_group)>1
    bar(X,'grouped','EdgeColor','none');
else
    bar(X,'EdgeColor','none','BarWidth',0.25);
end
legend(bar_group,'location','best')
hold on;
xticks(x)
xticklabels(xlabels)
xlim([0 numel(xlabels)+1]);
set(gca,'XTick',x,'XTickLabel',erase(xlabels,'_'),...
    'FontSize',8,'LineWidth',1);
ylabel(ylab);
box on;

% Sig stars
for i = 1:numel(xlabels)
    if numel(bar_group) > 1
        group_star_x = linspace(-0.2, 0.2, width(pvals));
        for iStar = 1:width(pvals)
            s = p2stars(pvals(i, iStar));
            if ~isempty(s)
                y = X(i,iStar);
                % lift star a bit above the bar
                lift = 0.002;
                if y >= 0
                    ystar = y + lift;
                    va = 'bottom';
                else
                    ystar = y - lift;
                    va = 'top';
                end
                text(i + group_star_x(iStar), ystar, s, 'HorizontalAlignment','center', ...
                    'VerticalAlignment',va, 'FontSize',10, 'FontWeight','bold');
            end
        end
    else
        s = p2stars(pvals(i,1));
        if ~isempty(s)
            y = X(i);
            % lift star a bit above the bar
            lift = 0.002;
            if y >= 0
                ystar = y + lift;
                va = 'bottom';
            else
                ystar = y - lift;
                va = 'top';
            end
            text(i, ystar, s, 'HorizontalAlignment','center', ...
                'VerticalAlignment',va, 'FontSize',10, 'FontWeight','bold');
        end
    end
end
ylim([min(X - 0.1,[],'all'), max(X + 0.1,[],'all')])
title( sprintf('%s', ttl));
end

function r2 = r2_from_preds(y, yhat)
y = y(:); yhat = yhat(:);
m = ~isnan(y) & ~isnan(yhat);
y = y(m); yhat = yhat(m);
if numel(y) < 2
    r2 = NaN; return
end
sse = sum((y - yhat).^2);
sst = sum((y - mean(y)).^2);
if sst == 0, r2 = NaN; else, r2 = 1 - sse/sst; end
end
