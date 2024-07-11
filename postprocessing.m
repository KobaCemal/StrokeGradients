clear
warning('off')
% Set global variables and parameters
derivatives_dir='/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives';
% cd(derivatives_dir)
subjects=readtable('/home/koba/Desktop/Stroke/scripts/subjects.csv');
behavioral=readtable('/home/koba/Desktop/Stroke/scripts/behavioral.csv');
parcellation_nii=load_nii('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/source/yeo400_resampled.nii.gz');
parcellation=reshape(parcellation_nii.img,[],1);
parcellation_mask_whole=parcellation>0;
parcellation_mask_left=parcellation>0 & parcellation <201; % first half is left hemisphere
parcellation_mask_right=parcellation>200;
parcellation_reduced=parcellation(parcellation_mask_whole>0);
[surf_lh, surf_rh] = load_conte69();
labeling = load_parcellation('schaefer',400);
[network_names, colormap] = fetch_yeo_networks_metadata(7);
schaefer_400 = fetch_parcellation('fsaverage5', 'schaefer', 400);

TR = 2;
samplingRate = 1 / TR;
lowFreq = 0.01; % Lower cutoff frequency in Hz
highFreq = 0.1; % Upper cutoff frequency in Hz
filterOrder = 5; % Filter order
[b, a] = butter(filterOrder, [lowFreq, highFreq] * 2 / samplingRate);
lags=-3:3; % TR, not seconds
num_control=sum(subjects.subj_type);
num_stroke=size(subjects,1)-num_control;

% gather the necessary files
func_list=char(table2array(readtable('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/source/func_files.txt',Delimiter=',',ReadVariableNames=false)));
mask_list=char(table2array(readtable('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/source/mask_files.txt',Delimiter=',',ReadVariableNames=false)));
confounds_list=char(table2array(readtable('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/source/confounds.txt',Delimiter=',',ReadVariableNames=false)));
seven_network_positions=table2array(readtable('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/source/7netpositions.txt',ReadVariableNames=false));
network_names={'Visual','Somatomotor','Dorsal Att.','Ventral Att.','Limbic','Control','Default Mode'};

load('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/corrmats_temp.mat')
% Loop starts here  -  no need to run if corrmats_temp.mat is loaded

%% Postprocessing: confound regression, band-pass filtering, lag calculation, saving niftis, lag correction, creating correlation matrices
corrmats=zeros(400,400,size(func_list,1));
corrmats_corrected=zeros(400,400,size(func_list,1));
corrmats_corrected_neg=zeros(400,400,size(func_list,1));

for i=1:size(func_list,1)
    disp(i)
    % load the files

    func_file=load_nii(strrep(func_list(i,:),' ', ''));
    confound_file=readtable(strrep(confounds_list(i,:),' ', ''),"FileType","text");
    %     mask_file=load_nii(strrep(mask_list(i,:),' ', ''));

    % mask the functional data and keep the confounds of interest
    func_data=reshape(func_file.img,[],size(func_file.img,4));
    %mask_data=reshape(mask_file.img,[],1);
    func_masked=func_data(logical(parcellation_mask_whole),:); % func data is converted to 2D and masked
    confounds_of_interest = {'framewise_displacement'  'csf' 'white_matter' 'trans_x' 'trans_y' 'trans_z' 'rot_x' 'rot_y' 'rot_z' ...
        'csf_derivative1' 'white_matter_derivative1' 'trans_x_derivative1' 'trans_y_derivative1' 'trans_z_derivative1' 'rot_x_derivative1' 'rot_y_derivative1' 'rot_z_derivative1' ...
        'csf_derivative1_power2' 'white_matter_derivative1_power2' 'trans_x_derivative1_power2' 'trans_y_derivative1_power2' 'trans_z_derivative1_power2' 'rot_x_derivative1_power2' 'rot_y_derivative1_power2' 'rot_z_derivative1_power2' ...
        'csf_power2' 'white_matter_power2' 'trans_x_power2' 'trans_y_power2' 'trans_z_power2' 'rot_x_power2' 'rot_y_power2' 'rot_z_power2'};

    confounds_data=[table2array(confound_file(:,confounds_of_interest)) (1:size(func_data,2))']; % adding linear trend
    confounds_data(isnan(confounds_data))=0; % change nan with zero

    % Regress te confounds our and apply band-pass filter
    func_clean=zeros(size(func_masked));
    for j=1:size(func_masked,1)
        [~,~,residuals] = regress(func_masked(j,:)',[ones(size(func_masked,2),1) confounds_data]);
        func_clean(j,:) = filter(b, a, residuals);
    end

    % Get the GS from GM (whole-brain if control, healthy hemisphere if stroke)
    sub_id=func_list(i,10:19);
    if strfind(func_list(i,:),'_space')-(strfind(func_list(i,:),'n-')+2)==1
        run=strfind(func_list(i,:),'n-')+2;
    else
        run=[strfind(func_list(i,:),'n-')+2 strfind(func_list(i,:),'n-')+3];
    end
    position=find(strcmp(sub_id,subjects.participant_id)==1);

    if subjects.subj_type(position)==1 % 1 is for control
        parcellation_in_brainmask=parcellation_mask_whole;
        gs=mean(func_data(parcellation_in_brainmask,:));
    elseif subjects.subj_type(position)==0
        if subjects.lesion_side(position) == 0 % 0 is for left
            parcellation_in_brainmask=parcellation_mask_right;
            gs=mean(func_data(parcellation_in_brainmask,:));
        elseif subjects.lesion_side(position) == 1
            parcellation_in_brainmask=parcellation_mask_left;
            gs=mean(func_data(parcellation_in_brainmask,:));
        end
    end

    % Calculate the lag
    corrs_lag=zeros(length(lags),size(func_clean,1));
    for j=1:length(lags)
        gsr_lag=circshift(gs,lags(j));
        corrs_lag(j,:)=corr(gsr_lag',func_clean');
    end
    corrs_lag(isnan(corrs_lag))=0;
    lag_pos=zeros(1,size(func_clean,1));
    for j=1:length(corrs_lag)
        vec=(corrs_lag(:,j));
        laginfo=(find(vec==max(vec)));
        if isempty(laginfo)
            lag_pos(j)=0;
        else
            lag_pos(j)=lags(laginfo(1));
        end
    end
    lags_all(i,:)=lag_pos;

    template=parcellation_nii;
    template_data=zeros(size(parcellation));
    template_data(logical(parcellation_mask_whole))=lag_pos*-1;
    template.img=reshape(template_data,size(parcellation_nii.img));
    delay_name='delays_custom/SUBID_run-RUN_delay_gm.nii.gz';
    delay_name=strrep(delay_name,'SUBID',sub_id);
    delay_name=strrep(delay_name,'RUN',func_list(i,run));
    save_nii(template,char(delay_name))

    % Fix the time series
    func_corrected=zeros(size(func_clean));
    for j=1:size(func_clean,1)
        func_corrected(j,:)=circshift(func_clean(j,:),lag_pos(j)*-1);
    end

    %     func_corrected_neg=zeros(size(func_clean));
    %     for j=1:size(func_clean,1)
    %         func_corrected_neg(j,:)=circshift(func_clean(j,:),lag_pos(j));
    %     end

    % Get the correlation matrices
    %     parcellation_in_brainmask=parcellation_mask_whole;
    mean_ts=zeros(max(parcellation_reduced),size(func_corrected,2));
    mean_ts_corrected=zeros(max(parcellation_reduced),size(func_corrected,2));
    %     mean_ts_corrected_neg=zeros(max(parcellation_reduced),size(func_corrected_neg,2));

    for j=1:max(parcellation_reduced)
        m=func_clean(parcellation_reduced==j,:);
        mean_ts(j,:)=mean(m(any(m,2),:));

        m=func_corrected(parcellation_reduced==j,:);
        mean_ts_corrected(j,:)=mean(m(any(m,2),:));

        %         m=func_corrected_neg(parcellation_reduced==j,:);
        %         mean_ts_corrected_neg(j,:)=mean(m(any(m,2),:));
    end

    corrmats(:,:,i)=corr(mean_ts');
    corrmats_corrected(:,:,i)=corr(mean_ts_corrected');
    %     corrmats_corrected_neg(:,:,i)=corr(mean_ts_corrected_neg');
end

%% Processing the matrices: taking the average for each subject and Z transformation
positions=zeros(size(corrmats,3),1);
corrmats_reduced=zeros(max(parcellation),max(parcellation),size(subjects,1));
corrmats_corrected_reduced=zeros(max(parcellation),max(parcellation),size(subjects,1));
for i=1:size(subjects,1)
    sub_id=subjects.participant_id(i);
    for j=1:size(corrmats,3)
        positions(j)=contains(func_list(j,:),sub_id);
    end
    total_runs(i,1)=sum(positions);
    corrmats_reduced(:,:,i)=atanh(mean(corrmats(:,:,logical(positions)),3));
    corrmats_corrected_reduced(:,:,i)=atanh(mean(corrmats_corrected(:,:,logical(positions)),3));
end
corrmats_reduced(isinf(corrmats_reduced))=1;
corrmats_corrected_reduced(isinf(corrmats_corrected_reduced))=1;


imagesc(mean(corrmats_reduced(:,:,1:num_control),3))
figure
imagesc(mean(corrmats_reduced_clean(:,:,1:num_control),3))

figure
imagesc(mean(corrmats_reduced(:,:,num_control+1:end),3))
figure
imagesc(mean(corrmats_corrected_reduced(:,:,1:num_control),3))
figure
imagesc(mean(corrmats_corrected_reduced(:,:,num_control+1:end),3))

%% Diagnostics of lags: mean - variance values and aberrancy calculation
% Variance and max values of lags

lags=load_nii('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/delays_custom/delays.nii.gz');
lags=reshape(squeeze((lags.img)),[],875);
lags=lags(parcellation_mask_whole,:);

% take the average of each subject
positions=zeros(size(corrmats,3),1);

% lags
for i=1:size(subjects,1)
    sub_id=subjects.participant_id(i);
    for j=1:size(corrmats,3)
        positions(j)=contains(func_list(j,:),sub_id);
    end
    lags_reduced(:,i)=mean(lags(:,logical(positions)),2);
end

% Mean lag 400 regions
for j=1:max(parcellation_reduced)
    for i=1:size(subjects,1)
        mean_lag(j,i)=mean(lags_reduced(parcellation_reduced==j,i));
    end
end
[h,p,c,stats]=ttest2(mean(mean_lag(:,1:num_control),2),mean(mean_lag(:,num_control+1:end),2))
[h,p,c,stats]=ttest2(var(mean_lag(:,1:num_control),[],2),var(mean_lag(:,num_control+1:end),[],2))
[h,p,c,stats]=ttest2(max((mean_lag(:,1:num_control)),[],2),max((mean_lag(:,num_control+1:end)),[],2))
[h,p,c,stats]=ttest2(min(abs(mean_lag(:,1:num_control)),[],2),min(abs(mean_lag(:,num_control+1:end)),[],2))
hold on;histogram(min((mean_lag(:,1:num_control)),[],2));histogram(min((mean_lag(:,num_control+1:end)),[],2));hold off
hold on;histogram(max((mean_lag(:,1:num_control)),[],2));histogram(max((mean_lag(:,num_control+1:end)),[],2));hold off


for i=1:400
    [h,p,c,stats]=ttest2(mean_lag(i,1:num_control),mean_lag(i,num_control+1:end));
    lag_ps(i)=p;
    lag_ts(i)=stats.tstat;
end
plot_hemispheres([mean(mean_lag(:,1:num_control),2) mean(mean_lag(:,num_control+1:end),2)], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'controls','Mean Lag-stroke'});
lagc=mean(mean_lag(:,1:num_control),2);
lags=mean(mean_lag(:,num_control+1:end),2);
plot_hemispheres([lag_ts'], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'lag diff roi',});


% Aberrancy
mean_control=mean(corrmats_reduced(:,:,1:num_control),3);
std_control=std(corrmats_corrected_reduced(:,:,1:num_control),[],3);

for i=1:size(subjects,1)
    abr(:,:,i)=(corrmats_reduced(:,:,i)-mean_control)./std_control;
    x=reshape(tril(abr(:,:,i),-1),[],1);
    mean_abr(i)=mean(x(x~=0));

    abr_cr(:,:,i)=(corrmats_corrected_reduced(:,:,i)-mean_control)./std_control;
    x=reshape(tril(abr_cr(:,:,i),-1),[],1);
    mean_abr_cr(i)=mean(x(x~=0));
end


mean_abr_control=mean_abr(1:num_control);
mean_abr_stroke=mean_abr(num_control+1:end);

mean_abr_control_corrected=mean_abr_cr(1:num_control);
mean_abr_stroke_corrected=mean_abr_cr(num_control+1:end);
hold on;histogram(mean_abr_stroke,'Normalization','probability');histogram(mean_abr_control,'Normalization','probability');hold off
figure
hold on;histogram(mean_abr_stroke_corrected,'Normalization','probability');histogram(mean_abr_control_corrected,'Normalization','probability');hold off

[h,p,c,stats]=ttest2(mean_abr_control,mean_abr_stroke)
[h,p,c,stats]=ttest2(mean_abr_control_corrected,mean_abr_stroke_corrected)

%% Calculating the gradients: scaling them and retrieving explained variances



% Get the gradients for each subject
[surf_lh, surf_rh] = load_conte69();
labeling = load_parcellation('schaefer',400);
labelin2=labeling.schaefer_400;
conn_matrices = load_group_fc('schaefer',400);
reference_gradient = GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
m = (conn_matrices.schaefer_400);
m=atanh(m);
m(isinf(m))=1;
reference_gradient = reference_gradient.fit(m);

for i=1:size(corrmats_reduced,3)
    disp(i)
    m=corrmats_reduced(:,:,i);
    %     m=atanh(m);
    %     m(isinf(m))=1;
    m(isnan(m))=1;
    before_correction=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
    before_correction=before_correction.fit(m,'reference',reference_gradient.gradients{1},'sparsity',90);
    gradients_before{i}=before_correction;
    m=corrmats_corrected_reduced(:,:,i);
    m(isnan(m))=1;
    %     m=atanh(m);
    %      m(isinf(m))=1;
    after_correction=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
    after_correction=before_correction.fit(m,'reference',reference_gradient.gradients{1});
    gradients_after{i}=after_correction;

    before_aligned(:,:,i)=gradients_before{i}.aligned{1};
    after_aligned(:,:,i)=gradients_after{i}.aligned{1};

    before_lambda(:,i)=gradients_before{i}.lambda{1}./sum(gradients_before{i}.lambda{1});
    after_lambda(:,i)=gradients_after{i}.lambda{1}./sum(gradients_after{i}.lambda{1});

end
%% Visualize the mean gradients, eigenvalues and network breakdowns
mean_before_control=mean(before_aligned(:,:,1:num_control),3);
mean_before_stroke=mean(before_aligned(:,:,num_control+1:end),3);

mean_after_control=mean(after_aligned(:,:,1:num_control),3);
mean_after_stroke=mean(after_aligned(:,:,num_control+1:end),3);

obj=plot_hemispheres(mean_before_control(:,[1 2 3]), {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
obj.colormaps(slanCM('viridis'))
pretty(obj)

obj=plot_hemispheres(mean_before_stroke(:,1:3), {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});

plot_hemispheres(mean_after_stroke(:,1:3), {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});

gradient_in_euclidean(mean_before_control(:,[1 2 3]),{surf_lh,surf_rh},labeling.schaefer_400);
obj.colormaps([.1 .1 .1; slanCM('bwr')])
gradient_in_euclidean(mean_before_stroke(:,[1 2]),{surf_lh,surf_rh},labeling.schaefer_400);

gradient_in_euclidean(mean_after_control(:,1:2),{surf_lh,surf_rh},labeling.schaefer_400);
gradient_in_euclidean(mean_after_stroke(:,1:2),{surf_lh,surf_rh},labeling.schaefer_400);

var_before_control=var(before_scaled(:,:,1:num_control),[],3);
var_before_stroke=var(before_scaled(:,:,num_control+1:end),[],3);

var_after_control=var(after_scaled(:,:,1:num_control),[],3);
var_after_stroke=var(after_scaled(:,:,num_control+1:end),[],3);

plot_hemispheres(var_before_control(:,1:3), {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
plot_hemispheres(var_before_stroke(:,1:3), {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
plot_hemispheres(var_after_control(:,1:3), {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
plot_hemispheres(var_after_stroke(:,1:3), {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
% Get the 7 networks breakdown
for i=1:7
    for j=1:10
        network_means_control_before(i,j)=mean(mean_before_control(seven_network_positions==i,j));
        network_means_control_after(i,j)=mean(mean_after_control(seven_network_positions==i,j));

        network_means_stroke_before(i,j)=mean(mean_before_stroke(seven_network_positions==i,j));
        network_means_stroke_after(i,j)=mean(mean_after_stroke(seven_network_positions==i,j));

        network_sd_control_before(i,j)=std(mean_before_control(seven_network_positions==i,j))/sqrt(sum(seven_network_positions==i));
        network_sd_control_after(i,j)=std(mean_after_control(seven_network_positions==i,j))/sqrt(sum(seven_network_positions==i));

        network_sd_stroke_before(i,j)=std(mean_before_stroke(seven_network_positions==i,j))/sqrt(sum(seven_network_positions==i));
        network_sd_stroke_after(i,j)=std(mean_after_stroke(seven_network_positions==i,j))/sqrt(sum(seven_network_positions==i));
    end
end

% draw the network breakdown of the mean gradients 
for j=1:3
    h=figure;axes('XTick', 1:7, 'XTickLabel', network_names, 'FontSize', 16);hold on;ylabel('Mean Gradient Score')
    bar(network_means_control_before(:,j), 'FaceColor', 'flat', 'CData', colormap);
    ylim([-0.1 0.1])
    errorbar(1:numel(network_means_control_before(:,j)), network_means_control_before(:,j), network_sd_control_before(:,j), network_sd_control_before(:,j), ...
        'Color', 'k', 'LineStyle', 'None');
    filename='/home/koba/Desktop/Stroke/figures/revision/GX_control_before.svg';
    filename=strrep(filename,'X',num2str(j));
    saveas(h,filename)

    h=figure;axes('XTick', 1:7, 'XTickLabel', network_names, 'FontSize', 16);hold on;ylabel('Mean Gradient Score')
    bar(network_means_stroke_before(:,j), 'FaceColor', 'flat', 'CData', colormap);
    ylim([-0.1 0.1])
    errorbar(1:numel(network_means_stroke_before(:,j)), network_means_stroke_before(:,j), network_sd_stroke_before(:,j), network_sd_stroke_before(:,j), ...
        'Color', 'k', 'LineStyle', 'None');
    filename='/home/koba/Desktop/Stroke/figures/revision/GX_stroke_before.svg';
    filename=strrep(filename,'X',num2str(j));
    saveas(h,filename)


    h=figure;axes('XTick', 1:7, 'XTickLabel', network_names, 'FontSize', 16);hold on;ylabel('Mean Gradient Score')
    bar(network_means_control_after(:,j), 'FaceColor', 'flat', 'CData', colormap);
    errorbar(1:numel(network_means_control_after(:,j)), network_means_control_after(:,j), network_sd_control_after(:,j), network_sd_control_after(:,j), ...
        'Color', 'k', 'LineStyle', 'None');
    ylim([-0.1 0.1])
    filename='/home/koba/Desktop/Stroke/figures/revision/GX_control_after.svg';
    filename=strrep(filename,'X',num2str(j));
    saveas(h,filename)

    h=figure;axes('XTick', 1:7, 'XTickLabel', network_names, 'FontSize', 16);hold on;ylabel('Mean Gradient Score')
    bar(network_means_stroke_after(:,j), 'FaceColor', 'flat', 'CData', colormap);
    errorbar(1:numel(network_means_stroke_after(:,j)), network_means_stroke_after(:,j), network_sd_stroke_after(:,j), network_sd_stroke_after(:,j), ...
        'Color', 'k', 'LineStyle', 'None');
    filename='/home/koba/Desktop/Stroke/figures/revision/GX_stroke_after.svg';
    filename=strrep(filename,'X',num2str(j));
    saveas(h,filename)

end

[correlation, feature] = meta_analytic_decoder(parcel2full(mean_after_control(:,1), schaefer_400), ...
    'template', 'fsaverage5');
wc = wordcloud(feature(correlation>0), correlation(correlation>0));
disp(correlation(1:3));

scree_plot(mean(before_lambda(:,1:num_control),2))
scree_plot(mean(before_lambda(:,num_control+1:end),2))
scree_plot(mean(after_lambda(:,1:num_control),2))
scree_plot(mean(after_lambda(:,num_control+1:end),2))

%% ROI-based comparison of the gradient scores and comparison of lambdas 
for i=1:max(parcellation)
    for j=1:10
        [b,bint,r] = regress(squeeze(before_aligned(i,j,:)),[ones(130,1) subjects.age subjects.gender subjects.handed]);
        before_aligned_res(i,j,:)=r;

        [b,bint,r] = regress(squeeze(after_aligned(i,j,:)),[ones(130,1) subjects.age subjects.gender subjects.handed]);
        after_aligned_res(i,j,:)=r;
    end
end

for i=1:size(corrmats,1)
    for j=1:10
        [h,p,c,stats]=ttest(squeeze(before_aligned(i,j,1:num_control)),squeeze(after_aligned(i,j,1:num_control)));
        control_p(i,j)=p;
        control_t(i,j)=stats.tstat;
        [h,p,c,stats]=ttest(squeeze(before_aligned(i,j,num_control+1:end)),squeeze(after_aligned(i,j,num_control+1:end)));
        stroke_p(i,j)=p;
        stroke_t(i,j)=stats.tstat;

        [h,p,c,stats]=ttest2(squeeze(before_aligned(i,j,1:num_control)),squeeze(before_aligned(i,j,num_control+1:end)));
        controlXstroke_before_p(i,j)=p;
        controlXstroke_before_t(i,j)=stats.tstat;
        [h,p,c,stats]=ttest2(squeeze(after_aligned(i,j,1:num_control)),squeeze(after_aligned(i,j,num_control+1:end)));
        controlXstroke_after_p(i,j)=p;
        controlXstroke_after_t(i,j)=stats.tstat;
    end
end
figure
imagesc(fdr_bh(control_p).*control_t)
figure
imagesc(fdr_bh(stroke_p).*stroke_t)
figure
imagesc(fdr_bh(controlXstroke_before_p).*controlXstroke_before_t)
figure
imagesc(fdr_bh(controlXstroke_after_p).*controlXstroke_after_t)

fdr_mask=(fdr_bh(control_p).*control_t);
plot_hemispheres(fdr_mask(:,1:3), {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
fdr_mask=(fdr_bh(stroke_p).*stroke_t);
plot_hemispheres(fdr_mask(:,1:3), {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
f1=fdr_bh(control_p).*control_t;
f2=fdr_bh(stroke_p).*stroke_t;

for i=1:10
    [h,p,c,stats]=ttest(before_lambda(i,1:num_control),after_lambda(i,1:num_control));
    lambdas_control(i,1)=p;
    lambdas_control(i,2)=stats.tstat;
    [h,p,c,stats]=ttest(before_lambda(i,num_control+1:end),after_lambda(i,num_control+1:end));
    lambdas_stroke(i,1)=p;
    lambdas_stroke(i,2)=stats.tstat;

    [h,p,c,stats]=ttest2(before_lambda(i,1:num_control),before_lambda(i,num_control+1:end));
    lambdas_controlXstroke_before(i,1)=p;
    lambdas_controlXstroke_before(i,2)=stats.tstat;
    [h,p,c,stats]=ttest2(after_lambda(i,1:num_control),after_lambda(i,num_control+1:end));
    lambdas_controlXstroke_after(i,1)=p;
    lambdas_controlXstroke_after(i,2)=stats.tstat;

end
fdr_bh(lambdas_stroke(1:3,1))
%% Max distance to zero, 3D variance 

% ED to control in 3D 
for i=1:max(parcellation)
    for j=1:size(subjects,1)
        EDtocontrol_before(i,j)=norm(mean(before_aligned(i,1:3,1:num_control),3)-before_aligned(i,1:3,j));
        EDtocontrol_after(i,j)=norm(mean(after_aligned(i,1:3,1:num_control),3)-after_aligned(i,1:3,j));
    end
end

for i=1:max(parcellation)
   [b,bint,r] = regress(EDtocontrol_before(i,:)',[ones(130,1) subjects.age subjects.gender subjects.handed]);
   edfunc3d_residuals_before(i,:)=r;

   [b,bint,r] = regress(EDtocontrol_after(i,:)',[ones(130,1) subjects.age subjects.gender subjects.handed]);
   edfunc3d_residuals_after(i,:)=r;
end

[h,p,c,stats]=ttest2(mean(EDtocontrol_before(:,1:num_control),2),mean(EDtocontrol_before(:,num_control+1:end),2))
[h,p,c,stats]=ttest2(mean(edfunc3d_residuals_before(:,1:num_control),2),mean(edfunc3d_residuals_before(:,num_control+1:end),2))

[h,p,c,stats]=ttest2(mean(EDtocontrol_after(:,1:num_control),2),mean(EDtocontrol_after(:,num_control+1:end),2))
[h,p,c,stats]=ttest2(mean(edfunc3d_residuals_after(:,1:num_control),2),mean(edfunc3d_residuals_after(:,num_control+1:end),2))


hold on;histogram(mean(EDtocontrol_before(:,1:num_control),2));histogram(mean(EDtocontrol_before(:,num_control+1:end),2));hold off
figure
hold on;histogram(mean(EDtocontrol_after(:,1:num_control),2));histogram(mean(EDtocontrol_after(:,num_control+1:end),2));hold off
xlabel('Mean ED_f_u_n_c across whole brain')
ylabel('Frequency')
legend('Control','Stroke')
dif1=mean(EDtocontrol_before(:,1:num_control),2)-mean(EDtocontrol_before(:,num_control+1:end),2);
dif2=mean(EDtocontrol_after(:,1:num_control),2)-mean(EDtocontrol_after(:,num_control+1:end),2);
[h,p,c,stats]=ttest2(dif1,dif2)




subjects.lesion_side==0
dif1=mean(EDtocontrol_before(:,1:num_control),2)-mean(EDtocontrol_before(:,num_control+1:end),2);
dif2=mean(EDtocontrol_after(:,1:num_control),2)-mean(EDtocontrol_after(:,num_control+1:end),2);

for i=1:10000
    ed_shuffled1=EDtocontrol_before(:,randperm(length(1:130)));
    ed_shuffled2=EDtocontrol_after(:,randperm(length(1:130)));
    shuffled_dif1(:,i)=mean(ed_shuffled1(:,1:num_control),2)-mean(ed_shuffled1(:,num_control+1:end),2);
    shuffled_dif2(:,i)=mean(ed_shuffled1(:,1:num_control),2)-mean(ed_shuffled2(:,num_control+1:end),2);

end
for i=1:400
    corrected1(i)=(abs(dif1(i)))>prctile(abs(shuffled_dif1(i,:)),95);
    corrected2(i)=(abs(dif2(i)))>prctile(abs(shuffled_dif2(i,:)),95);

end
plot_hemispheres([mean(EDtocontrol_before(:,1:num_control),2) mean(EDtocontrol_before(:,num_control+1:end),2) corrected1'.*dif1 corrected2'.*dif2], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400);

for i=1:400
    [h,p,c,stats]=ttest2(EDtocontrol_before(i,1:num_control),EDtocontrol_before(i,num_control+1:end));
    ps(i)=p;
    ts(i)=stats.tstat;
end

dd=(corrected1'.*dif1).*-1;
dd2=(corrected2'.*dif2).*-1;



for i=1:7
    network_means_EDtocontrol_before(i)=mean(dd(seven_network_positions==i));
    network_means_EDtocontrol_after(i)=mean(dd2(seven_network_positions==i));

    network_sd_EDtocontrol_before(i)=std(EDtocontrol_before(seven_network_positions==i))/sqrt(sum(seven_network_positions==i));
    network_sd_EDtocontrol_after(i)=std(EDtocontrol_after(seven_network_positions==i))/sqrt(sum(seven_network_positions==i));
end

figure;axes('XTick', 1:7, 'XTickLabel', network_names, 'FontSize', 16);hold on;ylabel('Mean ED_f_u_n_c difference')
bar(network_means_EDtocontrol_before, 'FaceColor', 'flat', 'CData', colormap);
errorbar(1:numel(network_means_EDtocontrol_before), network_means_EDtocontrol_before, network_sd_EDtocontrol_before, network_sd_EDtocontrol_before, ...
    'Color', 'k', 'LineStyle', 'None');

figure;axes('XTick', 1:7, 'XTickLabel', network_names, 'FontSize', 16);hold on;ylabel('Mean ED_f_u_n_c difference')
bar(network_means_EDtocontrol_after, 'FaceColor', 'flat', 'CData', colormap);
errorbar(1:numel(network_means_EDtocontrol_after), network_means_EDtocontrol_after, network_sd_EDtocontrol_after, network_sd_EDtocontrol_after, ...
    'Color', 'k', 'LineStyle', 'None');



% Do it for the gradients seperately

for i=1:max(parcellation)
    for j=1:size(subjects,1)
        for k=1:3
            EDtocontrol_before_sep(i,j,k)=norm(mean(before_aligned(i,k,1:num_control),3)-before_aligned(i,k,j));
            EDtocontrol_after_sep(i,j,k)=norm(mean(after_aligned(i,k,1:num_control),3)-after_aligned(i,k,j));
% 
%             EDtozero_before_sep(i,j,k)=norm(before_scaled(i,k,j));
%             EDtozero_after_sep(i,j,k)=norm(after_scaled(i,k,j));
        end
    end
end

[h,p,c,stats]=ttest2(mean(EDtocontrol_before_sep(:,1:num_control,1),2),mean(EDtocontrol_before_sep(:,num_control+1:end,1),2))
[h,p,c,stats]=ttest2(mean(EDtocontrol_before_sep(:,1:num_control,2),2),mean(EDtocontrol_before_sep(:,num_control+1:end,2),2))
[h,p,c,stats]=ttest2(mean(EDtocontrol_before_sep(:,1:num_control,3),2),mean(EDtocontrol_before_sep(:,num_control+1:end,3),2))


[h,p,c,stats]=ttest2(mean(EDtocontrol_after_sep(:,1:num_control,1),2),mean(EDtocontrol_after_sep(:,num_control+1:end,1),2))
[h,p,c,stats]=ttest2(mean(EDtocontrol_after_sep(:,1:num_control,2),2),mean(EDtocontrol_after_sep(:,num_control+1:end,2),2))
[h,p,c,stats]=ttest2(mean(EDtocontrol_after_sep(:,1:num_control,3),2),mean(EDtocontrol_after_sep(:,num_control+1:end,3),2))

dif11=mean(EDtocontrol_before_sep(:,1:num_control,1),2)-mean(EDtocontrol_before_sep(:,num_control+1:end,1),2)
dif12=mean(EDtocontrol_before_sep(:,1:num_control,2),2)-mean(EDtocontrol_before_sep(:,num_control+1:end,2),2)
dif13=mean(EDtocontrol_before_sep(:,1:num_control,3),2)-mean(EDtocontrol_before_sep(:,num_control+1:end,3),2)

dif21=mean(EDtocontrol_after_sep(:,1:num_control,1),2)-mean(EDtocontrol_after_sep(:,num_control+1:end,1),2);
dif22=mean(EDtocontrol_after_sep(:,1:num_control,2),2)-mean(EDtocontrol_after_sep(:,num_control+1:end,2),2);
dif23=mean(EDtocontrol_after_sep(:,1:num_control,3),2)-mean(EDtocontrol_after_sep(:,num_control+1:end,3),2);
plot_hemispheres([dif21 dif22 dif23], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Dif 1st ','Dif 2nd', 'Dif3rd'});

for i=1:1000
ed_shuffled=EDtocontrol_before(:,randperm(length(1:130)));

shuffled_dif(:,i)=mean(ed_shuffled(:,1:num_control),2)-mean(ed_shuffled(:,num_control+1:end),2);
end
for i=1:400
    corrected(i)=(abs(dif11(i)))>prctile(abs(shuffled_dif(i,:)),95);
end
 
plot_hemispheres([dif11 dif12 dif13], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Dif 1st ','Dif 2nd', 'Dif3rd'});


for i=1:10000
    ed_shuffled1_bef=EDtocontrol_before_sep(:,randperm(length(1:130)),1);
    shuffled_dif1_bef(:,i)=mean(ed_shuffled1_bef(:,1:num_control),2)-mean(ed_shuffled1_bef(:,num_control+1:end),2);

    ed_shuffled2_bef=EDtocontrol_before_sep(:,randperm(length(1:130)),2);
    shuffled_dif2_bef(:,i)=mean(ed_shuffled2_bef(:,1:num_control),2)-mean(ed_shuffled2_bef(:,num_control+1:end),2);

    ed_shuffled3_bef=EDtocontrol_before_sep(:,randperm(length(1:130)),3);
    shuffled_dif3_bef(:,i)=mean(ed_shuffled3_bef(:,1:num_control),2)-mean(ed_shuffled3_bef(:,num_control+1:end),2);

    ed_shuffled1_after=EDtocontrol_after_sep(:,randperm(length(1:130)),1);
    shuffled_dif1_after(:,i)=mean(ed_shuffled1_after(:,1:num_control),2)-mean(ed_shuffled1_after(:,num_control+1:end),2);

    ed_shuffled2_after=EDtocontrol_after_sep(:,randperm(length(1:130)),2);
    shuffled_dif2_after(:,i)=mean(ed_shuffled2_bef(:,1:num_control),2)-mean(ed_shuffled2_after(:,num_control+1:end),2);

    ed_shuffled3_after=EDtocontrol_after_sep(:,randperm(length(1:130)),3);
    shuffled_dif3_after(:,i)=mean(ed_shuffled3_after(:,1:num_control),2)-mean(ed_shuffled3_after(:,num_control+1:end),2);

end
for i=1:400
    corrected1_bef(i)=(abs(dif11(i)))>prctile(abs(shuffled_dif1_bef(i,:)),95);
    corrected2_bef(i)=(abs(dif12(i)))>prctile(abs(shuffled_dif2_bef(i,:)),95);
    corrected3_bef(i)=(abs(dif13(i)))>prctile(abs(shuffled_dif3_bef(i,:)),95);

    corrected1_after(i)=(abs(dif21(i)))>prctile(abs(shuffled_dif1_after(i,:)),95);
    corrected2_after(i)=(abs(dif22(i)))>prctile(abs(shuffled_dif2_after(i,:)),95);
    corrected3_after(i)=(abs(dif23(i)))>prctile(abs(shuffled_dif3_after(i,:)),95);
end
plot_hemispheres([corrected1_bef'.*dif11 corrected2_bef'.*dif12 corrected3_bef'.*dif13], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Dif 1st ','Dif 2nd', 'Dif3rd'});
plot_hemispheres([corrected1_after'.*dif21 corrected2_after'.*dif22 corrected3_after'.*dif23], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Dif 1st ','Dif 2nd', 'Dif3rd'});
g11=corrected1_bef'.*dif11;
g12=corrected2_bef'.*dif12;
g13=corrected3_bef'.*dif13;
g21=corrected1_after'.*dif21;
g22=corrected2_after'.*dif22;
g23=corrected3_after'.*dif23;

dd=g23*-1;
for i=1:7
        network_means_EDtocontrol_before(i)=mean(dd(seven_network_positions==i));
%         network_means_EDtocontrol_after(i)=mean(EDtocontrol_after(seven_network_positions==i));

        network_sd_EDtocontrol_before(i)=std(dd(seven_network_positions==i))/sqrt(sum(seven_network_positions==i));
%         network_sd_EDtocontrol_after(i)=std(EDtocontrol_after(seven_network_positions==i));
end

figure;axes('XTick', 1:7, 'XTickLabel', network_names, 'FontSize', 16);hold on;ylabel('Mean ED to control - before')
    bar(network_means_EDtocontrol_before, 'FaceColor', 'flat', 'CData', colormap);
     errorbar(1:numel(network_means_EDtocontrol_before), network_means_EDtocontrol_before, network_sd_EDtocontrol_before, network_sd_EDtocontrol_before, ...
         'Color', 'k', 'LineStyle', 'None');




[h,p,c,stats]=ttest2(dif11,dif21)
[h,p,c,stats]=ttest2(dif12,dif22)
[h,p,c,stats]=ttest2(dif13,dif23)

%% Correlation between difs, gradients, lags

% Corr with mean lag and gradients
[rho,p1] = corr(mean(mean_lag(:,1:num_control),2),mean(squeeze(before_aligned(:,1,:)),2)) % 1: 0.25 2: 0.7 3: 0.19
[rho,p2] = corr(mean(mean_lag(:,1:num_control),2),mean(squeeze(before_aligned(:,2,:)),2)) % 1: 0.25 2: 0.7 3: 0.19
[rho,p3] = corr(mean(mean_lag(:,1:num_control),2),mean(squeeze(before_aligned(:,3,:)),2)) % 1: 0.25 2: 0.7 3: 0.19



% Individual lag corrs with gradients or diffs from both side


lag_diff=mean(mean_lag(:,1:num_control),2)-mean_lag;
% corr(mean(mean_lag,2),mean(squeeze(after_aligned(:,2,1:num_control)),2))
% corr(mean(mean_lag,2),mean(squeeze(after_aligned(:,2,num_control+1:end)),2))
[rho,p1] = corr(mean(lag_diff(:,1:num_control),2),mean(EDtocontrol_before(:,1:num_control),2))
[rho,p1] = corr(mean(lag_diff(:,num_control+1:end),2),mean(EDtocontrol_before(:,num_control+1:end),2))
corr(mean(lag_diff(:,1:num_control),2),mean(EDtocontrol_after(:,1:num_control),2))
corr(mean(lag_diff(:,num_control+1:end),2),mean(EDtocontrol_after(:,num_control+1:end),2)) % lag diff is higher for stroke 
for i=1:130
    corrs(i,1)=corr(lag_diff(:,i),EDtocontrol_before(:,i))
    corrs(i,2)=corr(lag_diff(:,i),EDtocontrol_after(:,i))

    corrs(i,3)=corr(lag_diff(:,i),EDtocontrol_before(:,i))
    corrs(i,4)=corr(lag_diff(:,i),EDtocontrol_after(:,i))
end
[h,p,c,stats]=ttest2(corrs(1:num_control,4),corrs(num_control+1:end,4))
% for i=1:130
%     corrs(i,1)=corr(lag_diff(:,i),before_aligned(:,1,i))
%     corrs(i,2)=corr(lag_diff(:,i),after_aligned(:,1,i))
% 
%     corrs(i,3)=corr(lag_diff(:,i),before_aligned(:,2,i))
%     corrs(i,4)=corr(lag_diff(:,i),after_aligned(:,2,i))
% end
control_beh=behavioral.bigfactor3(1:num_control);
stroke_beh=behavioral.bigfactor3(num_control+1:end);
control_beh_nonan=control_beh(~isnan(control_beh));
stroke_beh_nonan=stroke_beh(~isnan(stroke_beh));


[h,p,c,stats]=ttest2(control_beh_nonan,stroke_beh_nonan)

[h,p,c,stats]=ttest((stroke_beh_nonan-mean(control_beh_nonan))/std(control_beh_nonan))

ix=~isnan(subjects.nihss_hospital);
[bb(:,1) bb(:,2)]=corr(EDtocontrol_before(:,ix)',subjects.nihss_hospital(ix));
[bc(:,1) bc(:,2)]=corr(EDtocontrol_after(:,ix)',subjects.nihss_hospital(ix));
% dd=corr((lag_diff(:,ix))',subjects.nihss_hospital(ix));
plot_hemispheres([bb(:,1).*(bb(:,2)<0.05) bc(:,1).*(bc(:,2)<0.05)], {surf_lh,surf_rh}, ...
'parcellation',labeling.schaefer_400, 'labeltext',{'2controlb','2controla'});

% Ed func with behaviorals -- calculate the difference between groups
ix=~isnan(behavioral.bigfactor2);
 ix=ix+[zeros(num_control,1);ones(num_stroke,1)];
%  ix=ix+[ones(num_control,1);zeros(num_stroke,1);];
%  ix=ix+(subjects.lesion_side==1);
ix=ix==2;
sum(ix)
[bb(:,1) bb(:,2)]=corr(EDtocontrol_before(:,ix)',behavioral.bigfactor2(ix));
[bc(:,1) bc(:,2)]=corr(EDtocontrol_after(:,ix)',behavioral.bigfactor2(ix));
% bbsig=bb(:,1).*((bb(:,2)<0.05));
% bcsig=bc(:,1).*((bc(:,2)<0.05));
bbsig=bb(:,1).*fdr_bh(bb(:,2));
bcsig=bc(:,1).*fdr_bh(bc(:,2));
plot_hemispheres([bbsig bcsig ], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400);
sum((abs(bbsig)>0))
sum((abs(bcsig)>0))
[h,p,c,stats]=ttest2(bbsig(abs(bbsig)>0),bcsig(abs(bcsig)>0))
% Difference between groups or hemispheres
[h,p,c,stats]=ttest2((atanh(bb(subjects.lesion_side(ix)==1))),(atanh(bb(subjects.lesion_side(ix)==0))))
% mean(abs(atanh(bb(subjects.subj_type(ix)==1))))
% mean(abs(atanh(bb(subjects.subj_type(ix)==0))))
beh2_before_fdr=bbsig;
beh2_after_fdr=bcsig;


dd=bbsig.*-1;
ee=bcsig.*-1;
for i=1:7
        network_means_EDtocontrol_before(i)=mean(dd(seven_network_positions==i));
      network_means_EDtocontrol_after(i)=mean(ee(seven_network_positions==i));

        network_sd_EDtocontrol_before(i)=std(dd(seven_network_positions==i))/sqrt(sum(seven_network_positions==i));
      network_sd_EDtocontrol_after(i)=std(ee(seven_network_positions==i))/sqrt(sum(seven_network_positions==i));
end

figure;axes('XTick', 1:7, 'XTickLabel', network_names, 'FontSize', 16);hold on;ylabel('Mean correlation')
    bar(network_means_EDtocontrol_before, 'FaceColor', 'flat', 'CData', colormap);
     errorbar(1:numel(network_means_EDtocontrol_before), network_means_EDtocontrol_before, network_sd_EDtocontrol_before, network_sd_EDtocontrol_before, ...
         'Color', 'k', 'LineStyle', 'None');
figure;axes('XTick', 1:7, 'XTickLabel', network_names, 'FontSize', 16);hold on;ylabel('Mean correlation')
    bar(network_means_EDtocontrol_after, 'FaceColor', 'flat', 'CData', colormap);
     errorbar(1:numel(network_means_EDtocontrol_after), network_means_EDtocontrol_after, network_sd_EDtocontrol_after, network_sd_EDtocontrol_after, ...
         'Color', 'k', 'LineStyle', 'None');


% Lag with behaviorals
ix=~isnan(behavioral.language_f);
ix=ix+[zeros(num_control,1);ones(num_stroke,1)];
% ix=ix+[ones(num_control,1);zeros(num_stroke,1);];
% ix=ix+(subjects.lesion_side==1);
ix=ix==2;
sum(ix)
[dd(:,1) dd(:,2)]=corr((lag_diff(:,ix))',behavioral.language_f(ix));
ix=~isnan(behavioral.motorl_f);
ix=ix+[zeros(num_control,1);ones(num_stroke,1)];
% ix=ix+[ones(num_control,1);zeros(num_stroke,1);];
% ix=ix+(subjects.lesion_side==1);
ix=ix==2;
sum(ix)
[dc(:,1) dc(:,2)]=corr((lag_diff(:,ix))',behavioral.motorl_f(ix));
ix=~isnan(behavioral.motorr_f);
ix=ix+[zeros(num_control,1);ones(num_stroke,1)];
% ix=ix+[ones(num_control,1);zeros(num_stroke,1);];
% ix=ix+(subjects.lesion_side==1);
ix=ix==2;
sum(ix)
[de(:,1) de(:,2)]=corr((lag_diff(:,ix))',behavioral.motorr_f(ix));

ix=~isnan(subjects.nihss_hospital);
[df(:,1) df(:,2)]=corr((lag_diff(:,ix))',subjects.nihss_hospital(ix));
plot_hemispheres([dd(:,1).*(dd(:,2)<0.05) dc(:,1).*(dc(:,2)<0.05) de(:,1).*(de(:,2)<0.05) df(:,1).*(df(:,2)<0.05)], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400);






[correlation, feature] = meta_analytic_decoder(parcel2full(double(bbsig<0), schaefer_400), ...
    'template', 'fsaverage5');
wc = wordcloud(feature(correlation>0), correlation(correlation>0));





%% Anatomical distance to lesion 

mask_coords=table2array(readtable('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/distance2lesion/empty_mask.txt')); % load the template


masklist=readtable('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/masklist.txt','ReadVariableNames', false); % load the subjects list that have masks 

    
% Load the distance 2 lesion info
distance2lesion=load_nii('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/distance2lesion/all_distances.nii.gz');
distance2lesion=reshape(squeeze((distance2lesion.img)),[],84);
distance2lesion=distance2lesion(parcellation_mask_whole,:);


for j=1:max(parcellation_reduced)
    for i=1:84
        ED_anat(j,i)=mean(distance2lesion(parcellation_reduced==j,i));
    end
end

% Get the subset of EDto zero and lag 
mask_index=zeros(130,1);
for i=1:size(masklist,1)
% 
   splt=split(masklist{i,2},'_');
   if strcmp(splt{1},'sub-PAT114' )|| strcmp(splt{1},'sub-PAT154' )
    fprintf('skip\n')
   else

   mask_index=mask_index+strcmp(subjects.participant_id,splt{1});
   end
end

ED_func_before=EDtocontrol_before(:,logical(mask_index));
ED_func_after=EDtocontrol_after(:,logical(mask_index));
lags_subset=mean_lag(:,logical(mask_index));
lag_diff_subset=mean(mean_lag(:,1:num_control),2)-lags_subset;
behavioral_subset=behavioral(logical(mask_index),:);
subjects_subset=subjects(logical(mask_index),:);


for i=1:size(ED_anat,2)
    within_sub_corrs_before(i)=corr(ED_func_before(:,i),ED_anat(:,i));
    within_sub_corrs_after(i)=corr(ED_func_after(:,i),ED_anat(:,i));
    within_sub_corrs_lag(i)=corr(lag_diff_subset(:,i),ED_anat(:,i));
end
[h,p,c,stats]=ttest(atanh(within_sub_corrs_lag))
[h,p,c,stats]=ttest(atanh(within_sub_corrs_after))
[h,p,c,stats]=ttest(atanh(within_sub_corrs_before)-atanh(within_sub_corrs_after))
mean(within_sub_corrs_after)
hold on;histogram(within_sub_corrs_before);histogram(within_sub_corrs_after);hold off
% ix=subjects_subset.lesion_site==6;
for i=1:400
    [between_sub_corrs_before(i,1) between_sub_corrs_before(i,2)]=corr(ED_func_before(i,num_control+1:end)',ED_anat(i,num_control+1:end)');
    [between_sub_corrs_after(i,1) between_sub_corrs_after(i,2)]=corr(ED_func_after(i,num_control+1:end)',ED_anat(i,num_control+1:end)');
    [between_sub_corrs_lag(i,1) between_sub_corrs_lag(i,2)]=corr(lag_diff_subset(i,num_control+1:end)',ED_anat(i,num_control+1:end)');
end
anat_before=between_sub_corrs_before(:,1).*(between_sub_corrs_before(:,2)<0.05)
anat_after=between_sub_corrs_after(:,1).*(between_sub_corrs_after(:,2)<0.05)

plot_hemispheres([ anat_before anat_after ((between_sub_corrs_after(:,2))) between_sub_corrs_lag(:,1).*(between_sub_corrs_lag(:,2)<0.05)], {surf_lh,surf_rh}, ...
'parcellation',labeling.schaefer_400);

bf4=between_sub_corrs_before(:,1).*((between_sub_corrs_before(:,2)<0.05))
dd=anat_before.*-1;
ee=anat_after.*-1;
for i=1:7
        network_means_EDtocontrol_before(i)=mean(dd(seven_network_positions==i));
         network_means_EDtocontrol_after(i)=mean(ee(seven_network_positions==i));

        network_sd_EDtocontrol_before(i)=std(dd(seven_network_positions==i))/sqrt(sum(seven_network_positions==i));
         network_sd_EDtocontrol_after(i)=std(ee(seven_network_positions==i))/sqrt(sum(seven_network_positions==i));
end

figure;axes('XTick', 1:7, 'XTickLabel', network_names, 'FontSize', 16);hold on;ylabel('Mean correlation')
    bar(network_means_EDtocontrol_before, 'FaceColor', 'flat', 'CData', colormap);
     errorbar(1:numel(network_means_EDtocontrol_before), network_means_EDtocontrol_before, network_sd_EDtocontrol_before, network_sd_EDtocontrol_before, ...
         'Color', 'k', 'LineStyle', 'None');

figure;axes('XTick', 1:7, 'XTickLabel', network_names, 'FontSize', 16);hold on;ylabel('Mean correlation')
    bar(network_means_EDtocontrol_after, 'FaceColor', 'flat', 'CData', colormap);
     errorbar(1:numel(network_means_EDtocontrol_after), network_means_EDtocontrol_after, network_sd_EDtocontrol_after, network_sd_EDtocontrol_after, ...
         'Color', 'k', 'LineStyle', 'None');

corr(mean(ED_func_before(:,1:num_control),2),mean(ED_anat(:,1:num_control),2))
corr(mean(ED_func_before(:,num_control+1:end),2),mean(ED_anat(:,num_control+1:end),2))
corr(mean(ED_func_after(:,1:num_control),2),mean(ED_anat(:,1:num_control),2))
corr(mean(ED_func_after(:,num_control+1:end),2),mean(ED_anat(:,num_control+1:end),2))



corr(mean(lag_diff(:,1:num_control),2),mean(EDtocontrol_after(:,1:num_control),2))
corr(mean(lag_diff(:,num_control+1:end),2),mean(EDtocontrol_after(:,num_control+1:end),2))
 for i=1:84
    corrs(i,1)=corr(lag_diff(:,i),ED_func_before(:,i))
    corrs(i,2)=corr(lag_diff(:,i),ED_func_after(:,i))

    corrs(i,3)=corr(lag_diff(:,i),ED_func_before(:,i))
    corrs(i,4)=corr(lag_diff(:,i),ED_func_after(:,i))
end

ix=~isnan(behavioral_subset.bigfactor3);
%  ix=ix+(subjects_subset.lesion_side==0);
% ix=ix==2;
sum(ix)
for i=1:400
    [temp_a(i,1) temp_a(i,2)]=corr(ED_func_before(i,ix)',behavioral_subset.bigfactor3(ix));
    [temp_b(i,1) temp_b(i,2)]=corr(ED_func_after(i,ix)',behavioral_subset.bigfactor3(ix));
    [temp_c(i,1) temp_c(i,2)]=corr(lag_diff_subset(i,ix)',behavioral_subset.bigfactor3(ix));
    [temp_d(i,1) temp_d(i,2)]=corr(ED_anat(i,ix)',behavioral_subset.bigfactor3(ix));
end

a=fdr_bh(temp_a(:,2)).*temp_a(:,1);
b=fdr_bh(temp_b(:,2)).*temp_b(:,1);
c=fdr_bh(temp_c(:,2)).*temp_c(:,1);
d=temp_d(:,1).*fdr_bh(temp_d(:,2));

plot_hemispheres([a b c d], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'ED_f_u_u_n_c before','ED_f_u_u_n_c after', 'Lag', 'ED_a_n_a_t'} );
% 
% 
% 
% 
% 
% 

% Lesion locations on surface
lesions_on_gm=load_nii('/media/koba/MULTIBOOT/LesionMaks_222/NiftiLesions/other/lesions_on_gm.nii.gz');
lesions_on_gm=reshape(squeeze((lesions_on_gm.img)),[],1);
lesions_on_gm=lesions_on_gm(parcellation_mask_whole,:);

for j=1:max(parcellation_reduced)
        mean_lesions_on_gm(j,1)=mean(lesions_on_gm(parcellation_reduced==j,:));
end

obj=plot_hemispheres([mean_lesions_on_gm], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400);
obj.colormaps(slanCM('magma'))

dd=mean_lesions_on_gm(:,1);
for i=1:7
        network_means_EDtocontrol_before(i)=mean(dd(seven_network_positions==i));
%         network_means_EDtocontrol_after(i)=mean(EDtocontrol_after(seven_network_positions==i));

        network_sd_EDtocontrol_before(i)=std(dd(seven_network_positions==i));
%         network_sd_EDtocontrol_after(i)=std(EDtocontrol_after(seven_network_positions==i));
end

figure;axes('XTick', 1:7, 'XTickLabel', network_names, 'FontSize', 16);hold on;ylabel('Mean ED to control - before')
    bar(network_means_EDtocontrol_before, 'FaceColor', 'flat', 'CData', colormap);
     errorbar(1:numel(network_means_EDtocontrol_before), network_means_EDtocontrol_before, network_sd_EDtocontrol_before, network_sd_EDtocontrol_before, ...
         'Color', 'k', 'LineStyle', 'None');

 for i=1:84
    [anats(i,1) anats(i,2)]=corr(ED_anat(:,i),mean_lesions_on_gm);
 end
 mean(anats)
[h,p,c,stats]=ttest(atanh(anats(:,1)))
histogram(atanh(anats(:,1)))

 for i=1:130
    [funcs(i,1) funcs(i,2)]=corr(EDtocontrol_before(:,i),mean_lesions_on_gm);
 end
 mean(funcs)
[h,p,c,stats]=ttest(atanh(funcs(:,1)))
histogram(atanh(funcs(:,1)))
%% Heterogeneity: Jaccard index 

% Load the resampled lesion masks 

mask_str='/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/NiftiLesions/resampled/SUB_lesion_resampled.nii.gz';
for i=1:84
    disp(i)
    for j=1:84

        mask_sub1=strrep(mask_str,'SUB',subjects_subset.participant_id{i});
        mask_sub1_img=load_nii(mask_sub1);
        mask_sub1_img=mask_sub1_img.img;

        mask_sub2=strrep(mask_str,'SUB',subjects_subset.participant_id{j});
        mask_sub2_img=load_nii(mask_sub2);
        mask_sub2_img=mask_sub2_img.img;

        % Find coordinates of lesions in matrix1 and matrix2
        [A_x, A_y, A_z] = ind2sub(size(mask_sub1_img), find(mask_sub1_img == 1));
        [B_x, B_y, B_z] = ind2sub(size(mask_sub2_img), find(mask_sub2_img == 1));

        % Calculate Dice Similarity
        union = sum(sum(sum(mask_sub1_img | mask_sub2_img)));
        inter = sum(sum(sum(mask_sub1_img & mask_sub2_img)));
        total_sum = sum(sum(sum(mask_sub1_img))) + sum(sum(sum(mask_sub2_img)));
        %intersection_dice = nnz(ismember([A_x, A_y, A_z], [B_x, B_y, B_z]));
        dice_similarity = 2 * inter / (total_sum);

        % Calculate Jaccard Similarity
        % intersection_jaccard = nnz(ismember([A_x, A_y, A_z], [B_x, B_y, B_z]));
        % union_jaccard = numel(unique([A_x, A_y, A_z; B_x, B_y, B_z], 'rows'));
        % jaccard_similarity = intersection_jaccard / union_jaccard;
        jaccard_similarity = inter/union;

        % Calculate centroid of lesions
        centroid_A = [mean(A_x), mean(A_y), mean(A_z)];
        centroid_B = [mean(B_x), mean(B_y), mean(B_z)];

        % Calculate Euclidean distance between centroids
        distance = norm(centroid_A - centroid_B);
        rho1=corr(ED_func_before(:,i),ED_func_before(:,j),'Type','Spearman');
        rho2=corr(ED_func_before(:,i),ED_func_after(:,j),'Type','Spearman');

        dicecoefs(i,j)=dice_similarity;
        jaccard(i,j)=jaccard_similarity;
        ED_between_lesions(i,j)=distance;
        ED_func_corr_before(i,j)=rho1;
        ED_func_corr_after(i,j)=rho2;
    end
end



dicecoefs1D=dicecoefs(tril(true(size(dicecoefs)), -1));
jaccard1D=jaccard(tril(true(size(jaccard)), -1));
ED_between_lesions1D=ED_between_lesions(tril(true(size(ED_between_lesions)), -1));
ED_func_corr_before1D=ED_func_corr_before(tril(true(size(ED_func_corr_before)), -1));
ED_func_corr_after1D=ED_func_corr_after(tril(true(size(ED_func_corr_after)), -1));
[r,p] = corr([dicecoefs1D jaccard1D ED_between_lesions1D ED_func_corr_before1D ED_func_corr_after1D ],'Type','Spearman');

dicecoefs_left=dicecoefs(subjects_subset.lesion_side==0,subjects_subset.lesion_side==0);
dicecoefs_right=dicecoefs(subjects_subset.lesion_side==1,subjects_subset.lesion_side==1);

jaccard_left=jaccard(subjects_subset.lesion_side==0,subjects_subset.lesion_side==0);
jaccard_right=jaccard(subjects_subset.lesion_side==1,subjects_subset.lesion_side==1);

ED_between_lesions_left=ED_between_lesions(subjects_subset.lesion_side==0,subjects_subset.lesion_side==0);
ED_between_lesions_right=ED_between_lesions(subjects_subset.lesion_side==1,subjects_subset.lesion_side==1);

ED_func_corr_before_left=ED_func_corr_before(subjects_subset.lesion_side==0,subjects_subset.lesion_side==0);
ED_func_corr_before_right=ED_func_corr_before(subjects_subset.lesion_side==1,subjects_subset.lesion_side==1);

ED_func_corr_after_left=ED_func_corr_after(subjects_subset.lesion_side==0,subjects_subset.lesion_side==0);
ED_func_corr_after_right=ED_func_corr_after(subjects_subset.lesion_side==1,subjects_subset.lesion_side==1);

dicecoefs_left1D=dicecoefs_left(tril(true(size(dicecoefs_left)), -1));
jaccard_left1D=jaccard_left(tril(true(size(jaccard_left)), -1));
ED_between_lesions_left1D=ED_between_lesions_left(tril(true(size(ED_between_lesions_left)), -1));
ED_func_corr_before_left1D=ED_func_corr_before_left(tril(true(size(ED_func_corr_before_left)), -1));
ED_func_corr_after_left1D=ED_func_corr_after_left(tril(true(size(ED_func_corr_after_left)), -1));

dicecoefs_right1D=dicecoefs_right(tril(true(size(dicecoefs_right)), -1));
jaccard_right1D=jaccard_right(tril(true(size(jaccard_right)), -1));
ED_between_lesions_right1D=ED_between_lesions_right(tril(true(size(ED_between_lesions_right)), -1));
ED_func_corr_before_right1D=ED_func_corr_before_right(tril(true(size(ED_func_corr_before_right)), -1));
ED_func_corr_after_right1D=ED_func_corr_after_right(tril(true(size(ED_func_corr_after_right)), -1));



[r1,p1] = corr([dicecoefs_left1D jaccard_left1D ED_between_lesions_left1D ED_func_corr_before_left1D ED_func_corr_after_left1D ]);
[r2,p2] = corr([dicecoefs_right1D(jaccard_right1D>0.2) jaccard_right1D(jaccard_right1D>0.2) ED_between_lesions_right1D(jaccard_right1D>0.2) ED_func_corr_before_right1D(jaccard_right1D>0.2) ED_func_corr_after_right1D(jaccard_right1D>0.2) ]);

