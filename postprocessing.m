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
[network_names, colormap] = fetch_yeo_networks_metadata(7);
schaefer_400 = fetch_parcellation('fsaverage5', 'schaefer', 400);
[surf_lh, surf_rh] = load_conte69();
labeling = load_parcellation('schaefer',400);


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

% Regress out the confounds - cancelled
% confounds=[subjects.age subjects.gender subjects.handed];
% for i=1:size(corrmats,1)
%     for j=1:size(corrmats,1)
%         [~,~,residuals] = regress(squeeze(corrmats_reduced(i,j,:)),[ones(size(confounds,1),1) confounds]);
%         corrmats_reduced_clean(i,j,:)=residuals;
%
%
%         [~,~,residuals] = regress(squeeze(corrmats_corrected_reduced(i,j,:)),[ones(size(confounds,1),1) confounds]);
%         corrmats_corrected_reduced_clean(i,j,:)=residuals;
%
%     end
% end
% for i=1:size(subjects,1)
%     for j=1:size(corrmats,1)
%         corrmats_reduced_clean(j,j,i)=1;
%         corrmats_corrected_reduced_clean(j,j,i)=1;
%     end
% end

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
    [h,p,c,stats]=ttest2(mean_lag(i,1:num_control),mean_lag(i,num_control+1:end))
    lag_ps(i)=p;
    lag_ts(i)=stats.tstat;
end
plot_hemispheres([mean(mean_lag(:,1:num_control),2) mean(mean_lag(:,num_control+1:end),2)], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'controls','Mean Lag-stroke'});

plot_hemispheres([lag_ts'], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'lag diff roi',});


labeling = load_parcellation('schaefer',400);

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
hold on;histogram(mean_abr_stroke);histogram(mean_abr_control);hold off
figure
hold on;histogram(mean_abr_stroke_corrected);histogram(mean_abr_control_corrected);hold off
[h,p,c,stats]=ttest2(mean_abr_control,mean_abr_stroke)
[h,p,c,stats]=ttest2(mean_abr_control_corrected,mean_abr_stroke_corrected)

%% Calculating the gradients: scaling them and retrieving explained variances
% Gradients from mean matrices -- not used
for i=2
    [surf_lh, surf_rh] = load_conte69();
    labeling = load_parcellation('schaefer',400);
    conn_matrices = load_group_fc('schaefer',400);
    reference_gradient = GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
    m = (conn_matrices.schaefer_400);
    m=atanh(m);
    m(isinf(m))=1;
    reference_gradient = reference_gradient.fit(m);

    plot_hemispheres(reference_gradient.gradients{1}(:,1:3), {surf_lh,surf_rh}, ...
        'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', ' Gradient 3'});
    gradient_in_euclidean([reference_gradient.gradients{1}(:,[1 2])] ,{surf_lh,surf_rh},labeling.schaefer_400);
    scree_plot(reference_gradient.lambda{1});

    m=mean(corrmats_reduced_clean(:,:,1:num_control),3);
    before_correction=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
    before_correction=before_correction.fit(m,'reference',reference_gradient.gradients{1});
    plot_hemispheres(before_correction.aligned{1}(:,1:3), {surf_lh,surf_rh}, ...
        'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
    gradient_in_euclidean(before_correction.aligned{1}(:,[1 2]),{surf_lh,surf_rh},labeling.schaefer_400);
    scree_plot(before_correction.lambda{1});


    m=mean(corrmats_corrected_reduced(:,:,1:num_control),3,"omitnan");
    after_correction=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
    after_correction=after_correction.fit(m,'reference',reference_gradient.gradients{1});
    plot_hemispheres(after_correction.aligned{1}(:,1:3), {surf_lh,surf_rh}, ...
        'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
    gradient_in_euclidean(after_correction.aligned{1}(:,[1 2]),{surf_lh,surf_rh},labeling.schaefer_400);
    scree_plot(after_correction.lambda{1});

    m=mean(corrmats_reduced(:,:,num_control+1:end),3,"omitnan");
    before_correction=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
    before_correction=before_correction.fit(m,'reference',reference_gradient.gradients{1});
    plot_hemispheres(before_correction.aligned{1}(:,1:3), {surf_lh,surf_rh}, ...
        'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
    gradient_in_euclidean(before_correction.aligned{1}(:,[1 2]),{surf_lh,surf_rh},labeling.schaefer_400);
    scree_plot(before_correction.lambda{1});


    m=mean(corrmats_corrected_reduced(:,:,num_control+1:end),3,"omitnan");
    after_correction=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
    after_correction=after_correction.fit(m,'reference',reference_gradient.gradients{1});
    plot_hemispheres(after_correction.aligned{1}(:,1:3), {surf_lh,surf_rh}, ...
        'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
    gradient_in_euclidean(after_correction.aligned{1}(:,[1 2]),{surf_lh,surf_rh},labeling.schaefer_400);
    scree_plot(after_correction.lambda{1});
end


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
    before_correction=before_correction.fit(m,'reference',reference_gradient.gradients{1});
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
% Rescale the gradients
% for i=1:10
%     for j=1:size(subjects,1)
%         centered=before_aligned(:,i,j)-mean(before_aligned(:,i,j));
%         before_scaled(:,i,j)=centered./std(centered);
% 
%         centered=after_aligned(:,i,j)-mean(after_aligned(:,i,j));
%         after_scaled(:,i,j)=centered./std(centered);
%     end
% end

% Clean the gradients
% for i=1:10
%     for j=1:400
%         [B,BINT,R,RINT,STATS] =regress(squeeze(before_scaled(j,i,:)),[ones(130,1) subjects.age subjects.gender ]);
%         before_residuals(j,i,:)=R;
% 
%         [B,BINT,R,RINT,STATS] =regress(squeeze(after_scaled(j,i,:)),[ones(130,1) subjects.age subjects.gender ]);
%         after_residuals(j,i,:)=R;
%     end
% end
% before_scaled=before_residuals;
% after_scaled=after_residuals;
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
obj.colormaps(slanCM('viridis'))

plot_hemispheres(mean_after_control(:,1:3), {S,S2}, ...
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

        network_sd_control_before(i,j)=std(mean_before_control(seven_network_positions==i,j));
        network_sd_control_after(i,j)=std(mean_after_control(seven_network_positions==i,j));

        network_sd_stroke_before(i,j)=std(mean_before_stroke(seven_network_positions==i,j));
        network_sd_stroke_after(i,j)=std(mean_after_stroke(seven_network_positions==i,j));
    end
end

% draw the network breakdown of the mean gradients
for j=1:3
    figure;axes('XTick', 1:7, 'XTickLabel', network_names, 'FontSize', 16);hold on;ylabel('Mean Gradient Score')
    bar(network_means_control_before(:,j), 'FaceColor', 'flat', 'CData', colormap);
     errorbar(1:numel(network_means_control_before(:,j)), network_means_control_before(:,j), network_sd_control_before(:,j), network_sd_control_before(:,j), ...
         'Color', 'k', 'LineStyle', 'None');

    figure;axes('XTick', 1:7, 'XTickLabel', network_names, 'FontSize', 16);hold on;ylabel('Mean Gradient Score')
    bar(network_means_stroke_before(:,j), 'FaceColor', 'flat', 'CData', colormap);
    errorbar(1:numel(network_means_stroke_before(:,j)), network_means_stroke_before(:,j), network_sd_stroke_before(:,j), network_sd_stroke_before(:,j), ...
         'Color', 'k', 'LineStyle', 'None');


    figure;axes('XTick', 1:7, 'XTickLabel', network_names, 'FontSize', 16);hold on;ylabel('Mean Gradient Score')
    bar(network_means_control_after(:,j), 'FaceColor', 'flat', 'CData', colormap);
    errorbar(1:numel(network_means_control_after(:,j)), network_means_control_after(:,j), network_sd_control_after(:,j), network_sd_control_after(:,j), ...
        'Color', 'k', 'LineStyle', 'None');

    figure;axes('XTick', 1:7, 'XTickLabel', network_names, 'FontSize', 16);hold on;ylabel('Mean Gradient Score')
    bar(network_means_stroke_after(:,j), 'FaceColor', 'flat', 'CData', colormap);
    errorbar(1:numel(network_means_stroke_after(:,j)), network_means_stroke_after(:,j), network_sd_stroke_after(:,j), network_sd_stroke_after(:,j), ...
        'Color', 'k', 'LineStyle', 'None');

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
 % With regressors of no interest
 for i=1:400
     for  j=1:10
         [B,BINT,R,RINT,STATS] = regress(squeeze(before_aligned(i,j,:)),[ones(130,1) subjects.subj_type subjects.age subjects.gender ]);
         ps(i,j)=STATS(3);
         bs(i,j)=B(2);

         [B,BINT,R,RINT,STATS] = regress(squeeze(after_aligned(i,j,:)),[ones(130,1) subjects.subj_type subjects.age subjects.gender ]);
         ps2(i,j)=STATS(3);
         bs2(i,j)=B(2);
     end
 end
imagesc(fdr_bh(ps(:,1:3)).*bs(:,1:3))
figure
imagesc(fdr_bh(ps2(:,1:3)).*bs2(:,1:3))
ax=(fdr_bh(ps(:,1:3)).*bs(:,1:3));
ay=(fdr_bh(ps2(:,1:3)).*bs2(:,1:3));
plot_hemispheres(ax(:,1:3), {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
plot_hemispheres(ay(:,1:3), {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});


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

%% Max distance to zero, 3D variance 

% % Euclidian distance to zero 
% % for i=1:max(parcellation)
% %     for j=1:size(subjects,1)
% %         ED_tozero_before(i,j)=norm(before_aligned(i,1:3,j));
% %         ED_tozero_after(i,j)=norm(after_aligned(i,1:3,j));
% %     end
% % end
% %  [h,p,c,stats]=ttest2(mean(ED_tozero_before(:,1:num_control),2),mean(ED_tozero_before(:,num_control+1:end),2))
% %  [h,p,c,stats]=ttest2(mean(ED_tozero_after(:,1:num_control),2),mean(ED_tozero_after(:,num_control+1:end),2))
% 
% % for i=1:400 
% %     [h,p,c,stats]=ttest2(ED_tozero_before(i,1:num_control),ED_tozero_before(i,num_control+1:end));
% %     ps(i)=p;
% %     ts(i)=stats.tstat;
% % end
% plot_hemispheres([mean(ED_tozero_before(:,1:num_control),2) mean(ED_tozero_before(:,num_control+1:end),2)], {surf_lh,surf_rh}, ...
%     'parcellation',labeling.schaefer_400, 'labeltext',{'control','ED2zero - stroke'});
% plot_hemispheres([mean(ED_tozero_after(:,1:num_control),2) mean(ED_tozero_after(:,num_control+1:end),2)], {surf_lh,surf_rh}, ...
%     'parcellation',labeling.schaefer_400, 'labeltext',{'control','ED2zero - stroke'});
% Distance to centroid
% mean_x = mean(before_scaled(:,1,1:num_control),3);
% mean_y = mean(before_scaled(:,2,1:num_control),3);
% mean_z = mean(before_scaled(:,3,1:num_control),3);
% variance_x = var(before_scaled(:,1,1:num_control),[],3);
% variance_y = var(before_scaled(:,2,1:num_control),[],3);
% variance_z = var(before_scaled(:,3,1:num_control),[],3);centroid = mean(mean(before_scaled(:,1:3,1:num_control),3));
% distances = sqrt(sum((before_scaled(:,1,1:num_control) - centroid).^2, 2));
% overall_variance = squeeze(mean(distances.^2));
% 
% mean_x = mean(before_scaled(:,1,num_control+1:end),3);
% mean_y = mean(before_scaled(:,2,num_control+1:end),3);
% mean_z = mean(before_scaled(:,3,num_control+1:end),3);
% variance_x = var(before_scaled(:,1,num_control+1:end),[],3);
% variance_y = var(before_scaled(:,2,num_control+1:end),[],3);
% variance_z = var(before_scaled(:,3,num_control+1:end),[],3);centroid = mean(mean(before_scaled(:,1:3,num_control+1:end),3));
% distances2 = sqrt(sum((before_scaled(:,1,num_control+1:end) - centroid).^2, 2));
% overall_variance2 = squeeze(mean(distances2.^2));
% [h,p,c,stats]=ttest2(overall_variance,overall_variance2)
% for i=1:400 
%     [h,p,c,stats]=ttest2(squeeze(distances(i,1,:)),squeeze(distances2(i,1,:)));
%     ps(i)=p;
%     ts(i)=stats.tstat;
% end
% 
% sum(fdr_bh(ps))

% ED to control in 3D 
for i=1:max(parcellation)
    for j=1:size(subjects,1)
        EDtocontrol_before(i,j)=norm(mean(before_aligned(i,1:3,1:num_control),3)-before_aligned(i,1:3,j));
        EDtocontrol_after(i,j)=norm(mean(after_aligned(i,1:3,1:num_control),3)-after_aligned(i,1:3,j));
    end
end

[h,p,c,stats]=ttest2(mean(EDtocontrol_before(:,1:num_control),2),mean(EDtocontrol_before(:,num_control+1:end),2))

[h,p,c,stats]=ttest2(mean(EDtocontrol_after(:,1:num_control),2),mean(EDtocontrol_after(:,num_control+1:end),2))
hold on;histogram(mean(EDtocontrol_before(:,1:num_control),2));histogram(mean(EDtocontrol_before(:,num_control+1:end),2));hold off
figure
hold on;histogram(mean(EDtocontrol_after(:,1:num_control),2));histogram(mean(EDtocontrol_after(:,num_control+1:end),2));hold off
xlabel('ED func')
ylabel('Frequency')
legend('Control','Stroke')
dif1=mean(EDtocontrol_before(:,1:num_control),2)-mean(EDtocontrol_before(:,num_control+1:end),2);
dif2=mean(EDtocontrol_after(:,1:num_control),2)-mean(EDtocontrol_after(:,num_control+1:end),2);
[h,p,c,stats]=ttest(dif1-dif2)

% dif1=mean(ED_tozero_before(:,1:num_control),2)-mean(ED_tozero_before(:,num_control+1:end),2);
% dif2=mean(ED_tozero_after(:,1:num_control),2)-mean(ED_tozero_after(:,num_control+1:end),2);
% [h,p,c,stats]=ttest(dif1-dif2)

plot_hemispheres([mean(EDtocontrol_before(:,1:num_control),2) mean(EDtocontrol_before(:,num_control+1:end),2) dif1.*(dif1<-0.01)], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'control','stroke', 'ED2control - diff.'});

% for i=1:400 not used-ROI-based
%     [h,p,c,stats]=ttest2(EDtocontrol_before(i,1:num_control),EDtocontrol_before(i,num_control+1:end));
%     ps(i)=p;
%     ts(i)=stats.tstat;
% end

dd=(dif1);
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

     figure;axes('XTick', 1:7, 'XTickLabel', network_names, 'FontSize', 16);hold on;ylabel('Mean ED to control - after')
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

plot_hemispheres([dif11 dif12.*(dif12<-0.006) dif13], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Dif 1st ','Dif 2nd', 'Dif3rd'});
gradient_in_euclidean(squeeze(mean(EDtocontrol_before_sep(:,1:num_control,1:3),2)) ,{surf_lh,surf_rh},labeling.schaefer_400)
gradient_in_euclidean(squeeze(mean(EDtocontrol_before_sep(:,num_control+1:end,1:3),2)) ,{surf_lh,surf_rh},labeling.schaefer_400)
gradient_in_euclidean(squeeze(mean(EDtocontrol_after_sep(:,1:num_control,1:2),2)) ,{surf_lh,surf_rh},labeling.schaefer_400)
gradient_in_euclidean(squeeze(mean(EDtocontrol_after_sep(:,num_control+1:end,1:2),2)) ,{surf_lh,surf_rh},labeling.schaefer_400)

% 
% [h,p,c,stats]=ttest2(mean(EDtozero_before_sep(:,1:num_control,1),2),mean(EDtozero_before_sep(:,num_control+1:end,1),2))
% [h,p,c,stats]=ttest2(mean(EDtozero_before_sep(:,1:num_control,2),2),mean(EDtozero_before_sep(:,num_control+1:end,2),2))
% [h,p,c,stats]=ttest2(mean(EDtozero_before_sep(:,1:num_control,3),2),mean(EDtozero_before_sep(:,num_control+1:end,3),2))
% 
% 
% [h,p,c,stats]=ttest2(mean(EDtozero_after_sep(:,1:num_control,1),2),mean(EDtozero_after_sep(:,num_control+1:end,1),2))
% [h,p,c,stats]=ttest2(mean(EDtozero_after_sep(:,1:num_control,2),2),mean(EDtozero_after_sep(:,num_control+1:end,2),2))
% [h,p,c,stats]=ttest2(mean(EDtozero_after_sep(:,1:num_control,3),2),mean(EDtozero_after_sep(:,num_control+1:end,3),2))
% 
% dif11=mean(EDtozero_before_sep(:,1:num_control,1),2)-mean(EDtozero_before_sep(:,num_control+1:end,1),2)
% dif12=mean(EDtozero_before_sep(:,1:num_control,2),2)-mean(EDtozero_before_sep(:,num_control+1:end,2),2)
% dif13=mean(EDtozero_before_sep(:,1:num_control,3),2)-mean(EDtozero_before_sep(:,num_control+1:end,3),2)
% % 
% plot_hemispheres([dif11 dif12 dif13], {surf_lh,surf_rh}, ...
%     'parcellation',labeling.schaefer_400, 'labeltext',{'Dif 1st ','Dif 2nd', 'Dif3rd'});
% gradient_in_euclidean(squeeze(mean(EDtozero_before_sep(:,1:num_control,1:2),2)) ,{surf_lh,surf_rh},labeling.schaefer_400)
% gradient_in_euclidean(squeeze(mean(EDtozero_before_sep(:,num_control+1:end,1:2),2)) ,{surf_lh,surf_rh},labeling.schaefer_400)
% gradient_in_euclidean(squeeze(mean(EDtozero_after_sep(:,1:num_control,1:2),2)) ,{surf_lh,surf_rh},labeling.schaefer_400)
% gradient_in_euclidean(squeeze(mean(EDtozero_after_sep(:,num_control+1:end,1:2),2)) ,{surf_lh,surf_rh},labeling.schaefer_400)


dif21=mean(EDtocontrol_after_sep(:,1:num_control,1),2)-mean(EDtocontrol_after_sep(:,num_control+1:end,1),2);
dif22=mean(EDtocontrol_after_sep(:,1:num_control,2),2)-mean(EDtocontrol_after_sep(:,num_control+1:end,2),2);
dif23=mean(EDtocontrol_after_sep(:,1:num_control,3),2)-mean(EDtocontrol_after_sep(:,num_control+1:end,3),2);
plot_hemispheres([dif21 dif22 dif23], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Dif 1st ','Dif 2nd', 'Dif3rd'});

[h,p,c,stats]=ttest(dif11-dif21)
[h,p,c,stats]=ttest(dif12-dif22)
[h,p,c,stats]=ttest(dif13-dif23)

% plot_hemispheres([dif11-dif21 dif12-dif22 dif13-dif23], {surf_lh,surf_rh}, ...
%     'parcellation',labeling.schaefer_400, 'labeltext',{'Dif 1st ','Dif 2nd', 'Dif3rd'});
% 
% dif21=mean(EDtozero_after_sep(:,1:num_control,1),2)-mean(EDtozero_after_sep(:,num_control+1:end,1),2);
% dif22=mean(EDtozero_after_sep(:,1:num_control,2),2)-mean(EDtozero_after_sep(:,num_control+1:end,2),2);
% dif23=mean(EDtozero_after_sep(:,1:num_control,3),2)-mean(EDtozero_after_sep(:,num_control+1:end,3),2);


%% Correlation between difs, gradients, lags

% Corr with mean lag and gradients
corr(mean(mean_lag(:,1:num_control),2),mean(squeeze(before_aligned(:,3,:)),2)) % 1: 0.25 2: 0.7 3: 0.19



% Individual lag corrs with gradients or diffs from both side


lag_diff=mean(mean_lag(:,1:num_control),2)-mean_lag;
% corr(mean(mean_lag,2),mean(squeeze(after_aligned(:,2,1:num_control)),2))
% corr(mean(mean_lag,2),mean(squeeze(after_aligned(:,2,num_control+1:end)),2))
corr(mean(lag_diff(:,1:num_control),2),mean(EDtocontrol_before(:,1:num_control),2))
corr((lag_diff(:,num_control+1:end))',(EDtocontrol_before(:,num_control+1:end))')
corr(mean(lag_diff(:,1:num_control),2),mean(EDtocontrol_after(:,1:num_control),2))
corr(mean(lag_diff(:,num_control+1:end),2),mean(EDtocontrol_after(:,num_control+1:end),2)) % lag diff is higher for stroke 
for i=1:130
    corrs(i,1)=corr(lag_diff(:,i),EDtocontrol_before(:,i))
    corrs(i,2)=corr(lag_diff(:,i),EDtocontrol_after(:,i))

    corrs(i,3)=corr(lag_diff(:,i),EDtocontrol_before(:,i))
    corrs(i,4)=corr(lag_diff(:,i),EDtocontrol_after(:,i))
end

% for i=1:130
%     corrs(i,1)=corr(lag_diff(:,i),before_aligned(:,1,i))
%     corrs(i,2)=corr(lag_diff(:,i),after_aligned(:,1,i))
% 
%     corrs(i,3)=corr(lag_diff(:,i),before_aligned(:,2,i))
%     corrs(i,4)=corr(lag_diff(:,i),after_aligned(:,2,i))
% end
% 
% ix=~isnan(subjects.nihss_hospital);
% [bb(:,1) bb(:,2)]=corr(EDtocontrol_before(:,ix)',subjects.nihss_hospital(ix));
% [bc(:,1) bc(:,2)]=corr(EDtocontrol_after(:,ix)',subjects.nihss_hospital(ix));
% % dd=corr((lag_diff(:,ix))',subjects.nihss_hospital(ix));
% plot_hemispheres([mean(lag_diff(:,num_control+1:end),2) mean(EDtocontrol_before(:,num_control+1:end),2)], {surf_lh,surf_rh}, ...
% 'parcellation',labeling.schaefer_400, 'labeltext',{'2controlb','2controla'});

% Ed func with behaviorals
ix=~isnan(behavioral.motorl_f);
ix=ix+[zeros(num_control,1);ones(num_stroke,1)];
% ix=ix+[ones(num_control,1);zeros(num_stroke,1);];

% ix=ix+(subjects.lesion_side==1);
ix=ix==2;
sum(ix)
[bb(:,1) bb(:,2)]=corr(EDtocontrol_before(:,ix)',behavioral.motorl_f(ix));
[bc(:,1) bc(:,2)]=corr(EDtocontrol_after(:,ix)',behavioral.motorl_f(ix));
plot_hemispheres([bb(:,1).*(bb(:,2)<0.01) bc(:,1).*(bc(:,2)<0.01) ], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400);


% Lag with behaviorals
ix=~isnan(behavioral.language_f);
ix=ix+[zeros(num_control,1);ones(num_stroke,1)];
% ix=ix+[ones(num_control,1);zeros(num_stroke,1);];
% ix=ix+(subjects.lesion_side==1);
ix=ix==2;
sum(ix)
[dd(:,1) dd(:,2)]=corr((mean_lag(:,ix))',behavioral.language_f(ix));
ix=~isnan(behavioral.motorl_f);
ix=ix+[zeros(num_control,1);ones(num_stroke,1)];
% ix=ix+[ones(num_control,1);zeros(num_stroke,1);];
% ix=ix+(subjects.lesion_side==1);
ix=ix==2;
sum(ix)
[dc(:,1) dc(:,2)]=corr((mean_lag(:,ix))',behavioral.motorl_f(ix));
ix=~isnan(behavioral.motorr_f);
ix=ix+[zeros(num_control,1);ones(num_stroke,1)];
% ix=ix+[ones(num_control,1);zeros(num_stroke,1);];
% ix=ix+(subjects.lesion_side==1);
ix=ix==2;
sum(ix)
[de(:,1) de(:,2)]=corr((mean_lag(:,ix))',behavioral.motorr_f(ix));

plot_hemispheres([dd(:,1).*(dd(:,2)<0.05) dc(:,1).*(dc(:,2)<0.05) de(:,1).*(de(:,2)<0.05)], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400);





% 
% [correlation, feature] = meta_analytic_decoder(parcel2full(abs(bb(:,1).*(bb(:,2)<0.05)), schaefer_400), ...
%     'template', 'fsaverage5');
% wc = wordcloud(feature(correlation>0), correlation(correlation>0));





%% Anatomical distance to lesion 

mask_coords=table2array(readtable('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/distance2lesion/empty_mask.txt')); % load the template


masklist=readtable('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/masklist.txt','ReadVariableNames', false); % load the subjects list that have masks 

% Create distance files in nifti
% for i=1:size(masklist,1)
% 
%    lesion_mask_str='distance2lesion/SUBID_lesion_coords.txt';
%    splt=split(masklist{i,2},'_');
%    lesion_mask_str=string(strrep(lesion_mask_str,'SUBID',splt(1)));
%    lesion_mask=table2array(readtable(lesion_mask_str),'ReadVariableNames', false); % load the lesion mask of the subject 
%    tosave_str='distance2lesion/SUBID_distance2lesion.txt';
% 
%    sub_distance=mask_coords;
% 
%    for j=1:size(sub_distance,1) 
%        current_coor=sub_distance(j,1:3); % go over each voxel of the template 
%        for k=1:size(lesion_mask,1)
%            distances2lesion(k)=norm(current_coor-lesion_mask(k,1:3));
%        end
%        sub_distance(j,4)=min(distances2lesion);
%        sub_distance(j,1)=sub_distance(j,1)*-1;
%        sub_distance(j,2)=sub_distance(j,2)*-1;
%        clear distances2lesion
%    end
%     writematrix(sub_distance,string(strrep(tosave_str,'SUBID',splt(1))), 'delimiter', ' ')
% end
%        
    
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

for i=1:400
    between_sub_corrs_before(i)=corr(ED_func_before(i,:)',ED_anat(i,:)');
    between_sub_corrs_after(i)=corr(ED_func_after(i,:)',ED_anat(i,:)');
    between_sub_corrs_lag(i)=corr(lag_diff_subset(i,:)',ED_anat(i,:)');
end
plot_hemispheres([between_sub_corrs_before' between_sub_corrs_after' between_sub_corrs_lag'], {surf_lh,surf_rh}, ...
'parcellation',labeling.schaefer_400);

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

% ix=~isnan(behavioral_subset.motorl_f);
%  ix=ix+(subjects_subset.lesion_side==0);
% ix=ix==2;
% sum(ix)
% for i=1:400
%     [temp_a(i,1) temp_a(i,2)]=corr(ED_func_before(i,ix)',behavioral_subset.motorl_f(ix));
%     [temp_b(i,1) temp_b(i,2)]=corr(ED_func_after(i,ix)',behavioral_subset.motorl_f(ix));
%     [temp_c(i,1) temp_c(i,2)]=corr(lag_diff(i,ix)',behavioral_subset.motorl_f(ix));
%     [temp_d(i,1) temp_d(i,2)]=corr(ED_anat(i,ix)',behavioral_subset.motorl_f(ix));
% end
% 
% a=(temp_a(:,2)<0.05).*temp_a(:,1);
% b=(temp_b(:,2)<0.05).*temp_b(:,1);
% c=(temp_c(:,2)<0.05).*temp_c(:,1);
% d=(temp_d(:,2)<0.05).*temp_d(:,1);
% 
% % a=fdr_bh(temp_a(:,2)).*temp_a(:,1);
% % b=fdr_bh(temp_b(:,2)).*temp_b(:,1);
% % c=fdr_bh(temp_c(:,2)).*temp_c(:,1);
% % d=fdr_bh(temp_d(:,2)).*temp_d(:,1);
% 
% 
% % a=temp_a(:,1);
% % a=a.*(abs(a)>0.3);
% % 
% % b=temp_b(:,1);
% % b=b.*(abs(b)>0.3);
% % 
% % c=temp_c(:,1);
% % c=c.*(abs(c)>0.3);
% % 
% % d=temp_d(:,1);
% % d=d.*(abs(d)>0.3);
% plot_hemispheres([a b c d], {surf_lh,surf_rh}, ...
%     'parcellation',labeling.schaefer_400);
% 
% 
% 
% 
% 
% 










%% Clustering with hctsa - canceled
for j=1
for i=1
    matrix = mean_lag';
    timeSeriesData = cell(size(matrix, 1), 1);
    for j = 1:size(matrix, 1)
        timeSeriesData{j} = matrix(j, :)';
    end

    labels=subjects.participant_id;
    controlEntries = repmat({'control'}, 24, 1);
    strokeEntries = repmat({'stroke'}, 106, 1);
    keywords = [controlEntries; strokeEntries];

    save('INP_test.mat','timeSeriesData','labels','keywords');
    TS_Init('INP_test.mat','catch22');close
    TS_Compute(false)
    TS_LabelGroups('raw',{'control','stroke'});
    TS_InspectQuality('summary')
    TS_Normalize('mixedSigmoid',[0.8,1.0]);
    distanceMetricRow = 'euclidean'; % time-series feature distance
    linkageMethodRow = 'average'; % linkage method
    distanceMetricCol = 'corr_fast'; % a (poor) approximation of correlations with NaNs
    accuracies(i,1)=TS_Classify('norm');

    matrix =before_scaled(:,ix)';
    timeSeriesData = cell(size(matrix, 1), 1);
    for j = 1:size(matrix, 1)

        timeSeriesData{j} = matrix(j, :)';
    end

    labels=(subjects.participant_id(ix));
    cellArray = cell(size(labels)); % Preallocate cell array

    % Logical indexing to assign 'left' or 'right' based on the values
    cellArray(subjects.lesion_side(ix) == 0) = {'a'};
    cellArray(subjects.lesion_site(ix) == 1) = {'b'};
    cellArray(subjects.lesion_site(ix) == 2) = {'c'};
    cellArray(subjects.lesion_site(ix) == 3) = {'d'};
    cellArray(subjects.lesion_site(ix) == 4) = {'e'};
    cellArray(subjects.lesion_site(ix) == 5) = {'f'};
    cellArray(subjects.lesion_site(ix) == 6) = {'g'};
    keywords = cellArray;

    save('INP_test.mat','timeSeriesData','labels','keywords');
    TS_Init('INP_test.mat','catch22');close
    TS_Compute(false)
    TS_LabelGroups('raw',{'a','b','c','d','e','f','g'});
    TS_InspectQuality('summary')
    TS_Normalize('mixedSigmoid',[0.8,1.0]);
    distanceMetricRow = 'euclidean'; % time-series feature distance
    linkageMethodRow = 'average'; % linkage method
    distanceMetricCol = 'corr_fast'; % a (poor) approximation of correlations with NaNs
    accuracies(i,2)=TS_Classify('norm');
end
TS_PlotDataMatrix('colorGroups',true)




ix=~isnan(subjects.lesion_site)

subjects.lesion_site(ix)

[idx,C] =kmeans((squeeze(before_scaled(:,1,:))),2)


heatmap(confusionmat(subjects.subj_type,idx-1))


% Cluster after correlation as in Striem-Amit 
stroke_subjects=subjects(num_control+1:end,:);
% tocluster=squeeze(EDtocontrol_after_sep(:,num_control+1:end,2))
tocluster=corr(EDtocontrol_after(:,num_control+1:end))'
tree = linkage(tocluster','average');
[~,T] = dendrogram(tree,10);

ix1=logical(sum(T==[4 8],2));
cl1=mean(tocluster(:,ix1),2);
ix2=logical(sum(T==[1 10],2));
cl2=mean(tocluster(:,ix2),2);
ix3=logical(sum(T==[7 9],2));
cl3=mean(tocluster(:,ix3),2);
cl4=mean(tocluster(:,find(T==4)),2)
cl5=mean(tocluster(:,find(T==5)),2)
cl6=mean(tocluster(:,find(T==6)),2)
cl7=mean(tocluster(:,find(T==7)),2)
cl8=mean(tocluster(:,find(T==8)),2)
cl9=mean(tocluster(:,find(T==9)),2)
cl10=mean(tocluster(:,find(T==10)),2)
plot_hemispheres([cl1 cl2 cl3], {surf_lh,surf_rh}, ...
'parcellation',labeling.schaefer_400, 'labeltext',{'1', '2', '3'});

for i=1:400
[B,BINT,R,RINT,STATS] = regress(EDtocontrol_before(i,:)',[ones(130,1) subjects.subj_type subjects.age subjects.gender]);
ps(i)=STATS(3);
bs(i)=B(2);

[B,BINT,R,RINT,STATS] = regress(EDtocontrol_after(i,:)',[ones(130,1) subjects.subj_type subjects.age subjects.gender]);
ps2(i)=STATS(3);
bs2(i)=B(2);
end
plot_hemispheres([(fdr_bh(ps).*bs)' (fdr_bh(ps2).*bs2)'], {surf_lh,surf_rh}, ...
'parcellation',labeling.schaefer_400, 'labeltext',{'1', '2'});
end









%% Group-level analysis -- The following code is just trial-errors, curated code for group-level is in other files
%

for i=1:size(corrmats,1)
    for j=1:10
        [h,p,c,stats]=ttest(squeeze(before_scaled(i,j,1:num_control)),squeeze(after_scaled(i,j,1:num_control)));
        control_p(i,j)=p;
        control_t(i,j)=stats.tstat;
        [h,p,c,stats]=ttest(squeeze(before_scaled(i,j,num_control+1:end)),squeeze(after_scaled(i,j,num_control+1:end)));
        stroke_p(i,j)=p;
        stroke_t(i,j)=stats.tstat;

        [h,p,c,stats]=ttest2(squeeze(before_scaled(i,j,1:num_control)),squeeze(before_scaled(i,j,num_control+1:end)));
        controlXstroke_before_p(i,j)=p;
        controlXstroke_before_t(i,j)=stats.tstat;
        [h,p,c,stats]=ttest2(squeeze(after_scaled(i,j,1:num_control)),squeeze(after_scaled(i,j,num_control+1:end)));
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

% Euclidian distance
for i=1:max(parcellation)
    for j=1:size(subjects,1)
        distances(i,j)=norm(mean(after_aligned(i,1:3,1:num_control),3)-after_aligned(i,1:3,j));
    end
end

plot_hemispheres([var(distances(:,1:num_control),[],2) var(distances(:,num_control+1:end),[],2)], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Distance','Distance'});
% Lags
lags=load_nii('delays_custom/delays.nii.gz');
lags=reshape(squeeze((lags.img)),[],875);
lags=lags(parcellation_mask_whole,:);

% take the average of each subject
positions=zeros(size(corrmats,3),1);
corrmats_reduced=zeros(max(parcellation),max(parcellation),size(subjects,1));
corrmats_corrected_reduced=zeros(max(parcellation),max(parcellation),size(subjects,1));
for i=1:size(subjects,1)
    sub_id=subjects.participant_id(i);
    for j=1:size(corrmats,3)
        positions(j)=contains(func_list(j,:),sub_id);
    end
    lags_reduced(:,i)=mean(lags(:,logical(positions)),2);
end

% Mean lag
for j=1:max(parcellation_reduced)
    for i=1:size(subjects,1)
        mean_lag(j,i)=mean(lags_reduced(parcellation_reduced==j,i));
    end
end
plot_hemispheres([mean(mean_lag(:,1:num_control),2)*-1 mean(mean_lag(:,num_control+1:end),2)*-1], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Mean Lag','Mean Lag'});
% Correlation with lag and gradient value
for i=1:size(subjects,1)
    for j=1:10
        lagXgra_before(i,j)=corr(squeeze(before_scaled(:,j,i)),mean_lag(:,i)*-1);
        lagXgra_after(i,j)=corr(squeeze(after_scaled(:,j,i)),mean_lag(:,i)*-1);
    end
end
hold on;histogram(lagXgra_before(27:end,4));histogram(lagXgra_before(1:26,4));hold off
figure
hold on;histogram(lagXgra_after(27:end,3));histogram(lagXgra_after(1:26,3));hold off
[h,p,c,stats]=ttest2(lagXgra_before(1:26,1),lagXgra_before(27:end,1))
[h,p,c,stats]=ttest2(lagXgra_after(1:26,3),lagXgra_after(27:end,3))

[h,p,c,stats]=ttest(lagXgra_before(1:26,2),lagXgra_after(1:26,2))
[h,p,c,stats]=ttest(lagXgra_before(27:end,3),lagXgra_after(27:end,3))

%
for i=1:size(subjects,1)
    for j=1:10
        dif=before_scaled(:,j,i)-after_scaled(:,j,i);
        lagXgra_before(i,j)=corr(dif,mean(mean_lag(:,1:26),2)*-1,"Type","Spearman");
        %         dif=mean(after_scaled(:,j,1:26),3)-after_scaled(:,j,i);
        %         lagXgra_after(i,j)=corr(dif,mean_lag(:,i)*-1);
    end
end
corr((mean(before_scaled(:,:,1:26),3)-mean(after_scaled(:,:,1:26),3)),mean(mean_lag(:,1:26),2))
corr((mean(before_scaled(:,:,27:end),3)-mean(after_scaled(:,:,27:end),3)),mean(mean_lag(:,27:end),2))
corr((mean(before_scaled,3)-mean(after_scaled,3)),mean(mean_lag,2))

% Node-wise correlation instead of whole-brain

for i=1:max(parcellation_reduced)
    for j=1:10
        tempcors_c(i,j)=corr(squeeze(before_scaled(i,j,1:26)),mean_lag(i,1:26)');
        tempcors_s(i,j)=corr(squeeze(after_scaled(i,j,27:end)),mean_lag(i,27:end)');
    end
end
%% TODO
% Count total nans in each mean matrix


% Intrahemispheric connectivity
load('corrmats_temp.mat')

% whole brain
m = (conn_matrices.schaefer_400);
corrmats_right=atanh(m(201:end,201:end));
corrmats_left=atanh(m(1:200,1:200));
ref_whole = GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
m(isinf(m))=1;

ref_whole = ref_whole.fit(m);
plot_hemispheres(ref_whole.gradients{1}(:,1:3), {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400);
gradient_in_euclidean([ref_whole.gradients{1}(:,[1 3])] ,{surf_lh,surf_rh},labeling.schaefer_400);


% seperate hemispheres

ref_lh = GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
corrmats_left(isinf(corrmats_left))=1;
ref_lh = ref_lh.fit(corrmats_left, 'reference', ref_whole.gradients{1}(1:200,:));
plot_hemispheres([ref_lh.aligned{1}(:,1:2);ref_lh.aligned{1}(:,1:2)] , {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400);
gradient_in_euclidean([ref_lh.aligned{1}(:,1:2);ref_lh.aligned{1}(:,1:2)],{surf_lh,surf_rh},labeling.schaefer_400);

ref_rh = GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
corrmats_right(isinf(corrmats_right))=1;
ref_rh = ref_rh.fit(corrmats_right, 'reference', ref_whole.gradients{1}(201:400,:));
plot_hemispheres([ref_rh.aligned{1}(:,1:3);ref_rh.aligned{1}(:,1:3)] , {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400);
gradient_in_euclidean([ref_rh.aligned{1}(:,[1 3]);ref_rh.aligned{1}(:,[1 3])],{surf_lh,surf_rh},labeling.schaefer_400);


% take the average of each subject
positions=zeros(size(corrmats,3),1);
corrmats_reduced=zeros(max(parcellation),max(parcellation),size(subjects,1));
for i=1:size(subjects,1)
    sub_id=subjects.participant_id(i);
    for j=1:size(corrmats,3)
        positions(j)=contains(func_list(j,:),sub_id);
    end
    corrmats_reduced(:,:,i)=atanh(mean(corrmats(:,:,logical(positions)),3));
end
corrmats_reduced(isinf(corrmats_reduced))=1;

% % Regress out the confounds
% confounds=[subjects.age subjects.gender subjects.handed];
% for i=1:size(corrmats,1)
%     for j=1:size(corrmats,1)
%         [~,~,residuals] = regress(squeeze(corrmats_reduced(i,j,:)),[ones(size(confounds,1),1) confounds]);
%         corrmats_reduced_clean(i,j,:)=residuals;
%     end
% end

% Get the gradients for each subject and do the group analysis
for i=1:size(corrmats_reduced,3)
    disp(i)
    m=corrmats_reduced(:,:,i);
    m(isnan(m))=1;
    whole_brain=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
    whole_brain=whole_brain.fit(m,'reference',ref_whole.gradients{1});
    whole_brain_gradients{i}=whole_brain;

    whole_brain_gradients_aligned(:,:,i)=whole_brain_gradients{i}.aligned{1};
    before_lambda(:,i)=whole_brain_gradients{i}.lambda{1}./sum(whole_brain_gradients{i}.lambda{1});
end
% Rescale the gradients
% for i=1:10
%     for j=1:size(subjects,1)
%         centered=whole_brain_gradients_aligned(:,i,j)-mean(whole_brain_gradients_aligned(:,i,j));
%         whole_brain_gradients_aligned_scaled(:,i,j)=centered./std(centered);
%     end
% end

wb_control=mean(whole_brain_gradients_aligned(:,:,1:num_control),3);
wb_stroke=mean(whole_brain_gradients_aligned(:,:,num_control+1:end),3);
plot_hemispheres(wb_control(:,[1 2 3]), {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
plot_hemispheres(wb_stroke(:,1:3), {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
gradient_in_euclidean(wb_control(:,[1 3]),{surf_lh,surf_rh},labeling.schaefer_400);
gradient_in_euclidean(wb_stroke(:,[1 3]),{surf_lh,surf_rh},labeling.schaefer_400);
scree_plot(mean(before_lambda(:,1:num_control),2))
scree_plot(mean(before_lambda(:,num_control+1:end),2))

% Do it for hemispheres
% Get the gradients for each subject and do the group analysis
for i=1:size(corrmats_reduced,3)
    disp(i)
    m=corrmats_reduced(1:200,1:200,i);
    m(isnan(m))=1;
    left_brain=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
    left_brain=left_brain.fit(m,'reference',ref_lh.gradients{1});
    left_brain_gradients{i}=left_brain;

    left_brain_gradients_aligned(:,:,i)=left_brain_gradients{i}.aligned{1};
    left_lambda(:,i)=left_brain_gradients{i}.lambda{1}./sum(left_brain_gradients{i}.lambda{1});

    m=corrmats_reduced(201:end,201:end,i);
    m(isnan(m))=1;
    right_brain=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
    right_brain=right_brain.fit(m,'reference',ref_rh.gradients{1});
    right_brain_gradients{i}=right_brain;

    right_brain_gradients_aligned(:,:,i)=right_brain_gradients{i}.aligned{1};
    right_lambda(:,i)=right_brain_gradients{i}.lambda{1}./sum(right_brain_gradients{i}.lambda{1});
end
% Rescale the gradients
% for i=1:10
%     for j=1:size(subjects,1)
%         centered=left_brain_gradients_aligned(:,i,j)-mean(left_brain_gradients_aligned(:,i,j));
%         left_brain_gradients_aligned_scaled(:,i,j)=centered./std(centered);
%         centered=right_brain_gradients_aligned(:,i,j)-mean(right_brain_gradients_aligned(:,i,j));
%         right_brain_gradients_aligned_scaled(:,i,j)=centered./std(centered);
%     end
% end
rh_control=mean(right_brain_gradients_aligned(:,:,1:num_control),3);
lh_control=mean(left_brain_gradients_aligned(:,:,1:num_control),3);
plot_hemispheres([rh_control(:,[1:3]);rh_control(:,[1:3])], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
plot_hemispheres([lh_control(:,[1:3]);lh_control(:,[1:3])], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
gradient_in_euclidean([rh_control(:,[1 3]);rh_control(:,[1 3])],{surf_lh,surf_rh},labeling.schaefer_400);
gradient_in_euclidean([lh_control(:,[1 3]);lh_control(:,[1 3])],{surf_lh,surf_rh},labeling.schaefer_400);

% Choose lesion site index
rh_stroke=mean(right_brain_gradients_aligned(:,:,subjects.lesion_side==0),3); % zero is for left
lh_stroke=mean(left_brain_gradients_aligned(:,:,subjects.lesion_side==1),3); % one is for right
plot_hemispheres([rh_stroke(:,[1:3]);rh_stroke(:,[1:3])], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
plot_hemispheres([lh_stroke(:,[1:3]);lh_stroke(:,[1:3])], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
gradient_in_euclidean([rh_stroke(:,[1 3]);rh_stroke(:,[1 3])],{surf_lh,surf_rh},labeling.schaefer_400);
gradient_in_euclidean([lh_stroke(:,[1 3]);lh_stroke(:,[1 3])],{surf_lh,surf_rh},labeling.schaefer_400);

% Group-level analysis

for i=1:size(corrmats,1)/2
    for j=1:10

        [h,p,c,stats]=ttest2(squeeze(right_brain_gradients_aligned_scaled(i,j,1:num_control)),squeeze(right_brain_gradients_aligned_scaled(i,j,subjects.lesion_side==0)));
        rh_p(i,j)=p;
        rh_t(i,j)=stats.tstat;
        [h,p,c,stats]=ttest2(squeeze(left_brain_gradients_aligned_scaled(i,j,1:num_control)),squeeze(left_brain_gradients_aligned_scaled(i,j,subjects.lesion_side==0)));
        lh_p(i,j)=p;
        lh_t(i,j)=stats.tstat;
    end
end
figure
imagesc(fdr_bh(rh_p).*rh_t)
figure
imagesc(fdr_bh(lh_p).*lh_t)
figure


fdr_mask=fdr_bh(stroke_p(:,1:3)).*stroke_t(:,1:3);
plot_hemispheres(fdr_mask(:,1:3), {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});



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

%% Work on one gradient

% take the average of each subject
positions=zeros(size(corrmats,3),1);
corrmats_reduced=zeros(max(parcellation),max(parcellation),size(subjects,1));
for i=1:size(subjects,1)
    sub_id=subjects.participant_id(i);
    for j=1:size(corrmats,3)
        positions(j)=contains(func_list(j,:),sub_id);
    end
    corrmats_reduced(:,:,i)=atanh(mean(corrmats(:,:,logical(positions)),3));
end
corrmats_reduced(isinf(corrmats_reduced))=1;

corrmats_control_whole=corrmats_reduced(:,:,1:num_control);
corrmats_stroke_whole=corrmats_reduced(:,:,num_control+1:end);

corrmats_control_right=(corrmats_reduced(201:end,201:end,1:num_control));
corrmats_control_left=(corrmats_reduced(1:200,1:200,num_control+1:end));

corrmats_stroke_right=(corrmats_reduced(201:end,201:end,1:num_control));
corrmats_stroke_left=(corrmats_reduced(1:200,1:200,num_control+1:end));

% whole brain
control_whole = GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
control_whole = control_whole.fit(mean(corrmats_control_whole,3),'reference',reference_gradient.gradients{1});
gradient_in_euclidean([control_whole.aligned{1}(:,[1 2])] ,{surf_lh,surf_rh},labeling.schaefer_400);

stroke_whole = GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
stroke_whole = stroke_whole.fit(mean(corrmats_stroke_whole,3,"omitnan"),'reference',reference_gradient.gradients{1});
gradient_in_euclidean([stroke_whole.aligned{1}(:,[1 2])] ,{surf_lh,surf_rh},labeling.schaefer_400);


% hemispheres
ref_lh = GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
corrmats_left(isinf(corrmats_left))=1;
ref_lh = ref_lh.fit(corrmats_left, 'reference', reference_gradient.gradients{1}(1:200,:));
plot_hemispheres([ref_lh.aligned{1}(:,1:2);ref_lh.aligned{1}(:,1:2)] , {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400);
gradient_in_euclidean([ref_lh.aligned{1}(:,1:2);ref_lh.aligned{1}(:,1:2)],{surf_lh,surf_rh},labeling.schaefer_400);

ref_rh = GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
corrmats_right(isinf(corrmats_right))=1;
ref_rh = ref_rh.fit(corrmats_right, 'reference', reference_gradient.gradients{1}(201:end,:));
gradient_in_euclidean([ref_rh.aligned{1}(:,[1 3]);ref_rh.aligned{1}(:,1:2)],{surf_lh,surf_rh},labeling.schaefer_400);

control_rh = GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
control_rh = control_rh.fit(mean(corrmats_control_whole(201:end,201:end,:),3),'reference',reference_gradient.gradients{1}(201:end,:));
gradient_in_euclidean([control_rh.aligned{1}(:,[1 3]);control_rh.aligned{1}(:,[1 2])] ,{surf_lh,surf_rh},labeling.schaefer_400);
stroke_rh = GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
stroke_rh = stroke_rh.fit(mean(corrmats_stroke_whole(201:end,201:end,:),3,"omitnan"),'reference',reference_gradient.gradients{1}(201:end,:));
gradient_in_euclidean([stroke_rh.aligned{1}(:,[1 3]);stroke_rh.aligned{1}(:,[1 2])] ,{surf_lh,surf_rh},labeling.schaefer_400);

stroke_whole = GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
stroke_whole = stroke_whole.fit(mean(corrmats_stroke_whole,3,"omitnan"),'reference',reference_gradient.gradients{1});
gradient_in_euclidean([stroke_whole.aligned{1}(:,[1 2])] ,{surf_lh,surf_rh},labeling.schaefer_400);




%% Make a Spin-permutation
% Use non-corrected data
% Check before and after
% Also check the intrahemispheric
%
