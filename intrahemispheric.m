

mean_control=mean(corrmats_corrected_reduced(:,:,1:num_control),3);%mean_control(isinf(mean_control))=1;
mean_control_right=mean(corrmats_corrected_reduced(201:end,201:end,1:num_control),3); %mean_control_right(isinf(mean_control_right))=1;
mean_control_left=mean(corrmats_corrected_reduced(1:200,1:200,1:num_control),3);%mean_control_left(isinf(mean_control_left))=1;
mean_stroke=mean(corrmats_corrected_reduced(:,:,num_control+1:end),3,"omitnan");%mean_stroke(isinf(mean_stroke))=1;
mean_stroke_right=mean(corrmats_corrected_reduced(201:end,201:end,subjects.lesion_side==0),3,"omitnan");%mean_stroke_right(isinf(mean_stroke_right))=1;
mean_stroke_left=mean(corrmats_corrected_reduced(1:200,1:200,subjects.lesion_side==1),3,"omitnan");%mean_stroke_left(isinf(mean_stroke_left))=1;

% references
[surf_lh, surf_rh] = load_conte69();
labeling = load_parcellation('schaefer',400);
conn_matrices = load_group_fc('schaefer',400);
m = (conn_matrices.schaefer_400);
corrmats_right=atanh(m(201:end,201:end));
corrmats_left=atanh(m(1:200,1:200));
ref_whole = GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
m(isinf(m))=1;

ref_whole = ref_whole.fit(m);
% plot_hemispheres(ref_whole.gradients{1}(:,1:3), {surf_lh,surf_rh}, ...
%     'parcellation',labeling.schaefer_400);
% gradient_in_euclidean([ref_whole.gradients{1}(:,[1:3])] ,{surf_lh,surf_rh},labeling.schaefer_400);

ref_lh = GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
corrmats_left(isinf(corrmats_left))=1;
ref_lh = ref_lh.fit(corrmats_left, 'reference', ref_whole.gradients{1}(1:200,:));
% plot_hemispheres([ref_lh.aligned{1}(:,1:2);ref_lh.aligned{1}(:,1:2)] , {surf_lh,surf_rh}, ...
%     'parcellation',labeling.schaefer_400);
% gradient_in_euclidean([ref_lh.aligned{1}(:,1:2);ref_lh.aligned{1}(:,1:2)],{surf_lh,surf_rh},labeling.schaefer_400);

ref_rh = GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
corrmats_right(isinf(corrmats_right))=1;
ref_rh = ref_rh.fit(corrmats_right, 'reference', ref_whole.gradients{1}(201:400,:));
% plot_hemispheres([ref_rh.aligned{1}(:,1:3);ref_rh.aligned{1}(:,1:3)] , {surf_lh,surf_rh}, ...
%     'parcellation',labeling.schaefer_400);
% gradient_in_euclidean([ref_rh.aligned{1}(:,[1 2]);ref_rh.aligned{1}(:,[1 2])],{surf_lh,surf_rh},labeling.schaefer_400);

% Calculate the gradients
controle_whole_gradient=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
controle_whole_gradient=controle_whole_gradient.fit(mean_control,'reference',ref_whole.gradients{1});
% gradient_in_euclidean([controle_whole_gradient.aligned{1}(:,[1 2])] ,{surf_lh,surf_rh},labeling.schaefer_400);

controle_right_gradient=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
controle_right_gradient=controle_right_gradient.fit(mean_control_right,'reference',ref_rh.aligned{1});
%gradient_in_euclidean([controle_right_gradient.aligned{1}(:,[1 2]);controle_right_gradient.aligned{1}(:,[1 2])] ,{surf_lh,surf_rh},labeling.schaefer_400);

controle_left_gradient=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
controle_left_gradient=controle_left_gradient.fit(mean_control_left,'reference',ref_lh.aligned{1});
% gradient_in_euclidean([controle_left_gradient.aligned{1}(:,[1 2]);controle_left_gradient.aligned{1}(:,[1 2])] ,{surf_lh,surf_rh},labeling.schaefer_400);

stroke_whole_gradient=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
stroke_whole_gradient=stroke_whole_gradient.fit(mean_stroke,'reference',ref_whole.gradients{1});
% gradient_in_euclidean([stroke_whole_gradient.aligned{1}(:,[1 3])] ,{surf_lh,surf_rh},labeling.schaefer_400);

stroke_right_gradient=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
stroke_right_gradient=stroke_right_gradient.fit(mean_stroke_right,'reference',ref_rh.aligned{1});
 %gradient_in_euclidean([stroke_right_gradient.aligned{1}(:,[1 2]);stroke_right_gradient.aligned{1}(:,[1 2])] ,{surf_lh,surf_rh},labeling.schaefer_400);

stroke_left_gradient=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
stroke_left_gradient=stroke_left_gradient.fit(mean_stroke_left,'reference',ref_lh.aligned{1});
% gradient_in_euclidean([stroke_left_gradient.aligned{1}(:,[1 2]);stroke_left_gradient.aligned{1}(:,[1 2])] ,{surf_lh,surf_rh},labeling.schaefer_400);

%
% Try to do with subject-wise data

for i=1:size(corrmats_reduced,3)
    disp(i)
    m=corrmats_reduced(:,:,i);
    m(isnan(m))=1;
    gradients_subjectwise=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
    gradients_subjectwise=gradients_subjectwise.fit(m,'reference',ref_whole.gradients{1});
    gradients_subjectwise_all{i}=gradients_subjectwise;
    gradients_subjectwise_aligned(:,:,i)=gradients_subjectwise.aligned{1};
    gradients_subjectwise_lambda(:,i)=gradients_subjectwise.lambda{1}./sum(gradients_subjectwise.lambda{1});

    m=corrmats_reduced(201:end,201:end,i);
    m(isnan(m))=1;
    gradients_subjectwise_right=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
    gradients_subjectwise_right=gradients_subjectwise.fit(m,'reference',ref_whole.gradients{1}(201:end,:));
    gradients_subjectwise_all_right{i}=gradients_subjectwise_right;
    gradients_subjectwise_right_aligned(:,:,i)=gradients_subjectwise_right.aligned{1};
    gradients_subjectwise_right_lambda(:,i)=gradients_subjectwise_right.lambda{1}./sum(gradients_subjectwise_right.lambda{1});

    m=corrmats_reduced(1:200,1:200,i);
    m(isnan(m))=1;
    gradients_subjectwise_left=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
    gradients_subjectwise_left=gradients_subjectwise.fit(m,'reference',ref_whole.gradients{1}(1:200,:));
    gradients_subjectwise_all_left{i}=gradients_subjectwise_left;
    gradients_subjectwise_left_aligned(:,:,i)=gradients_subjectwise_left.aligned{1};
    gradients_subjectwise_left_lambda(:,i)=gradients_subjectwise_left.lambda{1}./sum(gradients_subjectwise_left.lambda{1});

end
gradients_control_whole=gradients_subjectwise_aligned(:,:,1:num_control);cw=mean(gradients_control_whole,3);
gradients_control_left=gradients_subjectwise_left_aligned(:,:,1:num_control);cl=mean(gradients_control_left,3);
gradients_control_right=gradients_subjectwise_right_aligned(:,:,1:num_control);cr=mean(gradients_control_right,3);

gradients_stroke_whole=gradients_subjectwise_aligned(:,:,num_control+1:end);sw=mean(gradients_stroke_whole,3);
gradients_stroke_left=gradients_subjectwise_left_aligned(:,:,subjects.lesion_side==1);sl=mean(gradients_stroke_left,3);
gradients_stroke_right=gradients_subjectwise_right_aligned(:,:,subjects.lesion_side==0);sr=mean(gradients_stroke_right,3);

% gradient_in_euclidean([cl(:,1:3);cl(:,1:3)] ,{surf_lh,surf_rh},labeling.schaefer_400);
% gradient_in_euclidean([cr(:,1:3);cr(:,1:3)] ,{surf_lh,surf_rh},labeling.schaefer_400);
% gradient_in_euclidean([sl(:,1:3);sl(:,1:3)] ,{surf_lh,surf_rh},labeling.schaefer_400);
% gradient_in_euclidean([sr(:,1:3);sr(:,1:3)] ,{surf_lh,surf_rh},labeling.schaefer_400);


right_distance=zeros(size(gradients_stroke_right));
left_distance=zeros(size(gradients_stroke_left));

for i=1:200
    for j=1:10
        for k=1:size(gradients_stroke_right,3)
            right_distance(i,j,k)=norm(controle_right_gradient.aligned{1}(i,j)-gradients_stroke_right(i,j,k));
        end
    end
end
for i=1:200
    for j=1:10
        for k=1:size(gradients_stroke_left,3)
            left_distance(i,j,k)=norm(controle_left_gradient.aligned{1}(i,j)-gradients_stroke_left(i,j,k));
        end
    end
end

% Hemisphere-wise comparison
for j=1:10
        [h,p,c,stats]=ttest2(mean(squeeze(left_distance(:,j,:)),2),mean(squeeze(right_distance(:,j,:)),2));
        t_subs1(j)=stats.tstat;
        p_subs1(j)=p;
end

mean_left_distance=mean(left_distance,3);
mean_right_distance=mean(right_distance,3);



mean(squeeze(left_distance(:,j,:)),2)-mean(squeeze(left_distance(:,j,:)),2)



plot_hemispheres([mean(squeeze(left_distance(:,2,:)),2);mean(squeeze(right_distance(:,2,:)),2)] , {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400);



%%
right_distance_all=zeros(200,size(right_distance,3));
left_distance_all=zeros(200,size(left_distance,3));
for i=1:200
    for k=1:size(right_distance,3)
        right_distance_all(i,k)=norm(controle_right_gradient.aligned{1}(i,1:3)-gradients_stroke_right(i,1:3,k));
    end
end
for i=1:200
    for k=1:size(left_distance,3)
        left_distance_all(i,k)=norm(controle_left_gradient.aligned{1}(i,1:3)-gradients_stroke_left(i,1:3,k));
    end
end

intact3d_after=[mean(left_distance_all,2);mean(right_distance_all,2)]
intact1d_after=[mean(squeeze(left_distance(:,2,:)),2);mean(squeeze(right_distance(:,2,:)),2)]

[h,p,c,stats]=ttest2(mean(left_distance_all,2)',mean(right_distance_all,2)')
plot_hemispheres(intact3d, {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400);
dd=intact3d;
dd2=intact1d;

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

left_intact_before=left_distance_all;
right_intact_before=right_distance_all;

left_damaged_before=left_distance_all; % t test  1.4856   -1.7399    1.3237 p val 0.1382    0.0826    0.1864  3D -0.0680 0.9458
right_damaged_before=right_distance_all;

left_intact_after=left_distance_all; %           -1.8368   -3.4233   -2.6835 0.0670    0.0007    0.0076 -5.7298    1.9874e-08
right_intact_after=right_distance_all;

left_damaged_after=left_distance_all; %           0.1010   -2.9278   -0.0837 0.9196    0.0036    0.9333 -2.8189     0.0051
right_damaged_after=right_distance_all;


[h,p,c,stats]=ttest2(mean(left_intact_before,2)',mean(right_intact_before,2)')
[h,p,c,stats]=ttest2(mean(left_intact_after,2)',mean(right_intact_after,2)')


[h,p,c,stats]=ttest2(mean(left_damaged_before,2)',mean(right_damaged_before,2)')
[h,p,c,stats]=ttest2(mean(left_damaged_after,2)',mean(right_damaged_after,2)')

dif31=mean(left_intact_before,2)'-mean(right_intact_before,2)';
dif32=mean(left_intact_after,2)'-mean(right_intact_after,2)';
[h,p,c,stats]=ttest(dif31-dif32)

dif33=mean(left_damaged_before,2)'-mean(right_damaged_before,2)';
dif34=mean(left_damaged_after,2)'-mean(right_damaged_after,2)';
[h,p,c,stats]=ttest(dif33-dif34)
%%

for i=1:200
    for k=1:num_control
        right_control_all_after(i,k)=norm(controle_right_gradient.aligned{1}(i,1:3)-gradients_control_right(i,1:3,k));
    end
end
for i=1:200
    for k=1:num_control
        left_control_all_after(i,k)=norm(controle_left_gradient.aligned{1}(i,1:3)-gradients_control_left(i,1:3,k));
    end
end

[h,p,c,stats]=ttest2(mean(left_control_all_after,2)',mean(right_control_all_after,2)')


right_distance=zeros(size(gradients_control_right));
left_distance=zeros(size(gradients_control_left));

for i=1:200
    for j=1:10
        for k=1:size(gradients_control_right,3)
            right_distance_control(i,j,k)=norm(controle_right_gradient.aligned{1}(i,j)-gradients_control_right(i,j,k));
        end
    end
end
for i=1:200
    for j=1:10
        for k=1:size(gradients_control_left,3)
            left_distance_control(i,j,k)=norm(controle_left_gradient.aligned{1}(i,j)-gradients_control_left(i,j,k));
        end
    end
end

[~,p,~,stats]=ttest2(mean(squeeze((left_distance_control(:,1,:))),2),mean(squeeze((right_distance_control(:,1,:))),2))
[~,p,c,stats]=ttest2(mean(squeeze((left_distance_control(:,2,:))),2),mean(squeeze((right_distance_control(:,2,:))),2))
[~,p,~,stats]=ttest2(mean(squeeze((left_distance_control(:,3,:))),2),mean(squeeze((right_distance_control(:,3,:))),2))

dd=[mean(squeeze((left_distance_control(:,2,:))),2);mean(squeeze((right_distance_control(:,2,:))),2)]