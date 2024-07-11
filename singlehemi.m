ix=subjects.lesion_side==1;
dif1=mean(EDtocontrol_before(:,1:num_control),2)-mean(EDtocontrol_before(:,ix),2);
dif2=mean(EDtocontrol_after(:,1:num_control),2)-mean(EDtocontrol_after(:,ix),2);

for i=1:10000
    ed_shuffled1=EDtocontrol_before(:,randperm(length(1:(num_control+sum(ix)))));
    ed_shuffled2=EDtocontrol_after(:,randperm(length(1:(num_control+sum(ix)))));
    shuffled_dif1(:,i)=mean(ed_shuffled1(:,1:num_control),2)-mean(ed_shuffled1(:,num_control+1:end),2);
    shuffled_dif2(:,i)=mean(ed_shuffled1(:,1:num_control),2)-mean(ed_shuffled2(:,num_control+1:end),2);

end
for i=1:400
    corrected1(i)=(abs(dif1(i)))>prctile(abs(shuffled_dif1(i,:)),95);
    corrected2(i)=(abs(dif2(i)))>prctile(abs(shuffled_dif2(i,:)),95);

end
plot_hemispheres([mean(EDtocontrol_before(:,1:num_control),2)...
    mean(EDtocontrol_before(:,num_control+1:end),2) corrected1'.*dif1 corrected2'.*dif2], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400);

damage_on_left=corrected2'.*dif2;
damage_on_right=corrected2'.*dif2;


%%%%%
ix=subjects.lesion_side==1;
dif1=mean(EDtocontrol_before(:,1:num_control),2)-mean(EDtocontrol_before(:,ix),2);
dif2=mean(EDtocontrol_after(:,1:num_control),2)-mean(EDtocontrol_after(:,ix),2);

for i=1:10000
    ed_shuffled1=EDtocontrol_before(:,randperm(length(1:(num_control+sum(ix)))));
    ed_shuffled2=EDtocontrol_after(:,randperm(length(1:(num_control+sum(ix)))));
    shuffled_dif1(:,i)=mean(ed_shuffled1(:,1:num_control),2)-mean(ed_shuffled1(:,num_control+1:end),2);
    shuffled_dif2(:,i)=mean(ed_shuffled1(:,1:num_control),2)-mean(ed_shuffled2(:,num_control+1:end),2);

end
for i=1:400
    corrected1(i)=(abs(dif1(i)))>prctile(abs(shuffled_dif1(i,:)),95);
    corrected2(i)=(abs(dif2(i)))>prctile(abs(shuffled_dif2(i,:)),95);

end
plot_hemispheres([mean(EDtocontrol_before(:,1:num_control),2)...
    mean(EDtocontrol_before(:,num_control+1:end),2) corrected1'.*dif1 corrected2'.*dif2], {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400);

damage_on_left=corrected2'.*dif2;
damage_on_right=corrected2'.*dif2;