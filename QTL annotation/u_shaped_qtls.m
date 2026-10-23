%code is largely based on the following code, with some edits:
%/Users/shreyasg/Downloads/EF_joint_inference_Jun_2024/scripts_101envs/initial
%GO term analysis and gradient ordering/individual
%envs/u_shaped_qtls_no_extra_plots.m

A=readtable('final_model_cleaned_linear_effects_annotated_v7.txt');
A_size=size(A);
num_of_envs=102; %manual input, can change later

hub_pert_edges=[2,7,12,16,19,22,26,31,34,36,40,45,50,55,59,62,65,70,75,80,83,86,91,95,99,103];
base_env=[0,23,23,23,23,0,42,42,42,42,0,67,67,67,67,67,0,100,100,100,100,100,100,0,0];
hub_id=[0,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4]; %0=gelda; 1=gu; 2=li; 3=ynb; 4=ypd
pert_id=[1,2,3,4,5,1,2,3,4,5,1,6,2,3,4,5,1,2,3,4,5,7,8,9,1]; %1=temp; 2=caff; 3=cyclo; 4=glu; 5=gly; 6=4nqo; 7=kcl; 8=nacl; 9=ph

%I checked that all 4 of the variables above are equal to those in
%/Users/shreyasg/Downloads/EF_joint_inference_Jul_2025/9 - gradient
%ordering/gradient_ordered_qtls.m

u_shaped=zeros(1,25); %number of hub+perturb combos=25;

for j=1:25 %number of hub+perturb combos (loop to find QTL SNP #s that are U-shaped)
    count=1; 
    for i=2:A_size(1) %for each QTL
        clear vec
        if j~=22 %all but ypd+kcl envs which correspond to j=22; this is because kcl envs are out of order
        if base_env(j)~=0
        vec(1,1)=A{i,base_env(j)};
        vec(1,2:hub_pert_edges(1,j+1)-hub_pert_edges(1,j)+1)=A{i,hub_pert_edges(1,j):hub_pert_edges(1,j+1)-1};
        end
        if base_env(j)==0
        vec=A{i,hub_pert_edges(1,j):hub_pert_edges(1,j+1)-1};
        end
        end
        if j==22 %for kcl envs
        vec(1,1)=A{i,base_env(j)};
        vec(1,2:4)=A{i,87:89};
        vec(1,5)=A{i,86};
        vec(1,6)=A{i,90};
        end
        vec_size=size(vec);
        if vec_size(2)>4 %only looking at perturbation gradient with at least 5 points, so that we can detect two points going down and two points going up (stringent condition for U-shaped qtls)
        vec_diff=diff(vec);
        [n1,n2]=num_of_sign_changes(vec_diff);
        if n1==1 & n2>=2 & n2<=vec_size(2)-3 %checking that there is only 1 sign change and that at least 2 points are going up/down in the U-shape
        u_shaped(count,j)=i;
        count=count+1;
        
        %plotting these u-shaped things (this line to line 57 are new)
        % if count==2 %making a folder for each hub+pert combo
        %     str=split(A.Properties.VariableNames{1,hub_pert_edges(1,j)},'_');
        %     str2=strcat(str{1,1},{' '},str{2,1}); %name of the hub+pert combo
        %      if j~=1
        %      cd ..
        %      end
        %     mkdir(sprintf(str2{1,1}))
        %     cd(sprintf(str2{1,1}))
        % end
        
        % plot(1:vec_size(2),vec,'-*'); %effect size (y) vs pert level (x)
        % title(A{i,104}{:,:});
        % ylabel('Effect size');
        % xlabel('Pert. grad. level');
        % saveas(gcf,strcat(strrep(A{i,104}{:,:},'/','_'),'_',num2str(A{i,1}),'.png'));
        % close(gcf);
        end
        end
    end
end

u_shaped_size=size(u_shaped);
u_shaped_genes=strings(u_shaped_size(1)+1,u_shaped_size(2));

for i=2:u_shaped_size(1)+1 %finding genes that are u-shaped
    for j=1:u_shaped_size(2)
    if i==2
        str=split(A.Properties.VariableNames{1,hub_pert_edges(1,j)},'_');
        u_shaped_genes(1,j)=strcat(str{1,1},{' '},str{2,1});
    end

    if u_shaped(i-1,j)~=0 & strcmp(A.gene_name{u_shaped(i-1,j),1},'N/A')==0 %gene is not labelled N/A
    u_shaped_genes(i,j)=A.gene_name{u_shaped(i-1,j),1};
    end

    if u_shaped(i-1,j)~=0 & strcmp(A.gene_name{u_shaped(i-1,j),1},'N/A')==1 %gene is not labelled N/A
    u_shaped_genes(i,j)=A.intergenic_gene_name{u_shaped(i-1,j),1}; %don't put in gene name, put in intergenic gene name
    end
    end
end

%CALCULATING NUMBER OF U-SHAPED QTL AND SIGNIFICANCE TESTING

for j=1:25 %counting U-shaped QTLs and significance testing
    u_shaped_num(1,j)=nnz(u_shaped(:,j));
    str=split(A.Properties.VariableNames{1,hub_pert_edges(1,j)},'_');
    u_shaped_env_names(1,j)=strcat(str{1,1},{' '},str{2,1});

    if base_env(j)~=0
    pert_grad_size(1,j)=hub_pert_edges(1,j+1)-hub_pert_edges(1,j)+1; %number of points in perturbation gradient (including hub)
    end

    if base_env(j)==0
    pert_grad_size(1,j)=hub_pert_edges(1,j+1)-hub_pert_edges(1,j); %number of points in perturbation gradient (no real hub available)
    end

    u_shaped_env_sig_label(1,j)=u_shaped_env_names(1,j); %default label
    p(1,j)=1; %"default" p-value
    exp_u_shaped_num(1,j)=0; %expected number of u-shaped qtls

    if pert_grad_size(1,j)>4 %looking at gradients with at least 5 points, same as above
    %Calculate how many orderings (based on highest-lowest effect size) correspond to U-shapes and put that in the numerator instead of 2
        if pert_grad_size(1,j)==5
        p(1,j)=binocdf(u_shaped_num(1,j)-1,A_size(1)-1,12/factorial(pert_grad_size(1,j)),'upper');
        exp_u_shaped_num(1,j)=(A_size(1)-1)*12/factorial(pert_grad_size(1,j));
        end
        
        if pert_grad_size(1,j)==6
        p(1,j)=binocdf(u_shaped_num(1,j)-1,A_size(1)-1,40/factorial(pert_grad_size(1,j)),'upper');
        exp_u_shaped_num(1,j)=(A_size(1)-1)*40/factorial(pert_grad_size(1,j));
        end

        if p(1,j)<0.05
        u_shaped_env_sig_label(1,j)=strcat(u_shaped_env_names(1,j),'*');
        end

        if p(1,j)<10^(-4)
        u_shaped_env_sig_label(1,j)=strcat(u_shaped_env_names(1,j),'**');
        end

    end
end

%saving data for python
C1=cell(4,25);
C1(1,:)=u_shaped_env_names;
C1(2,:)=num2cell(u_shaped_num);
C1(3,:)=num2cell(exp_u_shaped_num);
C1(4,:)=num2cell(p);

T1=cell2table(C1,"RowNames",["grad_names","u_shaped_num","exp_u_shaped_num","p_value"]);
T1.Properties.VariableNames=u_shaped_env_names;
%uncomment if you need to save it again
writetable(T1,'u_shaped_additive_plotting.txt',Delimiter='\t',WriteRowNames=true);

scatter(1:25,u_shaped_num(1,:),"filled");
hold on
scatter(1:25,exp_u_shaped_num(1,:),"filled");
ax=gca;
ax.XTick=[1:25];
ax.XTickLabels=u_shaped_env_sig_label(1,:);
xtickangle(90);
xlim([0 26]);
legend('Additive QTLs with U-shaped ordering','Expected number of additive QTLs with U-shaped ordering');

%code to find out which envs each gene is U-shaped in

u_shaped_cell_array=cell(A_size(1),15); %can change 15 to 7, a qtl is max u-shaped in 5 diff envs

for i=1:A_size(1) %looping over all QTLs
    match=(u_shaped==i); %binary matrix showing if i can be found in the u_shaped matrix
    match=sum(match,1); %binary vector with 25 elements, 1 indicating that qtl is grad-ordered in that env
    match_size=size(match);
    
    u_shaped_cell_array{i,1}=i;
    if strcmp(A.gene_name{i,1},'N/A')==0
        u_shaped_cell_array{i,2}=A.gene_name{i,1};
    else
        u_shaped_cell_array{i,2}=A.intergenic_gene_name{i,1};
    end
    
    match_count=3;
    for j=1:match_size(2)
        if match(1,j)==1 %if the env matches a u_shaped env
            u_shaped_cell_array{i,match_count}=u_shaped_env_names{1,j};
            match_count=match_count+1;
        end
    end

end

%uncomment if you need to save this file again
%writecell(u_shaped_cell_array,'list_of_u_shaped_envs_for_each_qtl.txt','Delimiter','tab');

%continuing code from here on Sep 10, 2025. The main reason is to get the
%number of u-shaped qtls shared between different pairs of gradients (diff
%hub, diff pert; same hub, diff pert; diff hub, same pert.)
%Largely copying over code from
%/Users/shreyasg/Downloads/EF_joint_inference_Jul_2025/9 - gradient
%ordering/gradient_ordered_qtls.m

%TO MAKE FIGURES ON INTERSECTION OF U-SHAPED ORDERED QTLS
count2=1;
count3=1;
count4=1;

u_shaped_same_hub=zeros(1,2);
u_shaped_same_pert=zeros(1,2);
u_shaped_random=zeros(1,2); %counting just pairs of hub+perturb with neither a hub not a pert in common

%normalized metric(s) of similarity between pairs of envs' grad-ord qtls
%(jaccard index or overlap coeff) %add the code to calculate such things
u_shaped_same_hub_norm=zeros(1,1); %jaccard index
u_shaped_same_pert_norm=zeros(1,1);
u_shaped_random_norm=zeros(1,1);

u_shaped_same_hub_norm2=zeros(1,1); %overlap coeff
u_shaped_same_pert_norm2=zeros(1,1);
u_shaped_random_norm2=zeros(1,1);

u_shaped_interesect={}; %cell array that keeps track of every pair of "envs" (env_id_1 in col 1, env_id_2 in col 2) and the number of grad ord qtls in common
count5=1; %counter for the same

u_shaped_same_hub_qtls=zeros(0,0); %list of snp #s that are grad ordered across envs with the same hub
u_shaped_same_pert_qtls=zeros(0,0); %list of snp #s that are grad ordered across envs with the same pert
u_shaped_random_qtls=zeros(0,0); %list of snp #s that are grad ordered across random "envs" (not same hub nor same pert)

%edited the code below to add u_shaped_interesect{count5,6},
%u_shaped_interesect{count5,7}, and u_shaped_interesect{count5,8}

%after the first iteration of the code, I realized i have to remove all
%envs (11 of them) where i don't compute the number of u_shaped qtls as
%they have too few points in the gradient. Basically, need to remove all
%rows of u_shaped_interesect2 where exp_u_shaped_num of an env in that row
%(env names in col 1 and col 3) is 0.
%ACTUALLY, SIMPLEST WAY TO REMOVE THE ABOVE MENTIONED ENVS: REMOVE ENVS
%WITH OVERLAP COEFF = NaN. HAPPENS WHEN ONE OF THE ENVS HAS U_SHAPED_NUM=0
%AND THIS ONLY HAPPENS WHEN EXP_U_SHAPED_NUM=0 FOR OUR DATASET. CHECK THIS
%WHEN RUNNING.

for j=1:25 %number of hub+perturb combos (loops looking for pairs of hub+pertub combos with the same hub)
    for k=j+1:25 %number of hub+perturb combos (looking for the pairs of hub+perturb combos)
        if hub_id(j)==hub_id(k) & pert_id(j)~=pert_id(k) 
        u_shaped_same_hub(count2,:)=size(intersect(nonzeros(u_shaped(:,j)),nonzeros(u_shaped(:,k))));
        temp1=intersect(nonzeros(u_shaped(:,j)),nonzeros(u_shaped(:,k)));
        u_shaped_same_hub_qtls=cat(1,u_shaped_same_hub_qtls,temp1);
        
        u_shaped_same_hub_norm(count2,1)=u_shaped_same_hub(count2,1)/(u_shaped_num(1,j)+u_shaped_num(1,k)-u_shaped_same_hub(count2,1));
        u_shaped_same_hub_norm2(count2,1)=u_shaped_same_hub(count2,1)/min(u_shaped_num(1,j),u_shaped_num(1,k));
        count2=count2+1;

        u_shaped_interesect{count5,1}=u_shaped_env_names{1,j};
        u_shaped_interesect{count5,2}=u_shaped_num(1,j);
        u_shaped_interesect{count5,3}=u_shaped_env_names{1,k};
        u_shaped_interesect{count5,4}=u_shaped_num(1,k);
        u_shaped_interesect{count5,5}=u_shaped_same_hub(count2-1,1); %count2-1 because counter is incremented 4 lines above
        u_shaped_interesect{count5,6}='same hub, diff pert';
        u_shaped_interesect{count5,7}=u_shaped_same_hub_norm(count2-1,1);
        u_shaped_interesect{count5,8}=u_shaped_same_hub_norm2(count2-1,1);
        count5=count5+1;
        
        end

        if pert_id(j)==pert_id(k) & hub_id(j)~=hub_id(k)
        u_shaped_same_pert(count3,:)=size(intersect(nonzeros(u_shaped(:,j)),nonzeros(u_shaped(:,k))));
        temp2=intersect(nonzeros(u_shaped(:,j)),nonzeros(u_shaped(:,k)));
        u_shaped_same_pert_qtls=cat(1,u_shaped_same_pert_qtls,temp2);
        
        u_shaped_same_pert_norm(count3,1)=u_shaped_same_pert(count3,1)/(u_shaped_num(1,j)+u_shaped_num(1,k)-u_shaped_same_pert(count3,1));
        u_shaped_same_pert_norm2(count3,1)=u_shaped_same_pert(count3,1)/min(u_shaped_num(1,j),u_shaped_num(1,k));
        count3=count3+1;

        u_shaped_interesect{count5,1}=u_shaped_env_names{1,j};
        u_shaped_interesect{count5,2}=u_shaped_num(1,j);
        u_shaped_interesect{count5,3}=u_shaped_env_names{1,k};
        u_shaped_interesect{count5,4}=u_shaped_num(1,k);
        u_shaped_interesect{count5,5}=u_shaped_same_pert(count3-1,1); %count2-1 because counter is incremented 4 lines above
        u_shaped_interesect{count5,6}='diff hub, same pert';
        u_shaped_interesect{count5,7}=u_shaped_same_pert_norm(count3-1,1);
        u_shaped_interesect{count5,8}=u_shaped_same_pert_norm2(count3-1,1);
        count5=count5+1;
        end

        if hub_id(j)~=hub_id(k) & pert_id(j)~=pert_id(k)
        u_shaped_random(count4,:)=size(intersect(nonzeros(u_shaped(:,j)),nonzeros(u_shaped(:,k))));
        temp3=intersect(nonzeros(u_shaped(:,j)),nonzeros(u_shaped(:,k)));
        u_shaped_random_qtls=cat(1,u_shaped_random_qtls,temp3);
        
        u_shaped_random_norm(count4,1)=u_shaped_random(count4,1)/(u_shaped_num(1,j)+u_shaped_num(1,k)-u_shaped_random(count4,1));
        u_shaped_random_norm2(count4,1)=u_shaped_random(count4,1)/min(u_shaped_num(1,j),u_shaped_num(1,k));
        count4=count4+1;
        
        u_shaped_interesect{count5,1}=u_shaped_env_names{1,j};
        u_shaped_interesect{count5,2}=u_shaped_num(1,j);
        u_shaped_interesect{count5,3}=u_shaped_env_names{1,k};
        u_shaped_interesect{count5,4}=u_shaped_num(1,k);
        u_shaped_interesect{count5,5}=u_shaped_random(count4-1,1); %count2-1 because counter is incremented 4 lines above
        u_shaped_interesect{count5,6}='diff hub, diff pert';
        u_shaped_interesect{count5,7}=u_shaped_random_norm(count4-1,1);
        u_shaped_interesect{count5,8}=u_shaped_random_norm2(count4-1,1);
        count5=count5+1;
        end
    end
end

%as mentioned earlier, using the NaN entries of overlap coeff to filter out
%all cases when env1 or env2 is not considered for calculating number of
%u_shaped qtls
u_shaped_interesect = u_shaped_interesect(~isnan(cell2mat(u_shaped_interesect(:,8))), :);

%all entries of the following statement need to be 1, else the above NaN
%based filtering doesn't work
(exp_u_shaped_num==0)==(u_shaped_num==0) 

u_shaped_interesect2=[{'env1','u_shaped_num_env1','env2','u_shaped_num_env2','u_shaped_intersect','env1_env2_relationship','jaccard_index','overlap_coeff'};u_shaped_interesect];
writecell(u_shaped_interesect2,'number_of_shared_u_shaped_qtls_pairs_of_envs.txt','Delimiter','tab'); 

%plotting normalized metrics
%jaccard index

%%boxchart
% %the commented out code no longer works as the 3 arrays: u_shaped_...norm
% %contain elements where either env1 or env2 should have been filtered.
% %However, we can re-write the code here in some simple way to get it to
% %work. For now though, I'll do the plotting in python
% boxchart(ones(size(u_shaped_random_norm,1),1),u_shaped_random_norm)
% hold on
% boxchart(2*ones(size(u_shaped_same_hub_norm,1),1),u_shaped_same_hub_norm)
% hold on
% boxchart(3*ones(size(u_shaped_same_pert_norm,1),1),u_shaped_same_pert_norm)
% xticks([1,2,3]);
% xticklabels(["diff hub, diff pert","\bf{same hub}\rm, diff pert", "diff hub, \bf{same pert}"]);
% set(gca,'FontSize',15);
% xtickangle(0);
% 
% ylabel('Jaccard index of shared U-shaped ordered QTLs');
% p4=ranksum(u_shaped_random_norm,u_shaped_same_hub_norm);
% p5=ranksum(u_shaped_random_norm,u_shaped_same_pert_norm);
% p6=ranksum(u_shaped_same_hub_norm,u_shaped_same_pert_norm);

%continue copying over code from here