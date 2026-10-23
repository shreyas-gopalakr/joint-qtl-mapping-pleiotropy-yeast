%largely copied over from:
%/Users/shreyasg/Downloads/EF_joint_inference_Jun_2024/scripts_101envs/epistasis/gradient_ordering/grad_ord_epistasis2.m

opts = detectImportOptions('final_model_cleaned_epistatic_effects_annotated_v1.txt');
opts = setvartype(opts, 'Var1', 'string');
A = readtable('final_model_cleaned_epistatic_effects_annotated_v1.txt', opts);
A_size=size(A);
additive_QTL=383; %manual input of the number of additive QTL

%most of the following code was copied from gradient_ordered_qtls.m
num_of_envs=102; %manual input, can change later
hub_pert_edges=[2,7,12,16,19,22,26,31,34,36,40,45,50,55,59,62,65,70,75,80,83,86,91,95,99,103];
base_env=[0,23,23,23,23,0,42,42,42,42,0,67,67,67,67,67,0,100,100,100,100,100,100,0,0];
hub_id=[0,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4]; %0=gelda; 1=gu; 2=li; 3=ynb; 4=ypd
pert_id=[1,2,3,4,5,1,2,3,4,5,1,6,2,3,4,5,1,2,3,4,5,7,8,9,1]; %1=temp; 2=caff; 3=cyclo; 4=glu; 5=gly; 6=4nqo; 7=kcl; 8=nacl; 9=ph

grad_ord=zeros(1,25); %number of hub+perturb combos=25;

%only kcl appears to be in the wrong order, rewrite code for that part
%alone (vec=...)
for j=1:25 %number of hub+perturb combos (loop to find epistatic QTL SNP #s that are gradient-ordered)
    count=1;
    for i=additive_QTL+1:A_size(1) %for each epistatic QTL
        clear vec %ADDED THIS LINE MUCH LATER, PROBABLY NEED TO RERUN ALL ASSOCIATED CODE HERE
        %NEED TO DO THE SAME FOR EPISTASIS GRAD-ORDERED CODE TOO
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
        if all(diff(vec)>0) | all(diff(vec)<0)
        grad_ord(count,j)=i;
        count=count+1;
        end
    end
end

%plotting expected and observed number of epistatic gradient-ordered qtls
for j=1:25 %counting gradient-ordered QTLs and significance testing
    grad_ord_num(1,j)=nnz(grad_ord(:,j)); %counting gradient ordered QTL
    str=split(A.Properties.VariableNames{1,hub_pert_edges(1,j)},'_');
    grad_ord_env_names(1,j)=strcat(str{1,1},{' '},str{2,1}); %gradient ordered env/condition name

    if base_env(j)~=0
    pert_grad_size(1,j)=hub_pert_edges(1,j+1)-hub_pert_edges(1,j)+1; %number of points in perturbation gradient (including hub)
    end
    
    if base_env(j)==0
    pert_grad_size(1,j)=hub_pert_edges(1,j+1)-hub_pert_edges(1,j); %number of points in perturbation gradient (no real hub available)
    end

    p(1,j)=binocdf(grad_ord_num(1,j)-1,A_size(1)-additive_QTL,2/factorial(pert_grad_size(1,j)),'upper');
    exp_grad_ord_num(1,j)=(A_size(1)-additive_QTL)*2/factorial(pert_grad_size(1,j));
    
    grad_ord_env_sig_label(1,j)=grad_ord_env_names(1,j); %default label

    if p(1,j)<0.05
    grad_ord_env_sig_label(1,j)=strcat(grad_ord_env_names(1,j),'*');
    end
    
    if p(1,j)<10^(-8)
    grad_ord_env_sig_label(1,j)=strcat(grad_ord_env_names(1,j),'**');
    end

end

scatter(1:25,grad_ord_num(1,:),"filled");
hold on
scatter(1:25,exp_grad_ord_num(1,:),"filled");
ax=gca;
ax.XTick=[1:25];
ax.XTickLabels=grad_ord_env_sig_label(1,:);
xtickangle(90);
xlim([0 26]);
legend('Epistatic QTLs with gradient ordering','Expected number of epistatic QTLs with gradient ordering');

%only copied code from the original script till here, copy further if you 
% want to analyze things further

%saving data for python
C1=cell(4,25);
C1(1,:)=grad_ord_env_names;
C1(2,:)=num2cell(grad_ord_num);
C1(3,:)=num2cell(exp_grad_ord_num);
C1(4,:)=num2cell(p);

T1=cell2table(C1,"RowNames",["grad_names","grad_ord_num","exp_grad_ord_num","p_value"]);
T1.Properties.VariableNames=grad_ord_env_names;
writetable(T1,'grad_ord_epistatic_plotting.txt',Delimiter='\t',WriteRowNames=true);

%continued copying code from gradient_ordered_qtls.m

%TO MAKE FIGURES ON INTERSECTION OF GRADIENT ORDERED QTLS
count2=1;
count3=1;
count4=1;

grad_ord_same_hub=zeros(1,2);
grad_ord_same_pert=zeros(1,2);
grad_ord_random=zeros(1,2); %counting just pairs of hub+perturb with neither a hub not a pert in common

%normalized metric(s) of similarity between pairs of envs' grad-ord qtls
%(jaccard index or overlap coeff) %add the code to calculate such things
grad_ord_same_hub_norm=zeros(1,1); %jaccard index
grad_ord_same_pert_norm=zeros(1,1);
grad_ord_random_norm=zeros(1,1);

grad_ord_same_hub_norm2=zeros(1,1); %overlap coeff
grad_ord_same_pert_norm2=zeros(1,1);
grad_ord_random_norm2=zeros(1,1);

grad_ord_interesect={}; %cell array that keeps track of every pair of "envs" (env_id_1 in col 1, env_id_2 in col 2) and the number of grad ord qtls in common
count5=1; %counter for the same

grad_ord_same_hub_qtls=zeros(0,0); %list of snp #s that are grad ordered across envs with the same hub
grad_ord_same_pert_qtls=zeros(0,0); %list of snp #s that are grad ordered across envs with the same pert
grad_ord_random_qtls=zeros(0,0); %list of snp #s that are grad ordered across random "envs" (not same hub nor same pert)

%edited the code below to add grad_ord_interesect{count5,6},
%grad_ord_interesect{count5,7}, and grad_ord_interesect{count5,8}

for j=1:25 %number of hub+perturb combos (loops looking for pairs of hub+pertub combos with the same hub)
    for k=j+1:25 %number of hub+perturb combos (looking for the pairs of hub+perturb combos)
        if hub_id(j)==hub_id(k) & pert_id(j)~=pert_id(k) 
        grad_ord_same_hub(count2,:)=size(intersect(nonzeros(grad_ord(:,j)),nonzeros(grad_ord(:,k))));
        temp1=intersect(nonzeros(grad_ord(:,j)),nonzeros(grad_ord(:,k)));
        grad_ord_same_hub_qtls=cat(1,grad_ord_same_hub_qtls,temp1);
        
        grad_ord_same_hub_norm(count2,1)=grad_ord_same_hub(count2,1)/(grad_ord_num(1,j)+grad_ord_num(1,k)-grad_ord_same_hub(count2,1));
        grad_ord_same_hub_norm2(count2,1)=grad_ord_same_hub(count2,1)/min(grad_ord_num(1,j),grad_ord_num(1,k));
        count2=count2+1;

        grad_ord_interesect{count5,1}=grad_ord_env_names{1,j};
        grad_ord_interesect{count5,2}=grad_ord_num(1,j);
        grad_ord_interesect{count5,3}=grad_ord_env_names{1,k};
        grad_ord_interesect{count5,4}=grad_ord_num(1,k);
        grad_ord_interesect{count5,5}=grad_ord_same_hub(count2-1,1); %count2-1 because counter is incremented 4 lines above
        grad_ord_interesect{count5,6}='same hub, diff pert';
        grad_ord_interesect{count5,7}=grad_ord_same_hub_norm(count2-1,1);
        grad_ord_interesect{count5,8}=grad_ord_same_hub_norm2(count2-1,1);
        count5=count5+1;
        
        end

        if pert_id(j)==pert_id(k) & hub_id(j)~=hub_id(k)
        grad_ord_same_pert(count3,:)=size(intersect(nonzeros(grad_ord(:,j)),nonzeros(grad_ord(:,k))));
        temp2=intersect(nonzeros(grad_ord(:,j)),nonzeros(grad_ord(:,k)));
        grad_ord_same_pert_qtls=cat(1,grad_ord_same_pert_qtls,temp2);
        
        grad_ord_same_pert_norm(count3,1)=grad_ord_same_pert(count3,1)/(grad_ord_num(1,j)+grad_ord_num(1,k)-grad_ord_same_pert(count3,1));
        grad_ord_same_pert_norm2(count3,1)=grad_ord_same_pert(count3,1)/min(grad_ord_num(1,j),grad_ord_num(1,k));
        count3=count3+1;

        grad_ord_interesect{count5,1}=grad_ord_env_names{1,j};
        grad_ord_interesect{count5,2}=grad_ord_num(1,j);
        grad_ord_interesect{count5,3}=grad_ord_env_names{1,k};
        grad_ord_interesect{count5,4}=grad_ord_num(1,k);
        grad_ord_interesect{count5,5}=grad_ord_same_pert(count3-1,1); %count2-1 because counter is incremented 4 lines above
        grad_ord_interesect{count5,6}='diff hub, same pert';
        grad_ord_interesect{count5,7}=grad_ord_same_pert_norm(count3-1,1);
        grad_ord_interesect{count5,8}=grad_ord_same_pert_norm2(count3-1,1);
        count5=count5+1;
        end

        if hub_id(j)~=hub_id(k) & pert_id(j)~=pert_id(k)
        grad_ord_random(count4,:)=size(intersect(nonzeros(grad_ord(:,j)),nonzeros(grad_ord(:,k))));
        temp3=intersect(nonzeros(grad_ord(:,j)),nonzeros(grad_ord(:,k)));
        grad_ord_random_qtls=cat(1,grad_ord_random_qtls,temp3);
        
        grad_ord_random_norm(count4,1)=grad_ord_random(count4,1)/(grad_ord_num(1,j)+grad_ord_num(1,k)-grad_ord_random(count4,1));
        grad_ord_random_norm2(count4,1)=grad_ord_random(count4,1)/min(grad_ord_num(1,j),grad_ord_num(1,k));
        count4=count4+1;
        
        grad_ord_interesect{count5,1}=grad_ord_env_names{1,j};
        grad_ord_interesect{count5,2}=grad_ord_num(1,j);
        grad_ord_interesect{count5,3}=grad_ord_env_names{1,k};
        grad_ord_interesect{count5,4}=grad_ord_num(1,k);
        grad_ord_interesect{count5,5}=grad_ord_random(count4-1,1); %count2-1 because counter is incremented 4 lines above
        grad_ord_interesect{count5,6}='diff hub, diff pert';
        grad_ord_interesect{count5,7}=grad_ord_random_norm(count4-1,1);
        grad_ord_interesect{count5,8}=grad_ord_random_norm2(count4-1,1);
        count5=count5+1;
        end
    end
end

grad_ord_interesect2=[{'env1','grad_ord_num_env1','env2','grad_ord_num_env2','grad_ord_intersect','env1_env2_relationship','jaccard_index','overlap_coeff'};grad_ord_interesect];
writecell(grad_ord_interesect2,'number_of_shared_grad_ord_epistatic_qtls_pairs_of_envs.txt','Delimiter','tab'); 

%plotting normalized metrics
%jaccard index

%boxchart
boxchart(ones(size(grad_ord_random_norm,1),1),grad_ord_random_norm)
hold on
boxchart(2*ones(size(grad_ord_same_hub_norm,1),1),grad_ord_same_hub_norm)
hold on
boxchart(3*ones(size(grad_ord_same_pert_norm,1),1),grad_ord_same_pert_norm)
xticks([1,2,3]);
xticklabels(["diff hub, diff pert","\bf{same hub}\rm, diff pert", "diff hub, \bf{same pert}"]);
set(gca,'FontSize',15);
xtickangle(0);

ylabel('Jaccard index of shared gradient ordered QTLs');
p1=ranksum(grad_ord_random_norm,grad_ord_same_hub_norm);
p2=ranksum(grad_ord_random_norm,grad_ord_same_pert_norm);
p3=ranksum(grad_ord_same_hub_norm,grad_ord_same_pert_norm);

%overlap coefficient
%boxchart
boxchart(ones(size(grad_ord_random_norm2,1),1),grad_ord_random_norm2)
hold on
boxchart(2*ones(size(grad_ord_same_hub_norm2,1),1),grad_ord_same_hub_norm2)
hold on
boxchart(3*ones(size(grad_ord_same_pert_norm2,1),1),grad_ord_same_pert_norm2)
xticks([1,2,3]);
xticklabels(["diff hub, diff pert","\bf{same hub}\rm, diff pert", "diff hub, \bf{same pert}"]);
set(gca,'FontSize',15);
xtickangle(0);

ylabel('Overlap of shared gradient ordered QTLs');
p4=ranksum(grad_ord_random_norm2,grad_ord_same_hub_norm2);
p5=ranksum(grad_ord_random_norm2,grad_ord_same_pert_norm2);
p6=ranksum(grad_ord_same_hub_norm2,grad_ord_same_pert_norm2);


%code to find out which envs each gene is grad ordered in

grad_ord_cell_array=cell(A_size(1)-additive_QTL,30); %change 30 to the max number in the histogram above + 4 more cols for epistatic pair number, index and gene names (max num of envs a qtl is grad ordered in + 4)

for i=additive_QTL+1:A_size(1) %looping over all epistatic QTLs
    match=(grad_ord==i); %binary matrix showing if i can be found in the grad_ord matrix
    match=sum(match,1); %binary vector with 25 elements, 1 indicating that qtl is grad-ordered in that env
    match_size=size(match);
    
    grad_ord_cell_array{i-additive_QTL,1}=i;
    grad_ord_cell_array{i-additive_QTL,2}=A{i,1}; %value in first column of A 
    if strcmp(A.gene_name{i,1},'N/A')==0 %contains the gene name of the first qtl in the interaction
        grad_ord_cell_array{i-additive_QTL,3}=A.gene_name{i,1};
    else
        grad_ord_cell_array{i-additive_QTL,3}=A.intergenic_gene_name{i,1};
    end
    
    if strcmp(A.qtl2_gene_name{i,1},'N/A')==0
        grad_ord_cell_array{i-additive_QTL,4}=A.qtl2_gene_name{i,1};
    else
        grad_ord_cell_array{i-additive_QTL,4}=A.qtl2_intergenic_gene_name{i,1};
    end


    match_count=5;
    for j=1:match_size(2)
        if match(1,j)==1 %if the env matches a grad_ordered env
            grad_ord_cell_array{i-additive_QTL,match_count}=grad_ord_env_names{1,j};
            match_count=match_count+1;
        end
    end

end

writecell(grad_ord_cell_array,'list_of_gradient_ordered_envs_for_each_epistatic_qtl.txt','Delimiter','tab');

%analyzing results: do something similar to VLTE paper and see if any QTLs
%are significantly associated with a hub or perturbation