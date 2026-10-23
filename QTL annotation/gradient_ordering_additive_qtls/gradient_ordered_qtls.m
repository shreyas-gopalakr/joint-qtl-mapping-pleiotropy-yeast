%Copied over from
%/Users/shreyasg/Downloads/EF_joint_inference_May_2025/scripts/9 - gradient
%ordering

A=readtable("final_model_cleaned_linear_effects_annotated_v7.txt");

A_size=size(A);
num_of_envs=102; %manual input, can change later
hub_pert_edges=[2,7,12,16,19,22,26,31,34,36,40,45,50,55,59,62,65,70,75,80,83,86,91,95,99,103];
base_env=[0,23,23,23,23,0,42,42,42,42,0,67,67,67,67,67,0,100,100,100,100,100,100,0,0];
hub_id=[0,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4]; %0=gelda; 1=gu; 2=li; 3=ynb; 4=ypd
pert_id=[1,2,3,4,5,1,2,3,4,5,1,6,2,3,4,5,1,2,3,4,5,7,8,9,1]; %1=temp; 2=caff; 3=cyclo; 4=glu; 5=gly; 6=4nqo; 7=kcl; 8=nacl; 9=ph

grad_ord=zeros(1,25); %number of hub+perturb combos=25;

%only kcl appears to be in the wrong order, rewrote code for that part
%alone (vec=...)
for j=1:25 %number of hub+perturb combos (loop to find QTL SNP #s that are gradient-ordered)
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
        if all(diff(vec)>0) | all(diff(vec)<0)
        grad_ord(count,j)=i;
        count=count+1;
        end
    end
end

grad_ord_size=size(grad_ord);
grad_ord_genes=strings(grad_ord_size(1)+1,grad_ord_size(2));

for i=2:grad_ord_size(1)+1 %finding genes that are gradient-ordered
    for j=1:grad_ord_size(2)
    if i==2
        str=split(A.Properties.VariableNames{1,hub_pert_edges(1,j)},'_');
        grad_ord_genes(1,j)=strcat(str{1,1},{' '},str{2,1});
    end
    if grad_ord(i-1,j)~=0 & strcmp(A.gene_name{grad_ord(i-1,j),1},'N/A')==0 %gene is not labelled N/A
    grad_ord_genes(i,j)=A.gene_name{grad_ord(i-1,j),1};
    end

    if grad_ord(i-1,j)~=0 & strcmp(A.gene_name{grad_ord(i-1,j),1},'N/A')==1 %gene is labelled N/A
    grad_ord_genes(i,j)=A.intergenic_gene_name{grad_ord(i-1,j),1}; %don't put in gene name, put in intergenic gene name
    end
    end
end

%plotting expected and observed number of gradient-ordered qtls
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
    
    p(1,j)=binocdf(grad_ord_num(1,j)-1,A_size(1)-1,2/factorial(pert_grad_size(1,j)),'upper');
    exp_grad_ord_num(1,j)=(A_size(1)-1)*2/factorial(pert_grad_size(1,j));
    
    grad_ord_env_sig_label(1,j)=grad_ord_env_names(1,j); %default label
    
    if p(1,j)<0.05
    grad_ord_env_sig_label(1,j)=strcat(grad_ord_env_names(1,j),'\bf{*}');
    end
    
    if p(1,j)<10^(-8)
    grad_ord_env_sig_label(1,j)=strcat(grad_ord_env_names(1,j),'\bf{**}');
    end

end

%MAKE THE X-AXIS LABELS BOLD FOR THE SIGNIFICANT ENVS

%scatter plot
% scatter(1:25,grad_ord_num(1,:),50,"filled"); %15 is the marker size
% hold on
% scatter(1:25,exp_grad_ord_num(1,:),50,"filled"); %15 is the marker size
% ax=gca;
% ax.XTick=[1:25];
% ax.XTickLabels=grad_ord_env_sig_label(1,:);
% set(gca,'FontSize',18);
% xtickangle(90);
% xlim([0 26]);
% legend('Additive QTLs with gradient ordering','Expected number of additive QTLs with gradient ordering','FontSize',12);

%bar plot
grad_ord_num_bar_plot=[grad_ord_num.' exp_grad_ord_num.'];
bar(grad_ord_num_bar_plot);
set(gca,'FontSize',18);
ax2=gca;
ax2.XTick=[1:25];
ax2.XTickLabels=grad_ord_env_names;
xtickangle(90);
xlim([0 26]);
legend('Additive QTLs with gradient ordering','Expected number of additive QTLs with gradient ordering','FontSize',12);

%saving data for python
C1=cell(4,25);
C1(1,:)=grad_ord_env_names;
C1(2,:)=num2cell(grad_ord_num);
C1(3,:)=num2cell(exp_grad_ord_num);
C1(4,:)=num2cell(p);

T1=cell2table(C1,"RowNames",["grad_names","grad_ord_num","exp_grad_ord_num","p_value"]);
T1.Properties.VariableNames=grad_ord_env_names;
writetable(T1,'grad_ord_additive_plotting.txt',Delimiter='\t',WriteRowNames=true);


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
%writecell(grad_ord_interesect2,'number_of_shared_grad_ord_qtls_pairs_of_envs.txt','Delimiter','tab'); 

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
p4=ranksum(grad_ord_random_norm,grad_ord_same_hub_norm);
p5=ranksum(grad_ord_random_norm,grad_ord_same_pert_norm);
p6=ranksum(grad_ord_same_hub_norm,grad_ord_same_pert_norm);

%Finding "unique" gradient-ordered qtls (only shared between hubs but not
%randomly, only between pert but not randomly, etc).
%The next 25 lines of analysis (unique_same_pert_wrt_random=setdiff...) are
%exploratory, haven't yielded anything worth really analyzing further. The
%results also seem difficult/convoluted to analyze.
grad_ord_same_hub_qtls=unique(grad_ord_same_hub_qtls); %"unique" list of gradient-ordered qtls (no repetitions) shared between hubs
grad_ord_same_pert_qtls=unique(grad_ord_same_pert_qtls);
grad_ord_random_qtls=unique(grad_ord_random_qtls);

grad_ord_same_hub_qtls_size=size(grad_ord_same_hub_qtls);
for i=1:grad_ord_same_hub_qtls_size(1) %finding genes that are gradient-ordered shared b/w hub
    grad_ord_same_hub_genes(i,1)=A.gene_name(grad_ord_same_hub_qtls(i,1),1);
end

grad_ord_same_pert_qtls_size=size(grad_ord_same_pert_qtls);
for i=1:grad_ord_same_pert_qtls_size(1) %finding genes that are gradient-ordered shared b/w pert
    grad_ord_same_pert_genes(i,1)=A.gene_name(grad_ord_same_pert_qtls(i,1),1);
end

grad_ord_random_qtls_size=size(grad_ord_random_qtls);
for i=1:grad_ord_random_qtls_size(1) %finding genes that are gradient-ordered shared b/w random* envs
    grad_ord_random_genes(i,1)=A.gene_name(grad_ord_random_qtls(i,1),1);
end

intersection_same_hub_pert=intersect(grad_ord_same_hub_genes,grad_ord_same_pert_genes);
unique_same_hub_wrt_pert=setdiff(grad_ord_same_hub_genes,grad_ord_same_pert_genes);
unique_same_pert_wrt_hub=setdiff(grad_ord_same_pert_genes,grad_ord_same_hub_genes);

unique_same_hub_wrt_random=setdiff(grad_ord_same_hub_genes,grad_ord_random_genes); %very few of these. Does this mean something?
unique_same_pert_wrt_random=setdiff(grad_ord_same_pert_genes,grad_ord_random_genes); %very few of these. Does this mean something?

%Number of envs in which a QTL is gradient-ordered
temp=unique(grad_ord(:,1));
for i=2:grad_ord_size(2) %looping over envs
temp=cat(1,temp,unique(grad_ord(:,i))); %converting grad_ord into a single array (from a matrix)
end
grad_ord_qtl_occurences=temp;
blacklisted_qtls=[0]; %so that we don't count the 0s (they are not qtls)
grad_ord_qtl_occurences=grad_ord_qtl_occurences(ismember(grad_ord_qtl_occurences,blacklisted_qtls)==0); %removing entries from array that are 0
grad_ord_count=accumarray(grad_ord_qtl_occurences,1); %grad ord count is the number of times a qtl is gradient ordered, qtl id is just the row number of the array itself
figure;
H=histogram(grad_ord_count(2:end)); %not counting the first element of grad_ord_count as it corresponds to the intercept, the "fake" qtl #1

xlabel('Number of env. gradients in which QTL is gradient-ordered');
ylabel('Number of gradient-ordered QTLs');
%title('Histogram of number of envs in which QTL is gradient-ordered');
set(gca,"FontSize",16);

%did a lot of the editing from here to ~line 309 on the fly. If it ends up
%being important, redo it more carefully.
xlim([-0.5 15.5]);
xticks([0:1:15]);

y_pos_hist_text=[19;19;7;5;5;5;5]; %manual input of the y-positions of the gene names highlighted in the histogram
x_pos_hist_text=[9.6;10.6;11.6;12.6;13.6;14.6;15.6];

for i=9.5:1:15.5
hist_indices=find(grad_ord_count>i & grad_ord_count<i+1);
hist_text=A.gene_name(hist_indices,1);
    for j=1:size(hist_text,1)
        if strcmp(hist_text{j,1},'N/A')==1
            hist_text{j,1}=A.intergenic_gene_name{hist_indices(j,1),1};
        end
    end
text(x_pos_hist_text(i-8.5,1),y_pos_hist_text(i-8.5,1),hist_text,FontSize=7);
end

non_grad_ord_any_env=find(grad_ord_count==0); %finding out which qtls aren't
%gradient ordered in any env. 5/382 qtls not grad ord in any env. qtl #1 is
%a fake qtl, corresponds to the intercept, it doesn't count.

%code to find out which envs each gene is grad ordered in

grad_ord_cell_array=cell(A_size(1),17); %change 17 to the max number in the histogram above + 2 more cols for qtls pos and gene name (max num of envs a qtl is grad ordered in + 2)

for i=1:A_size(1) %looping over all QTLs
    match=(grad_ord==i); %binary matrix showing if i can be found in the grad_ord matrix
    match=sum(match,1); %binary vector with 25 elements, 1 indicating that qtl is grad-ordered in that env
    match_size=size(match);
    
    grad_ord_cell_array{i,1}=i;
    if strcmp(A.gene_name{i,1},'N/A')==0
        grad_ord_cell_array{i,2}=A.gene_name{i,1};
    else
        grad_ord_cell_array{i,2}=A.intergenic_gene_name{i,1};
    end
    
    match_count=3;
    for j=1:match_size(2)
        if match(1,j)==1 %if the env matches a grad_ordered env
            grad_ord_cell_array{i,match_count}=grad_ord_env_names{1,j};
            match_count=match_count+1;
        end
    end

end

writecell(grad_ord_cell_array,'list_of_gradient_ordered_envs_for_each_qtl.txt','Delimiter','tab');

%analyzing results: do something similar to VLTE paper and see if any QTLs
%are significantly associated with a hub or perturbation
