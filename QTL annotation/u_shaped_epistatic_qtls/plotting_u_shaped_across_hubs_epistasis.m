%Copied over from
%/Users/shreyasg/Downloads/EF_joint_inference_Jul_2025/scripts/12 -
%u_shaped/plotting_u_shaped_across_hubs/plotting_u_shaped_across_hubs.m

A=readtable("list_of_u_shaped_envs_for_each_epistatic_qtl.txt"); %got rid of the 382 empty rows at the end of the file, yay

%only considered the 4 hubs (not gel) and 5 perts that form a "complete
%set" (we have all combos of hubs and perts)
hubs=["gu ","li ","ynb ","ypd "]; %space after the name because it is a hub %hubs=["gel ","gu ","li ","ynb ","ypd "];
perts=[" temp"," caff"," cyclo"," glu"," gly"]; %space before the name because it is a pert %perts=[" temp"," caff"," cyclo"," glu"," gly"," 4nqo"," kcl"," nacl"," ph"]; 

num_hubs=size(hubs);
num_perts=size(perts);

A_size=size(A);

load("gradient_names.mat");
list_of_grads=grad_ord_env_names;

for i=1:num_perts(1,2) %looping over each pert
    for j=1:A_size(1,1) %looping over each QTL
        strtemp=erase(A{j,3:end},'gel temp'); %removing all occurences of gel temp (as mentioned above, only considered hubs and perts with "complete" hub x grad set
        m_temp=strfind(strtemp,perts(1,i));
        m=sum(~cellfun(@isempty,m_temp)); %cellfun(@isempty,m_temp) is a cell function that checks if each cell of m_temp is empty.
        %The ~ negates the output (marking non-empty cells as 1s) and the
        %sum counts the number of non-empty cells (i.e. matches with the
        %hub being considered)
        num_hubs_u_shaped(j,i)=m; %number of hubs in which pert i (column # is the idx here) is grad ord for qtl j (row # is the idx here)
    end
end

%saving data for python
C=[{"temp","caff","cyclo","glu","gly"}; num2cell(num_hubs_u_shaped)];
writecell(C,'number_of_hubs_with_u_shaped_epistatic_qtls.txt','Delimiter','\t');