function [y_1,y_2] = num_of_sign_changes(x)
x_size=size(x);
y_1=0; %number of sign changes
y_2=0; %position of sign changes
for i=1:x_size(2)-1
if (x(i)*x(i+1))<0 %indicating that x(i) and x(i+1) are of opposite sign, indicating a sign change
y_1=y_1+1;
y_2=i; %records the position of the LAST sign change
end
end
end