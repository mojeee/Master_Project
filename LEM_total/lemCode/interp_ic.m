%%%
%
% loads initial conditions of n numbers of cells
% and interpolates linearly to m numbers of cells
% 
% note: there is a ghost cell @ n+1 and m+1 
%%%

clear all;

load ic.csv;

dom = 1;%cm
col = 14;%number of columns of data
n = 1001;
m = 2021;

dx = 1/(m-1);%new dx

output = zeros((m+1),14);
output(1:m,1) = 0:dx:dom;

for i = 2:14,
  output(:,i) = interp1(ic(:,1),ic(:,i),output(:,1));
  output(m+1,i) = output(m,i);
end

output(m+1,1) = output(m,1)+dx;

f = fopen('newic.csv', 'w');
for row=1:(m+1)
  fprintf(f, '%18.16f %18.16f %18.14f %18.13f %18.16f %18.16f %18.16f %18.16f %18.16f %18.16f %18.16f %18.16f %18.16f %18.16f\n', output(row,:));
end
fclose(f);