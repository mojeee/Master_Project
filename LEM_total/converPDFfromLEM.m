%%%
%
%  Hong -- 2017/01/19
%
%  This code is for:
%  1) changing the output of the LEM PDF from bmp10.f from 51 bins to 50 bins
%  2) smoothing out any irregular behavior in the LEM PDF (say 0's)
%  3) Prints the output into Graham's OpenFOAM format (with his version of CSEFoam
%     back in November 2016)
%     --> I believe he has updated his format since then...
%
%  Several notes:
%  a) check that the PDF.csv file (the output file from bmp10.f) has no vertical 
%     line spacing, data is always one row after another for a total of 2500 rows
%  b) check that the ordering of the mean and variance is printed correctly
%     (though I think it will be)
%  c) check that the input file name is correct ( PDF.csv )
%  d) at the final step of the code, it automatically prints the some sample 
%     PDFs so you can check if everything worked out.
%     1 - yes, 0 - no
%  e) I can't seem to find my unconditional generation code, but that one is easy.
%     All you need to do:
%     i) take the conditional gradients from output of gradc.f and calculate the 
%        conditional SDR, which is:
%        conditional SDR = (conditional Diff coefficient)*(conditional gradients)^2
%        --> eg: conditional Diff coefficent of CO2 at each progress variable bin
%        --> eg: conditional thermal diffusivity at each temperature-based p.v. bin
%        need to add a line like: load gradients.csv (or whatever the gradient file is called)
%     ii) take the inner product of the conditional SDR with the PDF at each 
%         mean and variance
%     iii) make a table of 2500 doubles according ot Graham's format for OpenFOAM
%
%
%%%

clear all;
% clf;

load PDF.csv; %check this name, other than this, the code *SHOULD* be automatic...
x=zeros(1,50);
x(1:50) = 0.01:0.02:0.99;

% convert PDF from 51 to 50 bins
PDF = zeros(2500,50);

for i = 1:2500,
    PDF(i,1) = PDF1a(i,1) + 0.5*PDF1a(i,2);
    PDF(i,2:49) = 0.5*PDF1a(i,2:49) + 0.5*PDF1a(i,3:50);
    PDF(i,50) = 0.5*PDF1a(i,50) + PDF1a(i,51);
end


meancoor1=zeros(2500,1);
varcoor1 =zeros(2500,1);
for i = 1:2500,
   meancoor1(i)= sum( x.*PDF(i,:) );
   varcoor1(i) = sum(  PDF(i,:).*( meancoor1(i) - x ).^2  ) ...
                   /( meancoor1(i)*(1-meancoor1(i)) );
end

disp('Pdf data converted to 50 bins.');

points1 = zeros(2500,1);

meancoors=zeros(2500,1);
varcoors=zeros(2500,1);

for i=1:50,
    for j=1:50,
        meancoors( (i-1)*50 +j ) = 0.01+ 0.02*(i-1);
        varcoors( (i-1)*50 +j )  = 0.01+ 0.02*(j-1);
    end
end


for k=1:50,
    points1(:) = PDF(:,k);
    vq1 = griddata(meancoor1,varcoor1,points1,meancoors,varcoors);
    tempPDF1(:,k) = vq1;
end

% % modify the PDFs when there's a close-to-zero variance
% first of every 1 +50*n rows
for i = 1:50,
    tempPDF1((1+(i-1)*50),:) = 0.0;
    tempPDF1((1+(i-1)*50),i) = 1.0;
end
% modify the PDF for extremely low mean (first 50 rows)
for i = 1:50,
    tempPDF1(i,:) = 0.0;
    tempPDF1(i,1) = 1.0;
end
% modify the PDF for extremely high mean (last 50 rows)
for i = 2451:2500,
    tempPDF1(i,:) = 0.0;
    tempPDF1(i,50) = 1.0;
end

% check means / vars again
for index = 1:2500,
    tempmean1(index,1) = sum( x .* tempPDF1(index,:) );
    tempvar1(index,1) = sum(  tempPDF1(index,:) .*( tempmean1(index,1) - x ).^2  ) ...
                       /( tempmean1(index,1)*(1-tempmean1(index,1)) );
end

nanArray=zeros(2500,1);

% check for NaN's (if NaN, set PDF to generic max var @ correct mean value)
for index = 1:2500,
    if( isnan(tempmean1(index)) )        
        tempPDF1(index,1)  = (0.99 - meancoors(index) )/0.98;
        tempPDF1(index,50) = 1 - tempPDF1(index,1);
        tempPDF1(index,2:49) = 0.0;
    else

    end
    
end

% check means / vars again
for index = 1:2500,
    tempmean1(index,1) = sum( x .* tempPDF1(index,:) );
    tempvar1(index,1) = sum(  tempPDF1(index,:) .*( tempmean1(index,1) - x ).^2  ) ...
                       /( tempmean1(index,1)*(1-tempmean1(index,1)) );
end

binCentres = 0.01:0.02:0.99;
binCentres = binCentres;
cBins = 0.01:0.02:0.99;
segBins = 0.01:0.02:0.99;
% 
% 
% 
stringbin  = 'binCentres    50 (';
stringend  = ');';
stringmean = 'means         50 (';
stringseg  = 'segregations  50 (';

string1   = '                50 (';
string2   = ' )';
string3   = '         )';
string4   = '        50';
string5   = '        (';


% Now write into Graham's format
fileID = fopen('output.txt','w');
%binCentres
fprintf(fileID,'%s %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %s\n', stringbin, cBins, stringend);
% fprintf(fileID,'%s %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %s\n', stringbin, binCentres, stringend);
%means
fprintf(fileID,'%s %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %s\n', stringmean, cBins, stringend);
%segs
fprintf(fileID,'%s %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %s\n', stringseg, segBins, stringend);


for i = 1:50,
    for j = 1:50,
        fprintf(fileID,'%s %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %s\n', ...
                string1, tempPDF1(  ( (i-1)*50 +j ) , 1:50  ), string2);
    end
    fprintf(fileID,'%s\n %s\n %s\n', string3, string4, string5);
end
fclose(fileID);

i = 1;
k = 1;
while( i < 2501 && k > 0 )
    disp(i);
 
    subplot(5,1,1);
    plot(x,tempPDF1(i,:));
    axis([0 1 0 1]);
    strmean = ['mean = ',num2str(tempmean1(i))];
    strvar  = ['var = ',num2str(tempvar1(i))];
     
    text(0.7,0.7,strmean);
    text(0.7,0.6,strvar);
     
    subplot(5,1,2);
    plot(x,tempPDF1(i+1,:));
    axis([0 1 0 1]);
    strmean = ['mean = ',num2str(tempmean1(i+1))];
    strvar  = ['var = ',num2str(tempvar1(i+1))];
     
    text(0.7,0.7,strmean);
    text(0.7,0.6,strvar);
    
    subplot(5,1,3);
    plot(x,tempPDF1(i+2,:));
    axis([0 1 0 1]);
    strmean = ['mean = ',num2str(tempmean1(i+2))];
    strvar  = ['var = ',num2str(tempvar1(i+2))];
     
    text(0.7,0.7,strmean);
    text(0.7,0.6,strvar);
     
    subplot(5,1,4);
    plot(x,tempPDF1(i+3,:));
    axis([0 1 0 1]);
    strmean = ['mean = ',num2str(tempmean1(i+3))];
    strvar  = ['var = ',num2str(tempvar1(i+3))];
     
    text(0.7,0.7,strmean);
    text(0.7,0.6,strvar);
    
    subplot(5,1,5);
    plot(x,tempPDF1(i+4,:));
    axis([0 1 0 1]);
    strmean = ['mean = ',num2str(tempmean1(i+4))];
    strvar  = ['var = ',num2str(tempvar1(i+4))];
     
    text(0.7,0.7,strmean);
    text(0.7,0.6,strvar);    
     
    i = i + 5;
 
    prompt = 'continue?: 1-yes, 0-no';
    k = input(prompt);
end

