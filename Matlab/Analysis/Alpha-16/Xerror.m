function Xerror()



%               Manually change the variables below
%_________________________________________________________________%
dataPoints = 200;   % no. of data points (or No. of lines in file)
iso = 16;           % Number of isotopes in the network

i=0;                % initialize loop variable to 0
r = 3;              % # of calculations being used (For the mean of each X)
%_________________________________________________________________%
% initialize all arrays
            % Time arrays
Tlog = zeros(i);
T = zeros(i);
            %Fast arrays
Xflog = zeros(i);
Xf = zeros(i);
            % Ref arrays
Xrlog = zeros(i);
Xr = zeros(i);
            % Acc arrays
Xalog = zeros(i);
Xa = zeros(i);

            %EXTRA ARRAYS
IsoUnc = zeros(i);
Mean = zeros(i);
DevX = zeros(i);
DevXa = zeros(i);
DevXr = zeros(i);
%_________________________________________________________________%
% Open data file to be read in
fid = fopen('Alphafast1.data','r'); % FAST
fid2 = fopen('AlphaRef1.data','r'); % REF
fid3 = fopen('AlphaAcc1.data','r'); % ACC
%________________________________________________________________%

% read in lines from the files, 7 is the number of columns before the mass
% fractions are recorded in the output file from FERN. For structure of
% file, just look at typica plot output file
Columns = 7+iso; 

%________________________________________________________________%
%                   Read in files and create matrices 
%________________________________________________________________%

%If error occurs in while loop, check to see that the data files only
%contain numbers (First few lines of FERN output need to be deleted)
while (i < dataPoints)
    i = i+1;
    %Open %ASY-REF FILE
    U = fscanf(fid2,'%f %f %f',Columns);
    U = U';
    
    %Open PE-FAST FILE
    V = fscanf(fid,'%f %f %f',Columns);
    V = V';
    
    %Open PE-ACC FILE
    W = fscanf(fid3,'%f %f %f',Columns);
    W = W';
    
    %First column in each matrix will be Log time(s)
    Xf(i,1) = V(1);
    Xr(i,1) = U(1);
    Xa(i,1) = W(1);
    
    %Rest of columns (2-> iso) will be mass fractions
    for j=2:iso+1
        Xf(i,j) = V(6+j);
        Xr(i,j) = U(6+j);
        Xa(i,j) = W(6+j);
    end 
end % END WHILE-LOOP

%________________________________________________________________%
%                       Matrix  Manipulation
%________________________________________________________________%
% Mass Fraction array X, should now be a matrix [dataPoints x Isotopes]
% Example if using triple alpha w/ 200 plot outputs array size [200 x 3]

% Array is currently in Log value, so to get true values,
% get rid of logs: X = 10^(log(X)) for each array (X,Xr,Xa)
for i=1:dataPoints
    for j=2:iso
        Xf(i,j) = 10^(Xflog(i,j));
        Xa(i,j) = 10^(Xalog(i,j));
        Xr(i,j) = 10^(Xrlog(i,j));
    end
        T(i) = 10^(Tlog(i));    
    
end

%________________________________________________________________%
%           Find the Mean + Uncertainty of each matrix     
% Link; https://www.educba.com/uncertainty-formula/
%________________________________________________________________%

for i=1:dataPoints
    
    %First column is read in as time, instead of being 0
    Xf(i) = Xf(i);
    Xa(i) = Xa(i);
    Xr(i) = Xr(i);
    
    %First column of mean-matrix is set to be the log time
    % Doesn't matter which matrix is on the right, they should all have the
    % same time iteration, if not use the reference calculation
    Mean(i,1) = Xr(i);
    
    % Use mean and deviation from the mean for each isotope
    % Remember r is set to 3 (typically) for number of runs being used
    % Fill the second column of the Mean-matrix with mass fractions
    for j=2:iso+1
        Mean(i,j) = (Xf(i,j)+ Xa(i,j)+ Xr(i,j))/(r);
        
        %Sqaure of deviation
        DevXf(i,j) = (Xf(i,j) - Mean(i,j))^2; 
        DevXa(i,j) = (Xa(i,j) - Mean(i,j))^2;
        DevXr(i,j) = (Xr(i,j) - Mean(i,j))^2;
        
        % Isotopic Uncertainty
        IsoUnc(i,j) = ((DevXf(i,j) + DevXa(i,j) + DevXr(i,j))/(6))^(1/2);
    end
end
%___________________________________________________________________________%

% Combine the mean and Uncertainty to 1 matrix
% First column is log(t)
% Column 2 -> iso+1 is the mean of each iso
% Any column after iso+1 is uncertainty for indexed iso
C = cat(2,Mean,IsoUnc);

% Write out matrix with time, mean X, and Uncertainty to file for plotting
% Use Gnuplot script to plot the mean using uncertainty as error bars
writematrix(C,'MeanXE.txt','delimiter',' ');
type MeanXE.txt;



fclose(fid);
fclose(fid2);
fclose(fid3);
end
