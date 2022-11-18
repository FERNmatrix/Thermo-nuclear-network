function Xerror()


% Manually change the variables below

dataPoints = 1000; % no. of data points (or No. of lines in file)
iso = 16; % Number of isotopes in the network

i=0; % initialize loop variable to 0
r = 3; % # of calculations being used (For the mean of each X)

% initialize all arrays
            % Time arrays
Tlog = zeros(i);
T = zeros(i);
            %Fast arrays
Xlog = zeros(i);
X = zeros(i);
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

% Open data file to be read in
fid = fopen('AlphaHydroplot1.data','r'); % FAST
fid2 = fopen('AlphaHydroRef1.data','r'); % REF
fid3 = fopen('AlphaHydroAcc1.data','r'); % ACC

% read in lines from the file
Columns = 7+iso; 

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
end

for i=1:dataPoints
    
    %First column is read in as time, instead of being 0
    Xf(i) = Xf(i);
    Xa(i) = Xa(i);
    Xr(i) = Xr(i);
    
    %First column of mean is set to be the log time
    Mean(i,1) = Xr(i);
    
    %Use mean and deviation from the mean for each isotope
    for j=2:iso+1
        Mean(i,j) = (Xf(i,j)+ Xa(i,j)+ Xr(i,j))/(r);
        
        DevXf(i,j) = (Xf(i,j) - Mean(i,j))^2; %Sqaure of deviation
        DevXa(i,j) = (Xa(i,j) - Mean(i,j))^2;
        DevXr(i,j) = (Xr(i,j) - Mean(i,j))^2;
        
        %Uncertainty in each isotope
        IsoUnc(i,j) = ((DevXf(i,j) + DevXa(i,j) + DevXr(i,j))/(6))^(1/2);
    end
end
%Combine the mean and Uncertainty to 1 matrix
% First column is log t
%Column 2-iso+1 is the mean of each iso
%Any column after iso+1 is uncertainty for indexed iso
C = cat(2,Mean,IsoUnc);

%Write out matrix with time, mean X, and Uncertainty to file for plotting
writematrix(C,'MeanXE.txt','delimiter',' ');
type MeanXE.txt;

%for i=1:dataPoints
%Xax(i) = Mean(i,1);
%end
%Xax=Xax';
%Xax;
%}
%for i=1:dataPoints

 %   Y1 = Mean(i,2);
 %   Y2 = Mean(i,3);
 %   Y3 = Mean(i,4);
    
 %   E1 = IsoUnc(i,2);
 %   E2 = IsoUnc(i,3);
 %   E3 = IsoUnc(i,4);
%end

fclose(fid);
fclose(fid2);
fclose(fid3);
end
