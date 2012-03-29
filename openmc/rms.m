% Post Processing Script for CMFD runs

%% input information
path = '/media/Backup/opr_runs/';
hdfile = 'output.h5';
dim = 17;
Np = 177; % normalization constant (number of fissile zones)

%% read in 1 million case - openmc source

% other input info
run_start = 1;
run_end = 25;
cycle_start = 201;
cycle_end = 840;
var = 'openmc_src';

% preallocate a temp
src_tmp = zeros(dim,dim,run_end - run_start + 1);

% preallocate src
src_1 = cell(cycle_end - cycle_start + 1,1);

% begin loop over cycles
for j = cycle_start:cycle_end
      
    % begin loop over runs
    for i = run_start:run_end
       
        % read into a temp
        read_tmp = h5read(horzcat(path,'1mil/run',num2str(i),'/',hdfile),horzcat('/cycle',num2str(j),'/',var));
        
        % put into src temp
        src_tmp(:,:,i-run_start+1) = read_tmp(1,:,:,1);
        
    end
    
    % put into cell object
    src_1{j-cycle_start+1} = src_tmp;
    
    % display to user
    fprintf('Read in cycle: %d\n',j);
    
end

%% read in 1 million case - cmfd

% other input info
run_start = 1;
run_end = 25;
cycle_start = 201;
cycle_end = 840;
var = 'cmfd_source';

% preallocate a temp
src_tmp = zeros(dim,dim,run_end - run_start + 1);

% preallocate src
src_cmfd_1 = cell(cycle_end - cycle_start + 1,1);

% begin loop over cycles
for j = cycle_start:cycle_end
      
    % begin loop over runs
    for i = run_start:run_end
       
        % read into a temp
        read_tmp = h5read(horzcat(path,'1mil/run',num2str(i),'/',hdfile),horzcat('/cycle',num2str(j),'/',var));
        
        % put into src temp
        src_tmp(:,:,i-run_start+1) = read_tmp(1,:,:,1);
        
    end
    
    % put into cell object
    src_cmfd_1{j-cycle_start+1} = src_tmp;
    
    % display to user
    fprintf('Read in cycle: %d\n',j);
    
end

%% read in 64 million case


% other input info
run_start = 1;
run_end = 25;
cycle_start = 201;
cycle_end = 210;
var = 'openmc_src';

% preallocate a temp
src_tmp = zeros(dim,dim,run_end - run_start + 1);

% preallocate src
src_64 = cell(cycle_end - cycle_start + 1,1);

% begin loop over cycles
for j = cycle_start:cycle_end
      
    % begin loop over runs
    for i = run_start:run_end
       
        % read into a temp
        read_tmp = h5read(horzcat(path,'64mil/run',num2str(i),'/',hdfile),horzcat('/cycle',num2str(j),'/',var));
        
        % put into src temp
        src_tmp(:,:,i-run_start+1) = read_tmp(1,:,:,1);
        
    end
    
    % put into cell object
    src_64{j-cycle_start+1} = src_tmp;
    
    % display to user
    fprintf('Read in cycle: %d\n',j);
    
end

%% Compute Mean

% for last cycle of 64 mil
meanref = mean(src_64{10},3);

% make surface plot
surf(meanref)

%% Compute All RMS for 1 mil

% preallocate rms vector
rms_1_all = zeros(length(src_1),25);

% begin loop over cycles
for i = 1:length(src_1)
    
    for j = 1:25
        
        % compute rms and store
        rms_1_all(i,j) = sqrt(sum(sum((src_1{i}(:,:,j) - meanref).^2))*1/Np);
        
    end
    
end

%% Compute All RMS for 1 mil CMFD

% preallocate rms vector
rms_cmfd_1_all = zeros(length(src_cmfd_1),25);

% begin loop over cycles
for i = 1:length(src_cmfd_1)
    
    for j = 1:25
        
        % compute rms and store
        rms_cmfd_1_all(i,j) = sqrt(sum(sum((src_cmfd_1{i}(:,:,j) - meanref).^2))*1/Np);
        
    end
    
end

%% Compute All RMS for 64 mil

% preallocate rms vector
rms_64_all = zeros(length(src_64),25);

% begin loop over cycles
for i = 1:length(src_64)
    
    for j = 1:25
        
        % compute rms and store
        rms_64_all(i,j) = sqrt(sum(sum((src_64{i}(:,:,j) - meanref).^2))*1/Np);
        
    end
    
end

%% Plot All RMS

figure
loglog(linspace(1e6,640e6,640),rms_1_all*100,'b.');
hold on
loglog(linspace(1e6,640e6,640),rms_cmfd_1_all*100,'go');
loglog(linspace(64e6,640e6,10),rms_64_all*100,'r+');

%% Compute RMS for 1 mil

% preallocate rms vector
rms_1 = zeros(length(src_1),1);

% begin loop over cycles
for i = 1:length(src_1)
   
    % compute rms and store
    rms_1(i) = sqrt(sum(sum(((mean(src_1{i}(:,:,1:10),3) - meanref).^2)))*1/Np);
    
    % display
    fprintf('RMS computed for cycle: %d\n',i);
    
end

%% Compute RMS for 64 mil

% preallocate rms vector
rms_64 = zeros(length(src_64),1);

% begin loop over cycles
for i = 1:length(src_64)
   
    % compute rms and store
    rms_64(i) = sqrt((1/Np)*sum(sum((mean(src_64{i}(:,:,1:10),3) - meanref).^2)));
    
    % display
    fprintf('RMS computed for cycle: %d\n',i);
    
end