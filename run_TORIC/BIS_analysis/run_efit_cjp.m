%%%% This is the main script to run EFIT on the C-Mod data. %%%%
% Authored by: Lucas Spangher
% Things to keep in note here are the following:

% 1. The input file is the time_to_run_remaining_shots_12_04.csv file. This
% file contains the shot number, the start time, and the end time for each
% shot. The start time is the time (on the shot, in seconds) at which the EFIT will start running on the shot,
% and the end time is the time at which the EFIT will stop running. Feel free
% to modify this as you see fit. This table was created from flattop_times.py

% 2. You MUST create an efitname_model.{tree, characteristics, datafile} EFIT trees
% to run matlab: matlab -r -nodesktop -nosplash. Then, run "run_efit_cjp". matlab is dumb. 

filename = 'time_to_run_shots_all_1.csv';

% Load the CSV file into a matrix
data = readmatrix(filename, 'Delimiter', ',', "HeaderLines", 1);

% Extract the columns into separate arrays
shotlist = data(:, 1);

temp_nshots = length(shotlist);

% first_half = false; % Set to true if you want to flip the script to the other half of the shots

% note this is to split the Processing across two different machines. 
% if first_half
  % Select the first half of shotlist
%  shotlist = shotlist(1:floor(temp_nshots/2));
% else
  % Select the second half of shotlist
%  shotlist = shotlist(floor(temp_nshots/2)+1:end);
% end

nshots = length(shotlist);

start_time = data(:, 2);
end_time = data(:, 3);

ndigits = int2str(floor(log10(nshots)) + 1);
tree_dir = sprintf('%s/trees', pwd);
setenv('efit_cjp_path', tree_dir)

segfaults = [];

temp_seg_faults = cell(nshots, 1);

parfor ishot = 1:nshots
  shot = shotlist(ishot);
  start_time_i = start_time(ishot);
  end_time_i = end_time(ishot);

  fprintf(1,['Processing shot %10i  (%' ndigits 'i/%' ndigits ...
    'i  %6.2f%%)\n'], shot, ishot, nshots, ishot/nshots*100);

  unique_dir = sprintf('%s/dir_%d', pwd, shot);
  mkdir(unique_dir);
  original_dir = pwd;
  cd(unique_dir);

  % Time resolution to run EFITs
  dt = 1e-3; % [s]

  % number of timesteps is the floor of the difference between the start and
  % end times divided by the timestep

  ntimes = floor((end_time_i - start_time_i) / dt)

  [~, status] = mdsopen('efit_cjp', shot);
    disp(status)
    
  if (mod(status,2) == 0)
    mdsopen('efit_cjp', -1);
    mdstcl(['create pulse ' int2str(shot)]);
    mdsclose;

  else
    fprintf(1,'  Note: EFIT_CJP tree already exists for this shot... overwriting, bitches.\n');
    mdsopen('efit_cjp', -1);
    mdstcl(['create pulse ' int2str(shot)]);
    mdsclose;
  end

  % Create the script file that EFIT will read its inputs from.  Run EFIT in
  % burst mode (mode 10), with two bursts: a slow (standard) one, and a fast
  % (pre-disruption) one.  (By specifying a negative number for the number of
  % bursts, the C-Mod versions of the EFIT executables will sort the combined
  % burst timebases into numerical order so that the full timebase will be
  % monotonic.)

    % Use a unique filename for each iteration of the loop
    input_filename = sprintf('efit_script_input_%d.txt', shot);
    output_filename = sprintf('efit_script_output_%d.txt', shot);

    fileID = fopen(input_filename, 'w');
    fprintf(fileID, '9\n');
    fprintf(fileID, 'efit_cjp\n');
    fprintf(fileID, '\n');  
    fprintf(fileID, '%i\n', shot);
    fprintf(fileID, '%6.1f, %4.1f, %i\n', start_time_i*1000, dt*1e3, ntimes);
    fclose(fileID);

    % Call EFIT with the unique filenames

    % fprintf("Running '/usr/local/cmod/codes/efit/bin/fast_efitdd through segfault tolerant script\n")
    command = sprintf(['bash /home/cjperks/transport_world/run_TORIC/BIS_analysis/run_efit_segfault_tolerant.sh %s %s'], input_filename, output_filename);
    [status, cmdout] = unix(command);

    if status ~= 0
      fprintf('Error in processing shot %d: %s\n', shot, cmdout);
      fprintf('Segfault occurred on shot %d\n', shot);
      unix(sprintf('rm trees/efit_cjp_%d.datafile', shot));
      unix(sprintf('rm trees/efit_cjp_%d.tree', shot));
      unix(sprintf('rm trees/efit_cjp_%d.characteristics', shot));
      temp_seg_faults{ishot} = shot;

    end

    delete(input_filename);
    delete(output_filename);
    cd(original_dir);

    try
      rmdir(unique_dir, 's');  % Remove the directory and all its contents
    catch 
      fprintf('Could not remove directory %s\n', unique_dir);
    end

end

% Save segfaults as a CSV file
segfaults = [temp_seg_faults{:}];
csvwrite('segfaults_1.csv', segfaults);

