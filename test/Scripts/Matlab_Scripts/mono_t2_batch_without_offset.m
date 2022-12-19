% This m-file has been automatically generated using qMRgenBatch(mono_t2)
% Command Line Interface (CLI) is well-suited for automatization 
% purposes and Octave. 
%
% Please execute this m-file section by section to get familiar with batch
% processing for mono_t2 on CLI.
%
% Demo files are downloaded into mono_t2_data folder. 
%
% Written by: Agah Karakuzu, 2017
% =========================================================================
%% I- DESCRIPTION
qMRinfo('mono_t2'); % Describe the model
%% II- MODEL PARAMETERS
%%       a- create object    
Model = mono_t2; 
%%       b- modify options
%           |- This section will pop-up the options GUI. Close window to continue.
%           |- Octave is not GUI compatible. Modify Model.options directly.
#Model = Custom_OptionsGUI(Model); % You need to close GUI to move on. 
%% III- FIT EXPERIMENTAL DATASET
%%       a- set protocols via gui 
EchoTime  = [12.8000; 25.6000; 38.4000; 51.2000; 64.0000; 76.8000; 89.6000; 102.4000; 115.2000; 128.0000; 140.8000; 153.6000; 166.4000; 179.2000; 192.0000; 204.8000; 217.6000; 230.4000; 243.2000; 256.0000; 268.8000; 281.6000; 294.4000; 307.2000; 320.0000; 332.8000; 345.6000; 358.4000; 371.2000; 384.0000];
% EchoTime (ms) is a vector of [30X1]
Model.Prot.SEdata.Mat = [ EchoTime ];
Model.options.OffsetTerm=false
%%   
%%       b- load experimental data 
%          |- mono_t2 object needs 2 data input(s) to be assigned:
%          |-   SEdata
%          |-   Mask
data = struct();
% SEdata.nii.gz contains [260  320    1   30] data.
data.SEdata=double(load_nii_data('mono_t2_data/SEdata.nii.gz'));
% Mask.nii.gz contains [260  320] data.
data.Mask=double(load_nii_data('mono_t2_data/Mask.nii.gz'));
 
%%      b- fit dataset 
%             |- This section will fit data. 
startTime=clock()
FitResults = FitData(data,Model,0);
elapsed_time = etime (clock (), startTime)
%%       c- show fitting results 
%           |- Output map will be displayed.
%           |- If available, a graph will be displayed to show fitting in a voxel.
qMRshowOutput(FitResults,data,Model);
%%       d- Save results
%           |-  qMR maps are saved in NIFTI and in a structure FitResults.mat
%                that can be loaded in qMRLab graphical user interface
%           |-  Model object stores all the options and protocol.
%                It can be easily shared with collaborators to fit their 
%                own data or can be used for simulation.
FitResultsSave_nii(FitResults, 'mono_t2_data/SEdata.nii.gz');
Model.saveObj('mono_t2_Demo.qmrlab.mat');
%% V- SIMULATIONS
%     |- This section can be executed to run simulations for mono_t2.
%%       a- Single Voxel Curve
%           |- Simulates Single Voxel curves:
%                (1) use equation to generate synthetic MRI data
%                (2) add rician noise
%                (3) fit and plot curve
      x = struct;
      x.T2 = 100;
      x.M0 = 1000;
       Opt.SNR = 50;
      % run simulation
      figure('Name','Single Voxel Curve Simulation');
      FitResult = Model.Sim_Single_Voxel_Curve(x,Opt);
%%       b- Sensitivity Analysis 
%           |-    Simulates sensitivity to fitted parameters:
%                  (1) vary fitting parameters from lower (lb) to upper (ub) bound.
%                  (2) run Sim_Single_Voxel_Curve Nofruns times
%                  (3) Compute mean and std across runs
      %              T2            M0            
      OptTable.st = [1e+02         1e+03]; % nominal values
      OptTable.fx = [0             1]; %vary T2...
      OptTable.lb = [1             1]; %...from 1
      OptTable.ub = [3e+02         1e+04]; %...to 300
       Opt.SNR = 50;
       Opt.Nofrun = 5;
      % run simulation
      SimResults = Model.Sim_Sensitivity_Analysis(OptTable,Opt);
      figure('Name','Sensitivity Analysis');
      SimVaryPlot(SimResults, 'T2' ,'T2' );
      
      printf("Elapsed time=%f",elapsed_time)
%% VI- NOTES
% _No notes are available for this model._
% More information is available at https://qmrlab.readthedocs.io/en/master/mono_t2_batch.html
% Milford, D., et al. (2015). Mono-Exponential Fitting in T2-Relaxometry: Relevance of Offset and First Echo. PLoS One, 10(12), e0145255. 10.1371/journal.pone.0145255
