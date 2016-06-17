%script to simulate switching VAR model (SMS-VAR)

clc;
clear;

%generate same parameters
%rng(111);

%number of binary modes
Nm = 5; 

%number of phases
Nx = 3;

%dimensionality of VAR 
Ny = 4;

%length of time series
Lseq = 200;

%generate parameters
Dist = generateParamSVAR(Nx, Nm, Ny, Lseq);

%number of sequences to generate
Nseq = 50;

%generate data from SMS-VAR
[data, phase] = genSequenceSVAR(Nx, Nm, Ny, Nseq, Lseq, Dist);

%optional: anomalies can be injected into generated data for testing purposes
%currently only normal data is generated

disp('Data simulation done');

%number of processes to use: Nseq/Nproc should be an integer
Nproc = 5;

%number of EM iterations
NemIter = 5;

em = EM(data, Nx, 1, []);
em.run_parallel(Nproc, NemIter);

disp('Learning done');

%clean possible leftovers
system('rm -rf sim/res_*.mat');
system('rm -rf sim/data_*.mat');
system('rm -rf sim/output*.txt');

%evaluate sequences
em.evaluateSequences_parallel(Nproc);
    
disp('Evaluation done');

%the results are (res_1.mat, res_2.mat, ...) in folder sim/

%each file corresponds to each training sequence

%res.logLike - log likelihood for each time stamp in the sequence
%res.phaseErr - KL divergence for phase (see paper for details) for each
%time stamp in the sequence
    


