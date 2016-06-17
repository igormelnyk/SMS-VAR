# SMS-VAR
Semi-Markov Switching Vector Autoregressive (SMS-VAR) Model for Anomaly Detection in Aviation Systems

Demo version of the code to generate, train and evaluate SMS-VAR model. 

Run Matlab script runSim.m to execute the code. (Make sure the folder which contains the code has a subfolder named ‘sim’)

---Data---

The generated data consists of two parts: ModeData and SensorData. ModeData(1,:) is the time series of discrete modes (there are 2^Nm possible modes) and ModeData(2,:) is the time series of the mode durations. SensorData is Ny-dimensional time series of continuous data.
Optionally, anomalies can be injected into generated data for testing purposes. However, currently only normal data is generated.


---Training---

Parameter learning for SMS-VAR is done using EM. The current version of the code runs training (EM) and evaluation in parallel by starting several Matlab processes and dispatching a part of job to each of them. Please make sure that a command 'matlab' can be executed on the system which can start a new Matlab process. I tested it in Linux and it works fine. 


---Evaluation---

The results of execution are saved to folder sim. Each file corresponds to evaluating the built SMS-VAR model on each of the training sequence.

res.logLike - log likelihood for each time stamp in the sequence

res.phaseErr - KL divergence for phase (see paper for details) for each time stamp in the sequence.


The code is based on paper: I. Melnyk, A. Banerjee, B. Matthews, and N. Oza. Semi-Markov Switching Vector Autoregressive Model-based Anomaly Detection in Aviation Systems. KDD, 2016.

Please email melnyk@cs.umn.edu for questions or comments.
