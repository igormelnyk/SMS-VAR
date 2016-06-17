%function to simulate switching VAR model with independent binary modes

function [data, phase] = genSequenceSVAR(Nx, Nm, Ny, Nseq, Lseq, Dist)

data = cell(Nseq, 1);
phase = cell(Nseq, 1);

numModes = Dist.numModes;

for s=1:Nseq
    
    Dmodes = zeros(2, Lseq+1);
    Dsensors = zeros(Ny, Lseq+1);
    Phase = zeros(1, Lseq+1);
    
    
    %init
    Dmodes(1,1) = randi(numModes);
    Dmodes(2,1,:) = 1;
    Dsensors(:,1) = rand(Ny, 1);
    Phase(1) = randi(Nx);
    
    
    for i=2:Lseq+1
        
        %========== modes ============
        if Dmodes(2, i-1) == 1 %duration of prev mode expired
            
            %figure out what is the next mode state
            w = Dist.modeTransDist.transitions(:,Dmodes(1, i-1));
            Dmodes(1, i) = randsample(1:numModes, 1, true, full(w)); 
            
            %figure out what is the duration of next state (Poisson dist)
            mu = Dist.modeDurDist.durations(Dmodes(1, i));
            
            %make sure new state duration is not zero
            while(1)
                Dmodes(2, i) = poissrnd(mu);
                if Dmodes(2, i) ~= 0
                    break;
                end
            end
            
        else %same mode will continue since duration is not expired yet
                       
            Dmodes(1, i) = Dmodes(1, i-1);
            
            %decrement duration counter by one
            Dmodes(2, i) = Dmodes(2, i-1)-1;
        end
        
        %=========== phase =============
        %for specific combination of modes return p(xt|xt_1, mt, dt_1)
        mst = Dmodes(1, i);
        durs = Dmodes(2, i-1);
        
        %phase is originally in LOG scale, so get original scale
        T = exp(Dist.phaseTransDist.getValue(mst, durs));
        
        %new phase
        Phase(i) = randsample(1:Nx, 1, true, T(:, Phase(i-1)));
        
        %============ observations ============
        A = Dist.obsTransDist.As(:, :, Phase(i));
        mu = A*Dsensors(:, i-1);
        Sigma = eye(Ny);
        Dsensors(:, i) = mvnrnd(mu, Sigma);
    end
    
    Dmodes = Dmodes(:,2:end);
    Dsensors = Dsensors(:,2:end);
    Phase = Phase(2:end);
    
    d.ModeData = Dmodes;
    d.SensorData = Dsensors;    
    d.numModes = Nm;
    data{s} = d;
    
    phase{s} = Phase;    
end




















