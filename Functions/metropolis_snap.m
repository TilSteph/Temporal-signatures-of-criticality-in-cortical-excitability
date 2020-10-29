function [spin_snaps,Mmean, Mean_single_spin, lin_ix_save] = metropolis_snap(spin, kT, J, numIters, talk, flip_order)
% metropolis algorithm (from Mathworks)
% takes "snapshots" every numel(spin) iteration
% Adpated: TS, 12/2018

% Metropolis algorithm
spin_snaps = zeros(size(spin,1), size(spin,2), numIters/numel(spin)); % initialize snapshot variable
Mean_single_spin = zeros(numIters, 1);
lin_ix_save = zeros(numIters, 1);
Mmean = zeros(numIters/numel(spin), 1);
snap_count = 0;

for iter = 1 : numIters
    
    % Pick a random spin
    if exist('flip_order','var') % pseudo-random
        linearIndex = flip_order(iter);
    else % pick randomly
        linearIndex = randi(numel(spin));
    end
    lin_ix_save(iter) = linearIndex;
    [row, col]  = ind2sub(size(spin), linearIndex);
    
    % Find its nearest neighbors
    above = mod(row - 1 - 1, size(spin,1)) + 1;
    below = mod(row + 1 - 1, size(spin,1)) + 1;
    left  = mod(col - 1 - 1, size(spin,2)) + 1;
    right = mod(col + 1 - 1, size(spin,2)) + 1;
    
    neighbors = [      spin(above,col);
        spin(row,left);                spin(row,right);
                       spin(below,col)];
    
    % Calculate energy change if this spin is flipped
    dE = 2 * J * spin(row, col) * sum(neighbors);
    
    % Boltzmann probability of flipping
    prob = exp(-dE / kT);
    
    % Spin flip condition
    if dE <= 0 || rand() <= prob
        spin(row, col) = - spin(row, col);
    end
    
    Mean_single_spin(iter) = mean(spin(:));
    
    % take snapshot of current grid
    if mod(iter, numel(spin))==0 % after every loop through all spins
        snap_count = snap_count + 1;
        spin_snaps(:,:,snap_count) = spin;
        
        % mean magnetization
        Mmean(snap_count) = mean(spin(:));
    end
    
    if talk % logical
        % talk to the operator every 100th iteration
        if mod(iter, round(numIters/100)) == 0, disp(['Iteration #' num2str(iter) ' (' num2str(round(iter/(round(numIters/100)))) '%)']), end
    end
end
end

