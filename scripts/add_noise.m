function [noisy_projections, sigmaNoise] = ...
    add_noise(projections, sigmaNoiseFraction)
    
    ref = std(projections(:));
    sigmaNoise = sigmaNoiseFraction * ref;
    noise = normrnd(0, sigmaNoise, size(projections));
    noisy_projections = projections + noise;
end