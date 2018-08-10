function [edgeCount, pixelCount, probability, cloudProb]=fr_prob_calculation(sstFiles) 

% function [edgeCount, pixelCount, probability, cloudProb]=fr_prob(sstFiles) 
%
% Calculate frontal probabilities for a given profile list, no
% matter what is in this list (yearly, monthly, etc.)
% 
% Originally from prob_mois_ATLANTIC.m. Should check if all output
% are still needed.
%    
% F. Cyr - March 2013
    
% Load 1st image to get size
load(sstFiles(1,:)) % load 1st image
    
% Matrices initalization    
edgeCount = zeros(size(image_overlay)); % fronts detection
pixelCount  = zeros(size(image_overlay)); % cloud/ice free pixel
probability= zeros(size(image_overlay)); % probability
cloudProb = zeros(size(image_overlay)); % cloub prob.
imageCount = size(sstFiles, 1);

for i = 1:size(sstFiles, 1) 
    disp(sprintf('File %d of %d', i, imageCount))
    data = load(sstFiles(i,:));
    image_overlay=data.image_overlay;   
    prob_input=data.image_overlay;

    I=(prob_input==-2.5);
    edgeCount = edgeCount + I;
    
    P=((prob_input~=-4)&(prob_input~=-5)); % Ignore ICE (-4) and CONTINENT (-5)
    pixelCount = pixelCount + P; 
end

% Front frequency
I = find(pixelCount ~= 0); % avoid division by zero
probability(I) = edgeCount(I)./pixelCount(I)*100;

% Ice/cloud frequency (zeros if no pixel detected)
cloudProb(I) = ((imageCount-pixelCount(I))./imageCount)*100;  
                