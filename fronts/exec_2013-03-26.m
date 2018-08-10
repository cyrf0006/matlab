%During the night, the script should run the following...

fr_mean_prob('SST_atlantic', '~/data/front_data/AtlanticLatLon.mat', [1986:2010], [1:12], 'allMean')
fr_mean_prob('SST_atlantic', '~/data/front_data/AtlanticLatLon.mat', [1986:2010], [1:12], 'monthly')

fr_mean_prob('SST_hudson', '~/data/front_data/HudsonLatLon.mat', [1986:2010], [1:12], 'monthly')
fr_mean_prob('SST_baffin', '~/data/front_data/BaffinLatLon.mat', [1986:2010], [1:12], 'monthly')
fr_mean_prob('SST_pacific', '~/data/front_data/PacificLatLon.mat', [1986:2010], [1:12], 'monthly')

