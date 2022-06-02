# CEE362G_HBM
CEE362G final project

For spatial and temporal release of source reconstruction,

1. simulate data using create_2D_data.m. You can change number of detectors and locations (D_loc) easily, and this file will output your captured data at these detector locations in methane_data.txt. Currently, the velocity is 0 for simplicity. You can also change source locations (S_loc) and their release patterns, but to add sources, you'll have to modify a bit more, but pattern matching makes it easy to see what needs to be modified. We then add noise to the datapoints. The format of the data saved goes as follows: 1st column are the times of our measurements. Subsequent columns are the measurments at the detectors. Each column corresponds to a different detector.

2. after simulating the data, use air_pollution.m to reconstruct the time release signals at the 4 different sources. Just make sure the D_loc and S_loc matches with those in the create_2D_data.m file. 
