function theta = FlipWindDir(WindDir)
% -------------------------------------------------------------------------
% Function for flipping wind directions 180 derees.  
% 
% INPUT: - WindDir: wind direction where the wind is blowing towards (i.e. 
%          wind dir 315 is wind blowing towards NE)
% OUTPUT:- theta: windDir flipped 180 degrees.
% 
% Last edited: 04.May.2018, Torgeir
% -------------------------------------------------------------------------

% Indexes corresponding to angles > 180
idx = find(WindDir > 180);   
theta(idx) = WindDir(idx) - 180;

% Indexes corresponding to angles < 180 
idx = find(WindDir < 180);  
theta(idx) = WindDir(idx) + 180;

