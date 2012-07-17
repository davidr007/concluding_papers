function [Hypo_loc_out] = do_flips4inversion(Hypo_loc_in)
% 
% This routine flips the local coordinates for Hypo_loc_in events to ensure that
% that x2>0, y3>0 and z4>0 (i.e. we try to solve the inversion in this
% space)
%
% INPUTS:
% Hypo_loc_in   [double nEx3]
%               Events in a local coordinate system
%                 Hypo_loc_in=  [0      0       0;
%                                x2     0       0;      x2 = +/-
%                                x3     y3      0;      y3 = +/-
%                                x4     y4      z4;     z4 = +/-
%                                x5     y5      z5; 
%                                ...    ...     ...;
%                                xnE    ynE     znE];
%
% Hypo_loc_out  [double nEx3]
%               Events in a local coordinate system with x2>0, y3>0 and
%               z4>0
%                 Hypo_loc_in=  [0      0       0;
%                                x2     0       0;      x2 >0
%                                x3     y3      0;      y3 >0
%                                x4     y4      z4;     z4 >0
%                                x5     y5      z5; 
%                                ...    ...     ...;
%                                xnE    ynE     znE];
%
%
% David Robinson 
% 4 August 2008

Hypo_loc_out = Hypo_loc_in;

% Flip according to x2
if Hypo_loc_out(2,1) <0
    Hypo_loc_out(:,1) = -Hypo_loc_out(:,1);
end

% Flip according to y3
if Hypo_loc_out(3,2) <0
    Hypo_loc_out(:,2) = -Hypo_loc_out(:,2);
end
                    
% Flip according to z4
if Hypo_loc_out(4,3) <0
    Hypo_loc_out(:,3) = -Hypo_loc_out(:,3);
end
