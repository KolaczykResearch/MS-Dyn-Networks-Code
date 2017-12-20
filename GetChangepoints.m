%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Name: GetBreakpoints.m        Created: 07/28/00    Revised: 09/06/01
%
%% Usage: Extracts an optimal segmentation from probability arrays and
%         mask calculated with TriangleMDN.m  
%         NOTE: Function is defined in an implicit manner.
%
%% Inputs:  bp := vector of breakpoints (i.e., those calculated so far
%                 in moving through the array).
%            i := index of lefthand side of subinterval of x
%                 for which probabilities are to be compared
%            j := index of righthand side of subinterval of x
%                 for which probabilities are to be compared
%           t0 := triangular array of probability values under null
%                 hypothesis (i.e., no split)
%           t1 := triangular array of probability values under alternative
%                 hypothesis (i.e., split).
%           mm := mask indicating (with a '1') the decided split locations
%		          in the sequence (in relation to the arrays t0, t1).
%
%% Output:  bp := vector of breakpoints (i.e., including those just
%                calculated).
%
%% Calls:  Only internal Matlab functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bp = GetChangepoints(bp, i, j, t0, t1, mm)

if t0(i, j) > t1(i, j)  % ">" because we are doing minimization

   bp = [bp mm(i, j)];
   bp = GetChangepoints(bp, i, mm(i, j), t0, t1, mm);
   bp = GetChangepoints(bp, mm(i, j) + 1, j, t0, t1, mm);

end

