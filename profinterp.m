function [griddeddat] = profinterp(data, loc, grid, method)
%this fuction interpolates data from a profile or a set of profiles 
%(interpolating down each column separately) onto a grid
%unlike MATLAB's interp fuction, it will remove NaNs from the data and do a
%bunch of profiles at once. 

%Inputs: data, data locations, grid to interpolate onto, and method(optional)
%data and data locations must be the same size
%Outputs: data interpolated onto the grid

%allowable methods are those used in MATLAB's interp function
%'linear'| 'nearest' | 'next' | 'previous' | 'pchip' | 'cubic' | 'v5cubic' | 'makima' | 'spline'
%if no method is chosen, set it to be 'pchip' 

%Sophie Shapiro 2021

if ~exist('method','var')
    method = 'pchip';
end 

%make sure we have the right number of inputs and outputs
if ~(nargin == 3 | nargin == 4)
   error('profinterp:  Requires three or four inputs')
end 

if ~(nargout == 1)
   error('profinterp:  Requires one outputs')
end 

%make sure the inputs data and loc are both the same size
[md,nd] = size(data);
[ml,nl] = size(loc);

if (md ~= ml | nd ~= nl)
    error('profinterp: data and loc must have same dimensions')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Now we do the calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for j=1:size(data, 2) %number of profiles = num cols
        [~, ii, ~] = unique(loc(:, j)); %sort and remove duplicates
        dat=data(ii, j); lo=loc(ii,j);
        nvar = isfinite(dat); %get the indices of non-nan values in variable
        nangrids = isfinite(lo); %get the indices of non-nan values in grid var
        ns =and(nvar, nangrids)'; %get a nan index which is zero if either var or grid has a nan in that row
        if sum(ns) ~= 0 && sum(ns) ~= 1 %if the profile is all nans (nans has all 0s), skip it
            griddeddat(:, j)=interp1(lo(ns),dat(ns), grid, method, NaN);
            kk=find(grid>max(lo));%grab the max pr for mask edges
            griddeddat(kk,j)=NaN; %set NaNs below max pr
        else %set an all nan col to all nans 
            griddeddat(:, j) = NaN(size(grid, 2), 1);
        end
end