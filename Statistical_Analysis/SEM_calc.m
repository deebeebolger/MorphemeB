function sem=SEM_calc(vect, CI)
% SEM_calc - standard error of the mean, confidence interval
%
% function sem=SEM_calc(vect, CI)
%
%***********************************************************************

error(nargchk(1,2,nargin))

if isvector(vect)
    vect=vect(:);
end

if nargin==1
    stdCI = 1.96 ; % if the interval specified is 5%
elseif nargin==2  %if a different interval is specified
    CI = CI/2 ; %Convert to 2-tail
    stdCI = abs(norminv(CI,0,1)) ;  %uses inverse of the normal cdf.
end

sem = ( (std(vect)) ./ sqrt(sum(~isnan(vect))) ) * stdCI ;    % ignores nans
end