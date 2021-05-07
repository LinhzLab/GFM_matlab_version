function [measure] = measurefun(H, hH, type)
% function to evaluate the smallest nonzero canonical correlation between
% two set of variables.
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         May. 25, 2019
% Copyright (c) 2019, Liu Wei
% All rights reserved.

if(~exist('type', 'var'))
    type = 'canonical';
end
q = size(H, 2);
switch type 
    case 'canonical'
        [~, ~, r] = canoncorr(hH,H);
        measure = r(end);
    case 'Fnorm'
        measure = norm(H- hH, 'fro')^2/ numel(H);
end