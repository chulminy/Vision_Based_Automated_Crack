function [posMat, operMat] = ...
                            CompHaarPatch(position, wid_width, wid_height)

% position    : right-top position
% wid_width   : width of a window
% wid_height  : height of a window

% -------------------------------------------------------------------------
%                 |---------|                       |--------|--------|
%   |---------|   |||||||||||  |--------|-------|   |        ||||||||||
%   |||||||||||   |||||||||||  ||||||||||       |   |        ||||||||||
%   |||||||||||   |---------|  ||||||||||       |   |--------|--------|
%   |||||||||||   |         |  |--------|-------|   ||||||||||        |
%   |---------|   |         |                       ||||||||||        |
%                 |---------|                       |--------|--------|
%     Type 1         Type 2          Type 3               Type 4
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%   |----------|              
%   |||||||||||| 
%   ||||||||||||    |---------|---------|---------|
%   |----------|    |||||||||||         |||||||||||
%   |          |    |||||||||||         |||||||||||
%   |          |    |||||||||||         |||||||||||
%   |----------|    |---------|---------|---------|              
%   |||||||||||| 
%   ||||||||||||
%   |----------|              
%      Type 5                   Type 6                           
% -------------------------------------------------------------------------

% Example
% A = ones(20,20);
% A(1:10,1:10)   = 0;
% A(11:20,11:20) = 0;
% [numCorner operMat posMat] = CM_CompHaarPatch([1 1], 20, 20);
% intA = integralImage(A);
% intVal = 0;
% for ii=1:numCorner; 
%   intVal = intVal + intA(posMat{ii}(1), posMat{ii}(2))*operMat(ii);
% end
 
% set patch sizes (even number)
eps = 0.000001;
patch_w = fix(wid_width * (rand(1)-eps)); 
if rem(patch_w,2)==1; patch_w = patch_w + 1; end
if fix(patch_w/3)< 2; patch_w = 6; end 
% minimum patch_w is 6;

patch_h = fix(wid_height * (rand(1)-eps)); 
if rem(patch_h,2)==1; patch_h = patch_h + 1; end
if fix(patch_h/3)< 2; patch_h = 6; end 
% minimum patch_h is 6;

% set a patch position
width_range  = 0: (wid_width-patch_w);
height_range = 0: (wid_height-patch_h);

rt_pos = [datasample(width_range,1) datasample(height_range,1)];

% set a patch type
patchType = datasample(1:6,1);

numCorner   = 0;
operMat     = [];
posMat      = [];
if (patchType == 1)
    numCorner = 4;
    posMat = zeros(numCorner,2);
    A = position + rt_pos;
    B = A + [patch_w 0];
    C = A + [0 patch_h];
    D = A + [patch_w patch_h];
    posMat(1:numCorner,:) = [A;B;C;D];
    
    %(A+D-B-C)
    operMat = [1 -1 -1 1]'; 
elseif (patchType == 2)
    numCorner = 6;
    posMat = zeros(numCorner,2);
    A = position + rt_pos;
    B = A + [patch_w 0];
    C = A + [0 patch_h/2];
    D = A + [patch_w patch_h/2];
    E = A + [0 patch_h];
    F = A + [patch_w patch_h]; 
    posMat(1:numCorner,:) = [A;B;C;D;E;F];
    
    %(A+D-B-C)-(C+F-D-E)
    operMat = [1 -1 -2 2 1 -1]'; 
elseif (patchType == 3)
    numCorner = 6;
    posMat = zeros(numCorner,2);
    A = position + rt_pos;
    B = A + [patch_w/2 0];
    C = A + [patch_w 0];
    D = A + [0 patch_h];
    E = A + [patch_w/2 patch_h];
    F = A + [patch_w patch_h]; 
    posMat(1:numCorner,:) = [A;B;C;D;E;F];
    
    %(A+E-B-D)-(B+F-C-E)
    operMat = [1 -2 1 -1 2 -1]';     
elseif (patchType == 4)
    numCorner = 9;
    posMat = zeros(numCorner,2);
    A = position + rt_pos;
    B = A + [patch_w/2 0];
    C = A + [patch_w 0 ];
    D = A + [0 patch_h/2];
    E = A + [patch_w/2 patch_h/2];
    F = A + [patch_w patch_h/2]; 
    G = A + [0 patch_h]; 
    H = A + [patch_w/2 patch_h]; 
    I = A + [patch_w patch_h]; 
    posMat(1:numCorner,:) = [A;B;C;D;E;F;G;H;I];
    
    % -(A+E-B-D)+(B+F-C-E)+(D+H-G-E)-(E+I-F-H)
    operMat = [-1 2 -1 2 -4 2 -1 2 -1]';     
elseif (patchType == 5)
    numCorner = 8;
    posMat = zeros(numCorner,2);
    A = position + rt_pos;
    patch_h = patch_h - mod(patch_h,3);
    
    B = A + [patch_w 0];
    C = A + [0 patch_h/3];
    D = A + [patch_w patch_h/3];
    E = A + [0 2*patch_h/3];
    F = A + [patch_w 2*patch_h/3]; 
    G = A + [0 patch_h]; 
    H = A + [patch_w patch_h]; 
    
    posMat(1:numCorner,:) = [A;B;C;D;E;F;G;H];
    
    % (A+D-B-C)-(C+F-D-E)+(E+H-G-F)
    operMat = [1 -1 -2 2 2 -2 -1 1]';     
else
    numCorner = 8;
    posMat = zeros(numCorner,2);
    A = position + rt_pos;
    patch_w = patch_w - mod(patch_w,3);
    
    B = A + [1*patch_w/3 0];
    C = A + [2*patch_w/3 0];
    D = A + [3*patch_w/3 0];
    E = A + [0            patch_h];
    F = A + [1*patch_w/3  patch_h];
    G = A + [2*patch_w/3  patch_h];
    H = A + [3*patch_w/3  patch_h]; 
    
    posMat(1:numCorner,:) = [A;B;C;D;E;F;G;H];
    
    % (A+F-B-E)-(B+G-C-F)+(C+H-D-G)
    operMat = [1 -2 2 -1 -1 2 -2 1]';         
end

