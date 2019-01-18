function OUT = duplicate_row(IN)
% duplicate row of matrix 
% 
% rand(3)
% 
% ans =
% 
%     0.9706    0.8003    0.9157
%     0.9572    0.1419    0.7922
%     0.4854    0.4218    0.9595
% 
% duplicate_row(ans)
% 
% ans =
% 
%     0.9706    0.8003    0.9157
%     0.9706    0.8003    0.9157
%     0.9572    0.1419    0.7922
%     0.9572    0.1419    0.7922
%     0.4854    0.4218    0.9595
%     0.4854    0.4218    0.9595

SZ = size(IN);
OUT= zeros(SZ(1)*2,SZ(2));
OUT(1:2:SZ(1)*2-1,:) = IN;
OUT(2:2:SZ(1)*2,:) = IN;

end