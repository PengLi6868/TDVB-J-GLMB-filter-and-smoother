function [w_new,x_new,P_new,tzhi_new]= gaus_cap2(w,x,P,max_number,tzhi)

if length(w) > max_number
    [notused,idx]= sort(w,1,'descend');
    w_new= w(idx(1:max_number)); w_new = w_new * (sum(w)/sum(w_new));
    tzhi_new= tzhi(idx(1:max_number));
    x_new= x(:,idx(1:max_number));
    P_new= P(:,:,idx(1:max_number));
else
    x_new = x;
    tzhi_new = tzhi;
    P_new = P;
    w_new = w;
end
