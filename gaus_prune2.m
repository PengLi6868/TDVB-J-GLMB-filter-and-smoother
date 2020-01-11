function [w_new,x_new,P_new,tzhi_new]= gaus_prune2(w,x,P,elim_threshold,tzhi)

idx= find( w > elim_threshold );
w_new= w(idx);
tzhi_new = tzhi(idx);
x_new= x(:,idx);
P_new= P(:,:,idx);