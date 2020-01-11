function [The_track,J] = Analys_track(est,T)

% This model is used to rewrite the "est" into the trajectory form

N = 0;
for i =1 : T
    N = N + size(est.X{i,1},2);
end

n = 0;
J = 0;
while ( n < N )
    J = J + 1;
    track_flag = 0;
    L_now = [];
    The_track(J).m = zeros(4,T);
    The_track(J).K = [];
    The_track(J).Z = zeros(2,T);
    for k = 1 : T
        
        if size(est.L{k,1},2) == 0
            continue;
        end
        
        if track_flag == 0
            L_now = est.L{k,1}(:,1);
            The_track(J).m(:,k) = est.X{k,1}(:,1);
            The_track(J).K = k;
            The_track(J).Z(:,k) = est.Z{k,1}(:,1);
            est.L{k,1}(:,1) = [];
            est.Z{k,1}(:,1) = [];
            est.X{k,1}(:,1) = [];
            n = n + 1;
            track_flag = 1;
            continue;
        end
        
        for i = 1 : size(est.L{k,1},2)
            if L_now(1) == est.L{k,1}(1,i) && L_now(2) == est.L{k,1}(2,i)
                The_track(J).m(:,k) = est.X{k,1}(:,i);
                The_track(J).Z(:,k) = est.Z{k,1}(:,i);
                est.L{k,1}(:,i) = [];
                est.X{k,1}(:,i) = [];
                est.Z{k,1}(:,i) = [];
                n = n + 1;
                break;
            end
        end

    end
    
    for k = 1 : T
        if sum( The_track(J).m(:,k),1 ) ~= 0
            dela1 = k - 1;
            break;
        end
    end

    for k = T :-1: 1
         if sum( The_track(J).m(:,k),1 ) ~= 0
            dela2 = k + 1;
            break;
         end
    end
    if dela2 < T+1
        The_track(J).m(:,dela2:T ) = [];
        The_track(J).Z(:,dela2:T) = [];
    end
    if dela1 > 0
        The_track(J).m(:,1:dela1) = [];
        The_track(J).Z(:,1:dela1) = [];
    end
    The_track(J).K(2) = The_track(J).K(1) + size(The_track(J).m,2) - 1;
    

end
            
    