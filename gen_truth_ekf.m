function truth= gen_truth_ekf(model)

%variables
truth.K= 100;                   %length of data/number of scans
truth.X= cell(truth.K,1);             %ground truth for states of targets  
truth.N= zeros(truth.K,1);            %ground truth for number of targets
truth.L= cell(truth.K,1);             %ground truth for labels of targets (k,i)
truth.track_list= cell(truth.K,1);    %absolute index target identities (plotting)
truth.total_tracks= 0;          %total number of appearing tracks

%target initial states and birth/death times
nbirths= 10;
wturn = 2*pi/180;

xstart(:,1)  = [ 1000+1; -10; 1500-1; -10; wturn/8 ];        tbirth(1)  = 1;     tdeath(1)  = truth.K+1;
xstart(:,2)  = [ -250-1;  20; 1000+1; 3; -wturn/3 ];         tbirth(2)  = 10;    tdeath(2)  = truth.K+1;
xstart(:,3)  = [ -1500-1; 11; 250+1; 10; -wturn/2 ];          tbirth(3)  = 10;    tdeath(3)  = truth.K+1;
xstart(:,4)  = [ -1500; 20; 250; 0; 0 ];                                tbirth(4)  = 10;    tdeath(4)  = 66;
xstart(:,5)  = [ 250-1; 11; 750-1; 5; wturn/4 ];             tbirth(5)  = 20;    tdeath(5)  = 80;
xstart(:,6)  = [ -250+1; -12; 1000-1; -12; wturn/2 ];         tbirth(6)  = 40;    tdeath(6)  = truth.K+1;
xstart(:,7)  = [ 1000; 0; 1500; -10; wturn/4 ];                         tbirth(7)  = 40;    tdeath(7)  = truth.K+1;
xstart(:,8)  = [ 250; -25; 750; 0; -wturn/4 ];                          tbirth(8)  = 40;    tdeath(8)  = 80;
xstart(:,9)  = [ 1000; -25; 1500; 0; -wturn/4 ];                        tbirth(9)  = 60;     tdeath(9)  = truth.K+1;
xstart(:,10)  = [ 250; -20; 750; 15; wturn/4 ];                         tbirth(10)  = 60;    tdeath(10)  = truth.K+1;


%generate the tracks
for targetnum=1:nbirths
    targetstate = xstart(:,targetnum);
    for k=tbirth(targetnum):min(tdeath(targetnum),truth.K)
        truth.X{k}= [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} targetnum];
        truth.N(k) = truth.N(k) + 1;
        targetstate = gen_newstate_fn_ekf(model,targetstate,'noiseless');
     end
end
truth.total_tracks= nbirths;

