%%% This code replicates Arellano and Bond (1991, Table 1), Windmeijer (2005, Table
%%% 2)
%%% Dataset from STATA 'abdata', type .webuse abdata in stata command

clc;clear all;close all
tic
%%% Loading Data
dat = xlsread('abdata.xls',1,'B2:AN1032');  

year = dat(:,2);
y = dat(:,7);
w = dat(:,8);
k = dat(:,9);
ys = dat(:,10);
id = dat(:,13);
y_1 = dat(:,14);
y_2 = dat(:,15);
w_1 = dat(:,16);
k_1 = dat(:,17);
k_2 = dat(:,18);
ys_1 = dat(:,19);
ys_2 = dat(:,20);

yr_dummy =dat(:,21:29);
y_3 = dat(:,30);
w_2 = dat(:,31);
k_3 = dat(:,32);
ys_3 = dat(:,33);

n = size(y,1);

dy = y- y_1;
dy_1 =  y_1 - y_2;
dy_2 = y_2 -y_3;
dw = w - w_1;
dw_1 = w_1  - w_2;

dk = k - k_1;
dk_1 = k_1 -k_2;
dk_2 = k_2 - k_3;
dys = ys - ys_1;
dys_1 = ys_1- ys_2;
dys_2 = ys_2 - ys_3;
dyr_dummy =dat(:,34:39);


N = length(unique(id));
T = length(unique(year));

idx = repmat(id,1,N)==kron(ones(n,1),unique(id)');
ng = sum(idx)';

% Construct effective samples excluding first 3 rows in each individual
% due to lag specifications in model
idnx = idx;

for g=1:N
    
        idnx(sum(ng(1:g-1))+1:sum(ng(1:g-1))+3,g) = [0 0 0]';

end
nng = sum(idnx)';

% Construct Dependent Variables/Covariantes and Instruments

DY = [];
DX_AB = [];
DX_Wind = [];
mem = []; 
yearmem = [];
Z_AB =[];
Z_Wind = [];

for g=1:N
    
    
    Z_i = zeros(nng(g),(T-2)*(T-1)/2-1+T-3);
    Z_i_Wind = zeros(nng(g),(T-2)*(T-1)/2-1+T-3);
    
    DY = [DY; dy(idnx(:,g),:)];
    DX_AB = [DX_AB;  dy_1(idnx(:,g),:) dy_2(idnx(:,g),:) dw(idnx(:,g),:) dw_1(idnx(:,g),:) dk(idnx(:,g),:) dk_1(idnx(:,g),:) dk_2(idnx(:,g),:) dys(idnx(:,g),:) dys_1(idnx(:,g),:) dys_2(idnx(:,g),:) dyr_dummy(idnx(:,g),2:6) ones(nng(g),1)];
    DX_Wind = [DX_Wind;  dy_1(idnx(:,g),:) dy_2(idnx(:,g),:) dw(idnx(:,g),:) dw_1(idnx(:,g),:) dk(idnx(:,g),:) dys(idnx(:,g),:) dys_1(idnx(:,g),:) dyr_dummy(idnx(:,g),2:6) ones(nng(g),1)];
    mem = [mem; id(idnx(:,g),:)];
    yearmem = [yearmem; year(idnx(:,g),:)];
    
    t0 = year(idx(:,g),:);
    t0e = year(idnx(:,g),:);
    
   
% Constructions as in Arellano and Bond (1991), Windmeijer (2005)
   for s = 1:nng(g)

       t = t0e(s)-1976+1;
       z_i_temp = y(idx(:,g),:);

       if t0(1) == 1976

        Z_i(s,(t-3)*(t-2)/2+ t-4:(t-2)*(t-1)/2+t-4) = [1 z_i_temp(1:t-2)'];

       elseif t0(1) == 1977

       Z_i(s,(t-3)*(t-2)/2+ t-4:(t-2)*(t-1)/2+t-4) = [1 0 z_i_temp(1:t-3)'];

       else

       Z_i(s,(t-3)*(t-2)/2+ t-4:(t-2)*(t-1)/2+t-4) = [1 0 0 z_i_temp(1:t-4)'];

       end

   end
   
    Z_i_Wind = Z_i;
    Z_i = [Z_i dw(idnx(:,g),:) dw_1(idnx(:,g),:) dk(idnx(:,g),:) dk_1(idnx(:,g),:) dk_2(idnx(:,g),:) dys(idnx(:,g),:) dys_1(idnx(:,g),:) dys_2(idnx(:,g),:)];
    Z_AB = [Z_AB; Z_i];
    
    Z_i_Wind = [Z_i_Wind dw(idnx(:,g),:) dw_1(idnx(:,g),:) dk(idnx(:,g),:) dys(idnx(:,g),:) dys_1(idnx(:,g),:)];
    Z_Wind = [Z_Wind; Z_i_Wind];



% Different way of constructions    
%    
%   Z_i = zeros(nng(g),(T-2)*(T-1)/2-1);

%   for s = 1:nng(g);
%         
%        t = t0e(s)-1976+1;
%        z_i_temp = y(idx(:,g),:);
%        
%        if t0(1) == 1976
%            
%         Z_i(s,(t-3)*(t-2)/2:(t-2)*(t-1)/2-1) = [z_i_temp(1:t-2)'];
%         
%        elseif t0(1) == 1977
%            
%        Z_i(s,(t-3)*(t-2)/2+1:(t-2)*(t-1)/2-1) = [z_i_temp(1:t-3)'];
%            
%        else
%           
%        Z_i(s,(t-3)*(t-2)/2+2:(t-2)*(t-1)/2-1) = [z_i_temp(1:t-4)'];
%            
%        end
%     
%    end
%   
%     Z_i = [Z_i dw(idnx(:,g),:) dw_1(idnx(:,g),:) dk(idnx(:,g),:) dk_1(idnx(:,g),:) dk_2(idnx(:,g),:) dys(idnx(:,g),:) dys_1(idnx(:,g),:) dys_2(idnx(:,g),:) ones(nng(g),1)];
%     
%     Z = [Z; Z_i];
%    
    
    

end
toc

% Replicates Arellano and Bond (1991, Table 4)
[b2,sr,sw,s0,b2c,src,swc,s0c,b1,s1r,s10,J,pv] = twostep_gmm_wcrc_Dvector(DY,DX_AB,Z_AB,mem);
%[b2,sr,sw,s0,b2c,src,swc,s0c,b1,s1r,s10,J,pv] = twostep_gmm_wcrc_Dvector2(DY,DX_AB,Z_AB,mem);
[b_iter,s_iter,sc_iter,sw_iter,s0_iter,J,pv,iter] = iterated_gmm_wcrc(DY,DX_AB,Z_AB,mem);

Results_one_AB = [b1 s10 s1r  ]
Results_two_AB = [b2 s0 sr sw b2c s0c src swc]
Results_iter_AB = [b_iter s0_iter s_iter sw_iter]

% Replicates Windmeijer (2005, Table 2)
[b2,sr,sw,s0,b2c,src,swc,s0c,b1,s1r,s10,J,pv] = twostep_gmm_wcrc_Dvector(DY,DX_Wind,Z_Wind,mem);
%[b2,sr,sw,s0,b2c,src,swc,s0c,b1,s1r,s10,J,pv] = twostep_gmm_wcrc_Dvector2(DY,DX_Wind,Z_Wind,mem);
[b_iter,s_iter,sc_iter,sw_iter,s0_iter,J,pv,iter] = iterated_gmm_wcrc(DY,DX_Wind,Z_Wind,mem);

Results_one_Wind = [b1 s10 s1r  ]
Results_two_Wind = [b2 s0 sr sw b2c s0c src swc]
Results_iter_Wind = [b_iter s0_iter s_iter sw_iter]

Results_AB = [Results_one_AB(1:10,:) Results_two_AB(1:10,:) Results_iter_AB(1:10,:)];
Results_Wind = [Results_one_Wind(1:10,:) Results_two_Wind(1:10,:) Results_iter_Wind(1:10,:)];


% Exporting Tables

% columnLabels = {'Coeff','se','se-rc','Coeff','se', 'se-rc','se-wc','Coeff','se', 'se-rc','se-wc','Coeff','se', 'se-rc', 'se-rc-c', 'se-wc'};
% rowLabels = {'$n_{i,t-1}$','$n_{i,t-2}$','$w_{i,t}$','$w_{i,t-1}$','$k_{i,t}$','$k_{i,t-1}$','$k_{i,t-2}$','$ys_{i,t}$','$ys_{i,t-2}$','$n_{i,t-1}$'};
% 
% filename = (['Results_Arellano_Bond_Data_AB_table.tex']);   
% delete(filename);
% matrix2latex(Results_AB, filename, 'columnLabels', columnLabels,'rowLabels', rowLabels, 'alignment', 'c', 'format', '%-4.3f');
% 
% columnLabels = {'Coeff','se','se-rc','Coeff','se', 'se-rc','se-wc','Coeff','se', 'se-rc','se-wc','Coeff','se', 'se-rc', 'se-rc-c', 'se-wc'};
% rowLabels = {'$n_{i,t-1}$','$n_{i,t-2}$','$w_{i,t}$','$w_{i,t-1}$','$k_{i,t}$','$k_{i,t-1}$','$k_{i,t-2}$','$ys_{i,t}$','$ys_{i,t-2}$','$n_{i,t-1}$'};
% 
% filename = (['Results_Arellano_Bond_Data_Wind_table.tex']);   
% delete(filename);
% matrix2latex(Results_Wind, filename, 'columnLabels', columnLabels,'rowLabels', rowLabels, 'alignment', 'c', 'format', '%-4.3f');

% Results = zeros(40,8);
% for j = 1:10
%  Results(4*(j-1)+1:4*j-1,1) = Results_one_AB(j,:)';
%  Results(4*(j-1)+1:4*j,2) = Results_two_AB(j,1:4)';
%  Results(4*(j-1)+1:4*j,3) = Results_two_AB(j,5:end)';
%  Results(4*(j-1)+1:4*j,4) = Results_iter_AB(j,:)';
%  
%  Results(4*(j-1)+1:4*j-1,5) = Results_one_Wind(j,:)';
%  Results(4*(j-1)+1:4*j,6) = Results_two_Wind(j,1:4)';
%  Results(4*(j-1)+1:4*j,7) = Results_two_Wind(j,5:end)';
%  Results(4*(j-1)+1:4*j,8) = Results_iter_Wind(j,:)';
% end
% 
% columnLabels = {'one-step','two-step','two-step(centered)','iterated','one-step','two-step','two-step(centered)','iterated'};
% 
% filename = (['Results_table.tex']);   
% delete(filename);
% matrix2latex(Results, filename, 'alignment', 'c', 'format', '%-4.3f');
