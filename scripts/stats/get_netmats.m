clear all
%% compute modularity and graph metrics that can be used for predictive modelling

%% define paths

TR = 2.00;
addpath /home/ALDRECENTRUM/benjamin.garzon/Software/FSLnets_pack/FSLNets/
workdir = '/home/share/LeftHand/LHconnectivity/Berlin/melodic100.gica/';
bad_nets_file = 'bad_nets.txt';
ncomp = 100;

workdir = '/home/share/LeftHand/LHconnectivity/Berlin/melodic.gica/';
bad_nets_file = 'noninterestingnets.txt';
ncomp = 42;

cd(workdir)
group_maps=fullfile(workdir, 'groupmelodic.ica/melodic_IC');
ts_dir=fullfile(workdir, 'grot');


%% load data
ts=nets_load(ts_dir,TR,1);
bad_nets = load(bad_nets_file)';
ts.DD=setdiff([1:ncomp], bad_nets);

%% compute adjacencies and save
ts=nets_tsclean(ts,1);

netmats1 = nets_netmats(ts,1,'corr');       % full correlation (normalised covariances)
netmats5 = nets_netmats(ts,1,'ridgep');     % Ridge Regression partial, with rho=0.01

[Znet1,Mnet1]=nets_groupmean(netmats1,0); % full corr
[Znet5,Mnet5]=nets_groupmean(netmats5,0); % ridge

%mytable = readtable(fullfile(workdir, '../ICA_table.csv'));
netmats1_table = reshape_nets(netmats1);
dlmwrite(fullfile(workdir, 'netmats_full_reduced.csv'), netmats1_table);

netmats5_table = reshape_nets(netmats5);
dlmwrite(fullfile(workdir, 'netmats_ridge_reduced.csv'), netmats5_table);


%coefs = load('results/IQ_coefs.mat.txt');
%coefs_full = zeros(ncomp,ncomp);
%coefs_full(ts.DD, ts.DD) = coefs;
%nets_netweb(Znet1,coefs_full,ts.DD,group_maps,'netweb');

