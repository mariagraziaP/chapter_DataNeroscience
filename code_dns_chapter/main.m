
%% Load data

clc; clear; close all;

load('HCP_scfc_av.mat');

SC = sc_av;
FC = fc_av;
clear sc_av fc_av;

% Number of nodes
n = size(SC,1);

%% Compute network meausres

%%% Topological specialization

% Matching index
MI = matching_index_wei2(SC);
MI = MI + MI'; % Reflect to lower triangle

% Cosine similarity
CS = 1 - squareform(pdist(SC, 'cosine'));

%%% Modularity

gamma_sample = 1000;
in_gamma_low = 0.001;
in_gamma_high = 10;
[cons, CC, ~, ~, ciall, gam_range] = get_SCmodules(SC, gamma_sample, gamma_sample, n, in_gamma_low, in_gamma_high, 0);

%%% Communication

% Weight-to-length transformation
% This is necessary for the computation of shortest paths
L = 1./SC; 
L(L == inf) = 0;

% Shortest path length
SPL = distance_wei_floyd(L);

% Search information
SI = search_information(SC,L);
SI = (SI + SI')./2; % symmetrize

% Mean first passage time
S = diag(SC*ones(n,1)); % diagonal matrix of node strengths
SC_norm = S^(-1/2)*SC*S^(-1/2); % normalize SC to attenuate high-strength nodes
MFPT = mean_first_passage_time(SC_norm);

% MFPT = mean_first_passage_time(SC);
MFPT = (MFPT + MFPT')./2; % symmetrize

ind = find(triu(MFPT,1));
MFPTzsc = zscore(MFPT(ind));

newMFPT(ind) = MFPTzsc; 
newMFPT = newMFPT + newMFPT';

%% Matrix visualizations

figure;

subplot(2,4,1);
imagesc(log10(SC));
title('Structural connectivity (SC)');

subplot(2,4,2);
imagesc(FC);
title('Functional connectivity (FC)');

subplot(2,4,3);
imagesc(MI);
title('Macthing index (MI)');

subplot(2,4,4);
imagesc(CS);
title('Cosine similarity (CS)');

subplot(2,4,5);
imagesc(CC);
title('Modular co-classification (CC)');

subplot(2,4,6);
imagesc(SPL);
title('Shortest path length (SPL)');

subplot(2,4,7);
imagesc(SI);
title('Search information (SI)');

subplot(2,4,8);
imagesc(MFPT);
title('Mean first passage time (MFPT)');

%% Correlations with FC

corr_type = 'S';

% Upper triangle index
% idx = find(triu(ones(n),1));
idx = find(triu([[ones(n/2), zeros(n/2)]; [zeros(n/2), ones(n/2)]],1));

% Matching index
fc_corr_MI = corr(FC(idx), MI(idx), 'Type', corr_type);
% Co-classification matrix
fc_corr_CC = corr(FC(idx), CC(idx), 'Type', corr_type);
% Search information
fc_corr_SI = corr(FC(idx), SI(idx), 'Type', corr_type); 





