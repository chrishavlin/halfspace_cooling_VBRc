
% put VBR in the path
path_to_top_level_vbr='/home/chris/src/vbr';
addpath(path_to_top_level_vbr)
vbr_init

VBR.in.viscous.methods_list={'HZK2011'};
VBR.in.anelastic.methods_list={'xfit_premelt';};

VBR.in.SV.T_K=273 + linspace(800, 2000, 110)';% temperature [K]
sz = size(VBR.in.SV.T_K);
VBR.in.SV.phi = 1e-5 * ones(sz);
VBR.in.SV.dg_um=0.01*1e6* ones(sz); % grain size [um]
VBR.in.SV.P_GPa = 2 * ones(sz); % pressure [GPa]
VBR.in.SV.rho = 3300 * ones(sz); % density [kg m^-3]
VBR.in.SV.sig_MPa = .1 * ones(sz); % differential stress [MPa]
VBR.in.SV.chi = ones(sz); % composition factor
VBR.in.SV.f = [0.01]; % frequency
VBR.in.SV.Tsolidus_K=1200*ones(sz)+273; 

% with poroelastic, no small melt effect on viscosity 
VBR.in.elastic.methods_list={'anharmonic';'anh_poro'};
VBR = VBR_spine(VBR);


figure()
subplot(1,2,1)
plot(VBR.in.SV.T_K, VBR.out.anelastic.xfit_premelt.V/1e3)
hold on
plot(VBR.in.SV.Tsolidus_K, VBR.out.anelastic.xfit_premelt.V/1e3,'k')

subplot(1,2,2)
plot(VBR.in.SV.T_K./VBR.in.SV.Tsolidus_K,  VBR.out.anelastic.xfit_premelt.V/1e3)
  

figure()
subplot(1,2,1)
plot(VBR.in.SV.T_K, 1./VBR.out.anelastic.xfit_premelt.Q)
hold on
plot(VBR.in.SV.Tsolidus_K, 1./VBR.out.anelastic.xfit_premelt.Q,'k')

subplot(1,2,2)
plot(VBR.in.SV.T_K./VBR.in.SV.Tsolidus_K,  1./VBR.out.anelastic.xfit_premelt.Q)
  
