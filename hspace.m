%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CB_002_2D_HalfSpaceCooling.m
%
%  Calculate seismic properties for a half space cooling model, compares
%  results of two anelastic methods.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% put VBR in the path %%
  clear
  path_to_top_level_vbr='/home/chris/src/vbr';
  addpath(path_to_top_level_vbr)
  vbr_init

%% Set Thermodynamic State Using Half Space Cooling Model %%
%%   analytical solution:
%%     T(z,t)=T_surf + (T_asth - T_surf) * erf(z / (2sqrt(Kappa * t)))
%%   variables defined below

  % HF settings
  HF.Tsurf_C=0; % surface temperature [C]
  HF.Tasth_C=1400; % asthenosphere temperature [C]
  HF.V_cmyr=8; % half spreading rate [cm/yr]
  HF.Kappa=1e-6; % thermal diffusivity [m^2/s]
  HF.rho=3300; % density [kg/m3]
  HF.t_Myr=linspace(5,150,50)+1e-12; % seaflor age [Myrs]
  HF.z_km=linspace(0,350,50)'; % depth, opposite vector orientation [km]

  % HF calculations
  HF.s_in_yr=(3600*24*365); % seconds in a year [s]
  HF.t_s=HF.t_Myr*1e6*HF.s_in_yr; % plate age [s]
  HF.x_km=HF.t_s / (HF.V_cmyr / HF.s_in_yr / 100) / 1000; % distance from ridge [km]

  % calculate HF cooling model for each plate age
  HF.dT=HF.Tasth_C-HF.Tsurf_C;
  HF.T_C=zeros(numel(HF.z_km),numel(HF.x_km));
  for HFi_t = 1:numel(HF.t_s)
    HF.erf_arg=HF.z_km*1000/(2*sqrt(HF.Kappa*HF.t_s(HFi_t)));
    HF.T_C(:,HFi_t)=HF.Tsurf_C+HF.dT * erf(HF.erf_arg);
  end

% Load and set VBR parameters %
  VBR.in.elastic.methods_list={'anharmonic'};
  VBR.in.viscous.methods_list={'HK2003'};
  VBR.in.anelastic.methods_list={'andrade_psp'; 'xfit_mxw'; 'xfit_premelt'};
  VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); % unrelaxed elasticity
  VBR.in.elastic.anharmonic.Gu_0_ol = 75.5; % olivine reference shear modulus [GPa]
  VBR.in.SV.f = [0.01, 0.02, 0.04, 0.1];%  frequencies to calculate at

% copy halfspace model into VBR state variables, adjust units as needed
  dTdzad = 0.5; % adiabat degrees/km
  VBR.in.SV.T_K = HF.T_C+273 + dTdzad * HF.z_km; % set HF temperature, convert to K
  % construct pressure as a function of z, build matrix same size as T_K:
  HF.P_z=HF.rho*9.8*HF.z_km*1e3/1e9; % fix density, should account for thermal exp, compressibility
  VBR.in.SV.P_GPa = repmat(HF.P_z,1,numel(HF.t_s)); % pressure [GPa]

  sz=size(HF.T_C);
  
% calculate solidus 
h2o_ppm_s = 400 * ones(sz);
h2o_ppm_s(HF.z_km < 80, :) = 0;
wt_fraction_s = h2o_ppm_s * 1e-4;
wt_fraction_fluid = wt_fraction_s / 1e-2;
Tsol = SoLiquidus(VBR.in.SV.P_GPa*1e9, wt_fraction_fluid, 0, 'katz');

% set the other state variables as matrices of same size
sz=size(HF.T_C);
VBR.in.SV.rho = 3300 * ones(sz); % density [kg m^-3]
VBR.in.SV.sig_MPa = 10 * ones(sz); % differential stress [MPa]

VBR.in.SV.dg_um = 0.005 * 1e6 * ones(sz); % grain size [um]
VBR.in.SV.Tsolidus_K = Tsol.Tsol + 273; 

VBR.in.SV.phi = 0.01 * (VBR.in.SV.T_K > VBR.in.SV.Tsolidus_K) ; % melt fraction

% CALL THE VBR CALCULATOR %
  [VBR] = VBR_spine(VBR) ;

%% Build figures %%
  % input: T, Tsol, phi 
  figure('Position',[100 100 1200 400])
  ax1=subplot(1,3,1);
  contourf(HF.t_Myr,HF.z_km,VBR.in.SV.T_K - 273,20)
  colormap(ax1,hot)
  xlabel('Seaflor Age [Myr]')
  ylabel('Depth [km]')
  set(gca,'ydir','reverse')
  title('Temperature [C]')
  colorbar()
  
  ax2=subplot(1,3,2);
  contourf(HF.t_Myr,HF.z_km,VBR.in.SV.Tsolidus_K-273,20)
  colormap(ax2,hot)
  xlabel('Seaflor Age [Myr]')
  ylabel('Depth [km]')
  set(gca,'ydir','reverse')
  title('Solidus Temperature [C]')
  colorbar()
  
  ax3=subplot(1,3,3);
  contourf(HF.t_Myr,HF.z_km,VBR.in.SV.phi)  
  xlabel('Seaflor Age [Myr]')
  ylabel('Depth [km]')
  set(gca,'ydir','reverse')
  title('\phi')
  colorbar()



  % contour T(z,t)
  figure('Position',[100 100 1200 400])
  ax1=subplot(2,2,1);
  contourf(HF.t_Myr,HF.z_km,VBR.in.SV.T_K - 273,20)
  colormap(ax1,hot)
  xlabel('Seaflor Age [Myr]')
  ylabel('Depth [km]')
  set(gca,'ydir','reverse')
  title('Temperature [C]')
  colorbar()

  % contour shear wave velocity at different frequencies
  for i_f=1:3
     ax=subplot(2,2,i_f+1);
     % contourf(HF.t_Myr,HF.z_km,1./VBR.out.anelastic.xfit_premelt.Q(:,:,i_f), 'LineColor','none')
     contourf(HF.t_Myr,HF.z_km,VBR.out.anelastic.xfit_premelt.V(:,:,i_f), 'LineColor','none')
     colormap(ax,winter);
     xlabel('Seaflor Age [Myr]')
     ylabel('Depth [km]')
     set(gca,'ydir','reverse')
     title(['V [m/s] xfit premelt at ',num2str(VBR.in.SV.f(i_f)),' Hz'])
     colorbar()
  end
  
  figure('Position',[100 100 1200 400])
  ax1=subplot(2,2,1);
  contourf(HF.t_Myr,HF.z_km,VBR.in.SV.T_K - 273,20)
  colormap(ax1,hot)
  xlabel('Seaflor Age [Myr]')
  ylabel('Depth [km]')
  set(gca,'ydir','reverse')
  title('Temperature [C]')
  colorbar()

  % contour shear wave velocity at different frequencies
  for i_f=1:3
     ax=subplot(2,2,i_f+1);
     contourf(HF.t_Myr,HF.z_km,(1./VBR.out.anelastic.xfit_premelt.Q(:,:,i_f)), 30, 'LineColor','none')
     % contourf(HF.t_Myr,HF.z_km,VBR.out.anelastic.xfit_premelt.(:,:,i_f), 'LineColor','none')
     colormap(ax,winter);
     xlabel('Seaflor Age [Myr]')
     ylabel('Depth [km]')
     set(gca,'ydir','reverse')
     title(['1/Q xfit premelt at ',num2str(VBR.in.SV.f(i_f)),' Hz'])
     % set(gca, 'clim', [-10, -4])
     colorbar()
  end

figure('Position',[100 100 1200 400])
subplot(1,3,1)
plot(VBR.in.SV.Tsolidus_K(:,1)-273, HF.z_km(:,1),'--k', 'linewidth',2)
hold on
plot(VBR.in.SV.T_K(:,1:10:50)-273, HF.z_km(:,1))
xlabel('T [C]')
ylabel('z [km]')
set(gca,'ydir','reverse')

subplot(1,3,2)
i_f = 1;
plot(VBR.out.anelastic.xfit_premelt.V(:,1:10:50,i_f)/1e3, HF.z_km(:,1))
xlabel('V [km/s]')
set(gca,'ydir','reverse')

subplot(1,3,3)
i_f = 1;
plot(1./VBR.out.anelastic.xfit_premelt.Q(:,1:10:50,i_f), HF.z_km(:,1))
xlabel('1/Q')
set(gca,'ydir','reverse')
