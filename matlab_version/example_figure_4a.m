clc
clear all

calculate_boundary_reflection = true; %name is self-explanatory
dependent_scattering = true; %true= static structure factor, false = independent scattering assumption
effective_medium = true; %consider effective refractive index of the medium while calculating the boundary reflectivities
exact_scattering_phase = true;%true = exact scattering phase function, false= henyey greenstein approximation 

%problem definition
lambda=0.5; %freespace wavelength of incident ray in unit length 
radius=50*10^-3;
thickness=5;
% fv_list=[0.01,0.02,0.05,0.1,0.2];
fv_list=logspace(log10(0.01),log10(0.2),50)
polar_angle_deg=0; %incident angles. 0 = perpendicular to slab face. 90 parallel and should be avoided.


%optical properties
n_medium = 1; %real refractive index of substrate
k_medium = 0.0; %imaginary refractive index of substrate
n_pigment=1.5; %real refractive index of particle
k_pigment=0.; %imaginary refractive index of particle
n_substrat=1; %real refractive index of substrate
k_substrat=0; %imaginary refractive index of substrate

%numerical stability settings
photon_number=10^6; %number of rays that will be traced, higher the number more accurate the result
nang=10000; %discritization of scattering angle, for large size parameters (>1000) should be high (>50000) 
ref_nh_ind=zeros(length(fv_list),1);
tra_nh_ind=zeros(length(fv_list),1);
absorptance_ind=zeros(length(fv_list),1);
ref_nh_dep=zeros(length(fv_list),1);
tra_nh_dep=zeros(length(fv_list),1);
absorptance_dep=zeros(length(fv_list),1);

parfor i=1:length(fv_list)
    f_v=fv_list(i);
    dependent_scattering=false;
    [ref_nh_ind(i),tra_nh_ind(i),absorptance_ind(i)] = rte_monodisperse(calculate_boundary_reflection,dependent_scattering,effective_medium,exact_scattering_phase,lambda,thickness,radius,f_v,polar_angle_deg,n_medium,k_medium,n_pigment,k_pigment,n_substrat,k_substrat,photon_number, nang);
    dependent_scattering=true;
    [ref_nh_dep(i),tra_nh_dep(i),absorptance_dep(i)] = rte_monodisperse(calculate_boundary_reflection,dependent_scattering,effective_medium,exact_scattering_phase,lambda,thickness,radius,f_v,polar_angle_deg,n_medium,k_medium,n_pigment,k_pigment,n_substrat,k_substrat,photon_number, nang);
end


set(0, 'DefaultLineLineWidth', 2); %set thickness of all the lines = 2

figure('Renderer', 'painters', 'Position', [500 300 428 420]) % starting point and height - width of the frame

set(gca, 'ColorOrder', [0 0 0;0 0.5 0], 'NextPlot', 'replacechildren');% color of lines in the plot with the given order. remember it is periodic

hAx=gca;
semilogx(fv_list,tra_nh_ind,fv_list,tra_nh_dep)
hAx.XColor = [0 0 0];
hAx.YColor = [0 0 0];
hAx.LineWidth = 1.5;
axis square
hLg=legend('Independent','Dependent','Location','southeast');
hLg.LineWidth=1.5;
hLg.EdgeColor = [0 0 0];
xlabel('Particle volume fraction, f_v')
ylh=ylabel('Transmittance, T_n_h');
ylh.VerticalAlignment	= 'bottom'; %if it is not alligned well, try 'top' and 'bottom' too
xlim([0.01 0.2])
ylim([0 1])
hAx.XAxis.TickValues = [0.01 0.02 0.05 0.1 0.2];
% xticklabels({'10^{-4}','10^{-2}','0.5'})
set(gca,'FontSize',13)
set(gca,'XMinorTick','on','YMinorTick','on')
box on
% saveas(gcf,'fig_4a.png')
