function [ref_nh,tra_nh,absorptance] = fn_polydisperse()

%solution of radiative transfer equation in one layer pigmented plane parallel medium
%ray is incident from air to coating. coating is coated on a substrate
%substrate could be air or other material such as silver, glass etc.
%the code estimates spectral hemispherical reflectance, transmittance and absorptance
%can handle; independent and dependent scattering, boundary reflections, absorption in
%medium. can't handle; coherent backscattering, plascmon coupling, near field effects and polarized ray tracing.
%while calculating the scattering direction the code uses  henyey greenstein 
%function approximation

%phenomena settings
calculate_boundary_reflection = true; %name is self-explanatory
dependent_scattering = true; %true= static structure factor, false = independent scattering assumption
effective_medium = true; %consider effective refractive index of the medium while calculating the boundary reflectivities, this option is useless if calculate_boundary_reflection = false
exact_scattering_phase = false;%true = exact scattering phase function, false= henyey greenstein approximation 

%problem definition
lambda=0.5; %freespace wavelength of incident ray in unit length 
thickness=100;  %thickness of coating in unit length  
avg_radius=0.2; %radius of particle in unit length  
sigma = 0.0001;
f_v=0.10; %volume fraction. 0.01 corresponds to 1% 
polar_angle_deg=0; %incident angles. 0 = perpendicular to slab face. 90 is parallel and should be avoided.


%optical properties
n_medium = 1.33; %real refractive index of substrate
k_medium = 0; %imaginary refractive index of substrate
n_pigment=1.5; %real refractive index of particle
k_pigment=0; %imaginary refractive index of particle
n_substrat=1; %real refractive index of substrate
k_substrat=0; %imaginary refractive index of substrate

%numerical stability settings
photon_number=10^6; %number of rays that will be traced, higher the number more accurate the result
nang=20000; %discritization of scattering angle, for large size parameters (>1000) should be high (>50000) 

%%%%%%%%%%%%%%%%%%%%%%% end of inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rad_list=linspace(avg_radius-sigma*3,avg_radius+sigma*3,20);
rad_list(rad_list<=0)=[];
weight_vector = normpdf(rad_list,avg_radius,sigma); %the distribution.
normalize_weight_vector=trapz(rad_list,weight_vector);
weight_vector=weight_vector/normalize_weight_vector;

% figure
% plot(rad_list,weight_vector)

teta=linspace(eps,pi,nang)';%don't start with zero to avoid division by zero
u=cos(teta);
polar_angle_rad=polar_angle_deg*pi/180;

P11_r=zeros(nang,length(rad_list));
Csca_r=zeros(1,length(rad_list));
Cabs_r=zeros(1,length(rad_list));
g_r=zeros(1,length(rad_list));
S11 = [];
if calculate_boundary_reflection
    if effective_medium 
        [n_eff_m,k_eff_m] = maxwell_garnett(f_v,n_pigment,k_pigment,n_medium,k_medium);
    else
        n_eff_m = n_medium;
        k_eff_m = k_medium;
    end
    % calculate surface reflection from air to medium. 
    % medium to air and medium to substrate is calculated within snell.m.
    % air to medium is calculated here seperately since we need refraction angle
    teta_prime=F_fresnel_2(n_eff_m,k_eff_m,polar_angle_rad);
    cos_teta=cos(polar_angle_rad);
    sin_teta=sqrt(1-cos_teta*cos_teta);
    carpan2=1/(n_eff_m-1i*k_eff_m);
    sin_x2=sin_teta*carpan2;
    cos_x2=sqrt(1-sin_x2*sin_x2);
    carpan1=cos_teta/cos_x2;
    carpan3=cos_x2/cos_teta;
    E_parallel=(carpan1-carpan2)/(carpan1+carpan2);
    R_parallel=E_parallel*conj(E_parallel);
    E_orth=(carpan3-carpan2)/(carpan3+carpan2);
    R_orth=E_orth*conj(E_orth);
    reflectance=real(R_parallel+R_orth)*0.5;
    sur_reflection=reflectance;
else
    teta_prime=polar_angle_rad;
    sur_reflection=0;
    n_eff_m=1;
    k_eff_m=0;
    n_substrat=1;
    k_substrat=0;
end

for i=1:length(rad_list)
    radius=rad_list(i);
    V=4*pi*radius^3/3; %volume of a single sphere
    Area=pi*radius^2;
    x=2*pi*(n_medium + k_medium*1i)*radius/lambda; %size parameter
    x_reel=2*pi*(n_medium)*radius/lambda; %size parameter
    x0=2*pi*radius/lambda; %size parameter at free-space
    if abs(x)>0.2 && effective_medium
        warning('Effective medium theory is enabled but size parameter is larger than 0.2');
    end
    
    [Qext_ind,w0_ind,g_ind,S1,S2,P11_ind,QI,S_FACTOR,SCALE_FACTOR] = abs_mie_yang(radius,lambda,n_medium,k_medium,n_pigment,k_pigment,nang);
    % g_ind
    S11=S_FACTOR^2*0.5*(abs(S1).^2+abs(S2).^2)'/QI;
    Qsca_ind=Qext_ind*w0_ind;
    Qabs=Qext_ind*(1-w0_ind);
    if dependent_scattering
        S=SSF_correction(f_v, teta, lambda/abs(n_medium+k_medium*1i), radius);
        Qsca=2*trapz(teta,sin(teta) .* S11 .* S)/abs(x)^2;
        g=trapz(teta,sin(teta) .*cos(teta) .* S11 .* S)/trapz(teta,sin(teta) .* S11 .* S);
        if exact_scattering_phase
            P11_r(:,i) = (S11.*S)/trapz(teta,S11.*sin(teta).*S);
        end
    else
        Qsca=Qsca_ind;%scattering efficiency
        g=g_ind;
        if exact_scattering_phase
            P11_r(:,i) = S11/trapz(teta,sin(teta).*S11);
        end
    end
    
    
    Csca_r(i)=Area*Qsca;%scattering crossection[m^2]
    Cabs_r(i)=Area*Qabs;%absorption crossection[m^2]
    g_r(i)=g;%absorption crossection[m^2]
end

Csca = trapz(rad_list,Csca_r.*weight_vector);
Cabs = trapz(rad_list,Cabs_r.*weight_vector);
P11 = trapz(rad_list,(Csca_r.*P11_r.*weight_vector)')./trapz(rad_list,(Csca_r.*weight_vector)');
g = trapz(rad_list,g_r.*Csca_r.*weight_vector)/trapz(rad_list,Csca_r.*weight_vector);

if exact_scattering_phase
    flipped_u=flip(u);
    zero_to_one=linspace(0,1,nang)';
    cdf=zeros(nang,1);
    for i=2:nang
        cdf(i)=trapz(flipped_u(1:i),P11(1:i));%calculate cumulative distribution function
    end  
    inv_cdf_cos=cos(interp1(cdf,teta,zero_to_one,'linear','extrap'));%calculate inverse cumulative distribution function

    %%% when inv_cdf_cos(end) is not -1 it means nang is not high enough
    if ~isempty(S11)
        if inv_cdf_cos(end)>-0.99 
            warning('consider to increase nang value')
        end
    end
    %%% end of warning

else
    inv_cdf_cos=0;
end

r_3=trapz(rad_list,(weight_vector.*rad_list.^3)');
V_avg=(4/3)*pi*r_3;

alfa=f_v*Csca/V_avg;%scattering coefficient[1/m]
beta=f_v*Cabs/V_avg+(1-f_v)*4*pi*k_medium/lambda;%absorption coefficient[1/m] absorption of medium is implemented as a bulk property
mu_tot=alfa+beta;%extinction coefficient[1/m]
scat_prob=alfa/mu_tot;%scattering albedo also scattering probability

% tic
[ref_nh,tra_nh,absorptance]=monte_carlo(photon_number,sur_reflection,cos(teta_prime),thickness,scat_prob,mu_tot,n_eff_m,k_eff_m,n_substrat,k_substrat,inv_cdf_cos,g,exact_scattering_phase);
% toc