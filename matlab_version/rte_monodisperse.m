function [ref_nh,tra_nh,absorptance] = rte_monodisperse(calculate_boundary_reflection,dependent_scattering,effective_medium,exact_scattering_phase,lambda,thickness,radius,f_v,polar_angle_deg,n_medium,k_medium,n_pigment,k_pigment,n_substrat,k_substrat,photon_number, nang)

    
    %%%%%%%%%%%%%%%%%%%%%%% end of inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    teta=linspace(eps,pi,nang)';%don't start with zero to avoid division by zero
    u=cos(teta);
    polar_angle_rad=polar_angle_deg*pi/180;
    
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
        teta_prime=fn_fresnel(n_eff_m,k_eff_m,polar_angle_rad);
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
    
    V=4*pi*radius^3/3; %volume of a single sphere
    Area=pi*radius^2;
    x=2*pi*(n_medium + k_medium*1i)*radius/lambda; %size parameter
    x_reel=2*pi*(n_medium)*radius/lambda; %size parameter
    x0=2*pi*radius/lambda; %size parameter at free-space
    
    if abs(x)>0.2 && effective_medium
        warning('Effective medium theory is enabled but size parameter is larger than 0.2');
    end
    
    [Qext_ind,w0_ind,g_ind,S1,S2,P11_ind,QI,S_FACTOR,SCALE_FACTOR] = abs_mie(radius,lambda,n_medium,k_medium,n_pigment,k_pigment,nang);
    
    S11=S_FACTOR^2*0.5*(abs(S1).^2+abs(S2).^2)'/QI;
    Qsca_ind=Qext_ind*w0_ind;
    Qabs=Qext_ind*(1-w0_ind);
    if dependent_scattering
        S=SSF_correction(f_v, teta, lambda/abs(n_medium+k_medium*1i), radius);
        Qsca=2*trapz(teta,sin(teta) .* S11 .* S)/abs(x)^2;
        g=trapz(teta,sin(teta) .*cos(teta) .* S11 .* S)/trapz(teta,sin(teta) .* S11 .* S);
        if exact_scattering_phase
            P11 = (S11.*S)/trapz(teta,S11.*sin(teta).*S);
        end
    else
        Qsca=Qsca_ind;%scattering efficiency
        g=g_ind;
        if exact_scattering_phase
            P11 = S11/trapz(teta,sin(teta).*S11);
        end
    end
    
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
    
    Csca=Area*Qsca;%scattering crossection[m^2]
    Cabs=Area*Qabs;%absorption crossection[m^2]
    alfa=f_v*Csca/V;%scattering coefficient[1/m]
    beta=f_v*Cabs/V+(1-f_v)*4*pi*k_medium/lambda;%absorption coefficient[1/m] absorption of medium is implemented as a bulk property
    mu_tot=alfa+beta;%extinction coefficient[1/m]
    scat_prob=alfa/mu_tot;%scattering albedo also scattering probability
    
    [ref_nh,tra_nh,absorptance]=monte_carlo(photon_number,sur_reflection,cos(teta_prime),thickness,scat_prob,mu_tot,n_eff_m,k_eff_m,n_substrat,k_substrat,inv_cdf_cos,g,exact_scattering_phase)

end