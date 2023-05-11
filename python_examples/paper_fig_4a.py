import cyaxares
from numpy import zeros

import matplotlib.pyplot as plt
import matplotlib

fv_list=[0.01,0.02,0.05,0.1,0.2]
radius=50*10**-3;
thickness=5;

calculate_boundary_reflection=True
# dependent_scattering 
effective_medium=True
exact_scattering_phase=True
lambda_=0.5
polar_angle_deg=0
n_medium=1.0
k_medium=0.0
n_pigment=1.5
k_pigment=0
n_substrat=1
k_substrat=0
photon_number=10**5
nang=10000


r_dep=zeros((len(fv_list)))
t_dep=zeros((len(fv_list)))
a_dep=zeros((len(fv_list)))
r_ind=zeros((len(fv_list)))
t_ind=zeros((len(fv_list)))
a_ind=zeros((len(fv_list)))
for i, f_v in enumerate(fv_list):
    r_dep[i],t_dep[i],a_dep[i] = cyaxares.rte_monodisperse(calculate_boundary_reflection=calculate_boundary_reflection,dependent_scattering=True,effective_medium=effective_medium,exact_scattering_phase=exact_scattering_phase,lambda_=lambda_,thickness=thickness,radius=radius,f_v=f_v,polar_angle_deg=polar_angle_deg,n_medium=n_medium,k_medium=k_medium,n_pigment=n_pigment,k_pigment=k_pigment,n_substrat=n_substrat,k_substrat=k_substrat,photon_number=photon_number, nang=nang);
    r_ind[i],t_ind[i],a_ind[i] = cyaxares.rte_monodisperse(calculate_boundary_reflection=calculate_boundary_reflection,dependent_scattering=False,effective_medium=effective_medium,exact_scattering_phase=exact_scattering_phase,lambda_=lambda_,thickness=thickness,radius=radius,f_v=f_v,polar_angle_deg=polar_angle_deg,n_medium=n_medium,k_medium=k_medium,n_pigment=n_pigment,k_pigment=k_pigment,n_substrat=n_substrat,k_substrat=k_substrat,photon_number=photon_number, nang=nang);


plt.rc('lines', linewidth=2.5)

fig, ax = plt.subplots()
fig.set_size_inches(5.5, 5)

# change all spines
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)

# increase tick width
ax.tick_params(width=2)


line1, = ax.semilogx(fv_list, t_ind, color=[0,0,0], label='Independent')
line2, = ax.semilogx(fv_list, t_dep, color=[0,0.5,0], label='Dependent')

ax.set(xlabel='Particle volume fraction', ylabel='Transmittance',
       xlim =(0.01, 0.2), ylim =(0, 1),)
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
ax.legend(handlelength=4)
matplotlib.rcParams.update({'font.size': 14})
plt.legend(loc='center right')
legend = ax.legend()
legend.get_frame().set_linewidth(2)
legend.get_frame().set_edgecolor("black")
plt.xticks(fv_list, fv_list)
fig.savefig('figure_4a.png', dpi=300)
