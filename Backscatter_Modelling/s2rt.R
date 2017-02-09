# Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
# These MATLAB-based computer codes are made available to the remote
# sensing community with no restrictions. Users may download them and
# use them as they see fit. The codes are intended as educational tools
# with limited ranges of applicability, so no guarantees are attached to
# any of the codes. 
#
  #Code 11.1: S2RT/R Backscattering from Rayleigh Layer with Diffuse Upper
#Boundary

#Description: Code computes sigma_0_vv and sigma_0_hh for a weakly
#scattering Rayleigh layer with albedo a < 0.2. The layer is above a ground
#surface characterized by the PRISM model (code 10.5).

#Input Variables:
  #eps: complex dielectric constant of ground surface
#f: frequency (GHz)
#s: rms height of ground surface (m)
#a: Single-scattering albedo of Rayleigh layer (unitless)
#kappa_e: extinction coefficient of the Rayleigh layer (Np/m)
#d: thickness of the Rayleigh layer (m)
#theta: Incidence angle (degrees)

#Output Variables: 
  #sigma_0_vv (dB)
#sigma_0_hh (dB)

#Book reference: Section 11-2.2 and eq. 11.23 with n = 2

#Matlab code:
S2RTR_DiffuseUB <- function(eps,f,s,a,kappa_e,d,theta)  {

  theta_r <- theta*pi/180 #transform to radian

  kappa_s <- a * kappa_e  #scattering coefficient

  #-- call the PRISM-1 surface scattering model
  tmp <- PRISM1_ForwardModel(eps,theta,s,f)
  sig_s_vv <- tmp[1]
  sig_s_hh <- tmp[2]
  sig_s_hv <- tmp[3]
  sig_s_vv  <- 10^(sig_s_vv/10)
  sig_s_hh  <- 10^(sig_s_hh/10)

  #-- calculate transmissivity in layer
  tau <- kappa_e *d * sec(theta_r)
  T <- exp(-tau)

  #-- calculate reflectivity of lower interface
  tmp <- ReflTransm_PlanaryBoundary(1, eps, theta)
  t1 <- tmp[1]
  t2 <- tmp[2]
  gammah <- tmp[3]
  gammav <- tmp[4]
  t3 <- tmp[5]
  t4 <- tmp[6]
  t5 <- tmp[7]
  t6 <- tmp[8]

  #-- calculate the total backscattering coefficient according to eq 11.23
  sigma_0_vv <- T^2 * sig_s_vv + 0.75*a * cos(theta_r) *(1- T^2)*(1 + gammav^2 * T^2) + 3 * 2 * kappa_s * d * gammav * T^2

  sigma_0_hh = T^2 *sig_s_hh + 0.75*a * cos(theta_r) *(1- T.^2)*(1 + gammah^2 * T^2) + 3 * 2 * kappa_s * d * gammah * T^2


  sigma_0_vv = 10* log10(sigma_0_vv);
  sigma_0_hh = 10* log10(sigma_0_hh);


}


#Code 10.5: PRISM-1  - Forward Model

#Description: Code computes sigma_0 for all three polarization
#combinations, given the surface parameters.

#Input Variables:
  #eps = eps' - j eps'': Complex dielectric constant of the scattering
#medium
#theta: Incidence angle (degrees)
#s: rms height (m)
#f: Frequency (GHz)

#Output Products:
#sigma_0_vv (dB)
#sigma_0_hh (dB)
#sigma_0_hv (dB)

#Book Reference: Section 10-5

#Matlab Code: 

function [sig_0_vv sig_0_hh sig_0_hv] = PRISM1_ForwardModel(eps,theta,s,f)

ks = s .* (2*pi *f/0.3); # calculate roughness parameter

theta = theta .* (pi/180.0); # incidence angle in radian

gamma0 = Fresn_Refl0(eps); #reflectivity (normal incidence)

[gammav, gammah] = Fresn_Refl(eps, theta); 

p = ( 1- (2 .*theta/pi).^(1./(3*gamma0)) .* exp(-ks) ).^2; 
q = 0.23 .* sqrt(gamma0) .*( 1 - exp(-ks));

g = 0.70 .*(1 - exp(-0.65 .* ks.^1.8));

sigvv = g .* (cos(theta)).^3 ./sqrt(p) .* (gammav + gammah);

sig_0_vv = 10*log10(sigvv);
sig_0_hh = 10*log10(sigvv .* p);
sig_0_hv = 10*log10(sigvv .* q);

end

function gamma0 = Fresn_Refl0(eps)
# calculates Fresnel reflectivity at normal incidence.
gamma0 = (abs((1 - sqrt(eps)) ./ ( 1+ sqrt(eps)))).^2;
end

#-----------------------------------------------------------------
function [gammav, gammah] = Fresn_Refl(eps,theta)
# calculates Fresnel reflectivities of v and h-polarizations at given set
# of incidence angles.

[rho_v, rho_h] = refl_coef(theta, 1, eps);
gammav = (abs(rho_v)).^2;
gammah = (abs(rho_h)).^2;

end

#----------------------------------------------------------------
function [rho_v, rho_h] = refl_coef(the1, eps1, eps2)

# calculates the v and h-polarized reflection coefficients of a plane
# dielectric surface

n1 = sqrt(eps1);
n2 = sqrt(eps2);
costh2 = sqrt(1 - (n1.*sin(the1)./n2).^2);


rho_v = -(n2.*cos(the1) - n1.*costh2)./(n2.*cos(the1)+n1.*costh2);
rho_h = (n1.*cos(the1) - n2.*costh2)./(n1.*cos(the1)+n2.*costh2);
end


#Code 2.3: Oblique Reflection and Transmission @ Planar Boundry
#Description: Code computes the reflection coefficients, transmission
#coefficients, reflectivities and transmissivities for incidence in
#medium (medium 1) upon the planar boundary of a lossless or lossy
#medium (medium 2) at any incidence angle, for both h and v polarizations

#Input Variables:
#eps1: eps1r -j*eps1i: relative dielectric constant of medium 1
#eps2 = eps2r-j*eps2i: relative dielectric constant of medium 2
#theta1d: incidence angle in medium 1 in degrees

#Output Products:
#rhoh: reflection coefficient for h pol
#rhov: reflection coefficient for v pol
#gammah:reflectivity for h pol
#gammav: reflectivity for v pol
#tauh: transmission coefficient for h pol
#tauv: transmission coefficient for v pol
#Th: transmissivity for h pol
#Tv:transmissivity for v pol
#Book Reference: Sections 2-7 & 2-8

#Example call: [rhoh rhov gammah gammav tauh tauv Th Tv] = ReflTransm_PlanarBoundary(eps1, eps2, theta1d)
#Computes 

#MATLAB Code

function [rhoh rhov gammah gammav tauh tauv Th Tv] = ReflTransm_PlanarBoundary(eps1, eps2, theta1d)

theta1 = degtorad(theta1d);

sin_theta2 = sqrt(eps1)/sqrt(eps2).*sin(theta1);
cos_theta2 = sqrt(1 - sin_theta2.^2);

rhoh = (sqrt(eps1).*cos(theta1)-sqrt(eps2).*cos_theta2) ./ (sqrt(eps1).*cos(theta1) + sqrt(eps2).*cos_theta2);
rhov = (sqrt(eps1).*cos_theta2-sqrt(eps2).*cos(theta1)) ./ (sqrt(eps1).*cos_theta2 + sqrt(eps2).*cos(theta1));

tauh = 1 + rhoh;
tauv = (1 + rhov).*(cos(theta1)./cos_theta2);


gammah = abs(rhoh).^2;
gammav = abs(rhov).^2;

Th = 1-gammah;
Tv = 1-gammav;

end
