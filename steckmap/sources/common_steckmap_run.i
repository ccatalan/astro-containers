#include "STECKMAP/Pierre/POP/sfit.i"

// This is a yorick file demonstrating a basic run of STECKMAP
// it features the minimal steps required to interpret the spectrum, namely:
// 1 - read the spectrum
// 2 - create a SSP sequence
// 3 - run the fitting engine sfit (STECKMAP algorithm)

// NB: HELP CAN BE OBTAINED ABOUT THE FOLLOWING COMMANDS BY TYPING FOR INSTANCE:
// help, convert_all
// or
// info, convert_all


// 1 - reading the spectrum
file="PATH TO FITS FILE";
//file=fV; // the test case
a=convert_all(file,SNR0=100.,z0=0);
// SNR0=100 per pixel and z0=0 are the default values and can be adjusted
// should show the spectrum

// 2 - SSP sequence generation
// some parameters for model SSP sequence generation
ages=[1.e7,2.e10]; // ages domain in yr
wavel=[3000.,6000.]; // wavelength domain. If too wide, the output will be restricted to what is supplied by the chosen model.
basisfile="MILES"; // other common choices: "BC03", "PHR"
nbins=40; // number of age bins
dlambda=[]; // sampling of the basis. If void, use the original sampling of the chosen model.

b=bRbasis3(ages,wavel=[3000.,6000.],basisfile=basisfile,nbins=nbins,dlambda=dlambda); // generate basis
ws;plb,b.flux,b.wave; // pltos the basis


// 3 - Running the fitting engine

kin=1; // kin=0 -> no kinematics, kin=1 -> LOSVD search
epar=3; // with continuum matching (to deal with flux calibration errors)
nde=30; // number of control points for continuum matching (as a rule of thumb I try to have one control point per 100 Angstrom, so as not to change the depth/shape of the spectral line)
vlim=[-500.,500.]; // wavelength domain for the LOSVD
// NB: if no redshift is provided, by the user, the velocity domain should be wide enough to accomodate the large radial velocity of the observed galaxy, such as [-1500.,1500.] for instance
meval=500; // maximum number of evalutations during the optimisation
// MASK definition: a mask is defined as an array of dimension [2,n]  where n is the number of regions to be masked. Each region is delimited by its blue and red edge, as shown below:
RMASK=[[3100.,3120.],[6300.,6500.],[3900.,3910.]]; 

mux=1.e2; // smoothing parameter for the stellar age distribution
muz=1.e3; // smoothing parameter for the age-Z relation
muv=1.e1; // smoothing parameter for the LOSVD
// NB: it is usually good practice to play around with the smoothing parameters a bit to get a feel for it. values up to 10^5 larger or smaller than the default can be tried out.

L1="D2"; // smoothing kernel for stellar age distribution (here D2, i.e. square Laplacian)
L3="D1"; // smoothing kernel for age-metallicity relation (here D1, i.e. square gradient)

x=sfit(a,b,kin=kin,epar=epar,noskip=1,vlim=vlim,meval=meval,nde=30,RMASK=RMASK,mux=mux,muv=muv,muz=muz,L3=L3,L1=L1);
