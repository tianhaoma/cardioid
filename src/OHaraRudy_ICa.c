// ICaL  V,Cass:phiCaMK:PCaL :d,fFast,fSlow,fCaFast,fCaSlow,jCa,n,fCaMKFast,fCaCaMKFast
// ICaNa V,Nass:phiCaMK:PCaNa:d,fFast,fSlow,fCaFast,fCaSlow,jCa,n,fCaMKFast,fCaCaMKFast
// ICaK  V,Kss :phiCaMK:PCaK :d,fFast,fSlow,fCaFast,fCaSlow,jCa,n,fCaMKFast,fCaCaMKFast
#include <math.h>
#include <assert.h>
#include "OHaraRudy.h"
#include "OHaraRudy_ICa.h"
void OHaraRudy_ICaFunc(CELLPARMS *parmsPtr, double  *cell, int pOffset, DERIVED *derived, double dt )
{
   VOLTAGE *voltage = (VOLTAGE *)cell; 
   CONCENTRATIONS  *concentrations = (CONCENTRATIONS*) (cell + CONCENTRATIONS_OFFSET); 
   double V = voltage->Vm; 
   double Cass    =concentrations->Cass;
   double Nass    =concentrations->Nass;
   double Kss     =concentrations->Kss;

   double phiCaMK   = derived->phiCaMK; 

   PSTATE *pState = (PSTATE *)(cell+pOffset) ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   

   double d           =pState->d;
   double fFast       =pState->fFast;
   double fSlow       =pState->fSlow;
   double fCaFast     =pState->fCaFast;
   double fCaSlow     =pState->fCaSlow;
   double n           =pState->n;
   double jCa         =pState->jCa;
   double fCaMKFast   =pState->fCaMKFast;
   double fCaCaMKFast =pState->fCaCaMKFast;

   double AfSlow=1-AfFast;  
   double f = AfFast*fFast + AfSlow*fSlow; 

   double AfCaFast=0.3+0.6*sige(0.1*V-1.0); 
   double AfCaSlow=1-AfCaFast;  
   double fCa = AfCaFast*fCaFast + AfCaSlow*fCaSlow; 

   double AfCaMKFast=AfFast; 
   double AfCaMKSlow=1-AfCaMKFast;  
   double fCaMKSlow = fSlow; 
   double fCaMK = AfCaMKFast*fCaMKFast + AfCaMKSlow*fCaMKSlow; 

   double AfCaCaMKFast = AfCaFast; 
   double AfCaCaMKSlow = AfCaSlow; 
   double fCaCaMKSlow = fCaSlow; 
   double fCaCaMK = AfCaCaMKFast*fCaCaMKFast + AfCaCaMKSlow*fCaCaMKSlow; 


   double chi1 = d*(1-phiCaMK)*(f    *(1-n)+fCa    *n*jCa); 
   double chi2 = d*   phiCaMK *(fCaMK*(1-n)+fCaCaMK*n*jCa); 
   double chi  = chi1 + 1.1 * chi2; 
   double VFRT = V*FRT; 
   double VFFRT = F*VFRT; 
   double expVFRT = exp(VFRT); 
   double expCa =  SQ(expVFRT); // exp(zCa*VFRT);
   double expNa =     expVFRT;  // exp(zNa*VFRT);
   double expK  =     expVFRT;  // exp(zK *VFRT);
   double xCa = zCa*V*FRT; 
   double xNa = zNa*V*FRT; 
   double xK  = zK *V*FRT; 
   double xexpCa;  xexp(xexpCa,xCa); 
   double xexpNa;  xexp(xexpNa,xNa); 
   double xexpK ;  xexp(xexpK, xK ); 
   
   double psiCa   = zCa*F*(gammaCai*Cass*expCa-gammaCao*Cao)*xexpCa; 
   double psiCaNa = zNa*F*(gammaNai*Nass*expNa-gammaNao*Nao)*xexpNa; 
   double psiCaK  = zK *F*(gammaKi * Kss*expK -gammaKo *Ko )*xexpK ; 

   derived->I.CaL   = cP->PCaL *psiCa  *chi; 
   derived->I.CaNa  = cP->PCaNa*psiCaNa*chi; 
   derived->I.CaK   = cP->PCaK *psiCaK *chi; 

   double dMhu = sige(-(V+3.940)/4.230); 
   double dTau = 0.6 + 1/(exp(-0.05*(V+6.0)) + exp(0.09*(V+14.0))); 
   double dTauR = 1.0/dTau; 
   double dd = (dMhu-d)*dTauR;  // gate

   double fMhu = sige((V+19.58)/3.696); 
   double fFastTau = 7.0 + 1/(0.0045*exp(-(V+20)/10) + 0.0045*exp((V+20)/10)); 
   double fFastTauR = 1/fFastTau; 

   double fSlowTau  = 1000 + 1/(3.5e-5*exp(-(V+5)/4) + 3.5e-5*exp((V+5)/6)); 
   double fSlowTauR = 1/fSlowTau; 
   double dfSlow = (fMhu-fSlow)*fSlowTauR;  // gate
   double dfFast = (fMhu-fFast)*fFastTauR;  // gate

   double fCaMhu  =  fMhu; 
   double fCaFastTau = 7.0 + 1/(0.04*exp(-(V-4)/7) + 0.04*exp((V-4)/7)); 
   double fCaFastTauR = 1/fCaFastTau; 
   double fCaSlowTau = 100 + 1/(1.2e-4*exp(-(V)/3) + 1.2e-4*exp((V)/7)); 
   double fCaSlowTauR = 1/fCaSlowTau; 
   double dfCaSlow = (fCaMhu-fCaSlow)*fCaSlowTauR;  // gate
   double dfCaFast = (fCaMhu-fCaFast)*fCaFastTauR;  // gate

   double jCaMhu = fMhu; 
   double jCaTauR = 1/75.0; 
   double djCa = (jCaMhu-jCa)*jCaTauR;  // gate

   double fCaMKMhu = fMhu;
   double fCaMKFastTauR = fFastTauR/2.5; 
   double dfCaMKFast = (fCaMKMhu-fCaMKFast)*fCaMKFastTauR;  // gate

   double fCaCaMKFastMhu = fMhu; 
   double fCaCaMKFastTauR = fCaFastTauR/2.5; 
   double dfCaCaMKFast = (fCaCaMKFastMhu-fCaCaMKFast)*fCaCaMKFastTauR;  // gate

   double Kmn = 0.002; 
   double Kp2n=1000.0; 
   double Km2n=jCa*1.0; 
   double an = 1.0/(Kp2n/Km2n + SQ(SQ(1+Kmn/Cass)));
   double dn =  an*Kp2n -n * Km2n;  //gate

   ENDCODE()
   pState->d += dt*dd; 
   pState->fSlow += dt*dfSlow; 
   pState->fFast += dt*dfFast; 
   pState->fCaSlow += dt*dfCaSlow; 
   pState->fCaFast += dt*dfCaFast; 
   pState->jCa     += dt*djCa; 
   pState->fCaMKFast += dt*dfCaMKFast; 
   pState->fCaCaMKFast += dt*dfCaCaMKFast; 
   pState->n  += dt*dn; 

}