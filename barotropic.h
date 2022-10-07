! A macro for barotropic EOS

!
!  P = Cs**2 * rho + Kappa * rho ** Gamma
!
  ! ------------
  ! macro for cs
  ! ------------
#define GETCS(CS_,RHO_) \
  CS_ = sqrt(Cs**2 + Gamma * Kappa * (RHO_) ** (Gamma -1.d0) )
  ! -------------
  ! macro for cs2
  ! -------------
#define GETCS2(CS2_,RHO_) \
  CS2_ = Cs**2 + Gamma * Kappa * (RHO_) ** (Gamma -1.d0)
  ! ------------------------------------------------------------
  ! macro for csbar
  ! ------------------------------------------------------------
#define GETCSBAR(CSBAR_,RHOL_,RHOR_) \
  rhomax=max((RHOL_),(RHOR_)) ;\
  rhomin=min((RHOL_),(RHOR_)) ;\
  rho_bar=sqrt((RHOL_)*(RHOR_)) ;\
  delrho=max(5.0d-1*dlog(rhomax/rhomin),1.0d-5) ;\
  CSBAR_ = sqrt(Kappa*(rho_bar**(Gamma-1.0d0))*sinh(Gamma*delrho)/sinh(delrho) + Cs**2 )
  ! --------------------
  ! macro for Pressure
  ! --------------------
#define GETP(P_,RHO_) \
  P_ = Cs**2 * (RHO_) + Kappa*(RHO_)**Gamma

  ! --------------------
  ! macro for Temperature
  ! --------------------
#define GETTEMP(TEMP_,RHO_) \
  TEMP_ = ModelParam_temp *( 1.d0 + ((RHO_)/Rhocr) ** (Gamma -1.d0))

  
