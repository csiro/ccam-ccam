! --- physical
      real, parameter :: rearth=6.37122e6
      real, parameter :: r=287.,g=9.80616,cp=1004.64,cpv=1869.46
      real, parameter :: roncp=r/cp
! --- chemical
      real, parameter :: fAIR_MolM = 28.965,
     &           fH_MolM   =  1.00794,
     &           fC_MolM   = 12.011,
     &           fO_MolM   = 15.9994,
     &           fS_MolM   = 32.066,

     &           fH2_MolM  = 2*fH_MolM,			! molecular hydrogen
     &           fO2_MolM  = 2*fO_MolM,			! molecular oxygen
     &           fCO2_MolM =   fC_MolM   +   fO2_MolM,	! carbon dioxide
     &           fSO2_MolM =   fS_MolM   +   fO2_MolM,	! sulphur dioxide
     &           fSO4_MolM =   fSO2_MolM +   fO2_MolM,	! sulphate
     &           fCH4_MolM =   fC_MolM   + 2*fH2_MolM,	! methane
     &           fH2O_MolM =   fO_MolM   +   fH2_MolM	! water
