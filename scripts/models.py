#ODEs of different models

import numpy as np

#The simplest model of resource dependent model
def resource_dependent_model(molecules,t,alpha,beta,K,KTL,KTX,lamb1,lamb2):
    R,dT,dG,mT,mG,pT,pG,TX,TL = molecules
    return np.array(
        [
            -alpha*TX*dT*pT*R/(R+K)-alpha*TX*dG*pT*R/(R+K)-beta*TL*mT*R/(R+K)-beta*TL*mG*R/(R+K),
            0,
            0,
            alpha*TX*pT*R/(R+K),
            alpha*TX*pT*R/(R+K),
            beta*TL*(R/(R+K)),
            beta*TL*R/(R+K),
            -lamb1*TX,
            -lamb1*TL
        ]
    )
    
#resource dependent model that include saturation effect 
def resource_dependent_saturation_model(molecules,t,alpha,beta,K,KTL,KTX,lamb1,lamb2):
    R,dT,dG,mT,mG,pT,pG,TX,TL = molecules
    return np.array(
        [
            -alpha*TX*dT/(dT+dG+KTX)*pT*R/(R+K)-alpha*TX*dG/(dT+dG+KTX)*pT*R/(R+K)-beta*TL*mT/(mT+mG+KTL)*R/(R+K)-beta*TL*mG/(mT+mG+KTL)*R/(R+K),
            0,
            0,
            alpha*TX*dT/(dT+dG+KTX)*pT*R/(R+K),
            alpha*TX*dG/(dT+dG+KTX)*pT*R/(R+K),
            beta*TL*(mT/(mT+mG+KTL))*(R/(R+K)),
            beta*TL*mG/(mT+mG+KTL)*R/(R+K),
            -lamb1*TX,
            -lamb1*TL
        ]
    )
    

#define model
def Repressor_model(molecules,t,alpha,beta,K,Kr,KTL,KTX,lamb1,lamb2):
    T7_DNA,T7_RNA,T7,GFP_DNA,GFP_RNA,GFP,Repressor_DNA,Repressor_RNA,Repressor,R = molecules
    TX = 1
    TL = 1
    return np.array(
        [
            0,
            alpha*TX*T7_DNA*T7*R/(R+K)*(Kr/(Repressor+Kr)),
            beta*TL*T7_RNA*(R/(R+K)),
            0,
            alpha*TX*GFP_DNA*T7*R/(R+K),
            beta*TL*GFP_RNA*R/(R+K),
            0,
            alpha*TX*Repressor_DNA*T7*R/(R+K),
            beta*TL*Repressor_RNA*R/(R+K),
            -alpha*TX*T7_DNA*T7*R/(R+K)-alpha*TX*GFP_DNA*T7*R/(R+K)-alpha*TX*Repressor_DNA*T7*R/(R+K)-beta*TL*T7_RNA*R/(R+K)-beta*TL*GFP_RNA*R/(R+K)-beta*TL*Repressor_RNA*R/(R+K),
        ]
    )
    