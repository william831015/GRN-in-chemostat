def resource_dependent_model(molecules,t,alpha,beta,K,KTL,KTX,lamb1,lamb2):
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
    