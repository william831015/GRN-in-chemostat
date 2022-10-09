def dilution(molecule, DR=0.2):

    molecule_diluted = molecule*(1-DR)

    return molecule_diluted

def replenish(molecule_replenished, dilution_ratio, replenish_conc) : 

    molecule_replenished = replenish_conc * dilution_ratio

    return molecule_replenished

def dilute_out_species(molecules,DR=0.2):
    for i in (molecules):
        i = dilution(i,DR)
