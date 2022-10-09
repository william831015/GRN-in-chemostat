def dilution(molecule, DR=0.2):

    molecule_diluted = molecule*(1-DR)

    return molecule_diluted

def dilute_out_species(molecules):
    for i in (molecules):
        print(i)

def replenish(molecule_replenished, dilution_ratio, replenish_conc) : 

    molecule_replenished = replenish_conc * dilution_ratio

    return molecule_replenished