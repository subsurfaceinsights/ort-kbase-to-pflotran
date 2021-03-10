"""
Some useful tools for working with KBase outputs.

Includes:
    * Ion size dictionary
    * Rename KBase compounds to PFLOTRAN conventions
    * Sort PFLOTRAN outputs into useful groups for plotting
"""

import copy


def ion_size():
    """
    Library of ion size parameters.

    Mapped to formulas used in ModelSEED database

    From lookup table in excel spreadsheet released by Brian M. Tissue
    http://www.tissuegroup.chem.vt.edu/a-text/index.html

    This is only a partial list for now, can populate more as needed
    """
    sizes = {
            "Ag": 2.5,
            "Al": 9,
            "Ba": 5,
            "B": 8,
            "Br": 3,
            "BrO3": 3.5,
            "CHO3": 4.5,  # ModelSEED treats this as a synonym for CO3
            "Acetate": 4.5,
            "Ca": 6,
            "Cd": 5,
            "Cl": 3,
            "ClO2": 4,
            "ClO3": 3.5,
            "ClO4": 3.5,
            "CN": 3,
            "Co": 6,
            "Cr": 9,
            "CrO4": 4,
            "Cs": 2.5,
            "Cu": 6,
            "F": 3.5,
            "Fe2": 6,
            "Fe3": 9,
            "H": 9,
            "H3O+": 9,
            "I": 3,
            "IO3": 4,
            "IO4": 3.5,
            "K": 3,
            "Mg": 8,
            "Mn": 6,
            "MnO4": 3.5,
            "MoO4": 4.5,
            "NH4": 2.5,
            "NO2": 3,
            "NO3": 3,
            "Na": 4,
            "Nd": 9,
            "Ni": 6,
            "OH": 3.5,
            "PO4": 4,
            "HPO4": 4,
            "H2PO4": 4,
            "S": 5,
            "HS": 3.5,
            "SO3": 4.5,
            "HSO3": 4,
            "SO4": 4,
            "O4S2": 5,
            "Zn": 6,
            "H2O": 3
            }
    return sizes


def rename(compounds, f_index, name_index, charge_index, rmdisplay=1):
    """
    Rename modelseed compounds to match (approximately) PFLOTRAN conventions.

    Currently this mostly means adding charges (in symbol form) to
    the end of the name. The function also fixes some quirky compound names
    (e.g. H4N = NH4) and substitutes a generic formula in for biomass (C5H7O2N)
    If other compounds have 'null' formula this will need to be refined.

    **compounds**: ordered array from KBase FBA Objects (from the KBase API)5

    **f_index**: index of formula in input array

    **name_index**: index of ModelSEED name in input array

    **charge_index**: index of the compound charge in the input array

    **rmdisplay**: 0 to keep items at name_index, 1 to remove (*default 1*).
    For long formulas, the actual formula and the display name switch columns.
    We don't usually need the actual formula, but if we do, use this to keep it
    """
    newcpds = copy.deepcopy(compounds)

    for row in newcpds:
        row[charge_index] = int(row[charge_index])
    for row in newcpds:
        if row[f_index] == "H4N":
            row[f_index] = "NH4"
        elif row[f_index] == "O4S":
            row[f_index] = "SO4"
        elif row[f_index] == "HO4P":
            row[f_index] = "HPO4"
        elif row[f_index] == "CHO3":
            row[f_index] = "HCO3"
        elif row[f_index] == "CH4O":
            row[f_index] = "Methanol"
    # For short formulas, keep those and add charges.
    # For longer, use display names
    for row in newcpds:
        if len(row[f_index]) > 4:
            row[f_index] = (row[name_index].replace('-', '').replace('+', '')
                            .replace('/', '').replace('*', '').replace('(', '')
                            .replace(')', '').replace('[', '').replace(']', '')
                            .replace('{', '').replace('}', '').replace(',', '')
                            .replace('\'', ''))
        else:
            row[f_index] = row[f_index]

    # Append charges per pflotran typical style.
    pos_charges = [1, 2, 3]
    neg_charges = [-1, -2, -3]

    for new_cpd in newcpds:
        if new_cpd[charge_index] in pos_charges:
            new_cpd[f_index] += new_cpd[charge_index]*"+"
        elif new_cpd[charge_index] in neg_charges:
            new_cpd[f_index] += abs(new_cpd[charge_index])*"-"
        elif new_cpd[charge_index] == 0:
            new_cpd[f_index] += "(aq)"
        else:
            # Some ModelSEED compounds have crazy charges -> make strings
            new_cpd[f_index] += str(new_cpd[charge_index])
    # Biomass ends up with (aq) appended, but it doesn't go in the .dat or in
    # the main equation so we don't really care (for now)

    if rmdisplay == 1:
        for row in newcpds:
            row.pop(name_index)

    return newcpds


# Sorting function to make plotting cleaner and easier. Needs some work.
def speciessorter(df, timepoints=False):
    """
    Divides PFLOTRAN outputs in to useful groups.

    This function splits outputs from PFLOTRAN, with the .tec file read in as
    a dataframe (**df**) into logical groupings. This is useful for plotting as
    well as data management.

    **Inputs:**

    **df:** Dataframe containing the relevant columns

    **Timepoints:** Identify if the files being sorted are for distinct time
    points (no time column for individual timepoint files)

    **Outputs:**

    **Nspecies:** (Aqueous) nitrogen species involved in nitrogen cycling

    **Cspecies:** (Aqueous) carbon sources

    **Pspecies:** (Aqueous) phosphate-based phosphorous sources

    **Sspecies:** (Aqueous) sulfate-based sulfur species (just SO4 at the
    moment, but wanted to maintain flexibility)

    **ions:** aqueous ions (mostly metals)

    **solids:** Solid-phase compounds (e.g. CRF fertilizer)

    **physical:** Physical parameters included in the model (not grouped for
    plotting, just for reference)

    **biomass:** Biomass species
    """
    # Filter out fields that aren't useful to sort, depending on file type
    if not timepoints:
        newdf = df.drop(columns=['Time_[d]', 'pH', 'H+_[M]', 'O2(aq)_[M]'])
    else:
        newdf = df.drop(columns=['pH', 'H+_[M]', 'O2(aq)_[M]', 'X_[m]',
                                 'Y_[m]', 'Z_[m]'])
    # Now split out key categories as cleanly as possible. May need work.
    physical = [col for col in newdf.columns if ('Pressure' in col or
                                                 'Saturation' in col or
                                                 'Density' in col or
                                                 'X_g^l' in col or
                                                 'X_l^l' in col or
                                                 'X_g^g' in col or
                                                 'X_l^g' in col or
                                                 'Energy' in col or
                                                 'Temperature' in col or
                                                 'State' in col or
                                                 'Material' in col)]
    solids = [col for col in newdf.columns if '(s)' in col]
    biomass = [col for col in newdf.columns if ('(aq)' not in col and
                                                '-' not in col and
                                                '+' not in col and
                                                '(s)' not in col and
                                                '[M]' not in col and
                                                col not in physical)]
    nitrogen_species = [col for col in newdf.columns if ('N2' in col
                                                         or 'NO2-' in col
                                                         or 'NO3-' in col or
                                                         'NH4+' in col and
                                                         col not in solids)]
    phosphorus_species = [col for col in newdf.columns if ('PO4' in col and
                                                           col not in solids)]
    ions = [col for col in newdf.columns if ('H+' not in col and
                                             'ate' not in col and
                                             'ine' not in col and
                                             'ol' not in col and
                                             'O' not in col and
                                             '(aq)' not in col and
                                             '(s)' not in col and
                                             col not in nitrogen_species and
                                             col not in biomass and
                                             col not in physical and
                                             ('+' in col or '-' in col))]
    sulfur_species = [col for col in newdf.columns if ('SO4' in col and
                                                       col not in solids)]
    carbon_species = [col for col in newdf.columns if (col not in phosphorus_species and
                                                       col not in sulfur_species and
                                                       col not in solids and
                                                       col not in ions and
                                                       col not in nitrogen_species and
                                                       col not in biomass and
                                                       col not in physical)]

    return (nitrogen_species, carbon_species, phosphorus_species,
            sulfur_species, ions, solids, physical, biomass)
