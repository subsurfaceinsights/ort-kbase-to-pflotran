# -*- coding: utf-8 -*-##
"""
Script to download KBase FBA results and translate to PFLOTRAN stoichiometry.

In order for this script to work, objects in the KBase narrative must have the
following name structure:
models: <modelname>.model
fbas: <modelname>.<details>.fba

`modelname` should be a unique, human-readable identifier for the model so that
    outputs are easy to understand.
`details` are optional and can include media type/ constraints/ etc
"""

# Imports
import sys
import json
import os
import csv
import yaml
import numpy as np
import periodictable as ptab
from lib.installed_clients.WorkspaceClient import Workspace
from kbpf_tools import ion_size, rename


# Config file includes parameters for KBase API.
conf = "./config.yml"
with open(conf, 'r') as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)

# Details we need to find and access the narratives
token = cfg["KB_AUTH_TOKEN"]
kbase_endpoint = cfg["endpoint"]
ws_url = cfg["workspace-url"]
ws = Workspace(ws_url, token=token)
narrative = cfg["NAR_ID"]
nar_name = cfg["NAR_NAME"]

mediadir = f"./processed-outputs/{nar_name}/"

if not os.path.exists(mediadir):
    os.makedirs(mediadir)

pulldir = f"./pulled/{nar_name}/"

if not os.path.exists(pulldir):
    os.makedirs(pulldir)

# Import all Model Objects
Model_names = ws.list_objects(params={"workspaces": [narrative],
                                      "type": "KBaseFBA.FBAModel"})
Models = []
for name in Model_names:
    Models.append(name[1])

FBA_names = ws.list_objects(params={"workspaces": [narrative],
                                    "type": "KBaseFBA.FBA"})
FBAs = []
for name in FBA_names:
    FBAs.append(name[1])

# Print reaction strings to in a file.
all_reactions = []

# Import ion_size database before the loop, otherwise it throws an error
# Note this is manually maintained unless we can find a good database to poll
ion_size = ion_size()

for fba in FBAs:
    # Select the model based on fba name, and also use this to name the
    # resulting dat file. Filter out things classed as "fba" that aren't what
    # we're looking for (e.g. results from Predict Pathways app)
    if fba.split('.')[-1] == 'fba':
        ref_name = fba.split('.')[0]
        # Can't use startswith here, in case model names begin with the same
        # sequence (e.g. thaum_sbmlfixed, thaum_sbmlfixed_raw, etc)
        model_name = [model for model in Models if
                      model.split('.')[0] == ref_name]
        model_ref = f"{narrative}/{model_name[0]}"
        if len(model_name) > 1:
            print(str(len(model_name))+(" models with the same name, please "
                  "check your naming conventions."))
            print(("Reminder: In order for this script to work, objects in "
                   "the KBase narrative must have the following name "
                   "structure: \nmodels: <modelname>.model \n"
                   "fbas: <modelname>.<details>.fba \nmodelname should be a "
                   "unique, human-readable identifier for the model so that "
                   "outputs are easy to understand. If models and fba are "
                   "one-to-one, include details in modelname.\nDetails are"
                   " optional and can include media type, constraints, etc. if"
                   " the same model is to be grown on more than one media "
                   "condition (e.g. if draft models are generated using the "
                   "same genome but different gapfilling media."))
            sys.exit()

        # Use complex names, don't overwrite dat files
        if len(fba.split('.')) > 2:
            dot = '.'
            full_name = dot.join(fba.split('.')[0:-1])
        else:
            full_name = ref_name

        datfile = f"{pulldir}{full_name}.dat"

        # Call the models + fbas
        # The initial pull gets the target objects (dictionaries) wihin a list
        # within another dictionary
        model_data = ws.get_objects2({'objects':
                                      [{'ref': model_ref}]})['data'][0]['data']
        fba_ref = narrative+"/"+fba
        fba_data = ws.get_objects2({'objects':
                                    [{'ref': fba_ref}]})['data'][0]['data']

        # The keys we actually care about are the compounds. Exchange fluxes
        # come from the first, chemical formulas come from the second.
        exchanges = fba_data['FBACompoundVariables']
        compounds = model_data['modelcompounds']

        # Exchanges dictionary doesn't include formulas or names, so using
        # cpdid as the key. Splitting removes the full path from the id, as
        # those don't match the ones from the compounds array.
        sdict = {}
        exch = 0
        while exch < len(exchanges):
            sdict.update({exchanges[exch]['modelcompound_ref']
                          .split("/")[-1]:
                              [exchanges[exch]['modelcompound_ref']
                               .split("/")[-1], exchanges[exch]['value']]})
            exch += 1

        # Make a dictionary we can use to lookup formulas and charge by ID
        fdict = {}
        comp = 0
        while comp < len(compounds):
            # Some compounds have 'null' formula. To keep these unique and
            # identifiable and help with troubleshooting later, we should have
            # it populate those withsomething else
            if compounds[comp]['formula'] == 'null':
                fdict.update({compounds[comp]['id']:
                              [compounds[comp]['name'].split('_')[0],
                              compounds[comp]['name'],
                              compounds[comp]['charge']]})
            else:
                fdict.update({compounds[comp]['id']:
                              [compounds[comp]['name'].split('_')[0],
                              compounds[comp]['formula'],
                               compounds[comp]['charge']]})
            comp += 1

        # Write exchanges to a list
        Elist = list(sdict.values())
        for row in Elist:
            row.append(fdict[row[0]][0])
            row.append(fdict[row[0]][1])
            row.append(fdict[row[0]][2])

        # Now make an array from that list, since that's what the rest of the
        # code expects. This is kind of kludgy, but necessary to avoid errors
        Ea = np.array(Elist)
        Ea = list(map(tuple, Ea))

        # Previous versions used type <i8 for charge, but beginning Aug 2019
        # this throws "ValueError: invalid literal for int() with base 10:..."
        # Changing type to <f8 fixes the problem.
        Exchanges = np.core.records.fromrecords(Ea, dtype=[('id', 'U11'),
                                                           ('uptake', '<f8'),
                                                           ('name', 'U70'),
                                                           ('formula', 'U30'),
                                                           ('charge', '<f8')])

        # Create empty list to put active fluxes into
        nonZeros = []

        # Take only those with uptake != 0
        for row in Exchanges:
            if row['uptake'] != 0:
                nonZeros.append(row)

        # Don't stop even if the reaction doesn't proceed.
        all_reactions.append(full_name)
        if len(nonZeros) > 0:
            pass
        else:
            all_reactions.append("Reaction does not proceed")
            all_reactions.append("\n")
            continue

        dat_init = []

        # Biomass doesn't go in the .dat file
        # Will need to be carefull about other biomass formula styles, or other
        # compounds with null or blank formulas. We need to keep names too,
        # some compounds have the same formula.
        # This also removes weird characters from the names
        for row in nonZeros:
            row['formula'] = row['formula'].replace('*', '')
            if (row['name'].startswith('lps') or
                    row['name'].startswith('cpd11649')):
                row['formula'] = 'C175H317N5O101P6'
            if (row['formula'] != 'null' and row['formula'] != ''
                    and row['name'] != 'Biomass'):
                dat_init.append([row['formula'],
                                 row['name'],
                                 'tbd',
                                 row['charge'],
                                 'tbd'])

        # We want to know if more than just biomass gets omitted
        bioskips = f"{(len(nonZeros)-len(dat_init))} biomass compounds"
        " omitted from reaction string"

        i = 0
        full = len(dat_init)
        while i < full:
            form = ptab.formula(dat_init[i][0])
            mw = round(form.mass, 4)
            dat_init[i][4] = mw
            i += 1

        # Preprocessing to facilitate ion size lookup
        for row in dat_init:
            # If charge is 0, debye huckel goes to 0, so a shouldn't matter
            if row[3] == 0:
                row[2] = 0
            # Use ion size parameter from database wherever possible, fill in
            # "unknown" if not
            if row[2] == 'tbd':
                # All iron cpds have the same formula, we need to differentiate
                if row[0] == 'Fe':
                    row[2] = ion_size.get(row[1], "unknown")
                else:
                    row[2] = ion_size.get(row[0], "unknown")
            # Not available for all ions. Rough estimate based on Kielland 1937
            if row[2] == 'unknown':
                # single-element ions and small inorganics
                if len(row[0]) <= 3:
                    row[2] = 3
                # For unknown complex organics, which we assume will be long
                else:
                    row[2] = 6

                # Rename name and formula columns to PFLOTRAN conventions
        dat_final = rename(dat_init, 0, 1, 3)

        with open(datfile, 'w+') as csvfile:
            datwriter = csv.writer(csvfile, delimiter=' ',
                                   lineterminator='\n',
                                   quotechar="'",
                                   quoting=csv.QUOTE_NONNUMERIC)
            for row in dat_final:
                datwriter.writerow(row)

        # Create a dict that can be used to map renamed compounds to cpd ids
        # and formulas. Lookup by name (formulas are not unique).
        formula_lookup = {v[2]:
                          [k, v[2], v[3]] for k, v in sdict.items()}

        # Also compile a dict of biomasses
        biomass = {}
        # This should support community models, but has not been tested
        for x in model_data['biomasses']:
            bcpds = {}
            biomass.update({x['id']: {'ID': x['id'],
                                      "Ctot": 0,
                                      "Htot": 0,
                                      "Ntot": 0,
                                      "Otot": 0,
                                      "Ptot": 0,
                                      "Stot": 0}})
            for y in x['biomasscompounds']:
                bid = y['modelcompound_ref'].split("/")[-1]

                # We need to remove pseudoformulas/ unknowns, periodictable
                # formula doesn't understand them. However, we don't want to
                # remove valid formulas with Ra, Rn, Re, Rh, Rg, Rb, Ru, Rf.

                badform = ['R', 'Protein', 'DNA', 'RNA', 'synthesis',
                           'unknown', 'undefined', 'Biomass']
                Relements = ['Ra', 'Rn', 'Re', 'Rh', 'Rg', 'Rb', 'Ru', 'Rf']
                if (any(bad in fdict[bid][1] for bad in badform)
                        and not any(el in fdict[bid][1] for el in Relements)):
                    continue
                # bcpds.update({bid : {'coefficient' : y['coefficient']}})
                formula = ptab.formula(fdict[bid][1])
                components = formula.atoms
                compdict = {}
                for c, d in components.items():
                    compdict.update({str(c): d})
                bcpds.update({bid: {'coefficient': y['coefficient'],
                                    'formula': formula,
                                    'components': compdict}})
                complist = ['C', 'H', 'N', 'O', 'P', 'S']

                # Keep the output in mmol, as that's how we'd want to upload it
                # to KBase. If we decide we need mass (g) instead, we can get
                # it bymultiplying the += term by ptab.formula(thing).mass

                for thing in complist:
                    if thing in compdict:
                        mmol = (-1 * bcpds[bid]['coefficient'] *
                                compdict[thing])
                        biomass[x['id']][thing+"tot"] += mmol
        biofile = pulldir+full_name+"_biomass.json"
        biojson = json.dumps(biomass)
        bio_writer = open(biofile, "w")
        bio_writer.write(biojson)
        bio_writer.close()

        # Need an initialized table with a unique lookup in it
        # Formula is NOT unique (e.g. 3 compounds with formula 'Fe')
        media_init = []

        # This should exclude biomass also
        for row in nonZeros:
            row['formula'] = row['formula'].replace('*', '')
            if (row['name'].startswith('lps')
                    or row['name'].startswith('cpd11649')):
                row['formula'] = 'C175H317N5O101P6'
            if (row['formula'] != 'null' and row['formula'] != ''
                    and row['name'] != 'Biomass'):
                media_init.append([row['formula'],
                                   row['formula'],
                                   row['name'],
                                   row['charge']])

        media_final = rename(media_init, 0, 2, 3, 0)
        media_dict = {}
        for row in media_final:
            # Use "final display name" as key so we can use it for lookup later
            # Formula does NOT work as the lookup, not unique
            # (e.g. Fe = 3 different compounds)
            lookup_key = row[2]
            media_dict.update({row[0]:
                               [formula_lookup[lookup_key][0].split('_')[0],
                                formula_lookup[lookup_key][1], row[1]]})

        # Might as well print this to a file for later recall
        mediafile = mediadir+full_name+"_media.json"
        mediajson = json.dumps(media_dict)
        media_writer = open(mediafile, "w")
        media_writer.write(mediajson)
        media_writer.close()

        # Biomass term doesn't belong in the main reaction string
        fluxes = []
        bio_yield = []
        for row in nonZeros:
            if (row['formula'] == 'null' or row['formula'] == ''
                    or row['name'] == 'Biomass'):
                bio_yield.append([row['name'],
                                  row['uptake']])
            else:
                fluxes.append([row['formula'],
                               row['uptake'],
                               row['charge'],
                               row['name']])

        fluxnames = rename(fluxes, 0, 3, 2)

        # PFLOTRAN doesn't include H2O in reaction strings, so:

        rxn = []
        for row in fluxnames:
            if row[0] != 'H2O(aq)':
                rxn.append(row)

        # Some necessary operators for the final string
        space = " "
        plus = " + "
        becomes = " -> "
        reacts = []
        prods = []

        # Split up products and reactants so we can take the
        # all coefficients as positive
        for row in rxn:
            if row[1] < 0:
                prods.append(row)
            else:
                reacts.append(row)

        # Make a string of all reactants first. Iterate over all but the last
        reactstring = ''
        ireacts = 0
        totreacts = len(reacts)

        # Need to handle messy models with NO uptake
        if totreacts == 0:
            reactstring += ("No reactants for this reaction. Something is"
                            "probably wrong.")
        else:
            while ireacts < (totreacts - 1):
                reactstring += (str(abs(reacts[ireacts][1])) + space
                                + reacts[ireacts][0] + plus)
                ireacts += 1

            # Last one can't have + after it, so do this one separately
            reactstring += (str(abs(reacts[ireacts][1])) + space +
                            reacts[ireacts][0])

        # Handle products the same way
        prodstring = ''
        iprods = 0
        totprods = len(prods)

        # Need to handle messy models with ONLY uptake
        if totprods == 0:
            prodstring += ("No products for this reaction."
                           " Something is probably wrong.")
        else:
            while iprods < (totprods - 1):
                prodstring += (str(abs(prods[iprods][1])) + space +
                               prods[iprods][0] + plus)
                iprods += 1

            # Last one can't have + after it, so do this one separately
            prodstring += str(abs(prods[iprods][1])) + space + prods[iprods][0]

        overall = ''

        overall = reactstring + becomes + prodstring

        # Print biomass yield(s) separately.
        # Includes handling of multiple biomass terms (for community models)
        yield_terms = ''
        if len(bio_yield) == 0:
            bioyield = 'No biomass yield'
        elif len(bio_yield) == 1:
            bioyield = 'YIELD ' + str(abs(bio_yield[0][1]))
        else:
            ib = 0
            while ib < len(bio_yield) - 1:
                print(bio_yield[0][0] + ' YIELD ' + str(abs(bio_yield[0][1])))
                ib += 1

        # Print all reactions to a file at the end. Precede reactions with the
        # name of the fba, to make this file easier to read when there are many
        all_reactions.append(bioskips)
        all_reactions.append(overall)
        all_reactions.append(bioyield)
        all_reactions.append("\n")

# Now print the full reaction string
rxnfilename = pulldir+nar_name+"_OverallReactions.txt"
with open(rxnfilename, 'w+') as rxnfile:
    for row in all_reactions:
        rxnfile.write(row+"\n")
