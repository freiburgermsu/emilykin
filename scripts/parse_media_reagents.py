"""Parse the chemical composition of the synthetic media and feed reagents into JSON.

Source: 'Synthetic media and other reagents preparation_EM.xlsx'.

For every reagent the script records the as-PREPARED (stock / feed-bottle) and the
in-REACTOR (influent) concentration, in both mg/L and mM, then dissociates each
compound into its constituent ions and sums the ion concentrations.

The workbook is a hand-kept lab log, so concentrations are pulled from explicitly
identified cells (noted alongside each entry) rather than auto-detected. Conversions
follow the sheet's own conventions:
  * trace-element stocks are dosed at 0.3 mL/L (solution II) or 1 mL/L (solution I),
    so in-reactor = stock * dose_mL_per_L / 1000;
  * feed targets given "as COD" use the sheet's COD-equivalence factors;
  * feed targets given "as element-N" / "as element-P" are converted to the parent
    compound via molar mass (compound = target * MW_compound / atomic_mass_element).

Output: media_reagents_ions.json
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# Atomic / ionic molar masses (g/mol) used to weight the dissociated ions.
# ---------------------------------------------------------------------------
ION_MASS = {
    'Mg2+': 24.305, 'Ca2+': 40.078, 'Na+': 22.990, 'K+': 39.098, 'NH4+': 18.039,
    'Fe2+': 55.845, 'Fe3+': 55.845, 'Cu2+': 63.546, 'Mn2+': 54.938, 'Zn2+': 65.38,
    'Co2+': 58.933, 'Ni2+': 58.693, 'H+': 1.008,
    'Cl-': 35.453, 'OH-': 17.007, 'I-': 126.904, 'HCO3-': 61.017,
    'SO4_2-': 96.06, 'MoO4_2-': 159.94, 'NO2-': 46.005, 'NO3-': 62.004,
    'H2PO4-': 96.987, 'HPO4_2-': 95.979,
    'CH3COO-': 59.044, 'CH3CH2COO-': 73.071,
}
ATOMIC_N = 14.007
ATOMIC_P = 30.974

# ---------------------------------------------------------------------------
# Compound database: molar mass (g/mol) and ionic dissociation {ion: count}.
# `ions == {}`  -> treated as molecular (non/weakly-dissociating): no ion total.
# Molar masses match the workbook where it lists them, else standard values.
# ---------------------------------------------------------------------------
COMPOUNDS = {
    'MgSO4':         dict(mw=120.37, ions={'Mg2+': 1, 'SO4_2-': 1}),
    'MgSO4·7H2O':    dict(mw=246.48, ions={'Mg2+': 1, 'SO4_2-': 1}),
    'CaCl2·2H2O':    dict(mw=147.01, ions={'Ca2+': 1, 'Cl-': 2}),
    'NaHCO3':        dict(mw=84.01,  ions={'Na+': 1, 'HCO3-': 1}),
    'NH4Cl':         dict(mw=53.49,  ions={'NH4+': 1, 'Cl-': 1}),
    'FeCl3·6H2O':    dict(mw=270.30, ions={'Fe3+': 1, 'Cl-': 3}),
    'H3BO3':         dict(mw=61.83,  ions={}),                       # boric acid, molecular
    'CuSO4·5H2O':    dict(mw=249.68, ions={'Cu2+': 1, 'SO4_2-': 1}),
    'KI':            dict(mw=166.00, ions={'K+': 1, 'I-': 1}),
    'MnCl2·4H2O':    dict(mw=197.91, ions={'Mn2+': 1, 'Cl-': 2}),
    'Na2MoO4·2H2O':  dict(mw=241.95, ions={'Na+': 2, 'MoO4_2-': 1}),
    'ZnSO4·7H2O':    dict(mw=287.56, ions={'Zn2+': 1, 'SO4_2-': 1}),
    'CoCl2·6H2O':    dict(mw=237.93, ions={'Co2+': 1, 'Cl-': 2}),
    'NiCl2·6H2O':    dict(mw=237.69, ions={'Ni2+': 1, 'Cl-': 2}),
    'EDTA':          dict(mw=292.24, ions={}),                       # chelator, molecular
    'FeSO4·7H2O':    dict(mw=278.01, ions={'Fe2+': 1, 'SO4_2-': 1}),
    'HCl':           dict(mw=36.46,  ions={'H+': 1, 'Cl-': 1}),
    'NaOH':          dict(mw=40.00,  ions={'Na+': 1, 'OH-': 1}),
    'NaNO2':         dict(mw=69.00,  ions={'Na+': 1, 'NO2-': 1}),
    'NaNO3':         dict(mw=84.99,  ions={'Na+': 1, 'NO3-': 1}),
    'CH3COONa':      dict(mw=82.03,  ions={'Na+': 1, 'CH3COO-': 1}),
    'CH3CH2COONa':   dict(mw=96.06,  ions={'Na+': 1, 'CH3CH2COO-': 1}),
    'KH2PO4':        dict(mw=136.09, ions={'K+': 1, 'H2PO4-': 1}),
    'K2HPO4':        dict(mw=174.18, ions={'K+': 2, 'HPO4_2-': 1}),
    'Yeast extract': dict(mw=None,   ions={}),                       # undefined composition
    'ATU':           dict(mw=116.19, ions={}),                       # allylthiourea, molecular
}


def conc_block(prepared_mg_L, reactor_mg_L, formula):
    """Build a per-compound concentration block with mg/L + mM at both bases."""
    mw = COMPOUNDS[formula]['mw']

    def mM(mg_L):
        if mg_L is None or mw is None:
            return None
        return mg_L / mw

    block = {
        'prepared_mg_L': round(prepared_mg_L, 6) if prepared_mg_L is not None else None,
        'prepared_mM': round(mM(prepared_mg_L), 6) if mM(prepared_mg_L) is not None else None,
        'in_reactor_mg_L': round(reactor_mg_L, 6) if reactor_mg_L is not None else None,
        'in_reactor_mM': round(mM(reactor_mg_L), 6) if mM(reactor_mg_L) is not None else None,
    }
    return block


def compound_entry(name, formula, *, prepared_mg_L=None, reactor_mg_L=None, note=None):
    info = COMPOUNDS[formula]
    entry = {
        'name': name,
        'formula': formula,
        'molar_mass_g_per_mol': info['mw'],
        'concentration': conc_block(prepared_mg_L, reactor_mg_L, formula),
        'dissociation': info['ions'] or 'molecular (not dissociated)',
    }
    if note:
        entry['note'] = note
    return entry


def ion_totals(entries):
    """Sum ion concentrations (mM and mg/L) over a list of compound entries,
    separately for the prepared and in-reactor bases."""
    totals = {}
    for e in entries:
        ions = COMPOUNDS[e['formula']]['ions']
        if not ions:
            continue
        for basis, mM_key in (('prepared', 'prepared_mM'), ('in_reactor', 'in_reactor_mM')):
            cmpd_mM = e['concentration'][mM_key]
            if cmpd_mM is None:
                continue
            for ion, count in ions.items():
                slot = totals.setdefault(ion, {})
                slot[f'{basis}_mM'] = slot.get(f'{basis}_mM', 0.0) + cmpd_mM * count
    out = {}
    for ion, slot in sorted(totals.items()):
        rec = {}
        for basis in ('prepared', 'in_reactor'):
            mM = slot.get(f'{basis}_mM')
            if mM is not None:
                rec[f'{basis}_mM'] = round(mM, 6)
                rec[f'{basis}_mg_L'] = round(mM * ION_MASS[ion], 6)
        out[ion] = rec
    return out


def section(name, sheet, description, entries, *, notes=None):
    sec = {
        'sheet': sheet,
        'description': description,
        'compounds': entries,
        'ion_totals': ion_totals(entries),
    }
    if notes:
        sec['notes'] = notes
    return sec


# Convenience converters --------------------------------------------------
def from_cod(cod_mg_L, cod_factor):
    """compound mg/L from a target 'as COD' (cod_factor = g COD per g compound)."""
    return cod_mg_L / cod_factor


def from_element(target_mg_L, mw_compound, atomic_mass, n_atoms=1):
    """compound mg/L from a target expressed as element-N or element-P."""
    return target_mg_L * mw_compound / (atomic_mass * n_atoms)


def trace_reactor(stock_mg_L, dose_mL_per_L):
    return stock_mg_L * dose_mL_per_L / 1000.0


# ===========================================================================
# 1. Synthetic wastewater media — recipe 1 (without COD, P and nitrate)
#    Sheet 'Media ingredients_conc.' rows 2-10 (major: row6 influent, row9 stock g/L)
#    trace soln II dosed 0.3 mL/L (row6 = stock mg/L); trace soln I dosed 1 mL/L.
# ===========================================================================
media1 = []
# major components: prepared = stock g/L * 1000; in-reactor = influent mg/L (row 6)
media1 += [
    compound_entry('Magnesium sulfate (anhydrous)', 'MgSO4', prepared_mg_L=66.79 * 1000, reactor_mg_L=66.79,
                   note='Sheet also lists MgSO4·7H2O separately; both are summed in ion_totals (see notes).'),
    compound_entry('Magnesium sulfate heptahydrate', 'MgSO4·7H2O', prepared_mg_L=90 * 1000, reactor_mg_L=90),
    compound_entry('Calcium chloride dihydrate', 'CaCl2·2H2O', prepared_mg_L=56 * 1000, reactor_mg_L=14),
    compound_entry('Yeast extract', 'Yeast extract', prepared_mg_L=10 * 1000, reactor_mg_L=1),
    compound_entry('Allylthiourea (ATU)', 'ATU', prepared_mg_L=40 * 1000, reactor_mg_L=20),
    compound_entry('Sodium bicarbonate (buffer)', 'NaHCO3', prepared_mg_L=65 * 1000, reactor_mg_L=500),
]
# trace elements solution II (dosed 0.3 mL/L); prepared = stock mg/L (row 6)
TRACE_II = [
    ('Iron(III) chloride hexahydrate', 'FeCl3·6H2O', 1500),
    ('Boric acid', 'H3BO3', 150),
    ('Copper(II) sulfate pentahydrate', 'CuSO4·5H2O', 30),
    ('Potassium iodide', 'KI', 180),
    ('Manganese(II) chloride tetrahydrate', 'MnCl2·4H2O', 120),
    ('Sodium molybdate dihydrate', 'Na2MoO4·2H2O', 60),
    ('Zinc sulfate heptahydrate', 'ZnSO4·7H2O', 120),
    ('Cobalt(II) chloride hexahydrate', 'CoCl2·6H2O', 150),
    ('Nickel chloride hexahydrate', 'NiCl2·6H2O', 630),
    ('EDTA', 'EDTA', 15000),
]
for nm, fm, stock in TRACE_II:
    media1.append(compound_entry(nm, fm, prepared_mg_L=stock,
                                 reactor_mg_L=trace_reactor(stock, 0.3),
                                 note='trace elements solution II, dosed 0.3 mL/L'))
# trace elements solution I (dosed 1 mL/L)
for nm, fm, stock in [('Iron(II) sulfate heptahydrate', 'FeSO4·7H2O', 8340), ('EDTA', 'EDTA', 8340)]:
    media1.append(compound_entry(nm, fm, prepared_mg_L=stock,
                                 reactor_mg_L=trace_reactor(stock, 1.0),
                                 note='trace elements solution I, dosed 1 mL/L'))

# ===========================================================================
# 2. Synthetic wastewater media — recipe 2 (free of VFA, P, NO3-N and NO2-N)
#    Sheet 'Media ingredients_conc.' rows 32-35 (row35 = Conc. mg/L).
#    Major values are in-reactor mg/L; trace values are stock mg/L (dosed 0.3 mL/L).
# ===========================================================================
media2 = [
    compound_entry('Ammonium chloride', 'NH4Cl', reactor_mg_L=107),
    compound_entry('Magnesium sulfate heptahydrate', 'MgSO4·7H2O', reactor_mg_L=90),
    compound_entry('Calcium chloride dihydrate', 'CaCl2·2H2O', reactor_mg_L=14),
    compound_entry('Yeast extract', 'Yeast extract', reactor_mg_L=1),
    compound_entry('Allylthiourea (ATU)', 'ATU', reactor_mg_L=20),
    compound_entry('Sodium bicarbonate (buffer)', 'NaHCO3', reactor_mg_L=200),
]
TRACE_R2 = [
    ('Iron(III) chloride hexahydrate', 'FeCl3·6H2O', 1500),
    ('Boric acid', 'H3BO3', 150),
    ('Copper(II) sulfate pentahydrate', 'CuSO4·5H2O', 30),
    ('Potassium iodide', 'KI', 180),
    ('Manganese(II) chloride tetrahydrate', 'MnCl2·4H2O', 120),
    ('Sodium molybdate dihydrate', 'Na2MoO4·2H2O', 60),
    ('Zinc sulfate heptahydrate', 'ZnSO4·7H2O', 120),
    ('Cobalt(II) chloride hexahydrate', 'CoCl2·6H2O', 150),
    ('EDTA', 'EDTA', 10000),
]
for nm, fm, stock in TRACE_R2:
    media2.append(compound_entry(nm, fm, prepared_mg_L=stock,
                                 reactor_mg_L=trace_reactor(stock, 0.3),
                                 note='trace elements solution, dosed 0.3 mL/L'))

# ===========================================================================
# 3. Acid / base feed  (sheet 'acid and base feed')
#    HCl purchased 12.1 mol/L; feed bottle 0.75 mol/L. NaOH feed bottle 0.5 mol/L.
#    No fixed in-reactor target (dosed to a pH set-point) -> in_reactor = None.
# ===========================================================================
acid_base = [
    compound_entry('Hydrochloric acid (feed)', 'HCl',
                   prepared_mg_L=0.75 * COMPOUNDS['HCl']['mw'] * 1000,
                   note='feed-bottle 0.75 mol/L; purchased stock 12.1 mol/L (=441.2 g/L). '
                        'Dosed to a pH set-point, so no fixed in-reactor concentration.'),
    compound_entry('Sodium hydroxide (feed)', 'NaOH',
                   prepared_mg_L=0.5 * COMPOUNDS['NaOH']['mw'] * 1000,
                   note='feed-bottle 0.5 mol/L. Dosed to a pH set-point, so no fixed in-reactor concentration.'),
]

# ===========================================================================
# 4. COD + P + N feed  (sheet 'COD+P feed', representative recipe dated 2023-01-30)
#    prepared = "Chemical conc. in 1L bottle (g/L)" (row 15); in-reactor from row-11 targets.
# ===========================================================================
cod_feed = [
    compound_entry('Sodium acetate', 'CH3COONa',
                   prepared_mg_L=26.26795343283582 * 1000,
                   reactor_mg_L=from_cod(75, 0.7804878048780488),
                   note='carbon source; in-reactor from 75 mg/L-COD target (0.7805 gCOD/g).'),
    compound_entry('Sodium propionate', 'CH3CH2COONa',
                   prepared_mg_L=5.861933301492537 * 1000,
                   reactor_mg_L=from_cod(25, 1.1658165920682837),
                   note='carbon source; in-reactor from 25 mg/L-COD target (1.166 gCOD/g).'),
    compound_entry('Potassium phosphate monobasic', 'KH2PO4',
                   prepared_mg_L=12.001278866056811 * 1000,
                   reactor_mg_L=from_element(10, 136.09, ATOMIC_P),
                   note='P source; in-reactor from 10 mg/L PO4-P target.'),
    compound_entry('Sodium nitrite', 'NaNO2',
                   prepared_mg_L=5.9284487883365635 * 1000,
                   reactor_mg_L=from_element(5, 69.00, ATOMIC_N),
                   note='N source; in-reactor from 5 mg/L NO2-N target.'),
]

# ===========================================================================
# 5. Nitrate feed  (sheet 'nitrate feed', first documented recipe 2021-08-06)
# ===========================================================================
nitrate_feed = [
    compound_entry('Sodium nitrate', 'NaNO3',
                   prepared_mg_L=91.07142857142857 * 1000,
                   reactor_mg_L=from_element(100, 84.99, ATOMIC_N),
                   note='first documented recipe (2021-08-06): bottle 91.07 g/L, target 100 mg/L NO3-N. '
                        'Feed is time-varying; later recipes target 10-20 mg/L NO3-N (do not exceed 20).'),
]

# ===========================================================================
# 6. Nitrite feed test  (sheet 'nitrite feed test', representative 2021-12-06 recipe)
# ===========================================================================
nitrite_feed = [
    compound_entry('Sodium nitrite', 'NaNO2',
                   prepared_mg_L=201.6233766233766 * 1000,
                   reactor_mg_L=from_element(15, 69.00, ATOMIC_N),
                   note='representative recipe (2021-12-06): bottle 201.6 g/L, target 15 mg/L NO2-N. '
                        'Earlier recipes targeted 3 and 10 mg/L NO2-N (bottles 40.3 and 134.4 g/L).'),
]

# ===========================================================================
# 7. Biomass storage media  (sheet 'Biomass storage media')  — prepared only
# ===========================================================================
biomass_storage = [
    compound_entry('Sodium nitrite', 'NaNO2', prepared_mg_L=0.0345 * 1000,
                   note='0.5 mM NaNO2 in 1 L (34.5 mg).'),
    compound_entry('Ammonium chloride', 'NH4Cl', prepared_mg_L=0.02675 * 1000,
                   note='0.5 mM NH4Cl in 1 L (26.75 mg).'),
]

# ---------------------------------------------------------------------------
result = {
    '_meta': {
        'source_file': 'Synthetic media and other reagents preparation_EM.xlsx',
        'description': 'Molecules/ions and their concentrations for the synthetic '
                       'wastewater media and feed reagents of the EmilyKin bioreactor.',
        'concentration_bases': {
            'prepared': 'as-prepared stock / feed-bottle concentration',
            'in_reactor': 'concentration delivered to the reactor influent (after dilution / dosing)',
        },
        'units': 'mg_L = milligrams per litre; mM = millimoles per litre',
        'ion_naming': 'underscores denote charge magnitude, e.g. SO4_2- = sulfate (SO4^2-)',
        'assumptions': [
            'Trace-element stocks dosed at 0.3 mL/L (solution II) or 1 mL/L (solution I): '
            'in_reactor = prepared * dose / 1000.',
            'H3BO3 (boric acid) and EDTA are treated as molecular (not added to ion totals); '
            'EDTA in practice chelates the trace-metal cations.',
            'Yeast extract has no defined formula, so no molar or ion concentration is computed.',
            'Recipe 1 lists both anhydrous MgSO4 and MgSO4·7H2O; both are included as written and '
            'their Mg2+/SO4_2- contributions are summed, which may double-count magnesium.',
            'For the media, prepared concentrations live in separate stock bottles (each major '
            'component + one shared trace-element bottle), so a single prepared ion_total is not a '
            'physical solution; the in_reactor ion_total is the physically mixed reactor influent.',
            'Feed targets given as COD/element-N/element-P were converted to the parent compound '
            'using the molar masses and COD factors in the sheet.',
            'COD+P, nitrate and nitrite feeds are time-varying; one representative documented recipe '
            'is captured per feed (date noted).',
        ],
    },
    'media': {
        'synthetic_wastewater_no_COD_P_nitrate': section(
            'recipe 1', 'Media ingredients_conc.',
            'Synthetic wastewater media without COD, P and nitrate.', media1,
            notes='Major-component in-reactor values are the sheet "Influent Conc. (mg/L)" row; '
                  'prepared values are the stock "g/L" row x1000.'),
        'synthetic_wastewater_VFA_P_NO3_NO2_free': section(
            'recipe 2', 'Media ingredients_conc.',
            'Synthetic wastewater media free of VFA, P, NO3-N and NO2-N.', media2,
            notes='Major components give only an in-reactor "Conc. (mg/L)" (no separate stock).'),
    },
    'feeds': {
        'acid_base_feed': section('acid/base', 'acid and base feed',
                                  'Acid (HCl) and base (NaOH) pH-control feeds.', acid_base),
        'COD_P_N_feed': section('COD+P+N', 'COD+P feed',
                                'Carbon (acetate/propionate) + P (KH2PO4) + N (NaNO2) feed.', cod_feed),
        'nitrate_feed': section('nitrate', 'nitrate feed', 'Nitrate (NaNO3) N feed.', nitrate_feed),
        'nitrite_feed': section('nitrite', 'nitrite feed test', 'Nitrite (NaNO2) N feed.', nitrite_feed),
        'biomass_storage_media': section('biomass storage', 'Biomass storage media',
                                         'Biomass storage media (0.5 mM NaNO2 + 0.5 mM NH4Cl).',
                                         biomass_storage),
    },
}

out = REPO / 'media_reagents_ions.json'
out.write_text(json.dumps(result, indent=2, ensure_ascii=False))
print(f'wrote {out}')

# quick self-check printout
for path in [('media', 'synthetic_wastewater_no_COD_P_nitrate'),
             ('media', 'synthetic_wastewater_VFA_P_NO3_NO2_free')]:
    sec = result[path[0]][path[1]]
    print(f"\n{path[1]} in-reactor ions (mg/L):")
    for ion, rec in sec['ion_totals'].items():
        if 'in_reactor_mg_L' in rec:
            print(f"  {ion:>10}: {rec['in_reactor_mg_L']}")
