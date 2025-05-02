
import re
import pandas as pd
import json

# run
# python3 analysis_diet.py

# MATLAB code snippet provided by the user
matlab_code = """
function [modelOut]=useHighFiberDiet_AGORA(modelIn)
% assign a high fiber diet for the AGORA microbes
% used as constraints to predict the growth rates shown in Table S6 in
% Magnusdottir et al., Nat Biotech 2017
% defined similarly as in Heinken and Thiele, AEM 2015
% anaerobic medium

model=modelIn;
model=changeRxnBounds(model,model.rxns(strmatch('EX_',model.rxns)),0,'l');
% simple sugars and starch
model=changeRxnBounds(model,'EX_glc_D(e)',-0.03947,'l');
model=changeRxnBounds(model,'EX_gal(e)',-0.03947,'l');
model=changeRxnBounds(model,'EX_man(e)',-0.03947,'l');
model=changeRxnBounds(model,'EX_mnl(e)',-0.03947,'l');
model=changeRxnBounds(model,'EX_fuc_L(e)',-0.03947,'l');
model=changeRxnBounds(model,'EX_glcn(e)',-0.03947,'l');
model=changeRxnBounds(model,'EX_rmn(e)',-0.03947,'l');
model=changeRxnBounds(model,'EX_arab_L(e)',-0.04737,'l');
model=changeRxnBounds(model,'EX_drib(e)',-0.04737,'l');
model=changeRxnBounds(model,'EX_rib_D(e)',-0.04737,'l');
model=changeRxnBounds(model,'EX_xyl_D(e)',-0.04737,'l');
model=changeRxnBounds(model,'EX_oxa(e)',-0.11842,'l');
model=changeRxnBounds(model,'EX_lcts(e)',-0.01974,'l');
model=changeRxnBounds(model,'EX_malt(e)',-0.01974,'l');
model=changeRxnBounds(model,'EX_sucr(e)',-0.01974,'l');
model=changeRxnBounds(model,'EX_melib(e)',-0.01974,'l');
model=changeRxnBounds(model,'EX_cellb(e)',-0.01974,'l');
model=changeRxnBounds(model,'EX_tre(e)',-0.01974,'l');
model=changeRxnBounds(model,'EX_strch1(e)',-0.06818,'l');
% fiber
model=changeRxnBounds(model,'EX_amylopect900(e)',-0.0003472222,'l');
model=changeRxnBounds(model,'EX_amylose300(e)',-0.0010416667,'l');
model=changeRxnBounds(model,'EX_arabinan101(e)',-0.0036836935,'l');
model=changeRxnBounds(model,'EX_arabinogal(e)',-0.0004854997,'l');
model=changeRxnBounds(model,'EX_arabinoxyl(e)',-0.0067934783,'l');
model=changeRxnBounds(model,'EX_bglc(e)',-0.0000015625,'l');
model=changeRxnBounds(model,'EX_cellul(e)',-0.0006250000,'l');
model=changeRxnBounds(model,'EX_dextran40(e)',-0.0039062500,'l');
model=changeRxnBounds(model,'EX_galmannan(e)',-0.0003125000,'l');
model=changeRxnBounds(model,'EX_glcmannan(e)',-0.0007284382,'l');
model=changeRxnBounds(model,'EX_homogal(e)',-0.0028409091,'l');
model=changeRxnBounds(model,'EX_inulin(e)',-0.0104166667,'l');
model=changeRxnBounds(model,'EX_kestopt(e)',-0.0625000000,'l');
model=changeRxnBounds(model,'EX_levan1000(e)',-0.0003125000,'l');
model=changeRxnBounds(model,'EX_lmn30(e)',-0.0104166667,'l');
model=changeRxnBounds(model,'EX_lichn(e)',-0.0018382353,'l');
model=changeRxnBounds(model,'EX_pect(e)',-0.0007396450,'l');
model=changeRxnBounds(model,'EX_pullulan1200(e)',-0.0002604167,'l');
model=changeRxnBounds(model,'EX_raffin(e)',-0.1041666667,'l');
model=changeRxnBounds(model,'EX_rhamnogalurI(e)',-0.0003210616,'l');
model=changeRxnBounds(model,'EX_rhamnogalurII(e)',-0.0059148265,'l');
model=changeRxnBounds(model,'EX_starch1200(e)',-0.0002604167,'l');
model=changeRxnBounds(model,'EX_xylan(e)',-0.0007102273,'l');
model=changeRxnBounds(model,'EX_xyluglc(e)',-0.0002912395,'l');
% fat
model=changeRxnBounds(model,{'EX_arachd(e)'},-0.001664,'l');
model=changeRxnBounds(model,{'EX_chsterol(e)'},-0.002479,'l');
model=changeRxnBounds(model,{'EX_glyc(e)'},-0.899827,'l');
model=changeRxnBounds(model,{'EX_hdca(e)'},-0.198185,'l');
model=changeRxnBounds(model,{'EX_hdcea(e)'},-0.018258,'l');
model=changeRxnBounds(model,{'EX_lnlc(e)'},-0.179555,'l');
model=changeRxnBounds(model,{'EX_lnlnca(e)'},-0.008783,'l');
model=changeRxnBounds(model,{'EX_lnlncg(e)'},-0.008783,'l');
model=changeRxnBounds(model,{'EX_ocdca(e)'},-0.084641,'l');
model=changeRxnBounds(model,{'EX_ocdcea(e)'},-0.340722,'l');
model=changeRxnBounds(model,{'EX_octa(e)'},-0.006471,'l');
model=changeRxnBounds(model,{'EX_ttdca(e)'},-0.034338,'l');
% protein
model=changeRxnBounds(model,{'EX_ala_L(e)';'EX_ser_L(e)';'EX_cys_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_arg_L(e)';'EX_ile_L(e)';'EX_leu_L(e)';'EX_lys_L(e)';'EX_his_L(e)'},-0.15,'l');
model=changeRxnBounds(model,{'EX_asn_L(e)';'EX_asp_L(e)';'EX_thr_L(e)'},-0.225,'l');
model=changeRxnBounds(model,{'EX_glu_L(e)';'EX_met_L(e)';'EX_gln_L(e)';'EX_pro_L(e)';'EX_val_L(e)'},-0.18,'l');
model=changeRxnBounds(model,{'EX_phe_L(e)';'EX_tyr_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_gly(e)'},-0.45,'l');
model=changeRxnBounds(model,{'EX_trp_L(e)'},-0.08182,'l');
% minerals and vitamins
model=changeRxnBounds(model,{'EX_thm(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_ribflv(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_nac(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_btn(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_pnto_R(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_pydxn(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_pydx(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_pydam(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_fol(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_adpcbl(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_pheme(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_cbl1(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_chol(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_h2s(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_so4(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_spmd(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_na1(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_ca2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_cobalt2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_cl(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_k(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_fe2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_fe3(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_fe3dcit(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_mg2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_mn2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_mobd(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_cu2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_zn2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_sel(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_pime(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_mqn8(e)'},-1,'l');
model=changeRxnBounds(model,'EX_meoh(e)',-10,'l');
model=changeRxnBounds(model,'EX_h2o(e)',-10,'l');
modelOut=model;
end
"""

# -----------------------------------------------------------------------------
# useHighFiberDiet_AGORA: apply a “high‐fiber” diet medium to AGORA genome‐scale 
# metabolic models (Magnúsdóttir et al., Nat. Biotechnol. 2017; Heinken & Thiele)
#
# 1. Reset all exchange reaction lower bounds to zero:
#      model = changeRxnBounds(model, model.rxns(strmatch('EX_',model.rxns)), 0, 'l')
#
# 2. Simple sugars & starch (e.g., fru, glc, glc_D, gal, maltose) at –0.03947 to –0.06818 mmol/gDW/h
#
# 3. Complex fibers (inulin, pectin, cellulose, arabinoxylan, etc.) at much lower rates 
#    (–0.00026 to –0.10417 mmol/gDW/h)
#
# 4. Lipids & fatty acids (glycerol, cholesterol, LCFAs) at –0.00166 to –0.89983 mmol/gDW/h
#
# 5. Protein (amino acid) blocks grouped by uptake families (–0.15 to –1 mmol/gDW/h)
#
# 6. Micronutrients (vitamins, trace metals) all set to –1 mmol/gDW/h to ensure cofactor availability
#
# 7. Special cases:
#      • Methanogen‐specific: EX_meoh(e) → –10 mmol/gDW/h  
#      • Water exchange:     EX_h2o(e) → –10 mmol/gDW/h
#
# Returns the modified model object ready for high‐fiber diet simulations.
#
# References:
#   • useHighFiberDiet_AGORA.m on GitHub: 
#     https://github.com/ThieleLab/CodeBase/blob/master/Simulations_Magnusdottir_NatBiotech_2017/useHighFiberDiet_AGORA.m
#   • Magnúsdóttir et al., “Generation of genome‐scale metabolic reconstructions for 773 members 
#     of the human gut microbiota,” Nat. Biotechnol. 35(1):81–89 (2017).
# -----------------------------------------------------------------------------

# Parse lines
rows = []
current_type = None
for line in matlab_code.splitlines():
    stripped = line.strip()
    if stripped.startswith('%'):
        # Section header
        current_type = stripped.lstrip('%').strip()
    else:
        m = re.search(
            r"changeRxnBounds\s*\(\s*model\s*,\s*(\{.*?\}|'[^']*')\s*,\s*([-\d\.eE]+)\s*,",
            stripped
        )
        if m and current_type:
            rxn_part, value_str = m.group(1), m.group(2)
            value = float(value_str)
            # Clean and split reaction identifiers
            rxn_clean = re.sub(r"[\{\}'\"]", "", rxn_part)
            rxns = re.split(r"[;,\s]+", rxn_clean)
            for rxn in rxns:
                if rxn:
                    rows.append({
                        'reaction': rxn,
                        'type': current_type,
                        'value': value
                    })

# Create DataFrame
df = pd.DataFrame(rows)


# 1) Define your essential metabolites (as in the R code)
essential_metabolites = [
      "EX_12dgr180(e)", "EX_26dap_M(e)", "EX_2dmmq8(e)", "EX_2obut(e)", "EX_3mop(e)", "EX_4abz(e)", "EX_4hbz(e)",
    "EX_ac(e)", "EX_acgam(e)", "EX_acmana(e)", "EX_acnam(e)", "EX_ade(e)", "EX_adn(e)", "EX_adocbl(e)", 
    "EX_ala_D(e)", "EX_ala_L(e)", "EX_amet(e)", "EX_amp(e)", "EX_arab_D(e)", "EX_arab_L(e)", "EX_arg_L(e)",
    "EX_asn_L(e)", "EX_btn(e)", "EX_ca2(e)", "EX_cbl1(e)", "EX_cgly(e)", "EX_chor(e)", "EX_chsterol(e)",
    "EX_cit(e)", "EX_cl(e)", "EX_cobalt2(e)", "EX_csn(e)", "EX_cu2(e)", "EX_cys_L(e)", "EX_cytd(e)", 
    "EX_dad_2(e)", "EX_dcyt(e)", "EX_ddca(e)", "EX_dgsn(e)", "EX_fald(e)", "EX_fe2(e)", "EX_fe3(e)", "EX_fol(e)", 
    "EX_for(e)", "EX_gal(e)", "EX_glc_D(e)", "EX_gln_L(e)", "EX_glu_L(e)", "EX_gly(e)", "EX_glyc(e)", 
    "EX_glyc3p(e)", "EX_gsn(e)", "EX_gthox(e)", "EX_gthrd(e)", "EX_gua(e)", "EX_h(e)", "EX_h2o(e)", 
    "EX_h2s(e)", "EX_his_L(e)", "EX_hxan(e)", "EX_ile_L(e)", "EX_k(e)", "EX_lanost(e)", "EX_leu_L(e)", 
    "EX_lys_L(e)", "EX_malt(e)", "EX_met_L(e)", "EX_mg2(e)", "EX_mn2(e)", "EX_mqn7(e)", "EX_mqn8(e)", 
    "EX_nac(e)", "EX_ncam(e)", "EX_nmn(e)", "EX_no2(e)", "EX_ocdca(e)", "EX_ocdcea(e)", "EX_orn(e)", 
    "EX_phe_L(e)", "EX_pheme(e)", "EX_pi(e)", "EX_pnto_R(e)", "EX_pro_L(e)", "EX_ptrc(e)", "EX_pydx(e)", 
    "EX_pydxn(e)", "EX_q8(e)", "EX_rib_D(e)", "EX_ribflv(e)", "EX_ser_L(e)", "EX_sheme(e)", "EX_so4(e)", 
    "EX_spmd(e)", "EX_thm(e)", "EX_thr_L(e)", "EX_thymd(e)", "EX_trp_L(e)", "EX_ttdca(e)", "EX_tyr_L(e)", 
    "EX_ura(e)", "EX_val_L(e)", "EX_xan(e)", "EX_xyl_D(e)", "EX_zn2(e)", "EX_glu_D(e)", "EX_melib(e)", 
    "EX_chtbs(e)", "EX_metsox_S_L(e)", "EX_hdca(e)", "EX_gam(e)", "EX_indole(e)", "EX_glcn(e)", "EX_coa(e)", 
    "EX_man(e)", "EX_fum(e)", "EX_succ(e)", "EX_no3(e)", "EX_ins(e)", "EX_uri(e)", "EX_drib(e)", "EX_pime(e)", 
    "EX_lac_L(e)", "EX_glypro(e)", "EX_urea(e)", "EX_duri(e)", "EX_h2(e)", "EX_mal_L(e)", "EX_tre(e)", 
    "EX_orot(e)", "EX_glymet(e)", "EX_glyleu(e)", "EX_pydx5p(e)", "EX_so3(e)", "EX_nh4(e)"
]

# 2) Clean IDs to match your parsed `df['reaction']`
essential_clean = [m.replace('(', '_').replace(')', '') for m in essential_metabolites]

# 3) Identify which essentials are missing
existing = set(df['reaction'])
missing   = [m for m in essential_clean if m not in existing]

# 4) Append them with type='essential' and default value (-0.1)
if missing:
    df_ess = pd.DataFrame({
        'reaction': missing,
        'type':     ['essential'] * len(missing),
        'value':    [-0.1] * len(missing)
    })
    # concat back onto your original df
    df = pd.concat([df, df_ess], ignore_index=True)


# Step 4: Unmapped compounds
unmapped_compounds = [
    'EX_asn_L(e)', 'EX_gln_L(e)', 'EX_crn(e)', 'EX_elaid(e)', 'EX_hdcea(e)',
    'EX_dlnlcg(e)', 'EX_adrn(e)', 'EX_hco3(e)', 'EX_sprm(e)', 'EX_carn(e)',
    'EX_7thf(e)', 'EX_Lcystin(e)', 'EX_hista(e)', 'EX_orn(e)', 'EX_ptrc(e)',
    'EX_creat(e)', 'EX_cytd(e)', 'EX_so4(e)'
]
# Clean IDs
unmapped_clean = [m.replace('(', '_').replace(')', '') for m in unmapped_compounds]

# Remove any existing unmapped compounds from df
df = df[~df['reaction'].isin(unmapped_clean)].copy()

# Append unmapped compounds with default flux 50 and type 'unmapped'
df_unmapped = pd.DataFrame({
    'reaction': unmapped_clean,
    'type':     ['unmapped'] * len(unmapped_clean),
    'value':    [-50.0] * len(unmapped_clean)
})
df = pd.concat([df, df_unmapped], ignore_index=True)

# Step 5: Add choline if missing
ch = 'EX_chol_e'
if ch not in df['reaction'].values:
    df = pd.concat([df, pd.DataFrame([{
        'reaction': ch,
        'type':     'choline',
        'value':    -41.251
    }])], ignore_index=True)

# Step 6: Micronutrients adjustments
micronutrients = [
    'EX_adocbl(e)', 'EX_vitd2(e)', 'EX_vitd3(e)', 'EX_psyl(e)', 'EX_gum(e)',
    'EX_bglc(e)', 'EX_phyQ(e)', 'EX_fol(e)', 'EX_5mthf(e)', 'EX_q10(e)',
    'EX_retinol_9_cis(e)', 'EX_pydxn(e)', 'EX_pydam(e)', 'EX_pydx(e)',
    'EX_pheme(e)', 'EX_ribflv(e)', 'EX_thm(e)', 'EX_avite1(e)', 'EX_pnto_R(e)',
    'EX_na1(e)', 'EX_cl(e)', 'EX_k(e)', 'EX_pi(e)', 'EX_zn2(e)', 'EX_cu2(e)',
    'EX_btn(e)'
]
micron_clean = [m.replace('(', '_').replace(')', '') for m in micronutrients]

# Rule A: increase small fluxes by 100x
mask_A = df['reaction'].isin(micron_clean) & (df['value'].abs() <= 0.1)
df.loc[mask_A, 'value'] *= -100

# Rule B: ensure EX_pnto_R_e at least 0.1
mask_B = (df['reaction'] == 'EX_pnto_R_e') & (df['value'].abs() < 0.1)
df.loc[mask_B, 'value'] = -0.1

# Rule C: certain micronutrients >= 1
force1 = ['EX_fol_e', 'EX_arab_L_e', 'EX_xyl_D_e', 'EX_amp_e', 'EX_nh4_e', 'EX_cobalt2_e']
mask_C = df['reaction'].isin(force1) & (df['value'].abs() < 1)
df.loc[mask_C, 'value'] = -1

# Tag all micronutrient rows with type 'micronutrient'
df.loc[df['reaction'].isin(micron_clean), 'type'] = 'micronutrient'

def standardize_rxn(rxn: str) -> str:
    # 1) Turn “(e)” → “_e”
    s = rxn.replace('(e)', '_e')
    # 2) Remove any remaining brackets/parentheses
    s = re.sub(r'[\[\]\(\)]', '', s)
    # 3) Match EX_<core>_e
    m = re.match(r'^(EX)_(.+)_e$', s)
    if m:
        prefix, core = m.group(1), m.group(2)
        # double any underscores in the core
        core = core.replace('_', '__')
        return f'{prefix}_{core}_e'
    # fallback: return cleaned string
    return s

# apply as the very last step before saving
df['reaction'] = df['reaction'].apply(standardize_rxn)

df_unique = df.drop_duplicates(subset=['reaction'])

# Write to a new CSV
csv_path = 'exchange_bounds.csv'
df_unique.to_csv(csv_path, index=False)

print(f"Unique-reaction CSV written to {csv_path}")

# Assuming `df` contains your final DataFrame with columns ['reaction','type','value']
# and you’ve already standardized IDs and deduplicated:
exchange_list = df.drop_duplicates(subset=['reaction']).to_dict(orient='records')

# Build only the exchange_bounds part (you’ll add the globals manually)
config = {
    "exchange_bounds": exchange_list
}

# Write to JSON (or print to console)
with open('boundary_conditions.json', 'w') as fp:
    json.dump(config, fp, indent=2)

# Or to just print for manual copy‐paste:
print(json.dumps(config, indent=2))

# In this DataFrame, the `type` column categorizes each exchange reaction according to its dietary or modeling role:
# 
#   • 'starch'             – simple sugars and starch-derived carbohydrates (e.g., glucose, maltose, cellobiose),
#                            representing readily fermentable substrates in a high‐fiber diet.
#   • 'fiber'              – complex, non‐digested polysaccharides and oligosaccharides (e.g., inulin, cellulose,
#                            arabinoxylan), serving as slow‐release fermentable substrates for gut microbes.
#   • 'fat'                – lipid and fatty acid uptake reactions (e.g., glycerol, long‐chain fatty acids,
#                            cholesterol derivatives) reflecting dietary fats.
#   • 'protein'            – amino acid exchange reactions (e.g., EX_ala_L_e, EX_trp_L_e), modeling peptide
#                            digestion products.
#   • 'micronutrient'      – vitamins and trace minerals (e.g., EX_fol_e, EX_zn2_e) that are essential for
#                            cofactor requirements and enzyme function.
#   • 'essential'          – metabolites required by the metabolic models but not supplied by the diet
#                            (default uptake flux added at a minimal rate of -0.1 mmol/gDW/h).
#   • 'unmapped'           – compounds known to exist in the diet formulation but not yet mapped by the
#                            Diet Designer (added with a default flux of -50 mmol/gDW/h).
#   • 'choline'            – a special case nutrient (EX_chol_e) critical for certain pathways,
#                            enforced at -41.251 mmol/gDW/h if missing.
# 
# This classification enables downstream scripts to apply category‐specific constraints,
# defaults, and post‐processing rules before constructing the final exchange constraint set.
