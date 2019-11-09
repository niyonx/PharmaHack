import rdkit 
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd 
import numpy as numpy


file_function = '/Users/Alexis/Desktop/functions.tsv'
function_list = open(file_function, 'r').readlines()

for i in range(1,4):

  #create a DF for compound set
  file = '/Users/Alexis/Downloads/data/compound_set'+ str(i) +'.sdf'
  suppl = Chem.SDMolSupplier(file)
  
  '''
  list_smiles = []
  for mol in suppl:
    if mol is None: continue
    list_smiles.append(Chem.MolToSmiles(mol))

  data_frame = pd.DataFrame(list_smiles, columns=['smiles'])

  #iterate through the rdkit functions 
  for function in function_list:
    # create a csv to receive the inputs
    function_name = function.split('\n')[0]
    dict_mol = {}
'''
  descriptor_value_BalabanJ = []
  descriptor_value_BertzCT = []
  descriptor_value_ExactMolWt = []
  descriptor_value_HeavyAtomCount = [] 
  descriptor_value_HeavyAtomMolWt = []
  descriptor_value_MolLogP = []
  descriptor_value_MolMR = []

  descriptor_value_NHOHCount = []
  descriptor_value_NOCount =[]
  descriptor_value_NumAliphaticCarbocycles = []
  descriptor_value_NumAliphaticHeterocycles = []
  descriptor_value_NumAliphaticRings = []
  descriptor_value_NumAromaticCarbocycles = []
  descriptor_value_NumAromaticHeterocycles = []
  descriptor_value_NumAromaticRings = []
  descriptor_value_NumHAcceptors = []
  descriptor_value_NumHDonors = []
  descriptor_value_NumHeteroatoms = []
  descriptor_value_NumRadicalElectrons = []
  descriptor_value_NumRotatableBonds = []
  descriptor_value_NumSaturatedCarbocycles = []
  descriptor_value_NumSaturatedHeterocycles = []
  descriptor_value_NumSaturatedRings = []
  descriptor_value_NumValenceElectrons = []
  descriptor_value_RingCount = []

  descriptor_value_Asphericity = []
  descriptor_value_Eccentricity = []
  descriptor_value_InertialShapeFactor = []
  descriptor_value_SpherocityIndex = []

  descriptor_value_fr_Al_COO = []
  descriptor_value_fr_Al_OH = []
  descriptor_value_fr_Al_OH_noTert = []
  descriptor_value_fr_ArN = []
  descriptor_value_fr_Ar_COO = []
  descriptor_value_fr_Ar_N = []
  descriptor_value_fr_Ar_NH = []
  descriptor_value_fr_Ar_OH = []
  descriptor_value_fr_COO = []
  descriptor_value_fr_COO2 = []
  descriptor_value_fr_C_O = []
  descriptor_value_fr_C_O_noCOO = []
  descriptor_value_fr_C_S = []
  descriptor_value_fr_HOCCN = []
  descriptor_value_fr_Imine = []
  descriptor_value_fr_NH0 = []
  descriptor_value_fr_NH1 =[]
  descriptor_value_fr_NH2 = []
  descriptor_value_fr_N_O = []
  descriptor_value_fr_Ndealkylation1 = []
  descriptor_value_fr_Ndealkylation2 = []
  descriptor_value_fr_Nhpyrrole = []
  descriptor_value_fr_SH = []
  descriptor_value_fr_aldehyde = []
  descriptor_value_fr_alkyl_carbamate = []
  descriptor_value_fr_alkyl_halide = []
  descriptor_value_fr_allylic_oxid = []
  descriptor_value_fr_amide = []
  descriptor_value_fr_amidine = []
  descriptor_value_fr_aniline = []
  descriptor_value_fr_aryl_methy = []
  descriptor_value_fr_azide = []
  descriptor_value_fr_azo = []
  descriptor_value_fr_barbitur = []
  descriptor_value_fr_benzene = []
  descriptor_value_fr_benzodiazepine = []
  descriptor_value_fr_bicyclic = []
  descriptor_value_fr_diazo = []
  descriptor_value_fr_dihydropyridine =[]
  descriptor_value_fr_epoxide = []
  descriptor_value_fr_ester = []
  descriptor_value_fr_ether = []
  descriptor_value_fr_furan = []
  descriptor_value_fr_guanido = []
  descriptor_value_fr_halogen = []
  descriptor_value_fr_hdrzine = []
  descriptor_value_fr_hdrzone = []
  descriptor_value_fr_imidazole = []
  descriptor_value_fr_imide = []
  descriptor_value_fr_isocyan = []
  descriptor_value_fr_isothiocyan = []
  descriptor_value_fr_ketone = []
  descriptor_value_fr_lactam =[]
  descriptor_value_fr_lactone = []
  descriptor_value_fr_methoxy =[]
  descriptor_value_fr_morpholine = []
  descriptor_value_fr_nitrile = []
  descriptor_value_fr_nitro = []
  descriptor_value_fr_nitro_arom = []
  descriptor_value_fr_nitro_arom_nonortho = []
  descriptor_value_fr_nitroso = []
  descriptor_value_fr_oxazole =[]
  descriptor_value_fr_oxime = []
  descriptor_value_fr_para_hydroxylation = []
  descriptor_value_fr_phenol = []
  descriptor_value_fr_phenol_noOrthoHbond = []
  descriptor_value_fr_phos_acid = []
  descriptor_value_fr_phos_ester =[]
  descriptor_value_fr_piperdine = []
  descriptor_value_fr_piperzine = []
  descriptor_value_fr_priamide = []
  descriptor_value_fr_prisulfonamd = []
  descriptor_value_fr_pyridine =[]
  descriptor_value_fr_sulfide = []
  descriptor_value_fr_sulfonamd = []
  descriptor_value_fr_sulfone = []
  descriptor_value_fr_term_acetylene = []
  descriptor_value_fr_tetrazole = []
  descriptor_value_fr_thiazole = []
  descriptor_value_fr_thiocyan = []
  descriptor_value_fr_thiophene = []
  descriptor_value_fr_unbrch_alkane = []
  descriptor_value_fr_urea = []

  molecule_smile = []

  for mol in suppl:
    if mol is None: continue
    molecule_smile.append(Chem.MolToSmiles(mol))

    #bunch of curated descriptors based on diversity of variables
    #probed by the index and the information they give

    #poor descriptors
    descriptor_value_BalabanJ.append(Descriptors.BalabanJ(mol))
    descriptor_value_BertzCT.append(Descriptors.BertzCT(mol))
    descriptor_value_ExactMolWt.append(Descriptors.ExactMolWt(mol))
    descriptor_value_HeavyAtomCount.append(Descriptors.HeavyAtomCount(mol))
    descriptor_value_HeavyAtomMolWt.append(Descriptors.HeavyAtomMolWt(mol))
    descriptor_value_MolLogP.append(Descriptors.MolLogP(mol))
    descriptor_value_MolMR.append(Descriptors.MolMR(mol))

    #medium descriptors
    descriptor_value_NHOHCount.append(Descriptors.NHOHCount(mol))
    descriptor_value_NOCount.append(Descriptors.NOCount(mol))
    descriptor_value_NumAliphaticCarbocycles.append(Descriptors.NumAliphaticCarbocycles(mol))
    descriptor_value_NumAliphaticHeterocycles.append(Descriptors.NumAliphaticHeterocycles(mol))
    descriptor_value_NumAliphaticRings.append(Descriptors.NumAliphaticRings(mol))
    descriptor_value_NumAromaticCarbocycles.append(Descriptors.NumAromaticCarbocycles(mol))
    descriptor_value_NumAromaticHeterocycles.append(Descriptors.NumAromaticHeterocycles(mol))
    descriptor_value_NumAromaticRings.append(Descriptors.NumAromaticRings(mol))
    descriptor_value_NumHAcceptors.append(Descriptors.NumHAcceptors(mol))
    descriptor_value_NumHDonors.append(Descriptors.NumHDonors(mol))
    descriptor_value_NumHeteroatoms.append(Descriptors.NumHeteroatoms(mol))
    descriptor_value_NumRadicalElectrons.append(Descriptors.NumRadicalElectrons(mol))
    descriptor_value_NumRotatableBonds.append(Descriptors.NumRotatableBonds(mol))
    descriptor_value_NumSaturatedCarbocycles.append(Descriptors.NumSaturatedCarbocycles(mol))
    descriptor_value_NumSaturatedHeterocycles.append(Descriptors.NumSaturatedHeterocycles(mol))
    descriptor_value_NumSaturatedRings.append(Descriptors.NumSaturatedRings(mol))
    descriptor_value_NumValenceElectrons.append(Descriptors.NumValenceElectrons(mol))
    descriptor_value_RingCount.append(Descriptors.RingCount(mol))

    #subcategory of medium descriptors focussed on volume distribution
    #for more http://rdkit.org/docs/source/rdkit.Chem.Descriptors3D.html

    descriptor_value_Asphericity.append(Descriptors.Chem.Descriptors3D.Asphericity(mol))
    descriptor_value_Eccentricity.append(Descriptors.Chem.Descriptors3D.Eccentricity(mol))
    descriptor_value_InertialShapeFactor.append(Descriptors.Chem.Descriptors3D.InertialShapeFactor(mol))
    descriptor_value_SpherocityIndex.append(Descriptors.Chem.Descriptors3D.SpherocityIndex(mol))

    #strong descriptors
    descriptor_value_fr_Al_COO.append(Descriptors.fr_Al_COO(mol))
    descriptor_value_fr_Al_OH.append(Descriptors.fr_Al_OH(mol))
    descriptor_value_fr_Al_OH_noTert.append(Descriptors.fr_Al_OH_noTert(mol))
    descriptor_value_fr_ArN.append(Descriptors.fr_ArN(mol))
    descriptor_value_fr_Ar_COO.append(Descriptors.fr_Ar_COO(mol))
    descriptor_value_fr_Ar_N.append(Descriptors.fr_Ar_N(mol))
    descriptor_value_fr_Ar_NH.append(Descriptors.fr_Ar_NH(mol))
    descriptor_value_fr_Ar_OH.append(Descriptors.fr_Ar_OH(mol))
    descriptor_value_fr_COO.append(Descriptors.fr_COO(mol))
    descriptor_value_fr_COO2.append(Descriptors.fr_COO2(mol))
    descriptor_value_fr_C_O.append(Descriptors.fr_C_O(mol))
    descriptor_value_fr_C_O_noCOO.append(Descriptors.fr_C_O_noCOO(mol))
    descriptor_value_fr_C_S.append(Descriptors.fr_C_S(mol))
    descriptor_value_fr_HOCCN.append(Descriptors.fr_HOCCN(mol))
    descriptor_value_fr_Imine.append(Descriptors.fr_Imine(mol))
    descriptor_value_fr_NH0.append(Descriptors.fr_NH0(mol))
    descriptor_value_fr_NH1.append(Descriptors.fr_NH1(mol))
    descriptor_value_fr_NH2.append(Descriptors.fr_NH2(mol))
    descriptor_value_fr_N_O.append(Descriptors.fr_N_O(mol))
    descriptor_value_fr_Ndealkylation1.append(Descriptors.fr_Ndealkylation1(mol))
    descriptor_value_fr_Ndealkylation2.append(Descriptors.fr_Ndealkylation2(mol))
    descriptor_value_fr_Nhpyrrole.append(Descriptors.fr_Nhpyrrole(mol))
    descriptor_value_fr_SH.append(Descriptors.fr_SH(mol))
    descriptor_value_fr_aldehyde.append(Descriptors.fr_aldehyde(mol))
    descriptor_value_fr_alkyl_carbamate.append(Descriptors.fr_alkyl_carbamate(mol))
    descriptor_value_fr_alkyl_halide.append(Descriptors.fr_alkyl_halide(mol))
    descriptor_value_fr_allylic_oxid.append(Descriptors.fr_allylic_oxid(mol))
    descriptor_value_fr_amide.append(Descriptors.fr_amide(mol))
    descriptor_value_fr_amidine.append(Descriptors.fr_amidine(mol))
    descriptor_value_fr_aniline.append(Descriptors.fr_aniline(mol))
    descriptor_value_fr_aryl_methy.append(Descriptors.fr_aryl_methyl(mol))
    descriptor_value_fr_azide.append(Descriptors.fr_azide(mol))
    descriptor_value_fr_azo.append(Descriptors.fr_azo(mol))
    descriptor_value_fr_barbitur.append(Descriptors.fr_barbitur(mol))
    descriptor_value_fr_benzene.append(Descriptors.fr_benzene(mol))
    descriptor_value_fr_benzodiazepine.append(Descriptors.fr_benzodiazepine(mol))
    descriptor_value_fr_bicyclic.append(Descriptors.fr_bicyclic(mol))
    descriptor_value_fr_diazo.append(Descriptors.fr_diazo(mol))
    descriptor_value_fr_dihydropyridine.append(Descriptors.fr_dihydropyridine(mol))
    descriptor_value_fr_epoxide.append(Descriptors.fr_epoxide(mol))
    descriptor_value_fr_ester.append(Descriptors.fr_ester(mol))
    descriptor_value_fr_ether.append(Descriptors.fr_ether(mol))
    descriptor_value_fr_furan.append(Descriptors.fr_furan(mol))
    descriptor_value_fr_guanido.append(Descriptors.fr_guanido(mol))
    descriptor_value_fr_halogen.append(Descriptors.fr_halogen(mol))
    descriptor_value_fr_hdrzine.append(Descriptors.fr_hdrzine(mol))
    descriptor_value_fr_hdrzone.append(Descriptors.fr_hdrzone(mol))
    descriptor_value_fr_imidazole.append(Descriptors.fr_imidazole(mol))
    descriptor_value_fr_imide.append(Descriptors.fr_imide(mol))
    descriptor_value_fr_isocyan.append(Descriptors.fr_isocyan(mol))
    descriptor_value_fr_isothiocyan.append(Descriptors.fr_isothiocyan(mol))
    descriptor_value_fr_ketone.append(Descriptors.fr_ketone(mol))
    descriptor_value_fr_lactam.append(Descriptors.fr_lactam(mol))
    descriptor_value_fr_lactone.append(Descriptors.fr_lactone(mol))
    descriptor_value_fr_methoxy.append(Descriptors.fr_methoxy(mol))
    descriptor_value_fr_morpholine.append(Descriptors.fr_morpholine(mol))
    descriptor_value_fr_nitrile.append(Descriptors.fr_nitrile(mol))
    descriptor_value_fr_nitro.append(Descriptors.fr_nitro(mol))
    descriptor_value_fr_nitro_arom.append(Descriptors.fr_nitro_arom(mol))
    descriptor_value_fr_nitro_arom_nonortho.append(Descriptors.fr_nitro_arom_nonortho(mol))
    descriptor_value_fr_nitroso.append(Descriptors.fr_nitroso(mol))
    descriptor_value_fr_oxazole.append(Descriptors.fr_oxazole(mol))
    descriptor_value_fr_oxime.append(Descriptors.fr_oxime(mol))
    descriptor_value_fr_para_hydroxylation.append(Descriptors.fr_para_hydroxylation(mol))
    descriptor_value_fr_phenol.append(Descriptors.fr_phenol(mol))
    descriptor_value_fr_phenol_noOrthoHbond.append(Descriptors.fr_phenol_noOrthoHbond(mol))
    descriptor_value_fr_phos_acid.append(Descriptors.fr_phos_acid(mol))
    descriptor_value_fr_phos_ester.append(Descriptors.fr_phos_ester(mol))
    descriptor_value_fr_piperdine.append(Descriptors.fr_piperdine(mol))
    descriptor_value_fr_piperzine.append(Descriptors.fr_piperzine(mol))
    descriptor_value_fr_priamide.append(Descriptors.fr_priamide(mol))
    descriptor_value_fr_prisulfonamd.append(Descriptors.fr_prisulfonamd(mol))
    descriptor_value_fr_pyridine.append(Descriptors.fr_pyridine(mol))
    descriptor_value_fr_sulfide.append(Descriptors.fr_sulfide(mol))
    descriptor_value_fr_sulfonamd.append(Descriptors.fr_sulfonamd(mol))
    descriptor_value_fr_sulfone.append(Descriptors.fr_sulfone(mol))
    descriptor_value_fr_term_acetylene.append(Descriptors.fr_term_acetylene(mol))
    descriptor_value_fr_tetrazole.append(Descriptors.fr_tetrazole(mol))
    descriptor_value_fr_thiazole.append(Descriptors.fr_thiazole(mol))
    descriptor_value_fr_thiocyan.append(Descriptors.fr_thiocyan(mol))
    descriptor_value_fr_thiophene.append(Descriptors.fr_thiophene(mol))
    descriptor_value_fr_unbrch_alkane.append(Descriptors.fr_unbrch_alkane(mol))
    descriptor_value_fr_urea.append(Descriptors.fr_urea(mol))



    
    trans_dict = pd.DataFrame(list(zip(molecule_smile, descriptor_value_BertzCT, descriptor_value_ExactMolWt, descriptor_value_HeavyAtomCount, descriptor_value_HeavyAtomMolWt, descriptor_value_MolLogP, descriptor_value_MolMR, descriptor_value_NHOHCount, descriptor_value_NOCount, descriptor_value_NumAliphaticCarbocycles, descriptor_value_NumAliphaticHeterocycles, descriptor_value_NumAliphaticRings, descriptor_value_NumAromaticCarbocycles, descriptor_value_NumAromaticHeterocycles, descriptor_value_NumAromaticRings, descriptor_value_NumHAcceptors, descriptor_value_NumHDonors, descriptor_value_NumHeteroatoms, descriptor_value_NumRadicalElectrons, descriptor_value_NumRotatableBonds, descriptor_value_NumSaturatedCarbocycles, descriptor_value_NumSaturatedHeterocycles, descriptor_value_NumSaturatedRings, descriptor_value_NumValenceElectrons, descriptor_value_RingCount, descriptor_value_Asphericity, descriptor_value_Eccentricity, descriptor_value_InertialShapeFactor, descriptor_value_SpherocityIndex, descriptor_value_fr_Al_COO, descriptor_value_fr_Al_OH, descriptor_value_fr_Al_OH_noTert, descriptor_value_fr_ArN, descriptor_value_fr_Ar_COO, descriptor_value_fr_Ar_N, descriptor_value_fr_Ar_NH, descriptor_value_fr_Ar_OH, descriptor_value_fr_COO, descriptor_value_fr_COO2, descriptor_value_fr_C_O, descriptor_value_fr_C_O_noCOO, descriptor_value_fr_C_S, descriptor_value_fr_HOCCN, descriptor_value_fr_Imine, descriptor_value_fr_NH0, descriptor_value_fr_NH1, descriptor_value_fr_NH2, descriptor_value_fr_N_O, descriptor_value_fr_Ndealkylation1, descriptor_value_fr_Ndealkylation2, descriptor_value_fr_Nhpyrrole, descriptor_value_fr_SH, descriptor_value_fr_aldehyde, descriptor_value_fr_alkyl_carbamate, descriptor_value_fr_alkyl_halide, descriptor_value_fr_allylic_oxid, descriptor_value_fr_amide, descriptor_value_fr_amidine, descriptor_value_fr_aniline, descriptor_value_fr_aryl_methy, descriptor_value_fr_azide, descriptor_value_fr_azo, descriptor_value_fr_barbitur, descriptor_value_fr_benzene, descriptor_value_fr_benzodiazepine, descriptor_value_fr_bicyclic, descriptor_value_fr_diazo, descriptor_value_fr_dihydropyridine, descriptor_value_fr_epoxide, descriptor_value_fr_ester, descriptor_value_fr_ether, descriptor_value_fr_furan, descriptor_value_fr_guanido, descriptor_value_fr_halogen, descriptor_value_fr_hdrzine, descriptor_value_fr_hdrzone, descriptor_value_fr_imidazole, descriptor_value_fr_imide, descriptor_value_fr_isocyan, descriptor_value_fr_isothiocyan, descriptor_value_fr_ketone, descriptor_value_fr_lactam, descriptor_value_fr_lactone, descriptor_value_fr_methoxy, descriptor_value_fr_morpholine, descriptor_value_fr_nitrile, descriptor_value_fr_nitro, descriptor_value_fr_nitro_arom, descriptor_value_fr_nitro_arom_nonortho, descriptor_value_fr_nitroso, descriptor_value_fr_oxazole, descriptor_value_fr_oxime, descriptor_value_fr_para_hydroxylation, descriptor_value_fr_phenol, descriptor_value_fr_phenol_noOrthoHbond, descriptor_value_fr_phos_acid, descriptor_value_fr_phos_ester, descriptor_value_fr_piperdine, descriptor_value_fr_piperzine, descriptor_value_fr_priamide, descriptor_value_fr_prisulfonamd, descriptor_value_fr_pyridine, descriptor_value_fr_sulfide, descriptor_value_fr_sulfonamd, descriptor_value_fr_sulfone, descriptor_value_fr_term_acetylene, descriptor_value_fr_tetrazole, descriptor_value_fr_thiazole, descriptor_value_fr_thiocyan, descriptor_value_fr_thiophene, descriptor_value_fr_unbrch_alkane, descriptor_value_fr_urea)), columns = ['molecule_smile', 'descriptor_value_BertzCT', 'descriptor_value_ExactMolWt', 'descriptor_value_HeavyAtomCount', 'descriptor_value_HeavyAtomMolWt', 'descriptor_value_MolLogP', 'descriptor_value_MolMR', 'descriptor_value_NHOHCount', 'descriptor_value_NOCount', 'descriptor_value_NumAliphaticCarbocycles', 'descriptor_value_NumAliphaticHeterocycles', 'descriptor_value_NumAliphaticRings', 'descriptor_value_NumAromaticCarbocycles', 'descriptor_value_NumAromaticHeterocycles', 'descriptor_value_NumAromaticRings', 'descriptor_value_NumHAcceptors', 'descriptor_value_NumHDonors', 'descriptor_value_NumHeteroatoms', 'descriptor_value_NumRadicalElectrons', 'descriptor_value_NumRotatableBonds', 'descriptor_value_NumSaturatedCarbocycles', 'descriptor_value_NumSaturatedHeterocycles', 'descriptor_value_NumSaturatedRings', 'descriptor_value_NumValenceElectrons', 'descriptor_value_RingCount', 'descriptor_value_Asphericity', 'descriptor_value_Eccentricity', 'descriptor_value_InertialShapeFactor', 'descriptor_value_SpherocityIndex', 'descriptor_value_fr_Al_COO', 'descriptor_value_fr_Al_OH', 'descriptor_value_fr_Al_OH_noTert', 'descriptor_value_fr_ArN', 'descriptor_value_fr_Ar_COO', 'descriptor_value_fr_Ar_N', 'descriptor_value_fr_Ar_NH', 'descriptor_value_fr_Ar_OH', 'descriptor_value_fr_COO', 'descriptor_value_fr_COO2', 'descriptor_value_fr_C_O', 'descriptor_value_fr_C_O_noCOO', 'descriptor_value_fr_C_S', 'descriptor_value_fr_HOCCN', 'descriptor_value_fr_Imine', 'descriptor_value_fr_NH0', 'descriptor_value_fr_NH1', 'descriptor_value_fr_NH2', 'descriptor_value_fr_N_O', 'descriptor_value_fr_Ndealkylation1', 'descriptor_value_fr_Ndealkylation2', 'descriptor_value_fr_Nhpyrrole', 'descriptor_value_fr_SH', 'descriptor_value_fr_aldehyde', 'descriptor_value_fr_alkyl_carbamate', 'descriptor_value_fr_alkyl_halide', 'descriptor_value_fr_allylic_oxid', 'descriptor_value_fr_amide', 'descriptor_value_fr_amidine', 'descriptor_value_fr_aniline', 'descriptor_value_fr_aryl_methy', 'descriptor_value_fr_azide', 'descriptor_value_fr_azo', 'descriptor_value_fr_barbitur', 'descriptor_value_fr_benzene', 'descriptor_value_fr_benzodiazepine', 'descriptor_value_fr_bicyclic', 'descriptor_value_fr_diazo', 'descriptor_value_fr_dihydropyridine', 'descriptor_value_fr_epoxide', 'descriptor_value_fr_ester', 'descriptor_value_fr_ether', 'descriptor_value_fr_furan', 'descriptor_value_fr_guanido', 'descriptor_value_fr_halogen', 'descriptor_value_fr_hdrzine', 'descriptor_value_fr_hdrzone', 'descriptor_value_fr_imidazole', 'descriptor_value_fr_imide', 'descriptor_value_fr_isocyan', 'descriptor_value_fr_isothiocyan', 'descriptor_value_fr_ketone', 'descriptor_value_fr_lactam', 'descriptor_value_fr_lactone', 'descriptor_value_fr_methoxy', 'descriptor_value_fr_morpholine', 'descriptor_value_fr_nitrile', 'descriptor_value_fr_nitro', 'descriptor_value_fr_nitro_arom', 'descriptor_value_fr_nitro_arom_nonortho', 'descriptor_value_fr_nitroso', 'descriptor_value_fr_oxazole', 'descriptor_value_fr_oxime', 'descriptor_value_fr_para_hydroxylation', 'descriptor_value_fr_phenol', 'descriptor_value_fr_phenol_noOrthoHbond', 'descriptor_value_fr_phos_acid', 'descriptor_value_fr_phos_ester', 'descriptor_value_fr_piperdine', 'descriptor_value_fr_piperzine', 'descriptor_value_fr_priamide', 'descriptor_value_fr_prisulfonamd', 'descriptor_value_fr_pyridine', 'descriptor_value_fr_sulfide', 'descriptor_value_fr_sulfonamd', 'descriptor_value_fr_sulfone', 'descriptor_value_fr_term_acetylene', 'descriptor_value_fr_tetrazole', 'descriptor_value_fr_thiazole', 'descriptor_value_fr_thiocyan', 'descriptor_value_fr_thiophene', 'descriptor_value_fr_unbrch_alkane', 'descriptor_value_fr_urea'])
  #save the molecular set descriptors in a csv
  trans_dict.to_csv('./compound_set_'+ str(i) + '.csv' , sep=';')

