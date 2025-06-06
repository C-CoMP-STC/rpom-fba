from cobra.core.metabolite import Metabolite

ADDED_METABOLITES = [
    # DHPS Metabolism =============================================================================
    Metabolite(id="CPD-12693[e]",
               formula="C3H7O5S",
               charge=-1,
               name="DHPS",
               compartment="e"),
    Metabolite(id="CPD-12693[c]",
               formula="C3H7O5S",
               charge=-1,
               name="DHPS",
               compartment="c"),
    Metabolite(id="CPD-12694[c]",
               formula="C3H5O5S",
               charge=-1,
               name="1-hydroxy-2-oxo-3-sulfopropane",
               compartment="c"),
    Metabolite(id="CPD-12692[c]",
               formula="C3H7O5S",
               charge=-1,
               name="(R)-2,3-dihydroxypropane 1-sulfonate",
               compartment="c"),
    Metabolite(id="CPD-367[c]",
               formula="C3H4O6S",
               charge=-2,
               name="(2R)-3-sulfolactate",
               compartment="c"),
    Metabolite(id="TRIMETHYLAMINE-N-O[c]",
               formula="C3H9N1O1",
               charge=0,
               name="trimethylamine <i>N</i>-oxide",
               compartment="c"),
    Metabolite(id="SoxZY-L-Cysteine[c]",
               formula="",
               charge=0,
               name="",
               compartment="c"),
    Metabolite(id="Cytochromes-C-Oxidized[c]",
               formula="",
               name="an oxidized c-type cytochrome",
               compartment="c"),
    Metabolite(id="Cytochromes-C-Reduced[c]",
               formula="",
               name="a reduced c-type cytochrome",
               compartment="c"),
    Metabolite(id="SoxY-Thiocysteine-Sulfate[c]",
               formula="",
               name="a [SoxY protein]-S-sulfosulfanyl-L-cysteine",
               compartment="c"),
    Metabolite(id="SoxY-S-Thiocysteine[c]",
               formula="",
               name="a [SoxY protein]-S-sulfanyl-L-cysteine",
               compartment="c"),
    Metabolite(id="SoxY-Thiocysteine-S-Sulfate[c]",
               formula="",
               name="a [SoxY protein]-S-(2-sulfodisulfanyl)-L-cysteine",
               compartment="c"),
    Metabolite(id="SoxY-S-S-Thiocysteine[c]",
               formula="",
               name="a [SoxY protein]-S-disulfanyl-L-cysteine",
               compartment="c"),
    Metabolite(id="Pyruvate-dehydrogenase-lipoate[c]",
               formula="",
               name="a [pyruvate dehydrogenase E2 protein] N6-lipoyl-L-lysine",
               compartment="c"),
    Metabolite(id="Pyruvate-dehydrogenase-acetylDHlipoyl[c]",
               formula="",
               name="a [pyruvate dehydrogenase E2 protein] N6-S-acetyldihydrolipoyl-L-lysine",
               compartment="c"),
    Metabolite(id="Pyruvate-dehydrogenase-dihydrolipoate[c]",
               formula="",
               name="a [pyruvate dehydrogenase E2 protein] N6-dihydrolipoyl-L-lysine",
               compartment="c"),
    Metabolite(id="L-1-PHOSPHATIDYL-ETHANOLAMINE[c]",
               formula="",
               name="a phosphatidylethanolamine",
               compartment="c"),
    # From Transporter data ==================================================================
    Metabolite(id="CPD0-1265[e]",
               formula="C9H14O4",
               charge=-2,
               name="azelate",
               compartment="e"),
    Metabolite(id="CPD0-1265[p]",
               formula="C9H14O4",
               charge=-2,
               name="azelate",
               compartment="p"),
    Metabolite(id="CPD0-1265[c]",
               formula="C9H14O4",
               charge=-2,
               name="azelate",
               compartment="c"),
    Metabolite(id="CPD-335[e]",
               formula="C4H7O3",
               charge=-1,
               name="(R)-3-hydroxybutanoate",
               compartment="e"),
    Metabolite(id="CPD-335[p]",
               formula="C4H7O3",
               charge=-1,
               name="(R)-3-hydroxybutanoate",
               compartment="p"),
    Metabolite(id="CPD-335[c]",
               formula="C4H7O3",
               charge=-1,
               name="(R)-3-hydroxybutanoate",
               compartment="c"),
    Metabolite(id="CARNITINE[e]",
               formula="C7H15NO3",
               charge=0,
               name="carnitine",
               compartment="e"),
    Metabolite(id="CARNITINE[p]",
               formula="C7H15NO3",
               charge=0,
               name="carnitine",
               compartment="p"),
    Metabolite(id="CARNITINE[c]",
               formula="C7H15NO3",
               charge=0,
               name="carnitine",
               compartment="c"),
    Metabolite(id="CHOLINE[e]",
               formula="C5H14NO",
               charge=1,
               name="choline",
               compartment="e"),
    Metabolite(id="CHOLINE[p]",
               formula="C5H14NO",
               charge=1,
               name="choline",
               compartment="p"),
    Metabolite(id="CHOLINE[c]",
               formula="C5H14NO",
               charge=1,
               name="choline",
               compartment="c"),
    Metabolite(id="CIT[e]",
               formula="C6H5O7",
               charge=-3,
               name="citrate",
               compartment="e"),
    Metabolite(id="CIT[p]",
               formula="C6H5O7",
               charge=-3,
               name="citrate",
               compartment="p"),
    Metabolite(id="L-CYSTEATE[e]",
               formula="C3H6NO5S",
               charge=-1,
               name="L-cysteate",
               compartment="e"),
    Metabolite(id="L-CYSTEATE[p]",
               formula="C3H6NO5S",
               charge=-1,
               name="L-cysteate",
               compartment="p"),
    Metabolite(id="CPD-12693[p]",
               formula="C3H7O5S",
               charge=-1,
               name="DHPS",
               compartment="p"),
    Metabolite(id="FUM[e]",
               formula="C4H2O4",
               charge=-2,
               name="fumarate",
               compartment="e"),
    Metabolite(id="MAL[e]",
               formula="C4H4O5",
               charge=-2,
               name='malate',
               compartment="e"),
    Metabolite(id="GLC[e]",
               formula="C6H12O6",
               charge=0,
               name='glucose',
               compartment="e"),
    Metabolite(id="GLC[p]",
               formula="C6H12O6",
               charge=0,
               name='glucose',
               compartment="p"),
    Metabolite(id="BETA-D-XYLOSE[e]",
               formula="C5H10O5",
               charge=0,
               name="β-D-xylose",
               compartment="e"),
    Metabolite(id="BETA-D-XYLOSE[p]",
               formula="C5H10O5",
               charge=0,
               name="β-D-xylose",
               compartment="p"),
    Metabolite(id="GLYCEROL-3P[e]",
               formula="C3H7O6P",
               charge=-2,
               name="glycerol 3-phosphate",
               compartment="e"),
    Metabolite(id="N-ACETYL-D-GLUCOSAMINE[e]",
               formula="C8H15NO6",
               charge=0,
               name="N-acetyl-D-glucosamine",
               compartment="e"),
    Metabolite(id="N-ACETYL-D-GLUCOSAMINE[p]",
               formula="C8H15NO6",
               charge=0,
               name="N-acetyl-D-glucosamine",
               compartment="p"),
    Metabolite(id="N-ACETYL-D-GLUCOSAMINE[c]",
               formula="C8H15NO6",
               charge=0,
               name="N-acetyl-D-glucosamine",
               compartment="c"),
    Metabolite(id="CADAVERINE[e]",
               formula="C5H16N2",
               charge=2,
               name="cadaverine",
               compartment="e"),
    Metabolite(id="CADAVERINE[p]",
               formula="C5H16N2",
               charge=2,
               name="cadaverine",
               compartment="p"),
    Metabolite(id="PUTRESCINE[e]",
               formula="C4H14N2",
               charge=2,
               name="putrescine",
               compartment="e"),
    Metabolite(id="PUTRESCINE[p]",
               formula="C4H14N2",
               charge=2,
               name="putrescine",
               compartment="p"),
    Metabolite(id="PUTRESCINE[c]",
               formula="C4H14N2",
               charge=2,
               name="putrescine",
               compartment="c"),
    Metabolite(id="TAURINE[e]",
               formula="C2H7NO3S",
               charge=0,
               name="taurine",
               compartment="e"),
    Metabolite(id="TAURINE[p]",
               formula="C2H7NO3S",
               charge=0,
               name="taurine",
               compartment="p"),
    Metabolite(id="THYMIDINE[e]",
               formula="C10H14N2O5",
               charge=0,
               name="thymidine",
               compartment="e"),
    Metabolite(id="THYMIDINE[p]",
               formula="C10H14N2O5",
               charge=0,
               name="thymidine",
               compartment="p"),
    Metabolite(id="TRIMETHYLAMINE-N-O[e]",
               formula="C3H9NO",
               charge=0,
               name="trimethylamine <i>N</i>-oxide",
               compartment="e"),
    Metabolite(id="TRIMETHYLAMINE-N-O[p]",
               formula="C3H9NO",
               charge=0,
               name="trimethylamine <i>N</i>-oxide",
               compartment="p"),
    Metabolite(id="GLN[e]",
               formula = "C5H10N2O3",
               charge=0,
               name="L-glutamine",
               compartment="e"),
    Metabolite(id="L-ASPARTATE[e]",
               formula="C4H6NO4",
               charge=-1,
               name="L-aspartate",
               compartment="e"),
    Metabolite(id="ASN[e]",
               formula="C4H8N2O3",
               charge=0,
               name="L-asparagine",
               compartment="e"),
    Metabolite(id="VAL[e]",
               formula="C5H11NO2",
               charge=0,
               name="L-valine",
               compartment="e")
]
