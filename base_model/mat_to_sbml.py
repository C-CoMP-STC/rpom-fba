import os
from cobra.io import load_matlab_model, write_sbml_model

def main():
    models = ["base_model/Rpom_0.mat", "base_model/Rpom_02.mat", "base_model/Rpom_03.mat", "base_model/Rpom_04.mat", "base_model/Rpom_05.mat",
            "base_model/Rpom_06.mat", "base_model/Rpom_025.mat", "base_model/Rpom_035.mat", "base_model/Rpom_045.mat", "base_model/Rpom_055.mat"]

    for model in models:
        write_sbml_model(load_matlab_model(model), f"{os.path.splitext(model)[0]}.xml")


if __name__=="__main__":
    main()
