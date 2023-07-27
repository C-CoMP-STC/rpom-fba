import json
import os

import cobra
import matplotlib
import numpy as np
import pandas as pd
from cobra.io import read_sbml_model

matplotlib.use("Agg")
import matplotlib.pyplot as plt


MODEL_DIR = "clean_models/"
MODELS = ["Rpom_0.xml", "Rpom_02.xml", "Rpom_03.xml", "Rpom_04.xml", "Rpom_05.xml",
          "Rpom_06.xml", "Rpom_025.xml", "Rpom_035.xml", "Rpom_045.xml", "Rpom_055.xml"]


def main():
    rpom_models = []
    for model_file in MODELS:
        model = read_sbml_model(os.path.join(MODEL_DIR, model_file))
        rpom_models.append(model)
    
    
    for model in rpom_models:
        solution = model.optimize()
        print(f"Objective: {solution.objective_value}")

        print(f"{solution.reduced_costs[solution.reduced_costs != 0]}")

    with open("data/core_biomass.json", "r") as f:
        core_biomass = json.load(f)
        biomass_objective = core_biomass["core_biomass"]

    

if __name__ == "__main__":
    main()
