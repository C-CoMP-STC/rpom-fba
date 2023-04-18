from cobra.io import read_sbml_model
from memote.support.consistency import find_blocked_metabolites


def main():
    model = read_sbml_model("clean_models/Rpom_05.xml")
    solution = model.optimize()
    print (solution)


if __name__ == "__main__":
    main()
