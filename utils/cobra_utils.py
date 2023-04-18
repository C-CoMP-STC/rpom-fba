from cobra.io import read_sbml_model


def reactions_of_metabolite(model, metID):
    return model.reactions.query(lambda r: metID in [met.id for met in r.metabolites],
                                 attribute=None)


def path_to(model, met_from, met_to):
    pass


def diff_models(model1, model2):
    pass


def main():
    model1 = read_sbml_model("clean_models/Rpom_0.xml")
    model2 = read_sbml_model("clean_models/Rpom_02.xml")

    diff_models(model1, model2)


if __name__ == "__main__":
    main()
