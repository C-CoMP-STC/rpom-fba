from itertools import chain
import requests
import xmltodict
import urllib
import pandas as pd

WEBSVC = "https://websvc.biocyc.org/"


def get_session():
    s = requests.Session()
    s.post("https://websvc.biocyc.org/credentials/login/",
           data={"email": input("> Email: "),
                 "password": input("> Password: ")})
    return s


def pathways_of_organism(orgid, session=None):
    if session is None:
        session = get_session()

    query = urllib.parse.quote(f"[p:p<-{orgid}^^pathways]")
    r = session.get(f"{WEBSVC}xmlquery?{query}&detail=none")
    result = xmltodict.parse(r.content)
    result = result["ptools-xml"]["Pathway"]
    result = [pathway["@frameid"] for pathway in result]
    return result


def reactions_of_pathway(pathway, orgid, session=None):
    if session is None:
        session = get_session()

    query = urllib.parse.quote(f"[x:x<-{orgid}~{pathway}^reaction-list]")
    r = session.get(f"{WEBSVC}xmlquery?{query}&detail=full")
    result = xmltodict.parse(r.content)
    result = result["ptools-xml"]["Reaction"]

    def wrap_if_dict(r):
        return [r] if isinstance(r, dict) else r

    result = wrap_if_dict(result)

    def get_sign_left(direction):
        if direction in {"REVERSIBLE", "PHYSIOL-LEFT-TO-RIGHT", "IRREVERSIBLE-LEFT-TO-RIGHT", "LEFT-TO-RIGHT"}:
            return -1
        else:
            return 1

    df = pd.DataFrame([
        {"id": rxn["@frameid"],
         "reaction-direction": rxn["reaction-direction"],
         "metabolites": str({
             (cpd["Compound"]["@frameid"], cpd["compartment"]["cco"]["@frameid"] if "compartment" in cpd else "CCO-IN"):
             sign * int(cpd.get("coefficient", {"#text": 1})["#text"])

             for sign, cpd in chain([(get_sign_left(rxn["reaction-direction"]), cpd)
                                    for cpd in wrap_if_dict(rxn["left"])],
                                    [(-get_sign_left(rxn["reaction-direction"]), cpd)
                                    for cpd in wrap_if_dict(rxn["right"])])
         })
         }
        for rxn in result
    ])

    return df


def main():
    # print(pathways_of_organism("GCF_000011965"))
    # print(reactions_of_pathway("PWY-7980", "GCF_000011965"))

    # Entner-Doudoroff
    print(reactions_of_pathway("PWY-8004", "GCF_000011965"))


if __name__ == "__main__":
    main()
