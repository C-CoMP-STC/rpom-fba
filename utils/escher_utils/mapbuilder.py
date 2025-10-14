import json
from datetime import datetime

import cobra
from escher.validate import validate_map


class MapBuilder:
    def __init__(self,
                 name="new_map",
                 id=None,
                 description="",
                 template=None):

        # Load starting template
        if template is None:
            template = "utils/escher_utils/template.json"

        with open(template, 'r') as f:
            self.map = json.load(f)

        if not self.validate():
            raise ValueError(
                f"Template file {template} did not contain a valid Escher map.")

        # Set metadata accordingly
        if id is None:
            id = f'{name}__{datetime.now().strftime("%S_%M_%H__%d-%m-%Y")}'

        self.map[0]["map_name"] = name
        self.map[0]["map_id"] = id
        self.map[0]["map_description"] = description

        # id for adding reactions
        self.next_reaction_id = 0

    def add_reaction(self, reaction):
        new_reaction = self.serialize_reaction(reaction)

        new_reaction.update({
            # Visual attributes
            "label_x": self.next_reaction_id * 10,
            "label_y": self.next_reaction_id * 10,
            "segments": {}
        })

        self.map[1]["reactions"][str(self.next_reaction_id)] = new_reaction
        self.next_reaction_id += 1

    def serialize_reaction(self, reaction):
        return {
            # Biological attributes
            "name": reaction.name,  # reaction.notes["bigg.reaction"]?
            "bigg_id": reaction.id,
            "reversibility": reaction.reversibility,
            "gene_reaction_rule": reaction.gene_reaction_rule,
            "genes": [self.serialize_gene(g) for g in reaction.genes],
            "metabolites": [self.serialize_metabolite_coeff_pair(m) for m in reaction.metabolites.items()],
        }

    def serialize_gene(self, gene):
        return {
            "name": gene.name,
            "bigg_id": gene.id
        }

    def serialize_metabolite_coeff_pair(self, metabolite_coeff_pair):
        metabolite, coefficient = metabolite_coeff_pair
        return {
            "coefficient": coefficient,
            # metabolite.notes["bigg.metabolite"]?
            "bigg_id": metabolite.id
        }

    def validate(self, map=None, verbose=True):
        if map is None:
            map = self.map

        try:
            validate_map(map)
            return True
        except Exception as e:
            if verbose:
                print(e)
            return False

    def write(self, filepath):
        with open(filepath, "w") as f:
            f.write(json.dumps(self.map))


def test_map_builder():
    map = MapBuilder()

    assert map.map is not None
    assert map.validate()
    assert not map.validate(map=[{}], verbose=False)

    ecoli = cobra.io.load_model("iJO1366")
    map.add_reaction(ecoli.reactions[0])

    assert map.validate()


def main():
    test_map_builder()


if __name__ == "__main__":
    main()
