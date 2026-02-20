import sys
if ".." not in sys.path:
    sys.path.append("..")
    
import cobra
from cobra.io import read_sbml_model
# default model choice
model = read_sbml_model("../model/Rpom_05.xml")

# Some helper functions for exploring essential genes 
# would be good to source, sink, sources into a single function which takes arrays in future.

def check_essential_rxn(rxn_id, model=model):
    essential_rxn = model.reactions.get_by_id(rxn_id)
    with model:
        model.reactions.get_by_id("EX_glc").bounds = (-5.44, 0)

        sol = model.optimize()
        print(f"Growth rate with " + rxn_id + f": {sol.objective_value:.4f} 1/hr")
        baseline_flux = sol.fluxes[essential_rxn.id]
        print(f"Flux through " + rxn_id + f": {baseline_flux}\n")

        essential_rxn.knock_out()

        sol = model.optimize()
        print(f"Growth rate without " + rxn_id + f": {sol.objective_value:.4f} 1/hr")
    if sol.objective_value == 0:
        return True, baseline_flux
    else:
        return False, baseline_flux
    
def check_essential_metabolite(metabolite_id, model=model):
    essential_metabolite = model.metabolites.get_by_id(metabolite_id)
    with model:
        model.reactions.get_by_id("EX_glc").bounds = (-5.44, 0)

        sol = model.optimize()
        print(f"Growth rate with " + metabolite_id + f": {sol.objective_value:.4f} 1/hr")
        
        for r in essential_metabolite.reactions:
            r.knock_out()
        sol = model.optimize()
        print(f"Growth rate without " + metabolite_id + f": {sol.objective_value:.4f} 1/hr")
        
    if sol.objective_value == 0:
        return True
    else:
        return False

def get_flux(rxn_id, model=model):
    my_rxn = model.reactions.get_by_id(rxn_id)
    with model:
        model.reactions.get_by_id("EX_glc").bounds = (-5.44, 0)
        sol = model.optimize()
        print(f"Growth rate: {sol.objective_value:.4f} 1/hr")
        baseline_flux = sol.fluxes[my_rxn.id]
        print(f"Flux through {my_rxn.id}: {baseline_flux}")
        return baseline_flux
    
def add_source(rxn_id, to_add_id, model=model): 
    essential_rxn = model.reactions.get_by_id(rxn_id)
    with model:
        model.reactions.get_by_id('EX_glc').bounds = (-5.44, 0)

        # Baseline growth rate is 0, with tptase knocked out
        sol = model.optimize()
        is_essential, baseline_flux = check_essential_rxn(rxn_id)

        essential_rxn.knock_out()
        # Does adding free to_add fix things?
        source = model.add_boundary(
            model.metabolites.get_by_id(to_add_id),
            type="sink",
            reaction_id=to_add_id + "_source",
            lb=-1 * baseline_flux,  # Set to exactly the lost flux from above!
            ub=0
        )
        print(f"New reaction: {source.id}\n\t{source.reaction}\t\t{source.bounds}")

        sol = model.optimize()
        print(f"Growth rate with new reaction: {sol.objective_value:.2f} 1/hr")
        print()
        
def add_sink(rxn_id, to_sink_id, model=model): 
    essential_rxn = model.reactions.get_by_id(rxn_id)
    with model:
        model.reactions.get_by_id('EX_glc').bounds = (-5.44, 0)

        # Baseline growth rate is 0, with tptase knocked out
        sol = model.optimize()
        (is_essential, baseline_flux) = check_essential_rxn(rxn_id)

        essential_rxn.knock_out()
        # Does adding free to_add fix things?
        sink = model.add_boundary(
            model.metabolites.get_by_id(to_sink_id),
            type="sink",
            reaction_id=to_sink_id + "_sink",
            ub=baseline_flux,  # Set to exactly the lost flux from above!
            lb=0
        )
        print(f"New reaction: {sink.id}\n\t{sink.reaction}\t\t{sink.bounds}")

        sol = model.optimize()
        print(f"Growth rate with new reaction: {sol.objective_value:.2f} 1/hr")
        print()

def get_metabolite_reactions(met_id, model=model): 
    met = model.metabolites.get_by_id(met_id)

    rxns = list(met.reactions)          # list of reactions
    rxn_ids = [r.id for r in rxns]      # just IDs
    consuming = [r for r in met.reactions if r.metabolites[met] < 0]
    producing = [r for r in met.reactions if r.metabolites[met] > 0]

    print(f"Consuming reactions:")
    for rxn in consuming:
        print(f"\t{rxn.id}: {rxn.reaction}")

    print(f"Producing reactions:")
    for rxn in producing:
        print(f"\t{rxn.id}: {rxn.reaction}")
    print()
    return consuming, producing
        
        
def add_source_and_sink(rxn_id, to_add_id, to_sink_id, model=model): 
    essential_rxn = model.reactions.get_by_id(rxn_id)
    with model:
        model.reactions.get_by_id('EX_glc').bounds = (-5.44, 0)

        # Baseline growth rate is 0, with tptase knocked out
        sol = model.optimize()
        is_essential, baseline_flux = check_essential_rxn(rxn_id)

        essential_rxn.knock_out()
        # Does adding free to_add fix things?
        source = model.add_boundary(
            model.metabolites.get_by_id(to_add_id),
            type="sink",
            reaction_id=to_add_id + "_source",
            lb=-1 * baseline_flux,  # Set to exactly the lost flux from above!
            ub=0
        )
        print(f"New reaction: {source.id}\n\t{source.reaction}\t\t{source.bounds}")
        
        sink = model.add_boundary(
            model.metabolites.get_by_id(to_sink_id),
            type="sink",
            reaction_id=to_sink_id + "_sink",
            ub=baseline_flux,  # Set to exactly the lost flux from above!
            lb=0
        )
        print(f"New reaction: {sink.id}\n\t{sink.reaction}\t\t{sink.bounds}")

        sol = model.optimize()
        print(f"Growth rate with new reaction: {sol.objective_value:.2f} 1/hr")
        print()
        
def add_many_source_and_sink(rxn_id, to_add_ids, to_sink_ids, model=model): 
    essential_rxn = model.reactions.get_by_id(rxn_id)
    with model:
        model.reactions.get_by_id('EX_glc').bounds = (-5.44, 0)

        # Baseline growth rate is 0, with tptase knocked out
        sol = model.optimize()
        is_essential, baseline_flux = check_essential_rxn(rxn_id)

        essential_rxn.knock_out()
        # Does adding free to_add fix things?
        for to_add_id in to_add_ids:
            source = model.add_boundary(
                model.metabolites.get_by_id(to_add_id),
                type="sink",
                reaction_id=to_add_id + "_source",
                lb=-1 * baseline_flux,  # Set to exactly the lost flux from above!
                ub=0
            )
            print(f"New reaction: {source.id}\n\t{source.reaction}\t\t{source.bounds}")
        
        for to_sink_id in to_sink_ids:
            sink = model.add_boundary(
                model.metabolites.get_by_id(to_sink_id),
                type="sink",
                reaction_id=to_sink_id + "_sink",
                ub=baseline_flux,  # Set to exactly the lost flux from above!
                lb=0
            )
            print(f"New reaction: {sink.id}\n\t{sink.reaction}\t\t{sink.bounds}")

        sol = model.optimize()
        print(f"Growth rate with new reaction: {sol.objective_value:.2f} 1/hr")
        print()

def print_reaction_info(rxn_id, model=model):
    rxn = model.reactions.get_by_id(rxn_id)
    print(f"{rxn.id}:\n\t{rxn.reaction}\n")
    print(f"Bounds: {rxn.bounds}")
    print(f"Reversibility: {rxn.reversibility}")
    print(f"Genes: {rxn.genes}")
    print(f"Subsystem: {rxn.subsystem}")
    
def search_metabolites(words, model=model):
    cand = []
    for m in model.metabolites:
        s = (m.id + " " + m.name).lower()
        if any(k in s for k in words):
            cand.append(m)

    # Print results as aligned columns
    rows = [(m.id or "", m.name or "", m.formula or "") for m in cand]
    w_id = max([len("id"), *(len(r[0]) for r in rows)] or [2])
    w_name = max([len("name"), *(len(r[1]) for r in rows)] or [4])
    w_formula = max([len("formula"), *(len(r[2]) for r in rows)] or [7])

    print(f"{'id':<{w_id}} | {'name':<{w_name}} | {'formula':<{w_formula}}")
    print(f"{'-'*w_id}-+-{'-'*w_name}-+-{'-'*w_formula}")
    for mid, name, formula in rows:
        print(f"{mid:<{w_id}} | {name:<{w_name}} | {formula:<{w_formula}}")
        
def add_sources(rxn_id, metabolite_ids, model=model):
    essential_rxn = model.reactions.get_by_id(rxn_id)
    with model:
        model.reactions.get_by_id("EX_glc").bounds = (-5.44, 0)
        essential_rxn.knock_out()
        
        for metabolite_id in metabolite_ids:
            metabolite = model.metabolites.get_by_id(metabolite_id)
            source = model.add_boundary(
                model.metabolites.get_by_id(metabolite_id),
                type="sink",
                reaction_id=metabolite_id + "_source",
                lb=-1 * baseline_flux,  # Set to exactly the lost flux from above!
                ub=0
            )
            print(f"New reaction: {source.id}\n\t{source.reaction}\t\t{source.bounds}")
        sol = model.optimize()
        print(f"Growth rate with new reaction: {sol.objective_value:.2f} 1/hr")
        print()
            
            
            
        sol = model.optimize()
        print(f"Growth rate {sol.objective_value:.4f} 1/hr")
