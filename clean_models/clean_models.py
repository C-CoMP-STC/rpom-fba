"""
clean_models.py

Script for validating and cleaning the initial models (second step of the pipeline,
following model conversion from .mat to .sbml).
"""


def main():
    # TODO: EX_glc is reversible despite being "irreversible" in "rev"?
    # But exchange reactions should be reversible, so all the other EX reactions are wrong?
    pass


if __name__ == "__main__":
    main()
