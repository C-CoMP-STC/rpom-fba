import numpy as np
# import matplotlib
# matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

HIGHLIGHT_BLUE = "#1a66ff"
BACKGROUND_GRAY = "#cccccc"

C_BIOMASS = "#1a66ff"
C_GLUCOSE = "#ff953a"
C_ACETATE = "#ff5b82"

C_PROTEIN = "#312ab6"
C_LIPID = "#ffaf73"
C_RNA = "#ac5bc2"
C_PHB = "#ffd328"
C_DNA = "#c6184d"
C_MUREIN = "#aedda2"


def main():
    colors = [
        HIGHLIGHT_BLUE,
        BACKGROUND_GRAY,
        C_BIOMASS,
        C_GLUCOSE,
        C_ACETATE,
        C_PROTEIN,
        C_LIPID,
        C_RNA,
        C_PHB,
        C_DNA,
        C_MUREIN
    ]
    
    fig, ax = plt.subplots()
    ax.add_collection(PatchCollection([Rectangle((i, 0), 1, 1)
                                       for i in range(len(colors))],
                                      fc=colors))
    ax.set_xlim(0, len(colors))
    fig.savefig("colors.png")


if __name__ == "__main__":
    main()
