from argparse import ArgumentParser
from shutil import copyfile


GROWTH_RATES = ["0/h", "0.2/h", "0.25/h", "0.3/h", "0.35/h",
                "0.45/h", "0.4/h", "0.5/h", "0.55/h", "0.6/h"]


def process_args(args):
    growth_rate = args.growth

    # Set current model to use supplied growth rate
    copyfile(f"base_model/Rpom_{growth_rate.replace('/h', '').replace('.', '')}.xml",
             f"model.xml")


def cli():
    parser = ArgumentParser(
        description="Interface for FBA model of Ruegeria pomeroyi")

    parser.add_argument("-g", "--growth", choices=GROWTH_RATES, default="0.6/h",
                        help="Set model.xml to be the base model with the given growth rate.")

    args = parser.parse_args()
    process_args(args)


def main():
    cli()


if __name__ == "__main__":
    main()
