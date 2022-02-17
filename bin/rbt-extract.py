#!/usr/bin/env python

import click
from rxn_biocatalysis_tools import EnzymaticReaction


@click.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.argument("output_file", type=click.Path())
@click.option("--level", default=3)
@click.option(
    "--data",
    type=click.Choice(["reactants", "products"], case_sensitive=False),
    default="products",
)
def main(input_file: str, output_file: str, level: int, data: str):
    with open(input_file, "r") as f_in:
        with open(output_file, "w+") as f_out:
            if data == "products":
                for line in f_in:
                    rxn = EnzymaticReaction(line)
                    for smile in rxn.get_products_as_smiles():
                        f_out.write(f"{smile},{rxn.get_ec(level)},{rxn.get_ec()}\n")
            if data == "reactants":
                for line in f_in:
                    rxn = EnzymaticReaction(line)
                    for smile in rxn.get_reactants_as_smiles():
                        f_out.write(f"{smile},{rxn.get_ec(level)},{rxn.get_ec()}\n")


if __name__ == "__main__":
    main()
