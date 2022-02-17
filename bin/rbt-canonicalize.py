#!/usr/bin/env python

import click
from rxn_biocatalysis_tools import EnzymaticReaction, tokenize_smiles


@click.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.argument("output_file", type=click.Path())
@click.option("--pipe", is_flag=True)
def main(input_file: str, output_file: str, pipe: bool):
    pipe_char = ""
    if pipe:
        pipe_char = " |"

    smiles = []
    with open(input_file, "r") as f:
        for line in f:
            smiles.append(line.strip())

    with open(output_file, "w+") as f:
        for original_smiles in smiles:
            # In case of a predicted source which includes the EC number (always starting with "[v")
            # we need to take care of this
            smiles_part = original_smiles.replace(" ", "")
            ec_part = ""

            # Writing non-sense items that occur in the backward prediction as is
            try:
                ec_index = smiles_part.find("[v")
                if ec_index > -1:
                    ec_part = f" {smiles_part[ec_index:].strip()}".replace("][", "] [")
                    smiles_part = smiles_part[:ec_index].strip()

                # Using the EnzymaticReaction class here to get canonicalisation + ordering
                # ">>" is needed for it to be recoganised as a valid rxn smiles / smarts
                rxn = EnzymaticReaction(smiles_part + ">>")
                rxn.sort()
                sorted_canonicalised_smiles = ".".join(rxn.get_reactants_as_smiles())
                f.write(
                    f"{tokenize_smiles(sorted_canonicalised_smiles)}{pipe_char}{ec_part}\n"
                )
            except:
                f.write(f"{original_smiles}\n")


if __name__ == "__main__":
    main()
