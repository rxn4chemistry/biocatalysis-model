#!/usr/bin/env python

from random import shuffle, choice
import click
from rxn_biocatalysis_tools import EnzymaticReaction, tokenize_smiles


@click.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.argument("output_file", type=click.Path())
@click.option("--pipe", is_flag=True)
@click.option("--within-class", is_flag=True)
@click.option("--shuffle-only", is_flag=True)
def main(input_file: str, output_file: str, pipe: bool, within_class, shuffle_only):
    pipe_char = ""
    if pipe:
        pipe_char = " |"

    smiles = []
    ecs = []

    with open(input_file, "r") as f:
        for line in f:
            smiles_part = line.strip().replace(" ", "")
            ec_part = ""

            ec_index = smiles_part.find("[v")
            if ec_index > -1:
                ec_part = f" {smiles_part[ec_index:].strip()}".replace("][", "] [")
                smiles_part = smiles_part[:ec_index].strip()
            smiles.append(smiles_part)
            ecs.append(ec_part)

    ecs_shuffled = ecs.copy()
    shuffle(ecs_shuffled)

    if shuffle_only:
        with open(output_file, "w+") as f:
            for smi, ec in zip(smiles, ecs_shuffled):
                f.write(f"{tokenize_smiles(smi)}{pipe_char}{ec}\n")
    else:
        with open(output_file, "w+") as f:
            for smi, ec in zip(smiles, ecs):
                ecs_tmp = [item for item in ecs if item is not ec]
                if within_class:
                    ecs_tmp = [item for item in ecs if item[3] == ec[3]]

                f.write(f"{tokenize_smiles(smi)}{pipe_char}{choice(ecs_tmp)}\n")


if __name__ == "__main__":
    main()
