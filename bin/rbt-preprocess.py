#!/usr/bin/env python
"""Preprocess reactions for learning."""

from typing import List, Tuple
from pathlib import Path
from collections import Counter

import click
import pandas as pd
from tqdm import trange

from rdkit.Chem import AllChem as rdk

from rxn_biocatalysis_tools import (
    EnzymaticReaction,
    tokenize_enzymatic_reaction_smiles,
    disable_rdkit_logging,
)


def process_reaction(
    rxn: EnzymaticReaction, min_atom_count: int = 4
) -> EnzymaticReaction:
    """Removing precursors from products, remove molecules that couldn't be parsed by RDKit
    and remove single atoms from productsself.

    Args:
        rxn: An enzymatic reaction
        min_atom_count: Remove molecules with a heavy atom count smaller than this from the reaction

    Returns:
        The processed enzymatic reaction
    """
    rxn.remove_precursors_from_products()
    rxn.remove_none()
    rxn.products = [
        product
        for product in rxn.products
        if product.GetNumHeavyAtoms() >= min_atom_count
    ]
    return rxn


def remove_from_products(
    rxn: EnzymaticReaction, smart_patterns: List[str], smiles: List[str]
) -> EnzymaticReaction:
    """Remove the supplied cofactors from the products.

    Args:
        rxn: An enzymatic reaction
        cofactors: A list of SMILES of cofactors

    Returns:
        The enzymatic reaction with the cofactors removed from the products
    """

    if len(rxn.products) == 1:
        return rxn

    for smart_pattern in smart_patterns:
        mol_pattern = rdk.MolFromSmarts(smart_pattern)

        rxn.products = [
            mol for mol in rxn.products if not mol.HasSubstructMatch(mol_pattern)
        ]

        if len(rxn.products) == 1:
            return rxn

    rxn.products = [mol for mol in rxn.products if rxn.mol_to_smiles(mol) not in smiles]

    return rxn


def print_sources(rxns: List[EnzymaticReaction], title: str) -> None:
    """Print the sources and numbers of reactions per source as a table.

    Args:
        rxns: A list of enzymatic reactions
        title: The title of the table
    """
    sources = [rxn.source for rxn in rxns]

    print(f"\n\n{title}:")
    for key, value in dict(Counter(sources)).items():
        print(f"{key} -", str(value))

    print("Total:", len(rxns))


def write_splits(df: pd.DataFrame, ec_level: int, output_dir: Path) -> None:
    """Writing the training, validation, and testing src and tgt splits for a given EC level into a given output directory.

    Args:
        df: A pandas DataFrame containing the reaction strings as well as a field "ec_3" to grup the reactions by.
        ec_level: The EC level to be exported (0 to 4)
        output_dir: The path to the directory the outputs will be written to.
    """
    df_internal = df.copy()
    splits = {}
    splits["train"] = pd.DataFrame(columns=df_internal.columns)
    splits["valid"] = pd.DataFrame(columns=df_internal.columns)
    splits["test"] = pd.DataFrame(columns=df_internal.columns)

    df_internal[["rxn_str"]].to_csv(
        Path(output_dir, "combined.txt"), header=False, index=False
    )

    # Get the training set from reactions with unique products¨
    prod_val_cnts = df_internal.products.value_counts()
    df_unique_prods = df_internal[
        df_internal.products.isin(prod_val_cnts.index[prod_val_cnts.eq(1)])
    ]
    df_internal.drop(df_unique_prods.index, inplace=True)

    # Always group by x.x.x.- for splits to have a good coverage and
    # not miss sub-sub-classes
    for ec in df_internal.ec_1.unique():
        df_subset = df_internal[df_internal.ec_1 == ec]
        df_unique_prods_subset = df_unique_prods[df_unique_prods.ec_1 == ec]

        # Get 5% (of total) from rxns with unique products for test set
        n_total = len(df_subset) + len(df_unique_prods_subset)
        n_samples = round(n_total * 0.05)
        df_test = df_unique_prods_subset.sample(
            n=min(n_samples, len(df_unique_prods_subset)), random_state=42
        )
        df_unique_prods_subset.drop(df_test.index, inplace=True)

        # Join again to get training and validation sets
        df_subset = df_subset.append(df_unique_prods_subset)

        df_valid = df_subset.sample(n=n_samples, random_state=42)
        df_subset.drop(df_valid.index, inplace=True)

        splits["train"] = splits["train"].append(df_subset, ignore_index=True)
        splits["valid"] = splits["valid"].append(df_valid, ignore_index=True)
        splits["test"] = splits["test"].append(df_test, ignore_index=True)

    for key, value in splits.items():
        value["reactants"].to_csv(
            Path(output_dir, f"src-{key}.txt"), header=False, index=False
        )
        value["products"].to_csv(
            Path(output_dir, f"tgt-{key}.txt"), header=False, index=False
        )


@click.command()
@click.argument("input_files", type=click.Path(exists=True), nargs=-1)
@click.argument("output_path", type=click.Path(exists=True), nargs=1)
@click.option(
    "--remove-patterns",
    "remove_patterns_path",
    type=click.Path(exists=True),
    help="A file containing SMARTS patterns (one per line); molecules which match a pattern will be removed from the products",
)
@click.option(
    "--remove-molecules",
    "remove_molecules_path",
    type=click.Path(exists=True),
    help="A file containing molecules encoded as SMILES (one per line) to be removed from the products",
)
@click.option(
    "--ec-level",
    "ec_levels",
    multiple=True,
    type=int,
    default=[3],
    show_default=True,
    help="The EC level(s) to produce files for. This option can be repeated for multiple levels (--ec-level 1 --ec-level 2 ...)",
)
@click.option(
    "--max-products",
    "max_products",
    type=int,
    default=1,
    show_default=True,
    help="The maximum number of product molecules allowed in a reaction",
)
@click.option(
    "--min-atom-count",
    "min_atom_count",
    type=int,
    default=4,
    show_default=True,
    help="The minimum atom count of a molecule (molecules with less atoms are removed from the reactions)",
)
@click.option(
    "--bi-directional",
    "bi_directional",
    is_flag=True,
    help="Consider all reactions to be bi-directional",
)
@click.option(
    "--split-products",
    "split_products",
    is_flag=True,
    help="Wheter to split reactions with multiple products into multiple reactions with one product",
)
def main(
    input_files: str,
    output_path: str,
    remove_patterns_path: str,
    remove_molecules_path: str,
    ec_levels: Tuple[int, ...],
    max_products: int,
    min_atom_count: int,
    bi_directional: bool,
    split_products: bool,
):
    disable_rdkit_logging()
    remove_patterns = []
    remove_molecules = []

    if remove_patterns_path:
        for line in open(remove_patterns_path):
            if not line.startswith("//") and line.strip():
                remove_patterns.append(line.split("//")[0].strip())

    if remove_molecules_path:
        for line in open(remove_molecules_path):
            if not line.startswith("//") and line.strip():
                smiles = line.split("//")[0].strip()
                mol = rdk.MolFromSmiles(smiles)
                if mol:
                    remove_molecules.append(rdk.MolToSmiles(mol))

    enzymatic_reactions = []

    for input_file in input_files:
        source = Path(input_file).stem
        line_count = sum([1 for line in open(input_file, "r")])

        print(f"Parsing {source}...")
        with open(input_file, "r") as f:
            for line in f:
                rxn = EnzymaticReaction(line.strip(), canonical=True, source=source)
                enzymatic_reactions.append(rxn)
                if bi_directional:
                    enzymatic_reactions.append(rxn.reverse())
        print("Done.\n")

    print_sources(enzymatic_reactions, "Parsed Reactions")

    print(f"Processing reactions...")
    for i, rxn in enumerate(enzymatic_reactions):
        enzymatic_reactions[i] = process_reaction(rxn, min_atom_count)
    print("Done.\n")

    print("Removing patterns and molecules...")
    for i, rxn in enumerate(enzymatic_reactions):
        enzymatic_reactions[i] = remove_from_products(
            rxn, remove_patterns, remove_molecules
        )
    print("Done.\n")

    if split_products:
        print(
            "Splitting multi-product reactions into multiple one-product reactions..."
        )
        split_enzymatic_reactions = []
        for rxn in enzymatic_reactions:
            if len(rxn.products) == 1:
                split_enzymatic_reactions.append(rxn)
            else:
                for product in rxn.get_products_as_smiles():
                    split_enzymatic_reactions.append(
                        EnzymaticReaction.from_smarts_and_ec(
                            ".".join(rxn.get_reactants_as_smiles()) + ">>" + product,
                            rxn.get_ec(),
                            source=rxn.source,
                        )
                    )
        enzymatic_reactions = split_enzymatic_reactions
        print("Done.\n")

    print(
        f"Removing reactions with less than 2 reactans or more than {max_products} product(s)..."
    )
    enzymatic_reactions = [
        rxn
        for rxn in enzymatic_reactions
        if len(rxn.reactants) > 0
        and len(rxn.products) <= max_products
        and len(rxn.products) > 0
    ]
    print("Done.\n")

    print("Ordering molecules on both sides of the reactions...")
    for rxn in enzymatic_reactions:
        rxn.sort()
    print("Done.")

    print_sources(enzymatic_reactions, "Reactions after Filtering")

    df = pd.DataFrame({"rxn": enzymatic_reactions})

    print("Processing reactions for export...")

    cols = {}

    cols["rxn_str"] = []
    cols["ec"] = []
    cols["source"] = []

    for ec_level in range(1, 5):
        cols[f"ec_{ec_level}"] = []

    for ec_level in ec_levels:
        cols[f"rxn_str_ec{ec_level}"] = []
        cols[f"rxn_str_ec{ec_level}_tok"] = []

    for rxn in df.rxn:
        cols["rxn_str"].append(str(rxn))
        cols["ec"].append(rxn.get_ec())
        cols["source"].append(rxn.source)

        for ec_level in range(1, 5):
            cols[f"ec_{ec_level}"].append(rxn.get_ec(ec_level))

        for ec_level in ec_levels:
            rxn_str_level = rxn.to_string(ec_level)
            cols[f"rxn_str_ec{ec_level}"].append(rxn_str_level)

            cols[f"rxn_str_ec{ec_level}_tok"].append(
                tokenize_enzymatic_reaction_smiles(rxn_str_level, keep_pipe=True)
            )

    for key, value in cols.items():
        df[key] = value

    df = df[(df.ec != "") & (df.ec != " ")]
    df = df.drop_duplicates(subset=[f"rxn_str"], keep="first")

    print("Done.\n")

    print_sources(df.rxn, "Reactions after Deduplication")

    df[["rxn_str", "ec", "source"]].to_csv(
        Path(output_path, "combined_rxn_ec_sources.txt"), header=False, index=False
    )

    tasks = [f"exporting EC level {ec_level}" for ec_level in ec_levels]

    print("Exporting src and tgt files...")
    for ec_level in ec_levels:
        task = tasks.pop(0)

        df_tmp = df.copy()
        df_tmp = df_tmp.drop_duplicates(subset=[f"rxn_str_ec{ec_level}"], keep="first")

        # If someone finds a way on how to combne a split and assignment to two
        # columns with .loc, please change it. Otherwise surpress the warning
        # to get clean console output.
        pd.options.mode.chained_assignment = None
        df_tmp[["reactants", "products"]] = df_tmp[
            f"rxn_str_ec{ec_level}_tok"
        ].str.split(">>", expand=True)

        df_tmp.reactants = df_tmp.reactants.str.strip()
        df_tmp.products = df_tmp.products.str.strip()

        parent_dir = Path(output_path, f"experiments/{ec_level}")
        parent_dir.mkdir(parents=True, exist_ok=True)

        write_splits(df_tmp, ec_level, parent_dir)

        print(f"{task} complete ({len(df_tmp)} unique reactions).")


if __name__ == "__main__":
    main()
