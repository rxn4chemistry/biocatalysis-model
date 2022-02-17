"""Testing enzymatic reaction tokenizer."""
import pytest
from rxn_biocatalysis_tools import (
    tokenize_enzymatic_reaction_smiles,
    detokenize_enzymatic_reaction_smiles,
)


@pytest.fixture
def enzymatic_reaction_smiles():
    return "N[C@@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O|4.1.1.28>>NCCc1c[nH]c2ccc(F)cc12"


@pytest.fixture
def tok_enzymatic_reaction_smiles():
    return "N [C@@H] ( C c 1 c [nH] c 2 c c c ( F ) c c 1 2 ) C ( = O ) O [v4] [u1] [t1] [q28] >> N C C c 1 c [nH] c 2 c c c ( F ) c c 1 2"


@pytest.fixture
def enzymatic_reaction_smiles_part():
    return "N[C@@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O|4.1>>NCCc1c[nH]c2ccc(F)cc12"


@pytest.fixture
def tok_enzymatic_reaction_smiles_part():
    return "N [C@@H] ( C c 1 c [nH] c 2 c c c ( F ) c c 1 2 ) C ( = O ) O [v4] [u1] >> N C C c 1 c [nH] c 2 c c c ( F ) c c 1 2"


def test_tokenize_enzymatic_reaction_smiles(
    enzymatic_reaction_smiles: str, tok_enzymatic_reaction_smiles: str
):
    assert (
        tokenize_enzymatic_reaction_smiles(enzymatic_reaction_smiles)
        == tok_enzymatic_reaction_smiles
    )


def test_detokenize_enzymatic_reaction_smiles(
    enzymatic_reaction_smiles: str, tok_enzymatic_reaction_smiles: str
):
    assert (
        detokenize_enzymatic_reaction_smiles(tok_enzymatic_reaction_smiles)
        == enzymatic_reaction_smiles
    )


def test_tokenize_enzymatic_reaction_smiles_part(
    enzymatic_reaction_smiles_part: str, tok_enzymatic_reaction_smiles_part: str
):
    assert (
        tokenize_enzymatic_reaction_smiles(enzymatic_reaction_smiles_part)
        == tok_enzymatic_reaction_smiles_part
    )


def test_detokenize_enzymatic_reaction_smiles_part(
    enzymatic_reaction_smiles_part: str, tok_enzymatic_reaction_smiles_part: str
):
    assert (
        detokenize_enzymatic_reaction_smiles(tok_enzymatic_reaction_smiles_part)
        == enzymatic_reaction_smiles_part
    )
