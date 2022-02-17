"""Testing enzymatic reaction class."""
import pytest
from rxn_biocatalysis_tools import EnzymaticReaction


@pytest.fixture
def enzymatic_reaction():
    return EnzymaticReaction(
        "N[C@@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O|4.1.1.28>>NCCc1c[nH]c2ccc(F)cc12"
    )


def test_constructor(enzymatic_reaction: EnzymaticReaction):
    assert enzymatic_reaction.ec == ["4", "1", "1", "28"]
    assert str(enzymatic_reaction) == (
        "N[C@@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O|4.1.1.28>>NCCc1c[nH]c2ccc(F)cc12"
    )


def test_mol_to_smiles(enzymatic_reaction: EnzymaticReaction):
    assert (
        enzymatic_reaction.mol_to_smiles(enzymatic_reaction.products[0])
        == "NCCc1c[nH]c2ccc(F)cc12"
    )


def test_to_string(enzymatic_reaction: EnzymaticReaction):
    assert (
        enzymatic_reaction.to_string(0)
        == "N[C@@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O>>NCCc1c[nH]c2ccc(F)cc12"
    )
    assert (
        enzymatic_reaction.to_string(1)
        == "N[C@@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O|4>>NCCc1c[nH]c2ccc(F)cc12"
    )
    assert (
        enzymatic_reaction.to_string(2)
        == "N[C@@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O|4.1>>NCCc1c[nH]c2ccc(F)cc12"
    )
    assert (
        enzymatic_reaction.to_string(3)
        == "N[C@@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O|4.1.1>>NCCc1c[nH]c2ccc(F)cc12"
    )
    assert (
        enzymatic_reaction.to_string(4)
        == "N[C@@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O|4.1.1.28>>NCCc1c[nH]c2ccc(F)cc12"
    )
    assert (
        enzymatic_reaction.to_string()
        == "N[C@@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O|4.1.1.28>>NCCc1c[nH]c2ccc(F)cc12"
    )


def test_get_ec(enzymatic_reaction: EnzymaticReaction):
    assert enzymatic_reaction.get_ec(1) == "4"
    assert enzymatic_reaction.get_ec(2) == "4.1"
    assert enzymatic_reaction.get_ec(3) == "4.1.1"
    assert enzymatic_reaction.get_ec(4) == "4.1.1.28"


def test_reverse(enzymatic_reaction: EnzymaticReaction):
    assert (
        enzymatic_reaction.reverse().to_string()
        == "NCCc1c[nH]c2ccc(F)cc12|4.1.1.28>>N[C@@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O"
    )


def test_from_smarts_and_ec(enzymatic_reaction: EnzymaticReaction):
    assert enzymatic_reaction == EnzymaticReaction.from_smarts_and_ec(
        "N[C@@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O>>NCCc1c[nH]c2ccc(F)cc12", "4.1.1.28"
    )


def test_is_valid(enzymatic_reaction: EnzymaticReaction):
    assert EnzymaticReaction.is_valid(
        "NCCc1c[nH]c2ccc(F)cc12|4.1.1.28>>N[C@@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O"
    )
    assert not EnzymaticReaction.is_valid(
        "N[C@@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O>>NCCc1c[nH]c2ccc(F)cc12"
    )
    assert not EnzymaticReaction.is_valid(
        "N[C@@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O|>>NCCc1c[nH]c2ccc(F)cc12"
    )
    assert not EnzymaticReaction.is_valid(
        "N[C@@H](Cc1c[nH]c2ccc(F)cc12)C(=O)O>NCCc1c[nH]c2ccc(F)cc12"
    )
