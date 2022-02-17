"""Enzymatic reaction representation."""
import re
from typing import List, Any
from .chemical_reaction import ChemicalReaction
from rdkit.Chem import AllChem as rdk
from rdkit.Chem.rdchem import Mol

UNKNOWN_CHEMICAL_REGEX = re.compile(r"^(<.*>)$|^(<)|(>)$")


class EnzymaticReaction(ChemicalReaction):
    """Representation of an enzymatic reaction.

    Reactions containing enzyme are represented as reaction SMILES using a '|'
    to separate precursors from the EC number.
    """

    def __init__(
        self,
        enzymatic_reaction_smiles: str,
        remove_duplicates: bool = True,
        sanitize: bool = True,
        source: str = "unknown",
        **kwargs: Any,
    ):
        """Constructor for EnzymaticReaction.

        Args:
            enzymatic_reaction_smiles: an enzymatic reaction SMILES.
            remove_duplicates: duplicate removal. Defaults to True.
            sanitize: whether sanitization is enabled. Defaults to True.
            source: source for the enzymatic reaction. Defaults to "unknown".
        """
        vals = re.split(r">|\|", enzymatic_reaction_smiles)
        if len(vals) < 2:
            vals.append("")

        self.ec: List[str] = [level.strip() for level in vals[1].split(".")]
        self.source = source

        # hack
        self.kwargs = kwargs

        super().__init__(
            enzymatic_reaction_smiles.replace(f"|{vals[1]}", ""),
            remove_duplicates,
            sanitize,
            **kwargs,
        )

    def __str__(self) -> str:
        """Returns the extended reaction SMARTS of this instance (reactants|ec>agents>products).

        Returns:
            the extended reaction SMARTS representing this instance.
        """
        s = (
            ".".join(
                sorted([rdk.MolToSmiles(m, **self.kwargs) for m in self.reactants if m])
            )
            + ">"
            + ".".join(
                sorted([rdk.MolToSmiles(m, **self.kwargs) for m in self.agents if m])
            )
            + ">"
            + ".".join(
                sorted([rdk.MolToSmiles(m, **self.kwargs) for m in self.products if m])
            )
        )
        s_parts = s.split(">")

        if len(self.ec) > 0 and self.ec[0] != "":
            s_parts[0] += f'|{".".join(self.ec)}'
        return ">".join(s_parts).replace(" ", "")

    def __eq__(self, other: object) -> bool:
        """Compares the count, order, and SMILES string of each molecule in this reaction as well as the EC.

        Args:
            other: another EnzymaticReaction instance to be compared with this instance.

        Returns:
            whether this instance is equal to another.
        """
        if not isinstance(other, EnzymaticReaction):
            raise NotImplementedError(
                "EnzymaticReaction can be tested for equality with EnzymaticReaction objects"
            )
        return super().__eq__(other) and self.ec == other.ec

    def __hash__(self) -> int:
        """Get hash for the enzymatic reaction.

        Returns:
            enzymatic reaction hash.
        """
        return hash(str(self))

    def mol_to_smiles(self, mol: Mol) -> str:
        """Applies the kwargs supplied to the reaction to MolToSmiles for a given molecule.

        Args:
            mol: an RDKit molecule instance.

        Returns:
           the string representing the molecule.
        """
        return rdk.MolToSmiles(mol, **self.kwargs)

    def to_string(self, ec_depth: int = 4) -> str:
        """Get the string representing this reaction with a certain number of EC levels.

        Args:
            ec_depth: the number of EC classes to include (top-down). Defaults to 4.

        Returns:
           the string representing this reaction with the chosen levels of EC.
        """
        cpy = EnzymaticReaction(str(self))
        cpy.ec = cpy.ec[:ec_depth]
        return str(cpy).strip()

    def get_ec(self, ec_depth: int = 4) -> str:
        """Get the string representing the EC of this reaction.

        Args:
            ec_depth: the number of EC classes to include (top-down). Defaults to 4.

        Returns:
           the EC of the reaction as a string.
        """
        return ".".join(self.ec[:ec_depth]).strip()

    def reverse(self) -> "EnzymaticReaction":
        """Reverses the reaction (switching reactants and products).

        Returns:
            the reversed enzymatic reactions.
        """
        return EnzymaticReaction.from_smarts_and_ec(
            f"{'.'.join(self.get_products_as_smiles())}>>{'.'.join(self.get_reactants_as_smiles())}",
            self.get_ec(),
            self.source,
        )

    @staticmethod
    def from_smarts_and_ec(
        reaction_smiles: str, ec: str, source: str = "unknown"
    ) -> "EnzymaticReaction":
        """Creates an EnzymaticReaction instance from a reaction SMILES and an EC number.

         Args:
            reaction_smiles: a reaction SMILES.
            ec: EC number string representation.
            source: source for the enzymatic reaction. Defaults to "unknown".

        Returns:
            an EnzymaticReaction instance.
        """
        split = reaction_smiles.split(">>")
        return EnzymaticReaction(split[0] + "|" + ec + ">>" + split[1], source=source)

    @staticmethod
    def is_valid(enzymatic_reaction_smiles: str) -> bool:
        """Checks whether an enzymatic reaction SMILES (e.g. O.CO|1.2.3.4>>C(=O)O) is valid.

        Args:
            enzymatic_reaction_smiles: an enzymatic reaction SMILES.

        Returns:
            a bool indicating whether the supplied enzymatic reaction SMILES is valid.
        """
        if (
            "|" not in enzymatic_reaction_smiles
            or enzymatic_reaction_smiles.count(">") != 2
            or "|>>" in enzymatic_reaction_smiles
        ):
            return False

        return True
