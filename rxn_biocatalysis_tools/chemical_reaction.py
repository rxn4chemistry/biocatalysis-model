""" Contains the class Reaction representing unidirectional reactions. """
from enum import Enum
from typing import (
    Tuple,
    List,
    Any,
    overload,
    Iterable,
    Optional,
    Callable,
    Set,
    TypeVar,
    cast,
)
from rdkit.Chem import AllChem as rdk
from rdkit.Chem.rdchem import Mol

T = TypeVar("T")
V = TypeVar("V")


class ChemicalReactionPart(Enum):
    """An enum used to specify the three different parts of a chemical Reaction."""

    reactants = 1
    agents = 2
    products = 3


class ChemicalReaction:
    def __init__(
        self,
        reaction_smarts: str,
        remove_duplicates: bool = False,
        sanitize: bool = True,
        **kwargs: Any
    ):
        """Creates a new instance of type Reaction based on a reaction SMARTs.

        Args:
            reaction_smarts: A reaction smarts
            remove_duplicates: Whether to remove duplicates from within reactants, agents and products. Defaults to False.
            sanitize: Whether to pass the sanitize keyword argument to rdkit's MolFromSmiles.
            kwargs: Keyword arguments supplied to rdkit's MolToSmiles. Defaults: canonical=True
        """
        kwargs.setdefault("canonical", True)

        self.__reaction_smarts = reaction_smarts
        self.__remove_duplicates = remove_duplicates
        self.__smiles_to_mol_kwargs = kwargs
        self.reactants, self.agents, self.products = self.__reaction_to_mols(
            self.__reaction_smarts, sanitize
        )

    #
    # Overwrites / Virtuals
    #

    def __len__(self) -> int:
        """Returns the number of molecules participating in this reaction (reactants, agents, and products).

        Returns:
            The number of molecules participating in this reaction.
        """
        return len(self.reactants) + len(self.agents) + len(self.products)

    def __str__(self) -> str:
        """Returns the reaction SMARTS of this instance (reactants>agents>products).

        Returns:
            The reaction SMARTS representing this instance.
        """
        return (
            ".".join(
                [
                    rdk.MolToSmiles(m, **self.__smiles_to_mol_kwargs)
                    for m in self.reactants
                    if m
                ]
            )
            + ">"
            + ".".join(
                [
                    rdk.MolToSmiles(m, **self.__smiles_to_mol_kwargs)
                    for m in self.agents
                    if m
                ]
            )
            + ">"
            + ".".join(
                [
                    rdk.MolToSmiles(m, **self.__smiles_to_mol_kwargs)
                    for m in self.products
                    if m
                ]
            )
        )

    def __eq__(self, other) -> bool:
        """Compares the count, order, and SMILES string of each molecule in this reaction.

        Args:
            other: A reaction to be compared with this instance.

        Returns:
            Whether this instance is equal to another.
        """
        if len(self) != len(other):
            return False

        if len(self.reactants) != len(other.reactants):
            return False

        if len(self.agents) != len(other.agents):
            return False

        if len(self.products) != len(other.products):
            return False

        # We care what the output of MolToSmiles is, so check for equality
        # on this. Inefficient, but allows us to mutate the molecule at
        # any time
        molecules_self = (self.reactants, self.agents, self.products)
        molecules_other = (other.reactants, other.agents, other.products)
        for group_self, group_other in zip(molecules_self, molecules_other):
            for a, b in zip(group_self, group_other):
                if not self.__mols_equal(a, b):
                    return False

        return True

    #
    # Private Methods
    #

    def __reaction_to_mols(
        self, reaction_smarts: str, sanitize: bool
    ) -> Tuple[List[Mol], List[Mol], List[Mol]]:
        """Creates a tuple of lists of reactants, agents, and products as rdkit Mol instances from a reaction SMARTS.

        Args:
            reaction_smarts: A reaction SMARTS.
            sanitize: Whether to pass the sanitize keyword argument to rdkit's MolFromSmiles.

        Raises:
            ValueError: This error is raised if a non-valid reaction SMARTS is provided.

        Returns:
            A tuple of lists of reactants, agents, and products representing the reaction.
        """
        if reaction_smarts.count(">") != 2:
            raise ValueError('A valid SMARTS reaction must contain two ">".')

        raw_reactants, raw_agents, raw_products = tuple(
            [p.split(".") for p in reaction_smarts.split(">")]
        )

        if self.__remove_duplicates:
            raw_reactants = ChemicalReaction.remove_duplicates(raw_reactants)
            raw_agents = ChemicalReaction.remove_duplicates(raw_agents)
            raw_products = ChemicalReaction.remove_duplicates(raw_products)

        return (
            [
                rdk.MolFromSmiles(reactant, sanitize=sanitize)
                for reactant in raw_reactants
                if reactant != ""
            ],
            [
                rdk.MolFromSmiles(agent, sanitize=sanitize)
                for agent in raw_agents
                if agent != ""
            ],
            [
                rdk.MolFromSmiles(product, sanitize=sanitize)
                for product in raw_products
                if product != ""
            ],
        )

    def __mol_to_smiles(self, mol: Mol) -> str:
        """Encodes a molecule as a SMILES string by applying the rdkit MolToSmiles arguments supplied to this instantce.

        Args:
            mol: An rdkit Mol instance.

        Returns:
            The SMILES encoding of the input Mol.
        """
        return rdk.MolToSmiles(mol, **self.__smiles_to_mol_kwargs)

    def __mols_equal(self, first_mol: Mol, second_mol: Mol) -> bool:
        return self.__mol_to_smiles(first_mol) == self.__mol_to_smiles(second_mol)

    #
    # Pubilc Methods
    #

    def get_reactants_as_smiles(self) -> List[str]:
        """Returns the reactants of this reactions as a list of SMILES.

        Returns:
            A list of SMILES of the reactants.
        """
        return [
            rdk.MolToSmiles(reactant, **self.__smiles_to_mol_kwargs)
            for reactant in self.reactants
            if reactant
        ]

    def get_agents_as_smiles(self) -> List[str]:
        """Returns the agents of this reactions as a list of SMILES.

        Returns:
            A list of SMILES of the agents.
        """
        return [
            rdk.MolToSmiles(agent, **self.__smiles_to_mol_kwargs)
            for agent in self.agents
            if agent
        ]

    def get_products_as_smiles(self) -> List[str]:
        """Returns the products of this reactions as a list of SMILES.

        Returns:
            A list of SMILES of the products.
        """
        return [
            rdk.MolToSmiles(product, **self.__smiles_to_mol_kwargs)
            for product in self.products
            if product
        ]

    def find(self, pattern: str) -> Tuple[List[int], List[int], List[int]]:
        """Find the occurences of a SMARTS pattern within the reaction and returns a tuple
           of lists of indices in the reactants, agents, and products.

        Args:
            pattern: A SMARTS pattern.

        Returns:
            A tuple of lists of indices from the lists of reactants, agents, and products.
        """
        p = rdk.MolFromSmarts(pattern)

        # Avoid three method calls and do it directly
        return (
            [
                i
                for i, m in enumerate(self.reactants)
                if m and len(list(m.GetSubstructMatch(p))) > 0
            ],
            [
                i
                for i, m in enumerate(self.agents)
                if m and len(list(m.GetSubstructMatch(p))) > 0
            ],
            [
                i
                for i, m in enumerate(self.products)
                if m and len(list(m.GetSubstructMatch(p))) > 0
            ],
        )

    def find_in(self, pattern: str, reaction_part: ChemicalReactionPart) -> List[int]:
        """Finds a SMARTS pattern in a part (reactants, agents, products) of the reaction.

        Args:
            pattern: A SMARTS pattern.
            reaction_part (ChemicalReactionPart): The reaction part to search.

        Returns:
            A list of indices from the list of molecules representing the chosen reaction part.
        """
        p = rdk.MolFromSmarts(pattern)

        if reaction_part == ChemicalReactionPart.reactants:
            return [
                i
                for i, m in enumerate(self.reactants)
                if m and len(list(m.GetSubstructMatch(p))) > 0
            ]

        if reaction_part == ChemicalReactionPart.agents:
            return [
                i
                for i, m in enumerate(self.agents)
                if m and len(list(m.GetSubstructMatch(p))) > 0
            ]

        if reaction_part == ChemicalReactionPart.products:
            return [
                i
                for i, m in enumerate(self.products)
                if m and len(list(m.GetSubstructMatch(p))) > 0
            ]

        return []

    @overload
    def remove(self, indices: Tuple[List[int]]) -> None:
        ...

    @overload
    def remove(self, indices: Tuple[List[int], List[int]]) -> None:
        ...

    @overload
    def remove(self, indices: Tuple[List[int], List[int], List[int]]) -> None:
        ...

    def remove(self, indices):
        """Remove reactants, agents and products based on their index within the respective lists.

        Args:
            indices: The indices of the molecules to be removed from the reaction.
        """
        if len(indices) > 0:
            for idx in sorted(indices[0], reverse=True):
                del self.reactants[idx]

        if len(indices) > 1:
            for idx in sorted(indices[1], reverse=True):
                del self.agents[idx]

        if len(indices) > 2:
            for idx in sorted(indices[2], reverse=True):
                del self.products[idx]

    @overload
    def filter(self, indices: Tuple[List[int]]) -> None:
        ...

    @overload
    def filter(self, indices: Tuple[List[int], List[int]]) -> None:
        ...

    @overload
    def filter(self, indices: Tuple[List[int], List[int], List[int]]) -> None:
        ...

    def filter(self, indices):
        """Filter for reactants, agents and products based on their index within the respective lists. This is the complement to remove.

        Args:
            indices: The indices of the molecules to not be removed from the reaction.
        """
        if len(indices) > 0 and len(indices[0]) > 0:
            for idx in range(len(self.reactants) - 1, -1, -1):
                if idx not in indices[0]:
                    del self.reactants[idx]

        if len(indices) > 1 and len(indices[1]) > 0:
            for idx in range(len(self.agents) - 1, -1, -1):
                if idx not in indices[1]:
                    del self.agents[idx]

        if len(indices) > 2 and len(indices[2]) > 0:
            for idx in range(len(self.products) - 1, -1, -1):
                if idx not in indices[2]:
                    del self.products[idx]

    def sort(self, sort_reactants=True, sort_agents=True, sort_products=True) -> None:
        """Order the molecules participating in this reaction based on their SMILES strings.
           The rdkit MolToSmiles argument supplied to this instance will be applied.

        Args:
            sort_reactants: Whether to sort the reactants. Defaults to True.
            sort_agents: Whether to sort the agents. Defaults to True.
            sort_products: Whether to sort the products. Defaults to True.
        """
        # Mol to SMILES here is rather inefficient, but this allows for
        # changes to the Mol objects at any time in the lifecycle of the instance
        if sort_reactants:
            self.reactants = sorted(
                self.reactants,
                key=lambda m: rdk.MolToSmiles(m, **self.__smiles_to_mol_kwargs),
            )

        if sort_agents:
            self.agents = sorted(
                self.agents,
                key=lambda m: rdk.MolToSmiles(m, **self.__smiles_to_mol_kwargs),
            )

        if sort_products:
            self.products = sorted(
                self.products,
                key=lambda m: rdk.MolToSmiles(m, **self.__smiles_to_mol_kwargs),
            )

    def remove_precursors_from_products(self) -> None:
        """Removes prodcuts that are also found in reactants or agents."""

        reactants_smiles = self.get_reactants_as_smiles()
        agents_smiles = self.get_agents_as_smiles()
        products_smiles = self.get_products_as_smiles()

        for i, product in reversed(list(enumerate(products_smiles))):
            if product in reactants_smiles or product in agents_smiles:
                del self.products[i]

    def has_none(self) -> bool:
        """Checks whether the reactants, agents, or products contain None (usually due to failed rdkit MolFromSmiles).

        Returns:
            Whether the reactants, agents, or products contain None.
        """

        return None in self.reactants or None in self.agents or None in self.products

    def remove_none(self) -> None:
        """Removes all None values from the reactants, agents, and products."""
        self.reactants = [m for m in self.reactants if m is not None]
        self.agents = [m for m in self.agents if m is not None]
        self.products = [m for m in self.products if m is not None]

    #
    # Static Methods
    #
    @staticmethod
    def remove_duplicates(
        seq: Iterable[T], key: Optional[Callable[[T], V]] = None
    ) -> List[T]:
        """Remove duplicates and preserve order.
        Adapted from https://stackoverflow.com/a/480227
        Args:
            seq: sequence to remove duplicates from.
            key: what to base duplicates on, must be hashable.
                Defaults to the elements of seq.
        """
        if key is None:

            def key(x: T) -> V:
                return x  # type: ignore

        key = cast(Callable[[T], V], key)  # necessary for mypy

        seen: Set[V] = set()
        seen_add = seen.add
        return [x for x in seq if not (key(x) in seen or seen_add(key(x)))]
