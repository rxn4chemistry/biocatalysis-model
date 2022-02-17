"""Tokenizer for enzymatic reactions."""
import re

SMILES_TOKENIZER_PATTERN = r"(\%\([0-9]{3}\)|\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\||\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>>?|\*|\$|\%[0-9]{2}|[0-9])"
SMILES_REGEX = re.compile(SMILES_TOKENIZER_PATTERN)


def tokenize_enzymatic_reaction_smiles(rxn: str, keep_pipe=False) -> str:
    """Tokenize an enzymatic reaction SMILES in the form precursors|EC>>products.

    Args:
        rxn: an enzymatic reaction SMILES.
        keep_pipe: whether or not to keep the pipe separating the precursors from the EC as a token.
            Defaults to False.

    Returns:
        the tokenized enzymatic reaction SMILES.
    """
    parts = re.split(r">|\|", rxn)
    ec = parts[1].split(".")

    rxn = rxn.replace(f"|{parts[1]}", "")
    tokens = [token for token in SMILES_REGEX.findall(rxn)]
    arrow_index = tokens.index(">>")

    levels = ["v", "u", "t", "q"]

    if ec[0] != "":
        ec_tokens = [f"[{levels[i]}{e}]" for i, e in enumerate(ec)]
        if keep_pipe:
            ec_tokens.insert(0, "|")
        tokens[arrow_index:arrow_index] = ec_tokens

    return " ".join(tokens)


def detokenize_enzymatic_reaction_smiles(rxn: str) -> str:
    """Detokenize an enzymatic reaction SMILES in the form precursors|EC>>products.

    Args:
        rxn: a tokenized enzymatic reaction SMILES.

    Returns:
        the detokenized enzymatic reaction SMILES.
    """
    rxn = rxn.replace(" ", "")

    if "[v" in rxn and "|" not in rxn:
        pipe_index = rxn.index("[v")
        if pipe_index > -1:
            rxn = rxn[:pipe_index] + "|" + rxn[pipe_index:]

    if "[v" not in rxn and "|" in rxn:
        rxn = rxn.replace("|", "")

    if "|" not in rxn:
        return rxn

    precursor_split = rxn.split("|")

    if len(precursor_split) < 2:
        return ""

    reaction_split = precursor_split[1].split(">>")

    if len(reaction_split) < 2:
        return ""

    ec = (
        reaction_split[0]
        .replace("][", ".")
        .replace("[v", "")
        .replace("u", "")
        .replace("t", "")
        .replace("q", "")
        .replace("]", "")
    )

    return precursor_split[0] + "|" + ec + ">>" + reaction_split[1]


def tokenize_smiles(smiles: str) -> str:
    """
    Tokenize a SMILES molecule or reaction, and join the tokens with spaces.
    Args:
        smiles: SMILES string to tokenize, for instance 'CC(CO)=N>>CC(C=O)N'.
    Returns:
        SMILES string after tokenization, for instance 'C C ( C O ) = N >> C C ( C = O ) N'.
    """

    tokens = [token for token in SMILES_REGEX.findall(smiles)]
    return " ".join(tokens)
