"""Generic utilities."""


def disable_rdkit_logging() -> None:
    """Disables RDKit whiny logging."""
    import rdkit.rdBase as rkrb
    import rdkit.RDLogger as rkl

    logger = rkl.logger()
    logger.setLevel(rkl.ERROR)
    rkrb.DisableLog("rdApp.error")
