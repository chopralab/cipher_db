class InvalidRetrosyntheticTreeError(ValueError):
    pass


class InvalidChemicalNodeError(InvalidRetrosyntheticTreeError):
    pass


class InvalidReactionNodeError(InvalidRetrosyntheticTreeError):
    pass
