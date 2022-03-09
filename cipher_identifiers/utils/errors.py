class CompoundNotFoundError(Exception):
    pass


class SubstanceNotFoundError(Exception):
    pass


class ExperimentNotFoundError(Exception):
    pass


class ModelNotFoundError(Exception):
    pass


class BiomoleculeNotFoundError(Exception):
    pass


class ExternalDatabaseNotFoundError(Exception):
    pass


class BindingSiteNotFoundError(Exception):
    pass


class InvalidRequestError(Exception):
    pass


class InvalidSMILESError(Exception):
    pass


class CompoundAlreadyExistsError(Exception):
    pass
