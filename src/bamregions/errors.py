class BamregionsError(Exception):
    """Expected user-facing error."""


class AssemblyError(BamregionsError):
    """The BAM assembly could not be identified or validated."""


class GeneNotFoundError(BamregionsError):
    """A requested gene is absent from the selected annotation."""

