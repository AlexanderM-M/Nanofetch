class NanoFetchError(Exception):
    """Expected user-facing error."""


class AssemblyError(NanoFetchError):
    """The BAM assembly could not be identified or validated."""


class GeneNotFoundError(NanoFetchError):
    """A requested gene is absent from the selected annotation."""

