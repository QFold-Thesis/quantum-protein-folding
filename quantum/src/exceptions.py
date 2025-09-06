from logger import get_logger

logger = get_logger()


class ConformationEncodingError(Exception):
    """Exception raised for errors in the conformation encoding."""

    def __init__(
        self,
        message: str = "Invalid conformation encoding. Make sure that the encoding is of type ConformationEncoding.",
    ) -> None:
        logger.error(message)
        super().__init__(message)


class ChainLengthError(Exception):
    """Exception raised for errors in the chain length."""

    def __init__(
        self,
        message: str = "Invalid chain length. Make sure that the chain length is correct.",
    ) -> None:
        logger.error(message)
        super().__init__(message)


class InvalidAminoAcidError(Exception):
    """Exception raised for errors in the amino acid."""

    def __init__(
        self,
        message: str = "Invalid amino acid symbol. Make sure that the amino acid symbol is compatible with the solution.",
    ) -> None:
        logger.error(message)
        super().__init__(message)
