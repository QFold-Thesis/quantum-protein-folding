"""
Defines custom exceptions for protein folding simulations,
including errors for conformation encoding, chain length, and amino acid validity.
"""

from logger import get_logger

logger = get_logger()


class ConformationEncodingError(Exception):
    """Exception raised for errors in the conformation encoding."""

    def __init__(
        self,
        message: str = "Invalid conformation encoding. Make sure that the encoding is of type ConformationEncoding.",
    ) -> None:
        logger.exception(message)
        super().__init__(message)


class ChainLengthError(Exception):
    """Exception raised for errors in the chain length."""

    def __init__(
        self,
        message: str = "Invalid chain length. Make sure that the chain length is correct.",
    ) -> None:
        logger.exception(message)
        super().__init__(message)


class UnsupportedAminoAcidSymbolError(Exception):
    """Exception raised when amino acid symbol is not included in interaction."""

    def __init__(
        self,
        message: str = "Invalid amino acid symbol. Make sure that the amino acid symbol is compatible with the solution.",
    ) -> None:
        logger.exception(message)
        super().__init__(message)


class InvalidOperatorError(Exception):
    """Exception raised when an invalid quantum operator is encountered."""

    def __init__(
        self,
        message: str = "Invalid quantum operator encountered. Make sure that the operator is valid.",
    ) -> None:
        logger.exception(message)
        super().__init__(message)


class InvalidInteractionTypeError(Exception):
    """Exception raised when an invalid interaction type is encountered."""

    def __init__(
        self,
        message: str = "Invalid interaction type. Make sure that the interaction type is of type InteractionType.",
    ) -> None:
        logger.exception(message)
        super().__init__(message)
