class ConformationEncodingError(Exception):
    """Exception raised for errors in the conformation encoding."""

    def __init__(self, message: str = "Invalid conformation encoding.") -> None:
        super().__init__(message)
