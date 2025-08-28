class ConformationEncodingError(Exception):
    """Exception raised for errors in the conformation encoding."""

    def __init__(self, message: str = "Invalid conformation encoding. Make sure that the encoding is ConformationEncoding.") -> None:
        super().__init__(message)
