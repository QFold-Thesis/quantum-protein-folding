import logging
from datetime import datetime
from typing import TYPE_CHECKING

from constants import (
    DEFAULT_TIMEZONE,
    GLOBAL_LOGGER_NAME,
    LOGGER_DEFAULT_LEVEL,
    LOGS_DIRPATH,
)

if TYPE_CHECKING:
    from pathlib import Path


def get_logger(name: str = GLOBAL_LOGGER_NAME) -> logging.Logger:
    """
    Creates and configures a logger instance.

    Args:
        name (str): Name of the logger. Defaults to `GLOBAL_LOGGER_NAME`.

    Returns:
        logging.Logger: Configured logger instance.

    """
    logger = logging.getLogger(name)
    if not logger.hasHandlers():
        timestamp: str = datetime.now(tz=DEFAULT_TIMEZONE).strftime(
            "%Y_%m_%d__%H_%M_%S"
        )

        log_handlers: list[logging.Handler] = []

        stream_handler: logging.Handler = logging.StreamHandler()
        log_handlers.append(stream_handler)

        log_filepath: Path = LOGS_DIRPATH / f"{timestamp}.log"
        LOGS_DIRPATH.mkdir(parents=True, exist_ok=True)

        file_handler: logging.Handler = logging.FileHandler(log_filepath)
        log_handlers.append(file_handler)

        formatter = logging.Formatter(
            "%(asctime)-19s [%(levelname)-5s] | %(module)-20s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

        for handler in log_handlers:
            handler.setFormatter(formatter)
            logger.addHandler(handler)

        logger.setLevel(LOGGER_DEFAULT_LEVEL)
        logger.info("Logger and handlers initialized")
    return logger
