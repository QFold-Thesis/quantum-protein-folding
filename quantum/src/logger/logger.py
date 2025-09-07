import logging

from constants import GLOBAL_LOGGER_NAME


def get_logger(name: str = GLOBAL_LOGGER_NAME) -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.hasHandlers():
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            "%(asctime)-19s [%(levelname)-5s] | %(module)-20s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.DEBUG)
    return logger
