"""Logging configuration for GTEx pipeline."""

import logging
import sys
from pathlib import Path

# Format with timestamp, level, and message
LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def setup_logger(
    name: str,
    level: int = logging.INFO,
    log_file: str | None = None,
    format_string: str = LOG_FORMAT,
    date_format: str = DATE_FORMAT,
) -> logging.Logger:
    """
    Set up and configure a logger with console and optional file handlers.

    Creates a logger with standardized formatting for the GTEx pipeline.
    By default, logs to console. If log_file is provided, also logs to file.

    Args:
        name: Logger name (typically module or pipeline name)
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_file: Optional path to log file for persistent logging
        format_string: Custom format string for log messages
        date_format: Custom date format for timestamps

    Returns:
        Configured logger instance

    Examples:
        >>> logger = setup_logger("gtex_pipeline")
        >>> logger.info("Pipeline started")

        >>> logger = setup_logger("pipeline", log_file="pipeline.log", level=logging.DEBUG)
        >>> logger.debug("Detailed debug info")
    """
    # Get or create logger
    logger = logging.getLogger(name)

    # Clear existing handlers to avoid duplicates
    logger.handlers.clear()

    # Set log level
    logger.setLevel(level)

    # Create formatter
    formatter = logging.Formatter(format_string, datefmt=date_format)

    # Add console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # Add file handler if log_file specified
    if log_file:
        # Ensure directory exists
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # Prevent propagation to avoid duplicate logs
    logger.propagate = False

    return logger


def get_logger(name: str) -> logging.Logger:
    """
    Get an existing logger or create a default one.

    This is a convenience function for retrieving loggers without
    full configuration. If the logger doesn't exist, it creates
    one with default settings.

    Args:
        name: Logger name

    Returns:
        Logger instance (existing or newly created)

    Examples:
        >>> logger = get_logger("gtex_pipeline")
        >>> logger.info("Using existing or default logger")
    """
    logger = logging.getLogger(name)

    # If logger has no handlers, set up with defaults
    if not logger.handlers:
        return setup_logger(name)

    return logger
