"""Tests for logging configuration module."""

import logging

from gtex_pipeline.utils.logging_config import get_logger, setup_logger


class TestSetupLogger:
    """Tests for setup_logger function."""

    def test_setup_logger_creates_logger(self):
        """Test that setup_logger creates a valid logger."""
        logger = setup_logger("test_logger")
        assert logger is not None
        assert logger.name == "test_logger"
        assert isinstance(logger, logging.Logger)

    def test_setup_logger_with_file_handler(self, tmp_path):
        """Test that setup_logger can add file handler."""
        log_file = tmp_path / "test.log"
        logger = setup_logger("file_logger", log_file=str(log_file))
        assert logger is not None

        # Write a test log message
        logger.info("Test message")

        # Verify log file was created
        assert log_file.exists()
        content = log_file.read_text()
        assert "Test message" in content

    def test_setup_logger_with_level(self):
        """Test that setup_logger sets correct log level."""
        logger_debug = setup_logger("debug_logger", level=logging.DEBUG)
        assert logger_debug.level == logging.DEBUG

        logger_info = setup_logger("info_logger", level=logging.INFO)
        assert logger_info.level == logging.INFO

    def test_setup_logger_formats_messages(self, tmp_path, caplog):
        """Test that log messages are formatted with timestamp."""
        import re

        log_file = tmp_path / "formatted.log"
        logger = setup_logger("format_test", log_file=str(log_file))

        logger.info("Formatted message")

        content = log_file.read_text()
        # Check for timestamp pattern (YYYY-MM-DD HH:MM:SS)
        assert re.search(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}", content)
        assert "Formatted message" in content


class TestGetLogger:
    """Tests for get_logger function."""

    def test_get_logger_returns_existing_logger(self):
        """Test that get_logger returns previously created logger."""
        # Create logger first
        setup_logger("shared_logger")

        # Get same logger
        logger = get_logger("shared_logger")
        assert logger is not None
        assert logger.name == "shared_logger"

    def test_get_logger_creates_new_if_not_exists(self):
        """Test that get_logger creates new logger if needed."""
        logger = get_logger("new_logger")
        assert logger is not None
        assert logger.name == "new_logger"

    def test_get_logger_returns_consistent_instance(self):
        """Test that get_logger returns same instance for same name."""
        logger1 = get_logger("consistent_logger")
        logger2 = get_logger("consistent_logger")
        assert logger1 is logger2


class TestLoggerIntegration:
    """Integration tests for logging functionality."""

    def test_logger_outputs_to_console_by_default(self, capsys):
        """Test that logger outputs to console by default."""
        logger = setup_logger("console_test")

        logger.info("Console message")

        captured = capsys.readouterr()
        assert "Console message" in captured.out

    def test_logger_handles_different_log_levels(self, tmp_path):
        """Test that logger handles DEBUG, INFO, WARNING, ERROR levels."""
        log_file = tmp_path / "levels.log"
        logger = setup_logger("levels_test", log_file=str(log_file), level=logging.DEBUG)

        logger.debug("Debug message")
        logger.info("Info message")
        logger.warning("Warning message")
        logger.error("Error message")

        content = log_file.read_text()
        assert "Debug message" in content
        assert "Info message" in content
        assert "Warning message" in content
        assert "Error message" in content

    def test_logger_in_pipeline_context(self, tmp_path):
        """Test logger usage in simulated pipeline context."""
        log_file = tmp_path / "pipeline.log"
        logger = setup_logger("gtex_pipeline", log_file=str(log_file))

        logger.info("Pipeline started")
        logger.info("Loading GTEx data")
        logger.warning("Found 5 missing values")
        logger.info("Data processing complete")
        logger.info("Pipeline finished successfully")

        content = log_file.read_text()
        assert "Pipeline started" in content
        assert "Loading GTEx data" in content
        assert "Found 5 missing values" in content
        assert "Data processing complete" in content
        assert "Pipeline finished successfully" in content
