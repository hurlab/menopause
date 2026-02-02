"""Tests for custom exception classes."""

from gtex_pipeline.exceptions import (
    GtexFormatError,
    GtexIOError,
    GtexPipelineError,
    GtexProcessingError,
    GtexValidationError,
)


class TestGtexPipelineError:
    """Tests for base GtexPipelineError exception."""

    def test_base_exception_creation(self):
        """Test that base exception can be created with message."""
        error = GtexPipelineError("Test error message")
        assert str(error) == "Test error message"
        assert isinstance(error, Exception)

    def test_base_exception_inheritance(self):
        """Test that all custom exceptions inherit from base."""
        assert issubclass(GtexFormatError, GtexPipelineError)
        assert issubclass(GtexValidationError, GtexPipelineError)
        assert issubclass(GtexIOError, GtexPipelineError)
        assert issubclass(GtexProcessingError, GtexPipelineError)


class TestGtexFormatError:
    """Tests for GtexFormatError exception."""

    def test_format_error_creation(self):
        """Test that format error can be created with message."""
        error = GtexFormatError("Invalid file format")
        assert str(error) == "Invalid file format"
        assert isinstance(error, GtexPipelineError)

    def test_format_error_with_context(self):
        """Test format error with file path context."""
        error = GtexFormatError("Invalid format", file_path="/path/to/file.txt")
        assert "Invalid format" in str(error)
        assert "/path/to/file.txt" in str(error)


class TestGtexValidationError:
    """Tests for GtexValidationError exception."""

    def test_validation_error_creation(self):
        """Test that validation error can be created with message."""
        error = GtexValidationError("Data validation failed")
        assert str(error) == "Data validation failed"
        assert isinstance(error, GtexPipelineError)

    def test_validation_error_with_details(self):
        """Test validation error with validation details."""
        details = {"missing_values": 100, "invalid_genes": ["GENE1", "GENE2"]}
        error = GtexValidationError("Validation failed", details=details)
        assert "Validation failed" in str(error)
        assert "100" in str(error)


class TestGtexIOError:
    """Tests for GtexIOError exception."""

    def test_io_error_creation(self):
        """Test that IO error can be created with message."""
        error = GtexIOError("File not found")
        assert str(error) == "File not found"
        assert isinstance(error, GtexPipelineError)

    def test_io_error_with_file_path(self):
        """Test IO error with file path context."""
        error = GtexIOError("Cannot read file", file_path="/data/input.txt")
        assert "Cannot read file" in str(error)
        assert "/data/input.txt" in str(error)


class TestGtexProcessingError:
    """Tests for GtexProcessingError exception."""

    def test_processing_error_creation(self):
        """Test that processing error can be created with message."""
        error = GtexProcessingError("Processing failed")
        assert str(error) == "Processing failed"
        assert isinstance(error, GtexPipelineError)

    def test_processing_error_with_stage(self):
        """Test processing error with processing stage context."""
        error = GtexProcessingError("Error during processing", stage="tissue_filtering")
        assert "Error during processing" in str(error)
        assert "tissue_filtering" in str(error)
