"""Custom exception classes for GTEx pipeline."""


class GtexPipelineError(Exception):
    """
    Base exception class for all GTEx pipeline errors.

    All custom exceptions in the GTEx pipeline inherit from this base class,
    allowing for unified exception handling across the pipeline.

    Attributes:
        message: Human-readable error message
        context: Optional dictionary with additional error context
    """

    def __init__(self, message: str, **context):
        """
        Initialize base exception with message and optional context.

        Args:
            message: Human-readable error description
            **context: Additional context information (e.g., file_path, stage)
        """
        super().__init__(message)
        self.message = message
        self.context = context

    def __str__(self) -> str:
        """Return formatted error message with context."""
        if self.context:
            context_str = ", ".join(f"{k}={v}" for k, v in self.context.items())
            return f"{self.message} ({context_str})"
        return self.message


class GtexFormatError(GtexPipelineError):
    """
    Exception raised for GTEx file format validation failures.

    Raised when the input file does not conform to expected GTEx v10 format
    specifications, including missing headers, incorrect structure, or
    incompatible data types.
    """

    def __init__(self, message: str, file_path: str | None = None, **context):
        """
        Initialize format error with file path context.

        Args:
            message: Description of format issue
            file_path: Path to file that failed validation
            **context: Additional context information
        """
        if file_path:
            context["file_path"] = file_path
        super().__init__(message, **context)


class GtexValidationError(GtexPipelineError):
    """
    Exception raised for data quality validation failures.

    Raised when gene expression data fails quality checks, including
    excessive missing values, invalid TPM values, or duplicate identifiers.
    """

    def __init__(self, message: str, details: dict | None = None, **context):
        """
        Initialize validation error with details dictionary.

        Args:
            message: Description of validation failure
            details: Dictionary with validation details (e.g., missing_count)
            **context: Additional context information
        """
        if details:
            context["details"] = details
        super().__init__(message, **context)


class GtexIOError(GtexPipelineError):
    """
    Exception raised for file I/O operation failures.

    Raised when file read/write operations fail due to permissions,
    missing files, disk space issues, or other I/O-related problems.
    """

    def __init__(self, message: str, file_path: str | None = None, **context):
        """
        Initialize IO error with file path context.

        Args:
            message: Description of I/O error
            file_path: Path to file involved in I/O operation
            **context: Additional context information
        """
        if file_path:
            context["file_path"] = file_path
        super().__init__(message, **context)


class GtexProcessingError(GtexPipelineError):
    """
    Exception raised for general processing failures.

    Raised when data processing operations fail for reasons other than
    format validation or I/O errors, such as transformation failures,
    calculation errors, or unexpected processing states.
    """

    def __init__(self, message: str, stage: str | None = None, **context):
        """
        Initialize processing error with stage context.

        Args:
            message: Description of processing error
            stage: Processing stage where error occurred
            **context: Additional context information
        """
        if stage:
            context["stage"] = stage
        super().__init__(message, **context)
