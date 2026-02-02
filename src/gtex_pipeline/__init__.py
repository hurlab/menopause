"""
GTEx Data Processing Pipeline

A Python-based pipeline for processing GTEx v10 gene expression data
with focus on adipose tissue analysis for menopause research.
"""

__version__ = "1.0.0"
__author__ = "GoosLab"

from gtex_pipeline.exceptions import (
    GtexFormatError,
    GtexIOError,
    GtexPipelineError,
    GtexProcessingError,
    GtexValidationError,
)
from gtex_pipeline.exporter import export_processed_data
from gtex_pipeline.filter import filter_tissue
from gtex_pipeline.loader import load_gtex_data
from gtex_pipeline.validator import ValidationResult, validate_expression_data

__all__ = [
    "load_gtex_data",
    "filter_tissue",
    "validate_expression_data",
    "export_processed_data",
    "ValidationResult",
    "GtexPipelineError",
    "GtexFormatError",
    "GtexValidationError",
    "GtexIOError",
    "GtexProcessingError",
]
