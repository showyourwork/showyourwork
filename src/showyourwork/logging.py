import logging
from typing import Any, Dict, Optional

import snakemake

from showyourwork import paths


def get_logger(
    config: Optional[Dict[str, Any]] = None, **stream_kwargs: Any
) -> logging.Logger:
    logger = logging.getLogger("showyourwork")

    if not logger.handlers:
        # Root level
        logger.setLevel(logging.DEBUG)

        # All messages are emitted to the terminal
        stream_handler = snakemake.logging.ColorizingStreamHandler(**stream_kwargs)
        stream_handler.setLevel(logging.INFO)
        logger.addHandler(stream_handler)

        # If this is called from within a Snakemake workflow, then we'll get the
        # config object and we can write the logs to a file.
        if config is not None:
            log_file = paths.work(config).logs / "showyourwork.log"
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(logging.DEBUG)
            logger.addHandler(file_handler)

    return logger


def patch_snakemake_logging(config: Dict[str, Any]) -> None:
    verbose = config.get("verbose", False)
    logger = get_logger(config)
    snakemake_logger = snakemake.logging.logger

    for handler in snakemake_logger.logger.handlers:
        if isinstance(handler, logging.FileHandler):
            handler.setLevel(logging.DEBUG)
        else:
            if not verbose:
                handler.setLevel(logging.CRITICAL)

    # Custom Snakemake stdout handler
    if not hasattr(snakemake_logger, "custom_stream_handler"):
        snakemake_logger.custom_stream_handler = (
            snakemake.logging.ColorizingStreamHandler()
        )
        snakemake_logger.custom_stream_handler.setLevel(logging.ERROR)
        snakemake_logger.logger.addHandler(snakemake_logger.custom_stream_handler)

    # Custom Snakemake file handler
    if not hasattr(snakemake_logger, "custom_file_handler"):
        log_file = paths.work(config).logs / "snakemake.log"
        snakemake_logger.custom_file_handler = logging.FileHandler(log_file)
        snakemake_logger.custom_file_handler.setLevel(logging.DEBUG)
        snakemake_logger.logger.addHandler(snakemake_logger.custom_file_handler)

    # Emit any messages in "job_info" packets. These are logging strings that
    # are provided by the "message" entries in Snakemake rules.
    def custom_handler(msg: Dict[str, Any]) -> None:
        if snakemake_logger.is_quiet_about("rules"):
            return
        text = msg.get("msg", None)
        if text is None:
            return
        level = msg["level"]
        if level == "job_info":
            logger.info(text)
        elif level == "info" and text.startswith("Nothing to be done"):
            logger.info(text)

    if not verbose:
        snakemake_logger.log_handler.append(custom_handler)

    # Allow all conda messages to come through
    snakemake.deployment.conda.logger = logger
