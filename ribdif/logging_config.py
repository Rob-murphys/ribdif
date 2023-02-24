#!/usr/bin/env python3
import logging
from pathlib import Path

def cofigure_logging():
    
    class TracebackFilter(logging.Filter):
        def filter(self, record):
            return not record.exc_info  # exclude records with exc_info set
    
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    
    # create console handler with filter
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter('%(message)s'))
    console_handler.addFilter(TracebackFilter())
    logger.addHandler(console_handler)
    
    # create file handler without filter
    file_handler = logging.FileHandler("tmpfile.log", "w")
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(file_handler)
    return logger

def replace_log_file(outdir, logger):
    new_log_file = f"{outdir}/ribdif_log_file.log"
    Path.unlink(Path(new_log_file), missing_ok = True)
    file_handler = logging.FileHandler(new_log_file, "a")
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(file_handler)
    
    old_log_file = logger.handlers[1].baseFilename
    with open(logger.handlers[1].baseFilename, "r") as log_in, open(logger.handlers[2].baseFilename, "a") as log_out:
        log_out.write(log_in.read())
    logger.removeHandler(logger.handlers[1])
    Path.unlink(Path(old_log_file))
    return

    
        
    