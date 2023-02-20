#!/usr/bin/env python3
import logging

def cofigure_logging(outdir):
    logging.basicConfig(
        level=logging.INFO, 
        format='%(asctime)s [%(levelname)s] %(message)s', 
        handlers=[
            logging.FileHandler(f"{outdir}/ribdif_log_file.log"), 
            logging.StreamHandler()
        ]
    )


