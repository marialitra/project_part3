import os
import sys
import argparse
import re
import subprocess
import csv
import time

import torch
import esm
import numpy as np

# Import from Python existing libraries
from Bio import SeqIO
from typing import Optional, Dict, List, Tuple, Any
from collections import defaultdict
from pathlib import Path

# Import from user-defined files
from parse_files import parse_args_embed, parse_args_blast, parse_args_search, parse_neighbor_results
from parse_files import parse_blast_tsv, parse_ann_txt
from utils import Hit, filter_hits, write_top_hits, load_model, load_sequences, save_output
from utils import remap_output_ids, print_recall, compute_recall, print_QPS_dictionary, print_QPS, delete_files
from run_methods import run_method
from run_blast_methods import running_blast
from generate_reports import generate_per_query_report, generate_all_methods_report