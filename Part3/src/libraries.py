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

from Bio import SeqIO
from typing import Optional, Dict, List, Tuple
from collections import defaultdict
from pathlib import Path