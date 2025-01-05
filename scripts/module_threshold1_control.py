# Weier Guo Python
# Create: 01/03/2025
# INTRODUCTION: set all control files chireads to a chosen threshold.

## PACKAGES ##
import sys, math, os, time
import argparse 




## CALL ARGUMENTS ##
def parse_arguments():
    parser = argparse.ArgumentParser(description="Threshold1: set controls threshold.")
    parser.add_argument("--input", type=str, help="Input chiread file.")
    parser.add_argument("--control", type=str, nargs="*", help="Input control names.")
    parser.add_argument("--number", type=int, default=0, help="Choose control threshold (default=0).")
    parser.add_argument("--t1file", type=str, help="Output file with controls=chosen threshold.")
    return parser.parse_args()