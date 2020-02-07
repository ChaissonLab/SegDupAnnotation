#!/usr/bin/env bash
import sys

for line in sys.stdin:
    if line[0] == ">":
        name=line[1:]
    else:
        a=line.
