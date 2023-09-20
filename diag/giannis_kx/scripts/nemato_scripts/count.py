import os

"""Gives you the order of the poinc_t=****.bin files. Useful if you only want to process with Nemato a single file."""

path = os.getcwd()
poinc_files = [f for f in os.listdir(path) if f.startswith("poinc_t=")]
def timestamp(filestring):
    key = filestring.split(".")[0][8:]
    return int(key)

sfiles = sorted(poinc_files,key=timestamp)

