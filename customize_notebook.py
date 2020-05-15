#!/opt/hall-lab/python-3.7.0/bin/python

import json
import click

template = "VariantCalling.ipynb"

def file_paths(tsv):
    input_files = {}
    with open(tsv, 'r') as f:
        for line in f:
            (variable, file_path) = line.rstrip().split("\t")
            input_files[variable] = file_path
    return input_files

def format_file_path(variable, file_path):
    msg = '{} <- paste(dir, "{}", sep="/")\n'
    return msg.format(variable, file_path)

def update_notebook(template, tsv):
    with open(template, 'r') as f:
        data = json.load(f)
    input_files = file_paths(tsv)
    data['cells'][2]['source'] = [format_file_path(k, input_files[k]) for k in input_files.keys()]
    with open("new.ipynb", 'w') as f2:
        json.dump(data, f2, indent=4)

@click.command()
@click.option('--tsv', required=True, type=click.Path(exists=True), help='Map of input files for jupyter notebook')
def main(tsv):
    update_notebook(template, tsv)

if __name__=="__main__":
    main()
