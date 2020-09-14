#!/opt/hall-lab/python-3.7.0/bin/python

import json
import click

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

def update_notebook(template, tsv, output):
    with open(template, 'r') as f:
        data = json.load(f)
    input_files = file_paths(tsv)
    file_path_cells = [cell for cell in data['cells'] if 'metadata' in cell and 'tags' in cell['metadata'] and 'file_paths' in cell['metadata']['tags']]
    file_path_cells[0]['source'] = [format_file_path(k, input_files[k]) for k in input_files.keys()]
    with open(output, 'w') as f2:
        json.dump(data, f2, indent=4)

@click.command()
@click.option('--tsv', required=True, type=click.Path(exists=True), help='Map of input files for jupyter notebook')
@click.option('--template', required=True, type=click.Path(exists=True), help='Template for jupyter notebook')
@click.option('--output', required=True, type=click.Path(exists=False), help='Path to output notebook file')

def main(tsv, template, output):
    update_notebook(template, tsv, output)

if __name__=="__main__":
    main()
