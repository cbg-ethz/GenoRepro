import csv
import sys
import os
import typer
from pathlib import Path
from typing import List



def get_pos_diff(csv_path) -> list:
    pos_diff_list = []
    with open(csv_path, 'r') as csv_file:
        reader = csv.reader(csv_file, delimiter=',')
        header = next(reader)
        pos_diff = []
        for read in reader:
            read_name = read[0]
            get_pos = lambda x: [i.strip() for i in x[1:-1].split(",")]
            pos_g, pos_rep = get_pos(read[1]), get_pos(read[2])
            if "None" in pos_g + pos_rep:
                continue
            pos_g, pos_rep = list(map(int, pos_g)), list(map(int, pos_rep))
            is_zero = lambda x, y: (x == 0 and y != 0) or (y == 0 and x != 0)
            if is_zero(pos_g[0], pos_rep[0]):
                continue
            if is_zero(pos_g[1], pos_rep[1]):
                continue
            count_diff = lambda x,y: max(abs(x[0] - y[0]), abs(x[1] - y[1]))
            pos_diff_list.append(count_diff(pos_g, pos_rep))
        return pos_diff_list


app = typer.Typer(add_completion=False)
@app.command()
def main(input_csvs: List[Path] = typer.Argument(
            ...,
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            show_default=False,
            help="List of CSV files to compare; min: 2 files"
            ),
         output_path: Path = typer.Option(
            "./comparer_output.csv",
            "--output",
            "-o",
            help="Path to output directory or file."
            )
            ):
    general = []
    for f in input_csvs:
        pos_diff_list = get_pos_diff(f)
        general += pos_diff_list
        with open(output_path, mode='w') as output_file:
            writer = csv.writer(output_file, delimiter=',')
            writer.writerow(general)


if __name__ == "__main__":
    app()
