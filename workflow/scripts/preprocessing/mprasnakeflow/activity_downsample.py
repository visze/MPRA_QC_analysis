import click
import numpy as np
import pandas as pd
from pathlib import Path
from mpralib.mpradata import MPRABarcodeData, MPRAOligoData, CountSampling
from mpralib.utils.io import export_barcode_file


FRACTIONS = [i / 10 for i in range(1, 10)]


@click.command()
@click.option(
    "--reporter-experiment-barcode",
    "reporter_experiment_barcode_file",
    type=click.Path(exists=True),
    required=True,
    help="Input TSV file path",
)
@click.option("--output-folder", type=click.Path(dir_okay=True), required=True, help="Output folder path")
@click.option("--seed", type=int, default=42, help="Random seed for reproducibility (default: 42)")
def activity_downsample(reporter_experiment_barcode_file: str, output_folder: str, seed: int):
    """Create downsampled ratio_df for each fraction using reporter experiment barcode data and save as CSV files."""
    np.random.seed(seed)

    mpradata_bc = MPRABarcodeData.from_file(reporter_experiment_barcode_file)

    for fraction in FRACTIONS[:-1]:
        click.echo(f"Processing fraction {fraction:.1f}...")
        mpradata_bc.apply_count_sampling(CountSampling.RNA_AND_DNA, proportion=fraction)
        create_output_file(output_folder, fraction, mpradata_bc)
        mpradata_bc.drop_count_sampling()

    fraction = FRACTIONS[-1]  # Use the last fraction to initialize the sampling
    create_output_file(output_folder, fraction, mpradata_bc)


def create_output_file(output_folder: str, fraction: float, mpradata: MPRABarcodeData) -> None:
    output_path = Path(output_folder) / f"reporter_barcode_{fraction:.1f}.tsv.gz"
    export_barcode_file(mpradata, str(output_path))
    click.echo(f"Saved reporter barcode for fraction {fraction:.1f} to {output_path}")


if __name__ == "__main__":
    activity_downsample()
