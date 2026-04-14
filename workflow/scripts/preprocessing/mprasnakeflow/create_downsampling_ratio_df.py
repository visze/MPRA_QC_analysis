import click
import numpy as np
import pandas as pd
from pathlib import Path
from mpralib.mpradata import MPRABarcodeData, MPRAOligoData, CountSampling


FRACTIONS = [i / 10 for i in range(1, 11)]

@click.command()
@click.option("--reporter-experiment-barcode", "reporter_experiment_barcode_file", type=click.Path(exists=True), required=True, help="Input CSV file path")
@click.option("--output-folder", type=click.Path(dir_okay=True), required=True, help="Output folder path")
@click.option("--seed", type=int, default=42, help="Random seed for reproducibility (default: 42)")
def create_downsampling_ratio_df(reporter_experiment_barcode_file: str, output_folder: str, seed: int):
    """Create downsampled ratio_df for each fraction using reporter experiment barcode data and save as CSV files."""
    np.random.seed(seed)

    mpradata_bc = MPRABarcodeData.from_file(reporter_experiment_barcode_file)

    for fraction in FRACTIONS[:-1]:
        click.echo(f"Processing fraction {fraction:.1f}...")
        mpradata_bc.apply_count_sampling(CountSampling.RNA_AND_DNA, proportion=fraction)
        create_output_file(output_folder, fraction, mpradata_bc.oligo_data)
        mpradata_bc.drop_count_sampling()
    
    fraction = FRACTIONS[-1]  # Use the last fraction to initialize the sampling
    create_output_file(output_folder, fraction, mpradata_bc.oligo_data)
    

def create_output_file(output_folder: str, fraction: float, mpradata_oligo: MPRAOligoData) -> None:
    ratio_df = pd.DataFrame(mpradata_oligo.activity.T, columns=[f"RNA_DNA_ratio_log_rep{rep+1}" for rep in range(mpradata_oligo.n_obs)])
    rna_counts_df = pd.DataFrame(mpradata_oligo.rna_counts.T, columns=[f"RNA_rep{rep+1}" for rep in range(mpradata_oligo.n_obs)])
    dna_counts_df = pd.DataFrame(mpradata_oligo.dna_counts.T, columns=[f"DNA_rep{rep+1}" for rep in range(mpradata_oligo.n_obs)])
    ratio_df["cCRE"] = mpradata_oligo.oligos.values
    df = pd.concat([ratio_df, rna_counts_df, dna_counts_df], axis=1)
    df.set_index("cCRE", inplace=True)
    output_path = Path(output_folder) / f"ratio_df_{fraction:.1f}.csv.gz"
    df.to_csv(output_path, index=True, compression="gzip")
    click.echo(f"Saved activity per replicate for fraction {fraction:.1f} to {output_path}")


if __name__ == "__main__":
    create_downsampling_ratio_df()
