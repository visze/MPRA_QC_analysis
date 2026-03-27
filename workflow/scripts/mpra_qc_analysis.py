from .association_analysis import association
from .activity_analysis import activity
import click


@click.group(help="Command line interface of MPRA QC analysis.")
def main() -> None:
    pass


main.add_command(association)
main.add_command(activity)


if __name__ == "__main__":
    main()
