import argparse
import os

from capsid_net.pipeline import com
from capsid_net.pipeline import network
from capsid_net.pipeline import normals
from capsid_net.pipeline import pisa


def run_pipeline(args):
    os.makedirs(args.output, exist_ok=True)

    com_output = os.path.join(args.output, "center_of_mass.csv")

    com.run(argparse.Namespace(structure=args.structure, output=com_output))

    pisa.run_analysis(
        argparse.Namespace(
            interactions=args.pisa,
            classes=args.classes,
            output=args.output,
            rename_file=args.rename,
            grouping_file=args.grouping,
            group_capsomers=args.group_capsomers,
            group_zippers=args.group_zippers,
        )
    )

    normals.run_analysis(
        argparse.Namespace(
            com=com_output,
            prot_classes=args.classes,
            rename=args.rename,
            grouping=args.grouping,
            relatedness=os.path.join(args.output, "relatedness_matrix.csv"),
            custom_filter=args.custom_filter,
            output=args.output,
            node_colors=args.node_colors,
        )
    )

    network.run_analysis(
        argparse.Namespace(
            interaction=os.path.join(args.output, "interaction_matrix.csv"),
            com=com_output,
            grouping=args.grouping,
            rename=args.rename,
            classes=args.classes,
            output=args.output,
            use_saltbridges=args.use_saltbridges,
            pisa_mod_file=args.pisa_mod_file,
            tag_color_csv=args.tag_color_csv,
        )
    )


def build_parser():
    parser = argparse.ArgumentParser(prog="capsid-net", description="Capsid structure and PISA interaction analysis CLI.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    com.build_parser(subparsers.add_parser("com", help="Calculate chain centers of mass"))
    pisa.build_parser(subparsers.add_parser("pisa", help="Process PISA interaction tables"))
    normals.build_parser(subparsers.add_parser("normals", help="Compute capsomer centers, normals, and angle outputs"))
    network.build_parser(subparsers.add_parser("network", help="Render interaction network plots"))

    run_parser = subparsers.add_parser("run", help="Run the end-to-end pipeline")
    run_parser.add_argument("--structure", required=True, help="Path to the input structure file")
    run_parser.add_argument("--pisa", required=True, help="Path to the input PISA interactions CSV")
    run_parser.add_argument("--classes", required=True, help="Path to the protein classes CSV")
    run_parser.add_argument("--rename", required=True, help="Path to the rename CSV")
    run_parser.add_argument("--grouping", required=True, help="Path to the grouping CSV")
    run_parser.add_argument("--tag_color_csv", required=True, help="Path to the node tag/color CSV")
    run_parser.add_argument("--output", required=True, help="Output directory")
    run_parser.add_argument("--custom_filter", help="Optional group filter file for the normals stage")
    run_parser.add_argument("--node_colors", help="Optional node color CSV for the normals stage")
    run_parser.add_argument("--group_capsomers", action="store_true", help="Enable capsomer grouping in the PISA stage")
    run_parser.add_argument("--group_zippers", action="store_true", help="Enable zipper grouping in the PISA stage")
    run_parser.add_argument("--use_saltbridges", action="store_true", help="Use Dsb values for the network stage")
    run_parser.add_argument("--pisa_mod_file", help="Optional modified PISA CSV with Dsb values")
    run_parser.set_defaults(func=run_pipeline)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
