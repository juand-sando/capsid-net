import argparse
import os

from capsid_net.pipeline import com
from capsid_net.pipeline import normals as angles
from capsid_net.pipeline import network
from capsid_net.pipeline import preprocess


def run_pipeline(args):
    os.makedirs(args.output, exist_ok=True)

    com_output = os.path.join(args.output, "center_of_mass.csv")
    preprocess_config = os.path.join(args.output, "preprocess_config.json")

    com.run(argparse.Namespace(structure=args.structure, parser="auto", output=com_output))

    preprocess.run_analysis(
        argparse.Namespace(
            interactions=args.interactions,
            com=com_output,
            prot_classes=args.classes,
            output=args.output,
            rename=args.rename,
            grouping=args.grouping,
        )
    )

    angles.run_analysis(
        argparse.Namespace(
            config=preprocess_config,
            exclude=args.angles_exclude,
            output=args.output,
            graph_format=args.graph_format,
            tag_color_csv=args.tag_color_csv,
        )
    )

    network.run_analysis(
        argparse.Namespace(
            config=preprocess_config,
            exclude=args.network_exclude,
            layout_config=args.network_layout_config,
            output=args.output,
            graph_format=args.graph_format,
            tag_color_csv=args.tag_color_csv,
        )
    )


def build_parser():
    parser = argparse.ArgumentParser(
        prog="capsid-net",
        description="CLI for preprocessing a PISA interaction file and analyzing capsid interaction networks.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    com.build_parser(subparsers.add_parser("com", help="Calculate chain centers of mass"))
    preprocess.build_parser(
        subparsers.add_parser("preprocess", help="Prepare processed interactions, grouped COMs, angles, and a downstream config")
    )
    angles.build_parser(subparsers.add_parser("angles", help="Render capsid angle plots from preprocessed data"))
    network.build_parser(subparsers.add_parser("network", help="Render interaction network plots"))

    run_parser = subparsers.add_parser("run", help="Run the end-to-end pipeline")
    run_parser.add_argument("--structure", required=True, help="Path to the input structure file")
    run_parser.add_argument("--interactions", required=True, help="Path to the input interactions CSV")
    run_parser.add_argument("--classes", required=True, help="Path to the protein classes CSV")
    run_parser.add_argument("--rename", required=True, help="Path to the rename CSV")
    run_parser.add_argument("--grouping", required=True, help="Path to the grouping CSV")
    run_parser.add_argument("--tag_color_csv", required=True, help="Path to the node tag/color CSV")
    run_parser.add_argument("--output", required=True, help="Output directory")
    run_parser.add_argument("--graph-format", choices=["svg", "png"], default="svg", help="Output format for generated graph figures")
    run_parser.add_argument("--angles-exclude", help="Optional exclusion file for the angles plot")
    run_parser.add_argument("--network-exclude", help="Optional exclusion file for the network plot")
    run_parser.add_argument("--network-layout-config", help="Optional YAML file with network layout adjustment rules")
    run_parser.set_defaults(func=run_pipeline)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
