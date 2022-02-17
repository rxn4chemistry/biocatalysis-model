#!/usr/bin/env python
from typing import List, DefaultDict, Tuple, Dict
from collections import defaultdict
from pathlib import Path

import click
import pandas as pd
from tqdm import trange

from rxn_biocatalysis_tools import (
    EnzymaticReaction,
    detokenize_enzymatic_reaction_smiles,
)


def get_ec_prediction_accuracy(
    ground_truth: List[EnzymaticReaction],
    predictions: List[List[EnzymaticReaction]],
    pred_level: int = 3,
    group_by_level: int = 1,
    top_n: int = 1,
) -> Tuple[Dict[str, float], List[str], List[Tuple[str, str]]]:
    """Returns the prediction accuracy of EC classes at a given level.

    Args:
        ground_truth: The enzymatic reactions representing the ground truth.
        predictions: The enzymatic reactions representing the predictions.
        pred_level: The predicted EC level. Defaults to 3.
        group_by_level: The EC level by which to group and evaluate correctness (up to pred_level). Defaults to 1.
        top_n: The number of predictions per ground truth reaction to consider. Defaults to 1.

    Returns:
        Per-group and overall accuracies.
    """
    n_correct: DefaultDict[str, int] = defaultdict(int)
    n_total: DefaultDict[str, int] = defaultdict(int)

    n_correct_total = 0

    incorrect = []
    correct = []
    for gt, preds in zip(ground_truth, predictions):
        pred_ecs = [pred.get_ec(pred_level) for pred in preds[:top_n]]
        if gt.get_ec(pred_level) in pred_ecs:
            n_correct[gt.get_ec(group_by_level)] += 1
            n_correct_total += 1
            correct.append(gt.get_ec(pred_level))
        else:
            incorrect.append((gt.get_ec(pred_level), pred_ecs[0]))

        n_total[gt.get_ec(group_by_level)] += 1

    results = {}
    for key, value in n_total.items():
        results[key] = n_correct[key] / value

    results["overall"] = n_correct_total / len(ground_truth)

    return results, correct, incorrect


def get_accuracies(
    ground_truth: List[EnzymaticReaction],
    predictions: List[List[EnzymaticReaction]],
    group_by_level: int = 1,
    top_n: int = 1,
) -> Tuple[Dict[str, float], List[str], List[Tuple[str, str]]]:
    """Returns the reaction (plus EC number for backward) prediction accuracies.

    Args:
        ground_truth: The enzymatic reactions representing the ground truth.
        predictions: The enzymatic reactions representing the predictions.
        group_by_level: The EC level by which to group and evaluate correctness. Defaults to 1.
        top_n: The number of predictions per ground truth reaction to consider. Defaults to 1.

    Returns:
        Per-group and overall accuracies.
    """
    n_correct: DefaultDict[str, int] = defaultdict(int)
    n_total: DefaultDict[str, int] = defaultdict(int)
    n_correct_total = 0
    correct = []
    incorrect = []

    for gt, preds in zip(ground_truth, predictions):
        pred_rxns = [pred.to_string() for pred in preds[:top_n]]
        if gt.to_string() in pred_rxns:
            n_correct[gt.get_ec(group_by_level)] += 1
            n_correct_total += 1
            correct.append(gt.to_string())
        else:
            incorrect.append((gt.to_string(), pred_rxns[0]))

        n_total[gt.get_ec(group_by_level)] += 1

    results = {}
    for key, value in n_total.items():
        results[key] = n_correct[key] / value

    results["overall"] = n_correct_total / len(ground_truth)

    return results, correct, incorrect


@click.command()
@click.argument("project_path", type=click.Path(exists=True))
@click.option(
    "n_best_fw",
    "--n-best-fw",
    default=5,
    show_default=True,
    help="The number (n) of calculated tgt predictions per src.",
)
@click.option(
    "n_best_bw",
    "--n-best-bw",
    default=10,
    show_default=True,
    help="The number (n) of calculated src predictions per tgt.",
)
@click.option(
    "n_best_rtr",
    "--n-best-rtr",
    default=1,
    show_default=True,
    help="The number (n) of calculated (roundtrip) tgt predictions per predicted src.",
)
@click.option(
    "top_n_fw",
    "--top-n-fw",
    multiple=True,
    type=int,
    default=[1],
    show_default=True,
    help="The number (n) of forward predictions to consider in the evaluation.",
)
@click.option(
    "top_n_bw",
    "--top-n-bw",
    multiple=True,
    type=int,
    default=[1],
    show_default=True,
    help="The number (n) of backward predictions to consider in the evaluation.",
)
@click.option(
    "top_n_rtr",
    "--top-n-rtr",
    multiple=True,
    type=int,
    default=[1],
    show_default=True,
    help="The number (n) of roundtrip predictions to consider in the evaluation.",
)
@click.option(
    "name",
    "--name",
    help="Writes the result to a csv file with this name.",
)
@click.option(
    "top_n_range",
    "--top-n-range",
    is_flag=True,
    help="Whether to consider the forward, backward, and roundtrip predictions *up to* the respective top numbers (n) in the evaluation.",
)
@click.option(
    "--isomeric-smiles/--no-isomeric-smiles",
    "isomeric_smiles",
    default=True,
    help="Whether to evaluate stereochemistry.",
)
def main(
    project_path: str,
    n_best_fw: int,
    n_best_bw: int,
    n_best_rtr: int,
    top_n_fw: Tuple[int, ...],
    top_n_bw: Tuple[int, ...],
    top_n_rtr: Tuple[int, ...],
    name: str,
    top_n_range: bool,
    isomeric_smiles: bool,
):
    """The main function to evaluate biocatalysis model data."""

    #
    # Load the data supplied by options
    #

    base_path = Path(project_path)

    max_top_n_fw = max(top_n_fw)
    max_top_n_bw = max(top_n_bw)
    max_top_n_rtr = max(top_n_rtr)

    if top_n_range:
        top_n_fw = tuple(range(1, max_top_n_fw + 1))
        top_n_bw = tuple(range(1, max_top_n_bw + 1))
        top_n_rtr = tuple(range(1, max_top_n_rtr + 1))

    data = {}

    data["src-test"] = [line.strip() for line in open(Path(base_path, "src-test.txt"))]
    data["tgt-test"] = [line.strip() for line in open(Path(base_path, "tgt-test.txt"))]
    data["src-pred"] = [line.strip() for line in open(Path(base_path, "src-pred.txt"))]
    data["tgt-pred"] = [line.strip() for line in open(Path(base_path, "tgt-pred.txt"))]

    if Path(base_path, "tgt-pred-rtrp.txt").exists():
        data["tgt-pred-rtr"] = [
            line.strip() for line in open(Path(base_path, "tgt-pred-rtrp.txt"))
        ]

    rxns_ground_truth = []
    for i in trange(len(data["src-test"])):
        rxns_ground_truth.append(
            EnzymaticReaction(
                detokenize_enzymatic_reaction_smiles(
                    f"{data['src-test'][i]}>>{data['tgt-test'][i]}"
                ),
                isomericSmiles=isomeric_smiles,
            )
        )

    rxns_forward_pred = []
    for i in trange(len(data["src-test"])):
        offset = i * n_best_fw
        candidates = []
        for j in range(offset, offset + max_top_n_fw):
            try:
                candidates.append(
                    EnzymaticReaction(
                        detokenize_enzymatic_reaction_smiles(
                            f"{data['src-test'][i]}>>{data['tgt-pred'][j]}"
                        ),
                        isomericSmiles=isomeric_smiles,
                    )
                )
            except:
                candidates.append(EnzymaticReaction("C>>C"))

        rxns_forward_pred.append(candidates)

    rxns_backward_pred = []
    for i in trange(len(data["tgt-test"])):
        offset = i * n_best_bw
        candidates = []
        for j in range(offset, offset + max_top_n_bw):
            try:
                candidates.append(
                    EnzymaticReaction(
                        detokenize_enzymatic_reaction_smiles(
                            f"{data['src-pred'][j]}>>{data['tgt-test'][i]}"
                        ),
                        isomericSmiles=isomeric_smiles,
                    )
                )
            except:
                candidates.append(EnzymaticReaction("C>>C"))

        rxns_backward_pred.append(candidates)

    rxns_roundtrip_pred = []

    if Path(base_path, "tgt-pred-rtrp.txt").exists():
        for i in trange(len(data["tgt-test"])):
            offset = i * n_best_bw
            candidates = []
            for j in range(offset, offset + max_top_n_rtr):
                try:
                    candidates.append(
                        EnzymaticReaction(
                            detokenize_enzymatic_reaction_smiles(
                                f"{data['src-pred'][j]}>>{data['tgt-pred-rtr'][i * n_best_rtr * n_best_bw]}".replace(
                                    "\\\\", "\\"
                                )  # for some reason this is escaped here. Windows problem? Check on *nix
                            ),
                            isomericSmiles=isomeric_smiles,
                        )
                    )
                except:
                    candidates.append(EnzymaticReaction("C>>C"))

            rxns_roundtrip_pred.append(candidates)

    #
    # Calculate accuracies
    #

    results = []

    print("Assessing forward accuracy...")
    for top_n in top_n_fw:
        acc_forward, correct, incorrect = get_accuracies(
            rxns_ground_truth, rxns_forward_pred, 1, top_n
        )

        with open(Path(base_path, f"correct_fw_{top_n}.txt"), "w+") as f:
            f.write("\n".join(correct))

        with open(Path(base_path, f"incorrect_fw_{top_n}.txt"), "w+") as f:
            f.write("\n".join([gt_pred[0] + "," + gt_pred[1] for gt_pred in incorrect]))

        for key, value in acc_forward.items():
            results.append(
                {
                    "metric": "acc",
                    "type": "forward",
                    "top": top_n,
                    "ec": key,
                    "value": value,
                }
            )
    print("Done.\n")

    print("Assessing backward accuracy...")
    for top_n in top_n_bw:
        acc_backward, correct, incorrect = get_accuracies(
            rxns_ground_truth, rxns_backward_pred, 1, top_n
        )

        with open(Path(base_path, f"correct_bw_{top_n}.txt"), "w+") as f:
            f.write("\n".join(correct))

        with open(Path(base_path, f"incorrect_bw_{top_n}.txt"), "w+") as f:
            f.write("\n".join([gt_pred[0] + "," + gt_pred[1] for gt_pred in incorrect]))

        for key, value in acc_backward.items():
            results.append(
                {
                    "metric": "acc",
                    "type": "backward",
                    "top": top_n,
                    "ec": key,
                    "value": value,
                }
            )
    print("Done.\n")

    if Path(base_path, "tgt-pred-rtrp.txt").exists():
        print("Assessing roundtrip accuracy...")
        for top_n in top_n_rtr:
            acc_roundtrip, correct, incorrect = get_accuracies(
                rxns_ground_truth, rxns_roundtrip_pred, 1, top_n
            )

            with open(Path(base_path, f"correct_rtrp_{top_n}.txt"), "w+") as f:
                f.write("\n".join(correct))

            with open(Path(base_path, f"incorrect_rtrp_{top_n}.txt"), "w+") as f:
                f.write(
                    "\n".join([gt_pred[0] + "," + gt_pred[1] for gt_pred in incorrect])
                )

            for key, value in acc_roundtrip.items():
                results.append(
                    {
                        "metric": "acc",
                        "type": "roundtrip",
                        "top": top_n,
                        "ec": key,
                        "value": value,
                    }
                )
        print("Done.\n")

    print("Assessing EC accuracy...")
    for top_n in top_n_bw:
        acc_ec_pred, correct, incorrect = get_ec_prediction_accuracy(
            rxns_ground_truth, rxns_backward_pred, 3, 1, top_n
        )

        with open(Path(base_path, f"correct_ec_{top_n}.txt"), "w+") as f:
            f.write("\n".join(correct))

        with open(Path(base_path, f"incorrect_ec_{top_n}.txt"), "w+") as f:
            f.write("\n".join([gt_pred[0] + "," + gt_pred[1] for gt_pred in incorrect]))

        for key, value in acc_ec_pred.items():
            results.append(
                {
                    "metric": "acc",
                    "type": "ec",
                    "top": top_n,
                    "ec": key,
                    "value": value,
                }
            )
    print("Done.\n")

    df = pd.DataFrame(results)

    if name:
        df.to_csv(Path(base_path, f"{name}.csv"), index=False)


if __name__ == "__main__":
    main()
