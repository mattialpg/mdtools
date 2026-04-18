#!/home/mattia/.miniforge/envs/unimol/bin/python
import json
import sys
from pathlib import Path

from rdkit import Chem
from unipka import UnipKa

def build_payload(smiles, ph):
    model = UnipKa()
    df = model.get_distribution(smiles, pH=ph)
    acidic = model.get_acidic_macro_pka(smiles)
    basic = model.get_basic_macro_pka(smiles)

    if (
        df.empty
        or "smiles" not in df.columns
        or "population" not in df.columns
        or "charge" not in df.columns
    ):
        return {
            "smiles": smiles,
            "top_microstates": [],
            "acidic_macro_pka": acidic,
            "basic_macro_pka": basic,}

    ranked = df.sort_values(by="population", ascending=False)
    top = ranked.iloc[0]["smiles"]
    top_rows = ranked[["smiles", "charge", "population"]].head(5).to_dict("records")

    mol = Chem.MolFromSmiles(top)
    chosen = Chem.MolToSmiles(mol) if mol is not None else smiles
    return {
        "smiles": chosen,
        "top_microstates": top_rows,
        "acidic_macro_pka": acidic,
        "basic_macro_pka": basic,}


def write_report(payload, report_file, ph):
    rows = payload.get("top_microstates", []) or []
    acidic = payload.get("acidic_macro_pka")
    basic = payload.get("basic_macro_pka")

    lines = [f"Top microstates at pH {ph}:"]
    lines.append(f"{'smiles':<40} {'charge':>6} {'population':>14}")
    for row in rows:
        smi = str(row.get("smiles", ""))
        charge = row.get("charge", "")
        pop = row.get("population", "")
        try:
            pop_fmt = f"{float(pop):.6e}"
        except Exception:
            pop_fmt = str(pop)
        lines.append(f"{smi:<40} {str(charge):>6} {pop_fmt:>14}")

    lines.append("")
    lines.append(f"acidic macro pKa: {acidic}")
    lines.append(f"basic  macro pKa: {basic}")
    Path(report_file).write_text("\n".join(lines) + "\n")


def main():
    if len(sys.argv) < 2:
        raise SystemExit("Usage: unipka_ops.py <smiles> [ph] [report_file]")
    smiles = sys.argv[1]
    ph = float(sys.argv[2]) if len(sys.argv) > 2 else 7.4
    report_file = sys.argv[3] if len(sys.argv) > 3 else "unipka_report.txt"

    payload = build_payload(smiles, ph)
    write_report(payload, report_file, ph)
    print(json.dumps(payload))


if __name__ == "__main__":
    main()
