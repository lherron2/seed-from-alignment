"""CaCoFold prediction pipeline using Infernal/R-scape.

This expects local Rfam covariance models and Infernal binaries.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

from .parse_dotbracket import parse_dotbracket


def run(cmd: list[str], cwd: Path | None = None) -> None:
    subprocess.run(cmd, check=True, cwd=cwd)


def parse_fasta_id(fasta_path: Path) -> str:
    line = fasta_path.read_text().splitlines()[0]
    return line[1:].strip().split()[0] if line.startswith(">") else line.strip()


def parse_tblout_top_hit(tblout_path: Path) -> str | None:
    for line in tblout_path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        return parts[1]
    return None


def find_cacofold_sto(out_dir: Path) -> Path | None:
    for path in out_dir.glob("*.cacofold.sto"):
        return path
    for path in out_dir.glob("*.sto"):
        if "cacofold" in path.name:
            return path
    return None


def extract_ss_cons(sto_path: Path, seq_id: str) -> tuple[str, str]:
    seq = ""
    ss_cons = ""
    with sto_path.open() as fh:
        for line in fh:
            if line.startswith("#=GC SS_cons"):
                ss_cons += line.split()[2].strip()
            elif line.strip().startswith(seq_id):
                parts = line.split()
                if len(parts) >= 2:
                    seq += parts[1].strip()
    # Ungap
    ungapped_seq = []
    ungapped_ss = []
    for s, b in zip(seq, ss_cons):
        if s == "-":
            continue
        ungapped_seq.append(s)
        ungapped_ss.append(b)
    return "".join(ungapped_seq), "".join(ungapped_ss)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Infernal+CaCoFold and emit predictions.db")
    parser.add_argument("--cm-db", required=True, help="Path to Rfam covariance model database")
    parser.add_argument("--rscape", default="rscape", help="Path to R-scape binary")
    parser.add_argument("--infernal-bin", default=None, help="Path to Infernal bin dir")
    parser.add_argument("fasta", help="Input FASTA")
    parser.add_argument("outdir", help="Output dir")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fasta = Path(args.fasta)
    seq_id = parse_fasta_id(fasta)

    cmscan = "cmscan"
    cmfetch = "cmfetch"
    cmalign = "cmalign"
    if args.infernal_bin:
        cmscan = str(Path(args.infernal_bin) / "cmscan")
        cmfetch = str(Path(args.infernal_bin) / "cmfetch")
        cmalign = str(Path(args.infernal_bin) / "cmalign")

    tblout = outdir / "cmscan.tblout"
    run([cmscan, "--tblout", str(tblout), str(args.cm_db), str(fasta)])
    model = parse_tblout_top_hit(tblout)
    if model is None:
        # No hit: emit all dots.
        seq = fasta.read_text().splitlines()[1].strip()
        (outdir / "predictions.db").write_text(seq + "\n" + ("." * len(seq)) + "\n")
        return

    cm_path = outdir / f"{model}.cm"
    with cm_path.open("w") as fh:
        subprocess.run([cmfetch, str(args.cm_db), model], check=True, stdout=fh)

    aln_sto = outdir / "aln.sto"
    with aln_sto.open("w") as fh:
        subprocess.run([cmalign, str(cm_path), str(fasta)], check=True, stdout=fh)

    run([args.rscape, "--cacofold", "--nofigures", "--outdir", str(outdir), str(aln_sto)])

    cacofold_sto = find_cacofold_sto(outdir)
    if cacofold_sto is None:
        raise SystemExit("No CaCoFold .sto produced")

    seq, ss = extract_ss_cons(cacofold_sto, seq_id)
    if not seq or not ss:
        seq = fasta.read_text().splitlines()[1].strip()
        ss = "." * len(seq)

    # Write predictions.db (sequence + structure)
    (outdir / "predictions.db").write_text(seq + "\n" + ss + "\n")


if __name__ == "__main__":
    main()
