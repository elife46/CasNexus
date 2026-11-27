import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def compute_file_sha256(file_path: Path) -> str:
    hasher = hashlib.sha256()
    with file_path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def load_table(input_path: Path, sheet_name: str | None = None, sep: str = ",") -> pd.DataFrame:
    suffix = input_path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(input_path, sheet_name=sheet_name)
    if suffix == ".tsv":
        return pd.read_csv(input_path, sep="\t")
    return pd.read_csv(input_path, sep=sep)


def normalize_table(
    df: pd.DataFrame,
    column_map: dict[str, str],
    method: str = "library_sum",
    metadata: dict[str, str] | None = None,
) -> pd.DataFrame:
    # 1. Ham sütunları standart (canonical) isimlere eşle
    rename_dict = {src: dst for dst, src in column_map.items() if src in df.columns}
    work = df.rename(columns=rename_dict).copy()

    required = ["guide_id", "raw_signal"]
    missing = [col for col in required if col not in work.columns]
    if missing:
        raise ValueError(f"Eksik zorunlu sütunlar: {missing}")

    work["guide_id"] = work["guide_id"].astype(str).str.strip()
    work["raw_signal"] = pd.to_numeric(work["raw_signal"], errors="coerce")
    work = work.dropna(subset=["guide_id", "raw_signal"]).copy()

    if "mismatch_count" in work.columns:
        work["mismatch_count"] = pd.to_numeric(work["mismatch_count"], errors="coerce")

    # 2. Normalizasyon yöntemi
    if method == "library_sum":
        group_totals = work.groupby("guide_id")["raw_signal"].transform("sum")
        work["group_total_signal"] = group_totals
        work["norm_offtarget_load"] = np.where(
            group_totals > 0, work["raw_signal"] / group_totals, np.nan
        )
        work["norm_rpm"] = work["norm_offtarget_load"] * 1_000_000.0

    elif method == "ontarget_ratio":
        if "mismatch_count" not in work.columns:
            raise ValueError("ontarget_ratio yöntemi için mismatch_count sütunu gereklidir.")
        ontarget_df = (
            work[work["mismatch_count"] == 0]
            .groupby("guide_id")["raw_signal"]
            .first()
        )
        ref_signal = work["guide_id"].map(ontarget_df)
        work["ontarget_ref_signal"] = ref_signal
        work["norm_offtarget_load"] = np.where(
            (ref_signal.notna()) & (ref_signal > 0),
            work["raw_signal"] / ref_signal,
            np.nan,
        )
        work["norm_rpm"] = work["norm_offtarget_load"] * 1_000_000.0
    else:
        raise ValueError(f"Desteklenmeyen normalizasyon yöntemi: {method}")

    # 3. Metadata ekle
    if metadata:
        for k, v in metadata.items():
            work[k] = v

    return work


def main() -> None:
    parser = argparse.ArgumentParser(description="1. Adım: CHANGE-seq / Off-target Normalizasyonu")
    parser.add_argument("--input", required=True, help="Girdi dosyası (.xlsx, .csv, .tsv)")
    parser.add_argument("--output", default="step1_normalized.csv", help="Çıktı CSV yolu")
    parser.add_argument("--sheet", default=None, help="Excel sayfa adı (ör. CHANGE-seq_Supp_Table_3)")
    parser.add_argument("--sep", default=",", help="CSV ayırıcı (varsayılan ',')")
    parser.add_argument(
        "--method",
        default="library_sum",
        choices=["library_sum", "ontarget_ratio"],
        help="Normalizasyon yöntemi",
    )
    parser.add_argument(
        "--column-map",
        default=None,
        help='JSON sözlük: \'{"guide_id":"name","raw_signal":"CHANGEseq_reads","mismatch_count":"distance"}\'',
    )
    parser.add_argument("--source-paper", default="Lazzarotto_2020", help="Kaynak makale etiketi")
    parser.add_argument("--source-table", default="CHANGE-seq_Supp_Table_3", help="Kaynak tablo adı")
    parser.add_argument("--pub-year", default="2020", help="Yayın yılı")
    args = parser.parse_args()

    input_path = Path(args.input)
    default_map = {
        "guide_id": "name",
        "target_seq": "target",
        "offtarget_seq": "offtarget_sequence",
        "mismatch_count": "distance",
        "raw_signal": "CHANGEseq_reads",
    }
    column_map = json.loads(args.column_map) if args.column_map else default_map

    df_raw = load_table(input_path, sheet_name=args.sheet, sep=args.sep)
    file_hash = compute_file_sha256(input_path)

    metadata = {
        "Kaynak_Makale": args.source_paper,
        "Kaynak_Tablo": args.source_table,
        "Yayin_Yili": str(args.pub_year),
        "Kaynak_Dosya": input_path.name,
        "Dosya_SHA256": file_hash,
        "Normalizasyon_Yontemi": args.method,
    }

    df_norm = normalize_table(df_raw, column_map=column_map, method=args.method, metadata=metadata)
    output_path = Path(args.output)
    df_norm.to_csv(output_path, index=False)

    print(f"Tamamlandı: {len(df_norm)} satır işlendi -> {output_path}")
    print(f"Benzersiz gRNA sayısı: {df_norm['guide_id'].nunique()}")
    print(f"NaN oran sayısı (norm_offtarget_load): {df_norm['norm_offtarget_load'].isna().sum()}")


if __name__ == "__main__":
    main()