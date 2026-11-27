import argparse
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd


def dosya_sha256(dosya_yolu: Path) -> str:
    hasher = hashlib.sha256()
    with dosya_yolu.open("rb") as f:
        for parca in iter(lambda: f.read(8192), b""):
            hasher.update(parca)
    return hasher.hexdigest()


def normalize_et(
    df: pd.DataFrame,
    group_col: str = "name",
    signal_col: str = "CHANGEseq_reads",
    method: str = "library_sum",
    distance_col: str = "distance",
) -> pd.DataFrame:
    out = df.copy()
    out[signal_col] = pd.to_numeric(out[signal_col], errors="coerce")
    out = out.dropna(subset=[signal_col]).copy()

    if method == "library_sum":
        grup_toplam = out.groupby(group_col)[signal_col].transform("sum")
        out["gRNA_Toplam_Reads"] = grup_toplam
        out["CHANGEseq_Norm_Toplam"] = np.where(
            grup_toplam > 0, out[signal_col] / grup_toplam, np.nan
        )
        out["CHANGEseq_Norm_RPM"] = out["CHANGEseq_Norm_Toplam"] * 1_000_000.0
    elif method == "ontarget_ratio":
        if distance_col not in out.columns:
            raise ValueError(f"'{distance_col}' sütunu bulunamadı.")
        ontarget_map = (
            out[out[distance_col] == 0]
            .groupby(group_col)[signal_col]
            .first()
        )
        ref_val = out[group_col].map(ontarget_map)
        out["OnTarget_Ref_Reads"] = ref_val
        out["CHANGEseq_Norm_OnTarget"] = np.where(
            (ref_val.notna()) & (ref_val > 0), out[signal_col] / ref_val, np.nan
        )
    else:
        raise ValueError(f"Desteklenmeyen normalizasyon yöntemi: {method}")

    return out


def veri_setini_temizle(
    excel_dosya_yolu: str,
    norm_method: str = "library_sum",
    column_map: dict[str, str] | None = None,
) -> None:
    print("1. Veri temizleme ve normalizasyon işlemi başlatılıyor...")
    path_obj = Path(excel_dosya_yolu)
    dosya_hash = dosya_sha256(path_obj)

    default_map = {
        "guide_id": "name",
        "target": "target",
        "offtarget_seq": "offtarget_sequence",
        "distance": "distance",
        "raw_reads": "CHANGEseq_reads",
    }
    col_map = {**default_map, **(column_map or {})}

    # ---------------------------------------------------------
    # BÖLÜM 1: CHANGE-SEQ (IN VITRO) RİSK VERİSİNİN TEMİZLENMESİ
    # ---------------------------------------------------------
    print("Tablo 3 (CHANGE-seq) okunuyor...")
    df_change = pd.read_excel(path_obj, sheet_name="CHANGE-seq_Supp_Table_3")

    istenen_sutunlar = [
        col_map["guide_id"],
        col_map["target"],
        col_map["offtarget_seq"],
        col_map["distance"],
        col_map["raw_reads"],
    ]
    df_temiz = df_change[istenen_sutunlar].rename(
        columns={
            col_map["guide_id"]: "name",
            col_map["target"]: "target",
            col_map["offtarget_seq"]: "offtarget_sequence",
            col_map["distance"]: "distance",
            col_map["raw_reads"]: "CHANGEseq_reads",
        }
    ).copy()

    df_temiz = df_temiz.dropna(subset=["offtarget_sequence", "target"])
    df_temiz["offtarget_sequence"] = df_temiz["offtarget_sequence"].astype(str).str.upper().str.strip()
    df_temiz["target"] = df_temiz["target"].astype(str).str.upper().str.strip()
    df_temiz["distance"] = pd.to_numeric(df_temiz["distance"], errors="coerce")

    onceki_grna_sayisi = df_temiz["name"].nunique()

    # Gerçek PAM Dizisini Bulma (distance == 0 kuralı)
    gercek_pam_df = df_temiz[df_temiz["distance"] == 0][["name", "offtarget_sequence"]]
    pam_sozlugu = dict(zip(gercek_pam_df["name"], gercek_pam_df["offtarget_sequence"]))
    df_temiz["Gercek_23nt_Hedef"] = df_temiz["name"].map(pam_sozlugu)

    # distance == 0 kaydı olmayan gRNA'ları düşür ve kaybı raporla
    df_temiz = df_temiz.dropna(subset=["Gercek_23nt_Hedef"]).copy()
    sonraki_grna_sayisi = df_temiz["name"].nunique()
    print(
        f"-> distance==0 filtresi sonrası gRNA sayısı: {onceki_grna_sayisi} -> {sonraki_grna_sayisi} "
        f"(elenen: {onceki_grna_sayisi - sonraki_grna_sayisi})"
    )

    # NORMALİZASYON ENTEGRASYONU
    df_temiz = normalize_et(
        df_temiz,
        group_col="name",
        signal_col="CHANGEseq_reads",
        method=norm_method,
        distance_col="distance",
    )

    # Metadata sütunları (makale, tablo, yıl ve dosya hash)
    df_temiz["Kaynak_Makale"] = "Lazzarotto_2020"
    df_temiz["Kaynak_Tablo"] = "CHANGE-seq_Supp_Table_3"
    df_temiz["Yayin_Yili"] = 2020
    df_temiz["Dosya_SHA256"] = dosya_hash
    df_temiz["Normalizasyon_Yontemi"] = norm_method
    df_temiz["Deney_Turu"] = "in_vitro"
    df_temiz["Editor"] = "SpCas9"

    df_temiz.to_csv("casnexus_temiz_changeseq_egitim.csv", index=False)
    print("-> CHANGE-seq verisi temizlendi, normalize edildi ve CSV olarak kaydedildi.")

    # ---------------------------------------------------------
    # BÖLÜM 2: TARGETED SEQUENCING (RNP) DOĞRULAMA VERİSİNİN TEMİZLENMESİ
    # ---------------------------------------------------------
    print("\nTablo 7 (Off-target Doğrulama) okunuyor...")
    df_rnp = pd.read_excel(path_obj, sheet_name="CHANGE-seq_Supp_Table_7")

    rnp_sutunlar = [
        "Target Name",
        "Off-target Sequence",
        "Indel_RNP_1",
        "Total_RNP_1",
        "Indel_RNP_2",
        "Total_RNP_2",
        "Indel_RNP_3",
        "Total_RNP_3",
    ]
    df_rnp_temiz = df_rnp[rnp_sutunlar].dropna().copy()
    df_rnp_temiz["Off-target Sequence"] = (
        df_rnp_temiz["Off-target Sequence"].astype(str).str.upper().str.strip()
    )

    toplam_indel = (
        df_rnp_temiz["Indel_RNP_1"]
        + df_rnp_temiz["Indel_RNP_2"]
        + df_rnp_temiz["Indel_RNP_3"]
    )
    toplam_okuma = (
        df_rnp_temiz["Total_RNP_1"]
        + df_rnp_temiz["Total_RNP_2"]
        + df_rnp_temiz["Total_RNP_3"]
    )

    df_rnp_temiz["Gercek_Kesim_Yuzdesi"] = (toplam_indel / toplam_okuma) * 100
    df_rnp_temiz["Kaynak_Makale"] = "Lazzarotto_2020"
    df_rnp_temiz["Kaynak_Tablo"] = "CHANGE-seq_Supp_Table_7"
    df_rnp_temiz["Yayin_Yili"] = 2020
    df_rnp_temiz["Dosya_SHA256"] = dosya_hash

    nihai_test_sutunlari = [
        "Target Name",
        "Off-target Sequence",
        "Gercek_Kesim_Yuzdesi",
        "Kaynak_Makale",
        "Kaynak_Tablo",
        "Yayin_Yili",
        "Dosya_SHA256",
    ]
    df_rnp_temiz[nihai_test_sutunlari].to_csv(
        "casnexus_temiz_rnp_dogrulama.csv", index=False
    )
    print("-> RNP doğrulama verisi temizlendi ve CSV olarak kaydedildi.")
    print("\nİŞLEM BAŞARILI! Güncellenmiş 2 temiz CSV dosyası oluşturuldu.")


def main() -> None:
    varsayilan_yol = r"C:\Users\elife\OneDrive\Masaüstü\yeniden\CasNexus\41587_2020_555_MOESM3_ESM.xlsx"
    parser = argparse.ArgumentParser(description="CasNexus Veri Temizleme ve Normalizasyon")
    parser.add_argument("--input", default=varsayilan_yol, help="Girdi Excel dosya yolu")
    parser.add_argument(
        "--norm-method",
        default="library_sum",
        choices=["library_sum", "ontarget_ratio"],
        help="Normalizasyon yöntemi",
    )
    args = parser.parse_args()

    try:
        veri_setini_temizle(args.input, norm_method=args.norm_method)
    except FileNotFoundError:
        print(f"Hata: '{args.input}' dosyası bulunamadı. Lütfen --input yolunu kontrol edin.")


if __name__ == "__main__":
    main()
    
