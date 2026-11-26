import pandas as pd

def tablolari_birlestir():
    # 1. Elinizdeki son Cas-OFFinder + DeepHF tablosunu okuyoruz
    df_final = pd.read_csv("CRISPR_Final_Results (1).csv", sep=';')
    df_final['gRNA_sequence'] = df_final['gRNA_sequence'].astype(str).str.strip().str.upper()

    # 2. İlk başta oluşturduğumuz temiz laboratuvar (CHANGE-seq) verisini okuyoruz
    df_egitim = pd.read_csv("casnexus_temiz_changeseq_egitim.csv")
    df_egitim['Gercek_23nt_Hedef'] = df_egitim['Gercek_23nt_Hedef'].astype(str).str.strip().str.upper()

    # Sadece off-target (distance > 0) olan satırlardan gerçek deneysel okumaları ve 1/(distance+1) riskini özetliyoruz
    df_off = df_egitim[df_egitim['distance'] > 0].copy()
    df_off['Tekil_Risk'] = 1.0 / (df_off['distance'] + 1.0)

    lab_ozeti = df_off.groupby(['name', 'Gercek_23nt_Hedef']).agg(
        Deneysel_Risk_Skoru=('Tekil_Risk', 'sum'),
        Gercek_CHANGEseq_Reads=('CHANGEseq_reads', 'sum')
    ).reset_index()

    # 3. İki tabloyu 23 harfli dizi üzerinden birleştiriyoruz
    df_birlesik = pd.merge(
        lab_ozeti,
        df_final[['gRNA_sequence', 'PAM', 'Efficiency', 'Total_OffTargets']],
        left_on='Gercek_23nt_Hedef',
        right_on='gRNA_sequence',
        how='inner'
    )

    # Gereksiz tekrar eden sütunu silip kaydediyoruz
    df_birlesik.drop(columns=['Gercek_23nt_Hedef'], inplace=True)
    
    cikti_adi = "CasNexus_Optimizasyona_Hazir.csv"
    df_birlesik.to_csv(cikti_adi, index=False, sep=';')
    print(f"Başarılı! '{cikti_adi}' dosyası {len(df_birlesik)} satır olarak oluşturuldu.")
    print(df_birlesik.head())

if __name__ == "__main__":
    tablolari_birlestir()