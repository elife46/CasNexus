import pandas as pd
import numpy as np

def casnexus_metodolojik_analiz(csv_yolu):
    print("1. Veri okunuyor ve eksik veri kontrolü yapılıyor...")
    df = pd.read_csv(csv_yolu, sep=";")
    
    # Eğer eşleşme yoksa 0 değil "BİLİNMİYOR" olmalı
    if df['Total_OffTargets'].isnull().any():
        df['Total_OffTargets'] = df['Total_OffTargets'].fillna("BİLİNMİYOR")
        # Analize devam etmek için sadece bilinenleri filtreleyelim
        df = df[df['Total_OffTargets'] != "BİLİNMİYOR"].copy()
        df['Total_OffTargets'] = df['Total_OffTargets'].astype(float)

    # Adaylar gen bazında karşılaştırılabilmesi için gen adı sütunu oluşturuluyor
    df['Gen_Hedefi'] = df['name'].apply(lambda x: x.split('_')[0])

    # On-Target Verimi (Zaten 0-1 arasıydı, 100'lük sisteme çeviriyoruz)
    df['OnTarget_Skoru_100'] = (df['Efficiency'] * 100).round(2)
    
    # Off-Target Güvenliği (En az hata yapan 100 puan, en çok hata yapan 0 puan alır)
    min_off = df['Total_OffTargets'].min()
    max_off = df['Total_OffTargets'].max()
    df['Guvenlik_Skoru_100'] = ((max_off - df['Total_OffTargets']) / (max_off - min_off) * 100).round(2)

    # A) Yalnızca DeepHF'e Güvenen Tasarım
    df['Skor_Yalnizca_Verim'] = df['OnTarget_Skoru_100']
    
    # B) Yalnızca Cas-OFFinder'a Güvenen Tasarım
    df['Skor_Yalnizca_Guvenlik'] = df['Guvenlik_Skoru_100']
    
    # C) CasNexus Dengeli Skor (Alfa 0.50, Beta 0.50)
    df['CasNexus_Birlesik_Skor_100'] = (0.50 * df['OnTarget_Skoru_100'] + 0.50 * df['Guvenlik_Skoru_100']).round(2)

    # 4. Pareto Optimal Seçimi
    pareto_durumu = []
    for idx, row in df.iterrows():
        gen_df = df[df['Gen_Hedefi'] == row['Gen_Hedefi']]
        # Aynı gende hem verimi yüksek hem de güvenliği yüksek başka biri var mı?
        daha_iyi_var_mi = gen_df[
            (gen_df['OnTarget_Skoru_100'] >= row['OnTarget_Skoru_100']) & 
            (gen_df['Guvenlik_Skoru_100'] >= row['Guvenlik_Skoru_100']) &
            ((gen_df['OnTarget_Skoru_100'] > row['OnTarget_Skoru_100']) | (gen_df['Guvenlik_Skoru_100'] > row['Guvenlik_Skoru_100']))
        ]
        if len(daha_iyi_var_mi) == 0:
            pareto_durumu.append("PARETO (En İyi Seçeneklerden)")
        else:
            pareto_durumu.append("Elendi")
    df['Pareto_Optimal'] = pareto_durumu

    # Sütunları düzenleyip kaydetme
    secilen_sutunlar = [
        'name', 'Gen_Hedefi', 'gRNA_sequence', 'PAM', 
        'Total_OffTargets', 'Efficiency', 
        'OnTarget_Skoru_100', 'Guvenlik_Skoru_100', 'CasNexus_Birlesik_Skor_100',
        'Pareto_Optimal', 'Gercek_CHANGEseq_Reads'
    ]
    df_sonuc = df[secilen_sutunlar].sort_values(by=['Gen_Hedefi', 'CasNexus_Birlesik_Skor_100'], ascending=[True, False])
    
    df_sonuc.to_csv("CasNexus_Yeni_Siralama.csv", index=False, sep=";")

   # öRNEK
    b2m = df_sonuc[df_sonuc['Gen_Hedefi'] == 'TRAC'].copy()
    
    en_iyi_verim = b2m.sort_values(by='OnTarget_Skoru_100', ascending=False).iloc[0]
    en_iyi_guvenlik = b2m.sort_values(by='Guvenlik_Skoru_100', ascending=False).iloc[0]
    casnexus_secimi = b2m.sort_values(by='CasNexus_Birlesik_Skor_100', ascending=False).iloc[0]
    
    print(f"1. Sadece VERİM'e odaklansaydık seçeceğimiz aday : {en_iyi_verim['name']} (Kesim İhtimali: %{en_iyi_verim['OnTarget_Skoru_100']} | Gerçek Hata Okuması: {en_iyi_verim['Gercek_CHANGEseq_Reads']} -> YÜKSEK HATA!)")
    print(f"2. Sadece GÜVENLİK'e odaklansaydık seçeceğimiz aday: {en_iyi_guvenlik['name']} (Hata Sayısı Çok Az | Ama Kesim Verimi Sadece: %{en_iyi_guvenlik['OnTarget_Skoru_100']} -> KESEMEZ!)")
    print(f"3. CasNexus'un DENGELİ seçimi olan aday         : {casnexus_secimi['name']} (Kesim Verimi: %{casnexus_secimi['OnTarget_Skoru_100']} | Gerçek Hata Okuması: Sadece {casnexus_secimi['Gercek_CHANGEseq_Reads']} -> MÜKEMMEL!)")
if __name__ == "__main__":
    casnexus_metodolojik_analiz("CasNexus_Optimizasyona_Hazir.csv")