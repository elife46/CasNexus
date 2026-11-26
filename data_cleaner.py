import pandas as pd
import numpy as np

def veri_setini_temizle(excel_dosya_yolu):
    print("1. Veri temizleme işlemi başlatılıyor...")
    
    # ---------------------------------------------------------
    # BÖLÜM 1: CHANGE-SEQ (IN VITRO) RİSK VERİSİNİN TEMİZLENMESİ
    # (Grid Search / Meta-regresyon eğitim veriniz)
    # ---------------------------------------------------------
    print("Tablo 3 (CHANGE-seq) okunuyor...")
    df_change = pd.read_excel(excel_dosya_yolu, sheet_name='CHANGE-seq_Supp_Table_3')
    
    # 1. Gereksiz Sütunları Çöpe Atmak
    istenen_sutunlar = ['name', 'target', 'offtarget_sequence', 'distance', 'CHANGEseq_reads']
    df_temiz = df_change[istenen_sutunlar].copy()
    
    # 2. Eksik Verileri Silme ve Metin Standardizasyonu (Boşlukları alma, büyük harf)
    df_temiz = df_temiz.dropna(subset=['offtarget_sequence', 'target'])
    df_temiz['offtarget_sequence'] = df_temiz['offtarget_sequence'].astype(str).str.upper().str.strip()
    df_temiz['target'] = df_temiz['target'].astype(str).str.upper().str.strip()
    
    # 3. Gerçek PAM Dizisini Bulma (distance == 0 kuralı)
    # Hata sayısının 0 olduğu satırlar aslında gerçek on-target bölgesi ve gerçek PAM'ı içerir.
    gercek_pam_df = df_temiz[df_temiz['distance'] == 0][['name', 'offtarget_sequence']]
    
    # Her bir sgRNA (name) için o 23 harfli gerçek dizilimi bir sözlüğe kaydediyoruz
    pam_sozlugu = dict(zip(gercek_pam_df['name'], gercek_pam_df['offtarget_sequence']))
    
    # Bu sözlüğü ana tabloya yeni bir sütun olarak uyguluyoruz
    df_temiz['Gercek_23nt_Hedef'] = df_temiz['name'].map(pam_sozlugu)
    
    # 4. Kayıp Verileri "0" Yapmaktan Vazgeçmek
    # Eğer bir sgRNA'nın deneyde distance=0 olan gerçek okuması yoksa, uydurma veri girmiyoruz, siliyoruz.
    df_temiz = df_temiz.dropna(subset=['Gercek_23nt_Hedef'])
    
    # Metadata (Kimlik) sütunlarını ekleyelim
    df_temiz['Kaynak_Makale'] = 'Lazzarotto_2020'
    df_temiz['Deney_Turu'] = 'in_vitro'
    df_temiz['Editor'] = 'SpCas9'
    
    # Temiz CHANGE-seq verisini kaydet
    df_temiz.to_csv("casnexus_temiz_changeseq_egitim.csv", index=False)
    print("-> CHANGE-seq verisi temizlendi, gerçek PAM'lar eklendi ve CSV olarak kaydedildi.")

    
    # ---------------------------------------------------------
    # BÖLÜM 2: TARGETED SEQUENCING (RNP) DOĞRULAMA VERİSİNİN TEMİZLENMESİ
    # ---------------------------------------------------------
    print("\nTablo 7 (Off-target Doğrulama) okunuyor...")
    df_rnp = pd.read_excel(excel_dosya_yolu, sheet_name='CHANGE-seq_Supp_Table_7')
    
    # İhtiyacımız olan sütunları seçiyoruz (Kontrol ve dsODN sütunlarını çöpe atıyoruz)
    rnp_sutunlar = [
        'Target Name', 'Off-target Sequence',
        'Indel_RNP_1', 'Total_RNP_1',
        'Indel_RNP_2', 'Total_RNP_2',
        'Indel_RNP_3', 'Total_RNP_3'
    ]
    df_rnp_temiz = df_rnp[rnp_sutunlar].copy()
    
    # Eksik/Hatalı okumaları siliyoruz
    df_rnp_temiz = df_rnp_temiz.dropna()
    
    # Metin standardizasyonu
    df_rnp_temiz['Off-target Sequence'] = df_rnp_temiz['Off-target Sequence'].astype(str).str.upper().str.strip()
    
    # 5. 3 Tekrarlı Deneyleri Tek Bir Skora İndirgeme
    toplam_indel = df_rnp_temiz['Indel_RNP_1'] + df_rnp_temiz['Indel_RNP_2'] + df_rnp_temiz['Indel_RNP_3']
    toplam_okuma = df_rnp_temiz['Total_RNP_1'] + df_rnp_temiz['Total_RNP_2'] + df_rnp_temiz['Total_RNP_3']
    
    # (Indel / Total) * 100 formülüyle gerçek biyolojik risk/kesim yüzdesini hesaplıyoruz
    df_rnp_temiz['Gercek_Kesim_Yuzdesi'] = (toplam_indel / toplam_okuma) * 100
    
    # Modeli test ederken sadece bu 3 sütuna ihtiyacımız var
    nihai_test_sutunlari = ['Target Name', 'Off-target Sequence', 'Gercek_Kesim_Yuzdesi']
    df_rnp_temiz[nihai_test_sutunlari].to_csv("casnexus_temiz_rnp_dogrulama.csv", index=False)
    print("-> RNP doğrulama verisindeki 3 tekrar birleştirildi, yüzdeler hesaplandı ve CSV olarak kaydedildi.")
    
    print("\nİŞLEM BAŞARILI! Modeliniz için gereken 2 temiz CSV dosyası oluşturuldu.")

dosya_yolu = r"C:\Users\silah\OneDrive\Masaüstü\CasNexus\41587_2020_555_MOESM3_ESM.xlsx"

if __name__ == "__main__":
    try:
        veri_setini_temizle(dosya_yolu)
    except FileNotFoundError:
        print(f"Hata: '{dosya_yolu}' dosyası bulunamadı. Lütfen dosya adının ve yolunun doğru olduğundan emin olun.")