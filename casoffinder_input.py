import pandas as pd

# 3. SATIR: Parantez içinde sadece değişken adı (csv_dosya_yolu) olmalı
def casoffinder_girdi_hazirla(csv_dosya_yolu):
    # CSV dosyası noktalı virgül (;) ile ayrıldığı için sep=';' kullanıyoruz
    df = pd.read_csv(csv_dosya_yolu, sep=';')
    
    # 'gRNA_sequence' sütunundaki boşlukları temizleyip büyük harfe çeviriyoruz
    diziler_23nt = df['gRNA_sequence'].dropna().astype(str).str.strip().str.upper()
    
    # Son 3 harfi (PAM dizisini) kesip sadece ilk 20 harfi ([:20]) alıyoruz
    diziler_20nt = diziler_23nt.str[:20].drop_duplicates()
    
    cikti_dosyasi = "casoffinder_20nt_girdi.txt"
    with open(cikti_dosyasi, "w", encoding="utf-8") as f:
        for dizi in diziler_20nt:
            f.write(f"{dizi}\n")
            
    print(f"İşlem Başarılı! Toplam {len(diziler_20nt)} adet 20 harfli dizi '{cikti_dosyasi}' dosyasına kaydedildi.")
    print("\nKontrol amaçlı ilk 3 dizi:")
    for ornek in diziler_20nt.head(3):
        print(f" -> {ornek} (Uzunluk: {len(ornek)})")

if __name__ == "__main__":
    # 24. SATIR: Dosya yolunu başına 'r' koyarak buraya yazıyoruz
    dosya_adi = r"C:\Users\silah\Downloads\DeepHF_on-target_score(Sayfa1).csv"
    casoffinder_girdi_hazirla(dosya_adi)