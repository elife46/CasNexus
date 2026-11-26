import pandas as pd

def benzersiz_hedefleri_cikar(csv_dosya_yolu):
    print("Benzersiz hedefler tespit ediliyor...")
    
    # 1. Temizlenmiş eğitim verisini okuyoruz
    df = pd.read_csv(csv_dosya_yolu)
    
    # 2. 'name' ve 'Gercek_23nt_Hedef' sütunlarındaki tekrarları siliyoruz
    benzersiz_df = df[['name', 'Gercek_23nt_Hedef']].drop_duplicates()
    
    # 3. Sonuçları kolayca kopyalayabilmeniz için bir .txt dosyasına kaydediyoruz
    cikti_dosyasi = "benzersiz_deephf_girdileri.txt"
    
    with open(cikti_dosyasi, "w", encoding="utf-8") as f:
        f.write("DeepHF Girdileri (Gercek PAM Iceren 23nt Diziler)\n")
        f.write("-" * 55 + "\n")
        for index, row in benzersiz_df.iterrows():
            # Hem ismi hem de dizilimi araya boşluk koyarak yazdırıyoruz
            satir = f"{row['name']} : {row['Gercek_23nt_Hedef']}\n"
            f.write(satir)
            
    print(f"İşlem Başarılı! Toplam {len(benzersiz_df)} adet benzersiz hedef bulundu.")
    print(f"Diziler '{cikti_dosyasi}' dosyasına kaydedildi. Buradan kopyalayıp DeepHF'e verebilirsiniz.")

if __name__ == "__main__":
    # Bir önceki aşamada oluşturduğumuz temiz dosyanın adı
    dosya_yolu = "casnexus_temiz_changeseq_egitim.csv" 
    
    try:
        benzersiz_hedefleri_cikar(dosya_yolu)
    except FileNotFoundError:
        print(f"Hata: '{dosya_yolu}' bulunamadı. Lütfen temizleme kodunu çalıştırdığınız klasörde olduğunuza emin olun.")