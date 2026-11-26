import pandas as pd

# Önceki adımdaki temiz CSV dosyamızı okuyoruz
df = pd.read_csv("casnexus_temiz_changeseq_egitim.csv")

# Sadece 23 harfli dizileri alıp tekrarları siliyoruz
sadece_diziler = df['Gercek_23nt_Hedef'].drop_duplicates()

# DeepHF'e doğrudan yapıştırılacak formatta kaydediyoruz
with open("deephf_toplu_girdi.txt", "w", encoding="utf-8") as f:
    for dizi in sadece_diziler:
        f.write(f"{dizi}\n")

print(f"Toplu girdi hazır! 'deephf_toplu_girdi.txt' içindeki {len(sadece_diziler)} diziyi kopyalayıp araca yapıştırabilirsiniz.")