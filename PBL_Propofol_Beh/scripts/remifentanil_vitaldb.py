# Remifentanil nas 21 induções do VitalDB selecionadas no PBL_6 (coorte da Beatriz).
#
# Calcula, por caso:
#   - quando o remifentanil começou em relação ao propofol (s; negativo = antes);
#   - a taxa média de remifentanil (µg/kg/min) nos primeiros 5 e 20 min após o início do propofol;
#   - a Cp de remifentanil do modelo de Minto no início do propofol (ng/ml);
#   - se as trilhas do VitalDB estiverem na pasta (opcional): a diferença entre a Cp de Minto,
#     calculada com a taxa registrada (Orchestra/RFTN20_RATE), e a Cp que a própria bomba
#     mostrou (Orchestra/RFTN20_CP) — mediana do erro relativo, em %.
# E imprime as medianas, que são os valores iniciais do remifentanil na plataforma.
#
# Entradas:
#   pbl6   = pasta do PBL_6 (usa dados/vitaldb/processados/coorte.csv e
#            dados/vitaldb/processados/historico_remifentanil/case_*.csv)
#   trilhas = pasta com <caso>_RFTN20_RATE.csv e <caso>_RFTN20_CP.csv baixados de
#            https://api.vitaldb.net/<tid> (tid em https://api.vitaldb.net/trks), ou "-" para pular
#   saida  = arquivo CSV de saída (dados/remifentanil_vitaldb.csv)
# Uso:  python scripts/remifentanil_vitaldb.py ../PBL_6 trilhas dados/remifentanil_vitaldb.csv
# Precisa de numpy e scipy.
import sys, csv, os, glob, statistics
import numpy as np
from scipy.linalg import expm

pbl6, trilhas, saida = sys.argv[1], sys.argv[2], sys.argv[3]
proc = os.path.join(pbl6, "dados", "vitaldb", "processados")
with open(os.path.join(proc, "coorte.csv"), encoding="utf-8") as f:
    coorte = {int(r["caseid"]): r for r in csv.DictReader(f)}

def minto(idade, peso, altura, sexo):
    # Minto 1997, como em src/remifentanil.jl (massa magra de James)
    magra = 1.1 * peso - 128 * (peso / altura) ** 2 if sexo == "M" else 1.07 * peso - 148 * (peso / altura) ** 2
    da, dm = idade - 40, magra - 55
    return dict(V1=5.1 - 0.0201 * da + 0.072 * dm, V2=9.82 - 0.0811 * da + 0.108 * dm, V3=5.42,
                CL=2.6 - 0.0162 * da + 0.0191 * dm, Q2=2.05 - 0.0301 * da, Q3=0.076 - 0.00113 * da)

def matrizes(pk, dt):
    # solução exata do modelo de 3 compartimentos, em quantidade (µg); 1/min → 1/s
    k10, k12, k13 = pk["CL"] / pk["V1"], pk["Q2"] / pk["V1"], pk["Q3"] / pk["V1"]
    k21, k31 = pk["Q2"] / pk["V2"], pk["Q3"] / pk["V3"]
    A = np.array([[-(k10 + k12 + k13), k21, k31], [k12, -k21, 0], [k13, 0, -k31]]) / 60
    Phi = expm(A * dt)
    return Phi, np.linalg.solve(A, (Phi - np.eye(3)) @ np.array([1.0, 0, 0]))

def ler_trilha(caso, nome):
    t, v = [], []
    with open(os.path.join(trilhas, f"{caso}_RFTN20_{nome}.csv")) as f:
        next(f)
        for linha in f:
            a, b = linha.strip().split(",")
            if b != "":
                t.append(float(a)); v.append(float(b))
    return np.array(t), np.array(v)

linhas = []
for arq in sorted(glob.glob(os.path.join(proc, "historico_remifentanil", "case_*.csv"))):
    caso = int(os.path.basename(arq)[5:10]); c = coorte[caso]
    idade, peso, altura, sexo = float(c["age"]), float(c["weight_kg"]), float(c["height_cm"]), c["sex"]
    pk = minto(idade, peso, altura, sexo)
    with open(arq) as f:
        rows = list(csv.DictReader(f))
    inicio_propofol = int(float(rows[0]["propofol_onset_case_s"]))
    t = np.array([int(float(r["t_case_s"])) for r in rows])
    taxa = np.array([float(r["remi_rate_ml_h"]) if r["remi_rate_ml_h"] not in ("", "nan") else 0.0 for r in rows])  # mL/h
    # Cp de Minto no início do propofol (taxa constante em cada segundo; seringa de 20 µg/ml)
    Phi, Gam = matrizes(pk, 1.0); x = np.zeros(3); cp_inicio = None
    for ti, u in zip(t, taxa):
        if ti == inicio_propofol:
            cp_inicio = x[0] / pk["V1"]
        x = Phi @ x + Gam * (u * 20 / 3600)
    com = t[taxa > 0]
    def media(seg):
        j = (t >= inicio_propofol) & (t < inicio_propofol + seg)
        return float(taxa[j].mean() * 20 / 60 / peso)                # mL/h → µg/kg/min
    linha = dict(caso=caso, idade=idade, sexo=sexo, peso_kg=peso, altura_cm=altura,
                 inicio_remi_rel_propofol_s=int(com[0] - inicio_propofol) if len(com) else "",
                 taxa_media_5min_ugkgmin=round(media(300), 4), taxa_media_20min_ugkgmin=round(media(1200), 4),
                 cp_minto_inicio_propofol_ngml=round(float(cp_inicio), 3), erro_minto_vs_bomba_pct="")
    # comparação com a Cp da bomba, se as trilhas estiverem na pasta
    if trilhas != "-" and os.path.exists(os.path.join(trilhas, f"{caso}_RFTN20_RATE.csv")):
        tr, vr = ler_trilha(caso, "RATE"); tc, vc = ler_trilha(caso, "CP")
        dt = 0.1; Phi, Gam = matrizes(pk, dt)
        fim = min(tr[-1], tc[-1]); grade = np.arange(0, fim, dt); cp = np.zeros(len(grade))
        x = np.zeros(3); j = 0; u = 0.0
        for i, ti in enumerate(grade):
            while j < len(tr) and tr[j] <= ti:
                u = vr[j] * 20 / 3600; j += 1
            cp[i] = x[0] / pk["V1"]
            x = Phi @ x + Gam * u
        ok = (tc < fim) & (vc > 0.5)
        modelo = np.interp(tc[ok], grade, cp)
        linha["erro_minto_vs_bomba_pct"] = round(float(100 * np.median(np.abs(modelo - vc[ok]) / vc[ok])), 2)
    linhas.append(linha)

with open(saida, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(linhas[0])); w.writeheader(); w.writerows(linhas)
for k in ["inicio_remi_rel_propofol_s", "taxa_media_5min_ugkgmin", "taxa_media_20min_ugkgmin",
          "cp_minto_inicio_propofol_ngml", "erro_minto_vs_bomba_pct"]:
    v = [l[k] for l in linhas if l[k] != ""]
    if v:
        print(f"{k}: mediana {statistics.median(v)}  mínimo {min(v)}  máximo {max(v)}  (n = {len(v)})")
